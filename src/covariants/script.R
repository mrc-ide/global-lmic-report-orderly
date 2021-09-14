date_0 <- as.Date(date, "%Y-%m-%d")

#read in data
CoVariants_raw <- fromJSON(
  "https://raw.githubusercontent.com/hodcroftlab/covariants/master/cluster_tables/EUClusters_data.json"
)

start <- TRUE
#only interested in Delta so we need dates, total sequences and delta sequences
for(country in names(CoVariants_raw$countries)){
  country_df <- CoVariants_raw$countries[[country]]
  if("21A (Delta)" %in% names(country_df)){
    country_df <- as.data.frame(CoVariants_raw$countries[[country]][c("week", "21A (Delta)", "total_sequences")])
    country_df$delta_prop <- country_df[,2]/country_df[,3]
  } else{
    country_df <- as.data.frame(CoVariants_raw$countries[[country]][c("week")])
    country_df$delta_prop <- 0
  }
  country_df <- select(country_df, week, delta_prop) %>%
    mutate(iso3c = countrycode(country, "country.name", "iso3c"),
           iso3c = if_else(is.na(iso3c), country, iso3c))
  if(start){
    CoVariants_df <- country_df
    start <- FALSE
  } else{
    CoVariants_df <- rbind(CoVariants_df, country_df)
  }
}

#check countries whose iso3c cannot be matched
CoVariants_df %>%
  filter(nchar(iso3c) != 3) %>%
  pull(iso3c) %>%
  unique()
#an overseas territory, so we'll remove it
CoVariants_df <- filter(CoVariants_df, nchar(iso3c) == 3) %>%
  mutate(week = as.Date(week))

#examine
#ggplot2::ggplot(CoVariants_df) + ggplot2::geom_line(ggplot2::aes(x = week, y = delta_prop, colour = iso3c))

#data frame to hold shift characteristics
delta_characteristics <- data.frame(iso3c = unique(CoVariants_df$iso3c))

#calculate shift start (when delta >10% of cases) and end dates (when delta >90% of cases)
delta_characteristics <-
  CoVariants_df %>%
  group_by(iso3c) %>%
  arrange(week) %>%
  #only keep values where delta > 10% of sequences
  filter(delta_prop > 0.10) %>%
  filter(week == min(week)) %>%
  rename(start_date = week) %>%
  select(iso3c, start_date) %>%
  right_join(
    delta_characteristics
  )
delta_characteristics <-
  CoVariants_df %>%
  group_by(iso3c) %>%
  arrange(week) %>%
  #only keep values where delta > 10% of sequences
  filter(delta_prop > 0.90) %>%
  filter(week == min(week)) %>%
  rename(end_date = week) %>%
  select(iso3c, end_date) %>%
  right_join(
    delta_characteristics
  )

#only keep countries with start dates, and add un region/sub region for data
#imputation for countries with no data later on
delta_characteristics <- delta_characteristics %>%
  filter(!is.na(start_date)) %>%
  mutate(
    sub_region = countrycode(
      if_else(iso3c == "TWN", "CHN", iso3c),
      origin = "iso3c",
      destination = "un.regionsub.name"),
    continent = countrycode(iso3c, origin = "iso3c", destination = "continent")
  )

#check if any end dates are before the start dates
if(length(filter(delta_characteristics, start_date > end_date) %>%
          pull(iso3c)) > 0){
  #we just swap these around
  delta_characteristics <- delta_characteristics %>%
    mutate(
      start_date_old = start_date,
      end_date_old = end_date,
      start_date = min(start_date_old, end_date_old, na.rm = T),
      end_date = max(start_date_old, end_date_old, na.rm = T)
    ) %>%
    select(!c(start_date_old, end_date_old))
}
#check if the dates are equal
if(length(filter(delta_characteristics, start_date == end_date) %>%
          pull(iso3c)) > 0){
  #for these countries we'll assume the same average shift time as across the globe
  #we'll assume the start/end date is the middle
  median_days <- delta_characteristics %>%
    filter(start_date < end_date) %>%
    mutate(days = as.numeric(end_date - start_date)) %>%
    pull(days) %>%
    median()

  delta_characteristics <- delta_characteristics %>%
    mutate(
      old_start = start_date,
      old_end = end_date,
      start_date = if_else(
        identical(old_start, old_end),
        old_start - round(median_days/2),
        old_start
      ),
      end_date = if_else(
        identical(old_start, old_end),
        old_end + round(median_days/2),
        old_end
      )
    ) %>%
    select(!c(old_start,old_end))
}


#assume 25% immune escape (get source?)
delta_characteristics <- delta_characteristics %>%
  mutate(
    immune_escape = 0.25
  )

#calculate require dur_R for the shift period
delta_characteristics <- delta_characteristics %>%
  mutate(
    days_in_shift = as.numeric(end_date - start_date),
    required_dur_R = 1/((days_in_shift/360 - log(1-immune_escape))/days_in_shift)
  )

#add increased hospitalization
#(https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(21)00475-8/fulltext#seccestitle150)
delta_characteristics <- delta_characteristics %>%
  mutate(
    prob_hosp_multiplier = 1.45
  )

#print a plot to check our results
dir.create("calibration")
pdf("calibration/plot.pdf")
print(
  ggplot(delta_characteristics %>% arrange(start_date) %>%
           mutate(
             end_date = if_else(
               is.na(end_date),
               date_0,
               end_date
             )
           )) +
    geom_segment(aes(
      y = fct_reorder(iso3c, start_date),
      yend = fct_reorder(iso3c, start_date),
      x = start_date,
      xend = end_date
    ),
    alpha = 0.5,
    size = 2) +
    labs(
      x = "Date",
      y = NULL,
      title = "Duration of Delta Shift"
    )
)
dev.off()

#save data
saveRDS(delta_characteristics, "delta_characteristics.rds")
