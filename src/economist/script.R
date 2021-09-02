date_0 <- as.Date(date, "%Y-%m-%d")

##Get most recent fit from github repo
excess_deaths <- read.csv(
  "https://raw.githubusercontent.com/mrc-ide/covid-19-the-economist-global-excess-deaths-model/recent_fit/output-data/export_country.csv"
  ) %>%
  mutate(date = as.Date(date, origin = "1970-01-01")) %>%
  filter(date <= date_0)

if(max(excess_deaths$date) < date_0){
  warning(paste0("Excess mortality data ", date_0 - max(excess_deaths$date), " out of date."))
}

saveRDS(excess_deaths, "excess_deaths.rds")
