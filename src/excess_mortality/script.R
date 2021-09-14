RhpcBLASctl::blas_set_num_threads(1L)
RhpcBLASctl::omp_set_num_threads(1L)

system(paste0("echo Excess Fit for  ",iso3c))
date_0 <- as.Date(date)

version_min <- "0.6.9"
if(packageVersion("squire") < version_min) {
  stop("squire needs to be updated to at least v", version_min)
}

version_min <- "0.1.22"
if(packageVersion("nimue") < version_min) {
  stop("nimue needs to be updated to at least v", version_min)
}

## -----------------------------------------------------------------------------
## 1. GET INPUT DATA
## -----------------------------------------------------------------------------

## a. Get data from local files
## -----------------------------------------------------------------------------

# get data from economist script
data <- readRDS("excess_deaths.Rds")
data2 <- data[data$iso3c == iso3c, ]

## b. Sort out what is to be our death time series
## -----------------------------------------------------------------------------

# here I have just taken the maximum of either excess or covid on each day.
# but this is probably not the best idea as it likely overestimates covid deaths

# I'm not sure this does over estimate, if excess mortality is disconnected to
# covid deaths then we are essentially using reported deaths (an under-estimate).
# Only case I can think of for over-estimate is when excess deaths spikes due to
# lack of treatment etc, whilst covid deaths are well reported and tested for.
df <- data2
df$deaths <- as.integer(df$deaths)

#if no data or no deaths skip over and make empty file
if(nrow(df) == 0 | sum(df$deaths) == 0){
  #save empty files
  saveRDS(NULL, "res.rds")
  ggsave("fitting.pdf",width=12, height=12,NULL)

} else{
  #load par_inits
  pars_init_prev <- readRDS("pars_init.rds")

  #also check if start date is compatible with par_init since start date can
  #change with excess estimates
  if(!is.null(pars_init_prev[[iso3c]]$start_date)){
    data_start_date <- df %>% filter(deaths > 0) %>% filter(date == min(date)) %>% pull(date)
    par_start_date <- pars_init_prev[[iso3c]]$start_date %>%
      as.Date()
    if(
      par_start_date > data_start_date - 10 |
      par_start_date < data_start_date - 55
    ){
      #remove start date from par_init
      pars_init_prev[[iso3c]]$start_date <- NULL
    }
  }

  ## c. Any other parameters needed to be worked out
  ## -----------------------------------------------------------------------------

  # check that we have this iso3c in squire
  if(!(iso3c %in% squire::population$iso3c)) {
    stop("iso3c not found in squire")
  }
  country <- squire::population$country[match(iso3c, squire::population$iso3c)]
  pop <- squire::get_population(country)$n

  if(short_run) {
    replicates <- 2
    n_mcmc <- 20
    n_chains <- 1
  }


  ## -----------------------------------------------------------------------------
  ## 2. Delta Adjustments
  ## -----------------------------------------------------------------------------

  if(adjust_delta){
    #open data from covariants
    delta_characteristics <- readRDS("delta_characteristics.Rds") %>% ungroup()

    #get data or impute if not there
    if(iso3c %in% delta_characteristics$iso3c){
      this_iso3c <- iso3c
      delta_characteristics <- delta_characteristics %>%
        filter(iso3c == this_iso3c) %>%
        select(where(~is.numeric(.x)|is.Date(.x)))
    } else if(
      countrycode::countrycode(iso3c, origin = "iso3c", destination = "un.regionsub.name") %in%
      delta_characteristics$sub_region
    ){
      this_sub_region <- countrycode::countrycode(iso3c, origin = "iso3c", destination = "un.regionsub.name")
      #we then use the median values of all countries in that sub region
      delta_characteristics <- delta_characteristics %>%
        filter(
          sub_region == this_sub_region
        ) %>%
        summarise(across(
          where(~is.numeric(.x)|is.Date(.x)),
          ~median(.x, na.rm = T)
        ))
    } else if(
      countrycode::countrycode(iso3c, origin = "iso3c", destination = "continent") %in%
      delta_characteristics$continent
    ){
      this_continent <- countrycode::countrycode(iso3c, origin = "iso3c", destination = "continent")
      #we then use the median values of all countries in that contient
      delta_characteristics <- delta_characteristics %>%
        filter(
          continent == this_continent
        ) %>%
        summarise(across(
          where(~is.numeric(.x)|is.Date(.x)),
          ~median(.x, na.rm = T)
        ))
    } else{
      delta_characteristics <- delta_characteristics %>%
        summarise(across(
          where(~is.numeric(.x)|is.Date(.x)),
          ~median(.x, na.rm = T)
        ))
    }

    #we only use start date and immune escape atm
    delta_start_date <- as.Date(delta_characteristics$start_date)
    #assume shift takes 60 days
    days_in_shift <- 60
    #derive dur_R
    dur_R <- 1/((days_in_shift/360 - log(1-delta_characteristics$immune_escape))/days_in_shift)
    #get modifier on hospitalisation
    prob_hosp_multiplier <- delta_characteristics$prob_hosp_multiplier
  } else{
    #these settings should lead to no adjustment
    dur_R <- 365
    prob_hosp_multiplier <- 1
    delta_start_date <- NULL
    days_in_shift <- NULL
  }

  ## -----------------------------------------------------------------------------
  ## 3. Fit Model
  ## -----------------------------------------------------------------------------

  # fit model
  res <- fit_spline_rt(
    data = df,
    country = country,
    pop = pop,
    n_mcmc = as.numeric(n_mcmc),
    replicates = as.numeric(replicates),
    model = model,
    pars_obs_dur_R = dur_R,
    pars_obs_prob_hosp_multiplier = prob_hosp_multiplier,
    pars_obs_delta_start_date = delta_start_date,
    pars_obs_shift_duration = days_in_shift,
    n_chains = as.numeric(n_chains),
    pars_init_prev = pars_init_prev
  )

  ## -----------------------------------------------------------------------------
  ## 4. Summarise model for ease of viewing outputs and goodness of fit
  ## -----------------------------------------------------------------------------

  # remove the output for memory and ease
  output <- res$output
  res$output <- NULL

  # save output without output for memory
  saveRDS(res, "res.rds")
  res$output <- output

  # make a series of quick plots so we can check fits easily afterwards
  if (model == "SQUIRE") {
    rtp <- rt_plot_immunity(res)
  } else {
    rtp <- rt_plot_immunity_vaccine(res)
  }

  dp <- dp_plot(res)
  cdp <- cdp_plot(res)
  ar <- ar_plot(res)

  ggsave("fitting.pdf",width=12, height=12,
         cowplot::plot_grid(rtp$plot + ggtitle(country),
                            dp, cdp, ar, ncol = 1))
}
