
untibble <- function (tibble) {
  data.frame(unclass(tibble), check.names = FALSE, stringsAsFactors = FALSE)
}  ## escape the nonsense


get_NRSTrips <- function(){
  NRSTrips <- read_csv(paste0(rawD,.Platform$file.sep,"nrs_trips.csv"), na = "(null)") %>% 
    rename(Station = STATION_NAME, Latitude = Y_COORD, Longitude = X_COORD, SampleDateLocal = TRIP_START_DATETIME_LOCAL, 
           NRScode = NRS_CODE, StationDepth_m = STATION_DEPTH, SampleDepth_m = SAMPLEDEPTH_M, stateCode = STATE_CODE,
           DaylightSavings = DAYLIGHT_SAVINGS_Y_N, UTCoffsetH = UTC_OFFSET_H) %>%
    # select(-SAMPLEDEPTH_M) %>%
    mutate(Year = year(SampleDateLocal),
           Month = month(SampleDateLocal),
           Day = day(SampleDateLocal),
           Time_24hr = str_sub(SampleDateLocal, -8, -1), # hms doesn"t seem to work on 00:00:00 times
           tz = tz_lookup_coords(Latitude, Longitude, method = "fast"),
           SampleDateUTC = with_tz(force_tzs(SampleDateLocal, tz, roll = TRUE), "UTC"),
           SampleDateLocal = as.character(SampleDateLocal)) %>% 
    select(-c(tz, DaylightSavings, UTCoffsetH, stateCode)) %>%
    distinct()
  return(NRSTrips)
}

# Bring in copepod information table with sizes etc.
get_ZooInfo <- function(){
  ZInfo <- read_csv(paste0(rawD,.Platform$file.sep,"taxon_info.csv"), na = "(null)") %>% 
    dplyr::rename( "TaxonName" = "TAXON_NAME") %>% 
    untibble()
}
