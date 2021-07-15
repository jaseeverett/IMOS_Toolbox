suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(lutz)
  library(data.table)
})

## Functions for operating
untibble <- function (tibble) {
  data.frame(unclass(tibble), check.names = FALSE, stringsAsFactors = FALSE)
}  ## escape the nonsense

raw <- "https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/master/Plankton/RawData/"
output <- "https://raw.githubusercontent.com/PlanktonTeam/IMOS_Toolbox/master/Plankton/Output/"


# Bring in all sample data
getLFTrips <- function(){
  LFSamp <- read_csv(paste0(raw, "BGC_LFish_Samples.csv"), na = "",
                     col_types = cols(FLAG_COMMENT = col_character())) %>% 
    rename(i_Sample = I_SAMPLE_ID, TripCode = TRIP_CODE, Station = STATIONNAME, Latitude = LATITUDE, Longitude = LONGITUDE, SampleDateLocal = SAMPLEDATE_LOCAL,  
           ProjectName = PROJECTNAME, Volume_m3 = VOLUME_M3, Vessel = VESSEL, TowType = TOWTYPE, GearDepth_m = GEARDEPTH_M,
           GearMesh_um = GEARMESH_UM, WaterDepth_m = BATHYM_M, Temp_DegC = TEMPERATURE_C, Salinity = SALINITY,
           Comments = COMMENTS, QC_Flag = QC_FLAG, FlagComments = FLAG_COMMENT) %>%
    mutate(Year = year(SampleDateLocal),
           Month = month(SampleDateLocal),
           Day = day(SampleDateLocal),
           Time_24hr = str_sub(SampleDateLocal, -8, -1), # hms doesn"t seem to work on 00:00:00 times
           tz = tz_lookup_coords(Latitude, Longitude, method = "fast"),
           SampleDateUTC = with_tz(force_tzs(SampleDateLocal, tz, roll = TRUE), "UTC")) %>%
    select(i_Sample:SampleDateLocal, Year:SampleDateLocal, Latitude:FlagComments)
  return(LFSamp)
}

# Bring in larval fish data
getLFData <- function(){
  LFData <- fread(paste0(raw, "BGC_LFish_CountRaw.csv"), na = "") %>% 
    rename(i_Sample = I_SAMPLE_ID, TripCode = TRIP_CODE, ScientificName = SCIENTIFICNAME, SPCode = SPCODE,
           Taxon_Count = TAXON_COUNT, Comments = COMMENTS) 
  return(LFData)
}

# Make count data of all larval fish
getLFCountAll <- function(){
  LFCount <- LFSamp %>% left_join(LFData %>% select(-Comments), by = c("i_Sample", "TripCode")) %>%
  mutate(Header = paste(ScientificName, SPCode, sep = " ")) %>%
  select(-ScientificName, -SPCode) %>%
  arrange(Header) %>%
  pivot_wider(names_from = Header, values_from = Taxon_Count, values_fill = 0) %>%
  arrange(SampleDateLocal)
  return(LFCount)
}

LFSamp <- getLFTrips()
LFData <- getLFData()
LFCount <- getLFCountAll()

fwrite(LFCount, "Output/LFCount_All.csv")

# Make count data of all larval fish
getLFCountBGC <- function(){
  LFCountBGC <- LFSamp %>% 
    filter(grepl('IMOS', ProjectName)) %>%
    left_join(LFData %>% select(-Comments), by = c("i_Sample", "TripCode")) %>%
    mutate(Header = paste(ScientificName, SPCode, sep = " ")) %>%
    select(-c(ScientificName, SPCode, Temp_DegC, Salinity)) %>%
    arrange(Header) %>%
    pivot_wider(names_from = Header, values_from = Taxon_Count, values_fill = 0) %>%
    arrange(SampleDateLocal)
  return(LFCountBGC)
}

LFCountBGC <- getLFCountBGC()

fwrite(LFCountBGC, "Output/LFCount_BGC.csv") 
