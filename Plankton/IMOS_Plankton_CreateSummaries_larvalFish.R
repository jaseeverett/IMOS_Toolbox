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
