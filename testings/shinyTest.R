setwd("/home/nicomic/Projects/Rdevel/rcell2/R/")
devtools::load_all(recompile = T)

library(rcell2)
library(tidyverse)
library(profvis)
library(bench)

setwd("/home/nicomic/Projects/Colman/HD/reentry/analisis_Far1NG_2019Dai")

pdata.file <- "./data/20191212_Timecourse_FAR1-NG_RtCC/TC_Far1.csv"
pdata <- read_tsv(file = pdata.file,
                  col_types = cols(
                    col_integer(),
                    col_factor(),
                    col_double(),
                    col_integer(),
                    col_integer()
                  ),
                  locale = locale(decimal_mark = "."))

# cell.data <- rcell2::cell.load.alt(path = "./data/20191212_Timecourse_FAR1-NG_RtCC/split/")
cell.data <- rcell2::cell.load.alt(path = "./data/20191212_Timecourse_FAR1-NG_RtCC/split_vcellid/")
# cell.data <- rcell2::cell.load.alt(path = "./data/20191212_Timecourse_FAR1-NG_RtCC/split.o/")

paths.magick <- rcell2::magickPaths(cell.data)

saved_data2 <-readRDS("/tmp/saverdata2.RDS")

if(F){
  source("/home/nicomic/Projects/Rdevel/rcell2/testings/shinyTest.R")
  # debugonce(polyFilterApply)
  # debugonce(applyFilter)
  # debugonce(calculateTruth)
  # debugonce(applyFilter)
  # debugonce(magickCell)
  # debug(hover_closest)
  # profvis({
    cdata <- cell.data$d #%>% filter(ucid == 2000000281)
    saved_data3 <- cdata %>%
      shinyCell(pdata, paths.magick,
                plotType = "Dots",
                # initial_facet = "pos ~ .",
                # filters = saved_data$filters,
                facet_grid_option = TRUE,
                n_max = 4^2, max_size = 720)
    # saved_data4 <- saved_data3$cdata %>% filter(filter) %>% 
    #   shinyCell(pdata, paths.magick,
    #             plotType = "Dots",
    #             # initial_facet = "pos ~ .",
    #             # filters = saved_data3$filters,
    #             facet_grid_option = TRUE,
    #             n_max = 7^2, max_size = 720)
  # })
  cdata <- cell.data$d %>% 
    filter(ucid %in% c(
      1000000152
    ))
  
  magickCell(cdata, paths.magick)
}

