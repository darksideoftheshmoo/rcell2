# A script to run rcell2's shinyApp for filtering cell data
# Set the appropriate paths and execute from the commandline (using Rscript) or from RStudio

library(rcell2)
options(shiny.reactlog = TRUE)

# Load "cell data"
cell_data <- rcell2::load_cell_data("~/Projects/Programación/Rdevel/rcell2/data/image_samples/")
cdata <- cell_data$data

# Load "pdata"
pdata <- readr::read_csv("~/Projects/Programación/Rdevel/rcell2/data/image_samples/pdata.csv")

# extract image paths and channels from cell.data
paths <- cell_data$images

# run shinyCell. Filtered output would be saved in saved_data
# A usage example follows:
saved_data <- rcell2::shinyCell(cdata, pdata, paths,
                                # filters = saved_data$filters,
                                plotType = "Dots", facet_grid_option = T, facets_scale_free = "free")
