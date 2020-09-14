options(shiny.reactlog = TRUE)

# this will not work as is
d <- readRDS("data/celldata.RDS")
pdata <- readr::read_csv("data/pdata.csv")

# extract cellid dataframe
cdata <- d$data
cdata <- merge(cdata, pdata, by = "pos")

# extract image paths and channels from cell.data
paths <- d$images
paths$path <- "data/20191121_experiment/"
paths$file <- sub("//", "/", paste(paste(paths$path, "/", paths$image, sep = "")))  # make complete (not full) paths

# not really necessary I think.
channels <- d$channels

# run shinyCell. Filtered output would be saved in saved_data
# A usage example follows:
saved_data <- cellMagick::shinyCell(cdata, pdata, paths,
                                    # filters = saved_data$filters,
                                    plotType = "Hex", facet_grid_option = T, facets_scale_free = "free")

# saveRDS(saved_data, "data/saved_magick.RDS")
saved_data <- readRDS("data/saved_magick.RDS")
saved_data <- cellMagick::shinyCell(cdata, pdata, paths, #initial_facet = "strain~alphaF",
                                    filters = saved_data$filters,
                                    plotType = "Dots", facet_grid_option = T, facets_scale_free = "free")
