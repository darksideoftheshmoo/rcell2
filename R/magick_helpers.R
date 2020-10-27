#' Extaer paths para cellMagick del objeto cell.data
#'
#' @param cell.data directory where images are stored, full path.
#' @return A dataframe with paths.
# @examples
# cell.args <- cellArgs(path = path)
#' @export
magickPaths <- function(cell.data){
  cdata.map <- cell.data$d.map
  paths <- cdata.map %>%
    mutate(channel = paste0(toupper(channel), "FP"))
  
  paths <- bind_rows( # Bind it with itself, to get entries for BF as well.
    paths %>%
      select(pos, t.frame, channel, fluor) %>%
      rename(file = fluor),
    paths %>%
      filter(flag == 0) %>% #  Since there can be duplicates in the "bright" column, keep just one channel
      select(pos, t.frame, channel, bright) %>%
      rename(file = bright) %>%
      mutate(channel = "BF")
  ) %>%
    mutate(path = dirname(file),
           is.out = FALSE)
  
  paths <- bind_rows(paths,
                     paths %>%  # bind it with itself, out files are named exactly the same, but with an extra ".out.tif"
                       mutate(file = paste0(file, ".out.tif"),
                              channel = paste0(channel, ".out"),
                              is.out = TRUE)) %>%
    mutate(image = basename(file))
  
  return(paths)
}

#' Display an image in rmarkdown with knitr
#'
#' @param imgs a magick image or cellMagick output
#' @param .prefix A string prepended to the file name, "" by default.
#' @param .resize a cellMagick image resize string as "200x200" (default NULL, for no resizing).
#' @param .path Directory where the output should be saved.
#' @param .include Set to true to use knitr::include_graphics directly.
#' @return An path to a temporary image file.
# @examples
# cell.args <- cellArgs(path = path)
#' @export
magickForKnitr <- function(imgs, .prefix = "tile", .resize = NULL, .path = tempdir(), .include = F){
  
  dir.create(.path, recursive = T)
  
  temp <- tempfile(tmpdir = .path, fileext = ".png", pattern = .prefix)
  
  if(class(imgs) == "list") imgs <- imgs$img
  
  imgs %>% 
    {if(is.null(.resize)) . else magick::image_resize(., .resize)} %>% 
    magick::image_write(path = temp, format = "png")
  
  print(temp)
  
  if(.include) knitr::include_graphics(temp)
  
  return(temp)
}
