library(tidyverse)
masks_tsv_path <- "~/Projects/Colman/gitlabs_acl/rtcc/far1/analisis_Far1NG_2019Dai/data/20191212_Timecourse_FAR1-NG_RtCC/full_dataset/Position01/out_all_masks.tsv"
shape_pixtype = "b"
shape_flagtype = 0

masks_coords <- readr::read_tsv(masks_tsv_path,
                                col_types = cols(cellID = readr::col_integer(),
                                                 t.frame = readr::col_integer(),
                                                 flag = readr::col_integer(),
                                                 x = readr::col_integer(),
                                                 y = readr::col_integer(),
                                                 pixtype = readr::col_factor(levels = c("b", "i")))
) %>% 
  {if(!is.null(shape_pixtype)) filter(., pixtype %in% shape_pixtype) else .} %>% 
  {if(!is.null(shape_flagtype)) filter(., flag %in% shape_flagtype) else .} %>% 
  mutate(id = factor(paste(cellID, t.frame, flag, pixtype,sep = "_")))

masks_coords2 <- readr::read_tsv(masks_tsv_path) %>% 
  {if(!is.null(shape_pixtype)) filter(., pixtype %in% shape_pixtype) else .} %>% 
  {if(!is.null(shape_flagtype)) filter(., flag %in% shape_flagtype) else .} %>% 
  mutate(id = factor(paste(cellID, t.frame, flag, pixtype,sep = "_")))

object.size(masks_coords) %>% print(units = "MiB")
object.size(masks_coords2) %>% print(units = "MiB")
