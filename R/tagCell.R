# Added to remove NOTES from devtools:check()
# https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
# globalVariables(strsplit("channel choice counts h i level pos saved treatment w x y", " ")[[1]])
# "A horrible hack" (?)

# I have found the cause of re-rendering plots.
# clientData$output_scatterplot_width changes 10 PIXELS for no reason and triggers a re-render.
# Maybe its the scrollbar or something...

#' Filtrar cdata usando gráficos y dibujando regiones
#'
#' @param cdata dataframe of "cell data".
#' @param pdata dataframe "position data".
#' @param paths dataframe of image paths.
#' @param cell_tags list of named vectors corresponding to tag groups and tags: list(named_item1 = c(option1, option2, ...), named_item2 ...).
#' @param randomize_ucids Randomize ucid order.
#' @param tag_box_size size of the image crop in pixels (integer).
#' @param cell_resize resize of the image crop in pixels (integer).
#' @param tag_channels_select a vector giving names for the image channels: c("BF", "YFP.out", etc....).
#' @param equalize_images Use magick's function to "equalize" the images.
#' @param normalize_images Use magick's function to "normalize" the images.
#' @param seed seed for random sampling of images.
#' @param tmp_output_file File path into which tagging information will be dumped by user request. NULL by default, to automatically create and append to a tmp file.
#' @param tag_ggplot a ggplot object to display in the second tab, may be used for something someday.
#' @param max.frames Max number of t.frames to render in the cell strip. Set to 0 to disable.
#' @param tags.df Previous tag dataframe, used to restore or view previous tags in the app (restores tags that are named in the cell_tags list).
#' @param verbose Print debugging messages (with levels at either 0, 1 or 2).
# @param n_max max number of boxes in the image.
# @param ... extra arguments, not used.
#' @return Lots of stuff.
#' @examples
#' path <- "/mac/apesta/trololololol/"
#' 
#' cell.data <- rcell2::cell.load.alt(path = path)
#' 
#' image.paths <- cell.data$d.paths  # Si usaste load_cell es: image.paths <- rcell2::magickPaths(cell.data)
#' 
#' pdata <- read_tsv(paste0(path, "pdata.csv"))
#' 
#' cdata <- left_join(cell.data$d, pdata)
#' 
#' p <- ggplot() + 
#'   geom_line(aes(x=t.frame, y=cf.y, group=ucid))
#' 
#' tag_channels_select <- c("BF", "BF.out", "YFP", "YFP.out")
#' 
#' saved <- rcell2::tagCell(cdata,
#'                          pdata, 
#'                          image.paths,
#'                          cell_tags = list(far1_drop = c(TRUE,
#'                                                         FALSE),
#'                                           budding =   c("emergence",
#'                                                         "division", 
#'                                                         "shmoo_o_algo"),
#'                                           artifact =  c("segmentation",
#'                                                         "crowding",
#'                                                         "out_of_focus",
#'                                                         "interesante",
#'                                                         "death",
#'                                                         "flown_away",
#'                                                         "not_a_cell")
#'                          ),
#'                          tag_channels_select = tag_channels_select,
#'                          equalize_images = T,
#'                          normalize_images = F,
#'                          n_max = 50,
#'                          tag_box_size = 75,
#'                          cell_resize = 300,
#'                          tag_ggplot = p,
#'                          tmp_output_file = "../output/annotations/progress.csv", 
#'                          debug_messages = F
#'                          )
#'                          
#' @import shiny ggplot2 magick keys
#' @importFrom grDevices rgb
#' @importFrom utils head
#' @export
tagCell <- function(cdata,
                    pdata,
                    paths,
                    cell_tags,
                    randomize_ucids = FALSE,
                    tag_box_size = 50,
                    cell_resize=NULL,
                    tag_channels_select=c("BF", "BF.out"),
                    # n_max=10,
                    seed = 1,
                    tmp_output_file=NULL,
                    tag_ggplot = NULL,
                    equalize_images = F,
                    normalize_images = F,
                    # prev.annot.df=NULL,  # TO-DO: implement resume annotations
                    max.frames=10,
                    tags.df=NULL,
                    verbose=0
                    ){
  
  # To-do
  # Invalid input$facet generates warnings and errors, this should be handled. Also, only "~", "." and "+" are handled in forumlas.
  # Implement more-than-2 variable faceting. The third and ith faceting variables of the brush are stored in "panelvar3" and so on (?)
  # Integrate polygon filter functionality, currently the drawn polygons do nothing (except show up).
  
  # Hotkeys (for the keys package)
  hotkeys <- c(
    "left", 
    "right",
    "shift+left", 
    "shift+right"
  )
  
  # Debug message level
  stopifnot(verbose %in% 0:2)
  switch (verbose+1,
          {
            runtime_messages <- F
            debug_messages <- F
          },
          {
            runtime_messages <- T
            debug_messages <- F
          },
          {
            runtime_messages <- T
            debug_messages <- T
          }
  )
  
  # Check NAs in ucid variable
  if(any(is.na(cdata[["ucid"]]))) stop("\ntagCell: ucid variable contains NA values")
  
  # Check ucid type and convert to integer
  ucid_class_check <- class(cdata[["ucid"]])
  if(ucid_class_check != "integer"){
    
    if(ucid_class_check == "factor"){
      warning(paste("\ntagCell: cohercing factor ucid to integer type"))
      cdata <- dplyr::mutate(cdata, ucid = as.integer(as.character.factor(ucid)))
      
    } else {
      warning(paste("\ntagCell: cohercing", ucid_class_check, "ucid to integer type"))
      cdata <- mutate(cdata, ucid = as.integer(ucid))
    }
  }
  
  # Progress file
  if(is.null(tmp_output_file)){
    tmp_output_file <- tempfile(tmpdir = tempdir(), fileext = ".txt", pattern = "tagCell_progress")
  } else {
    dir.create(dirname(normalizePath(tmp_output_file)), recursive = T, showWarnings = F)
  }
  if(debug_messages) print(paste("Will append tagging progress to file:", tmp_output_file))

    
  # Setup environments for the shiny app, from this environment
  environment(tagCellServer) <- environment()
  environment(tagCellUi) <- environment()
  
  #### RUN APP ####
  saved <- shiny::runApp(list(ui = tagCellUi(), server = tagCellServer))
  
  #### RETURN RESULT ####
  # Devolver una lista con los objetos cdata cfilter y los stringFilters
  return(saved)
}
            
#' Pivot cell tags to a cdata-joinable dataframe, with one hot encoding
#' 
#' @param tags.df Output from tagCell.
#' @param exclude.cols Character vector with names of columns which should be removed from input.
#' 
#' @details 
#' Ver:
#' * `~/Projects/Academia/Doctorado/gitlabs_acl/rtcc/far1/analisis_Far1_arresto-lavado/R/analisis_pos_2_a_7_v7_tags_analysis.Rmd`
#' * https://stackoverflow.com/questions/55288338/r-multi-hot-encoding-among-multiple-columns
#' * https://stackoverflow.com/a/63454411/11524079
#' @export
tags.to.onehot <- function(tags.df, exclude.cols = c("pos", "cellID", "viewed")){
  annotations.dt <- tags.df %>% filter(!is.na(t.frame)) %>% 
    # Remove some redundant ID columns
    {.[,!names(.) %in% exclude.cols]}
  
  # Replace boolean values with strings
  # annotations.dt <- annotations #%>%
    # mutate(far1_drop = ifelse(far1_drop, "t_drop", "false_drop")) %>%
    # mutate(far1_deloc = ifelse(far1_deloc, "far1_deloc", "far1_no_deloc"))
  
  # Replace missing values with a string
  # annotations.dt[is.na(annotations.dt)] <- "not_tagged"
  
  # Convert all columns to factor type
  annotations.dt <- mutate_all(annotations.dt, factor)
  
  # Set data table
  data.table::setDT(annotations.dt)
  
  # Melt data table by ucid and t.frame
  annotations.dt <- annotations.dt %>% data.table::melt(id.vars = c("ucid","t.frame"))
  
  # Paste variable names with their values.
  # These end up as column names after casting (see below).
  annotations.dt <- annotations.dt %>% 
    unite(value, variable, value, sep = ".") %>% 
    mutate_all(factor)
  
  # Esta funcion va a calcular el valor de una celda solo cuando
  # hay múltiples valores para ella en la tabla original (esos vienen de columna "ind", ver más abajo).
  # Se usa entonces cuando las filas no son identificadas únicamente por las columnas de ID (aunque quizás se use siempre :shrug:).
  aggr.fun <- function(x) length(x) > 0
  
  annotations.dt <- data.table::dcast(
    # Add "ind" column, filled with ones.
    data.table::setDT(annotations.dt)[,ind:=1],  # [1,],  # [c(1,1,2,3),],
    # ?
    fun.aggregate	= aggr.fun,
    # Specify what columns are identifiers (LHS),
    # and which columns have values that will be cast to columns (RHS).
    ucid+t.frame~value,
    # The column holding the values that should be cast
    value.var = "ind",
    # The fill value for missing values for a combination if ID columns.
    fill=FALSE
  ) %>%
    # Convert ucid and t.frame factors to characters
    mutate(ucid = as.character.factor(ucid),
           t.frame = as.character.factor(t.frame)) %>%
    # Then convert everything to numeric type
    mutate_all(as.numeric)
  
  return(annotations.dt)
}
