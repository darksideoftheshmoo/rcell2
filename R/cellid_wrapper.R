#' Correr Cell-ID desde R usando .C()
#'
#' @param args el comando de cellid entero, tal como se ejecutaria en bash "cell -p ..."
#' @param debug_flag Set to 0 to disable CellID printf messages.
# @useDynLib rcell2 CellID
#' @export
#' @return Exit code from CellID
cellid <- function(args, debug_flag=0){

  # args <- "~/Software/cellID-linux/cell -p /home/nicomic/Projects/Colman/HD/scripts/cellMagick/data/images/parameters.txt -b /tmp/Rtmp7fjlFo/file2b401093d715 -f /tmp/Rtmp7fjlFo/file2b402742f6ef -o /home/nicomic/Projects/Colman/HD/uscope/20200130_Nico_screen_act1_yfp/1/Position001/out"
  argv <- strsplit(args, " ")[[1]]
  argc <- length(argv)
  
  if(debug_flag != 0) print("Printing argv and argc before .C() call to CellID.")
  if(debug_flag != 0) print(argv)
  if(debug_flag != 0) print(argc)

  exit_code <- 0
  
  if(exit_code != 1) stop(paste("CellID is not bundled in this branch, see master_cellid", exit_code))
  
  return(exit_code)
}

#' Function to run CellID
#'
#' @param arguments An argument data.frame, as built by rcell2::cellArgs2.
#' @param cell.command Path to the binary executable (get if from https://github.com/darksideoftheshmoo/cellID-linux). Or "cellBUILTIN" for the builtin binary.
#' @param no_cores Position-wise parallelization,internally capped to number of positions in cell.args.
#' @param dry Do everything without actually running CellID, print the commands that would have been issued.
#' @param debug_flag Set to 0 to disable CellID printf messages (builtin CellID only).
#' @param label_cells_in_bf Set to TRUE to enable labeling of cells with their CellID in the BF output image (CellID option '-l').
#' @param fill_interior_pixels Set to TRUE to fill each cell interior area in the output image file with intensity-labeled pixels (CellID option '-i').
#' @param output_coords_to_tsv Set to TRUE to write cell interior and boundary pixels data to a .tsv file in the output directory (CellID option '-m').
#' @param encode_cellID_in_pixels Set to TRUE to write cell interior and boundary pixels with intensity-encoded CellIDs and blank the rest of the image (CellID option '-s').
#' @param ignore.stdout Set to FALSE to see CellID output from a system call.
#' @param intern Set to TRUE to save CellID output from a system call to a file in the "out" directories (one per position) and the commands to a file at the first "path" in the arguments data.frame.
#' @return A dataframe with one column indicating the issued commands. Use rcell2::load_cell_data to get the results from the CellID output, typically located at the images path.
# @examples
# cell(cell.args, path = path)
#' @import purrr dplyr stringr tidyr doParallel readr parallel
#' @rawNamespace import(foreach, except = c("when", "accumulate"))
#' @importFrom purrr map
#' @export
cell2 <- function(arguments,
                  cell.command = NULL,
                  # cell.command = "cellBUILTIN",
                  # cell.command = "~/Software/cellID-linux/cell",
                  # cell.command = "~/Projects/Rdevel/rcell2/bin/cell",
                  no_cores = NULL, 
                  debug_flag=0,
                  dry = F,
                  label_cells_in_bf = F,
                  fill_interior_pixels = F,
                  output_coords_to_tsv = F,
                  encode_cellID_in_pixels = F,
                  ignore.stdout = T, intern = F){
  
  positions <- arguments$pos %>% unique()
  n_positions <- positions %>% length()
  n_times <- arguments$t.frame %>% unique() %>% length()
  
  # Create output directories
  for(d in unique(arguments$output)) dir.create(d, showWarnings = F)
  
  # Run CellID
  if(is.null(no_cores)) no_cores <- parallel::detectCores() - 1  # Problema rarísimo: se repiten rows cada "no_cores" posiciones
  cl <- parallel::makeCluster(
    min(n_positions,
        no_cores), 
    outfile = tempfile(pattern = "dopar", tmpdir = "/tmp/", fileext = ".log")
  )
  parallel::clusterExport(cl, "arguments")
  doParallel::registerDoParallel(cl)
  
  sent_commands <- foreach::foreach(pos=positions) %dopar% {
  # sent_commands <- list()
  # for(pos in positions){
    arguments_pos <- arguments[arguments$pos == pos,]
    print(arguments_pos)
    
    bf_rcell2 <- tempfile(tmpdir = arguments_pos$output[1],
                          fileext = ".txt",
                          pattern = "bf_rcell2.")
    fl_rcell2 <- tempfile(tmpdir = arguments_pos$output[1],
                          fileext = ".txt",
                          pattern = "fl_rcell2")
    
    base::write(x = paste0(arguments_pos$path, "/", arguments_pos$bf), file = bf_rcell2)
    base::write(x = paste0(arguments_pos$path, "/", arguments_pos$image), file = fl_rcell2)

    # if(is.null(cell.command)) cell.command <- system.file("cell", package = "rcell2", mustWork = T)
    if(is.null(cell.command)) stop("\nError: cell.command must point to an existing CellID binary on your system.")
    
    command <- paste0(normalizePath(cell.command),
                      " -b ", bf_rcell2,
                      " -f ", fl_rcell2,
                      " -o ", paste0(normalizePath(arguments_pos$output[1]), "/out"),
                      " -p ", arguments_pos$parameters[1],
                      {if(label_cells_in_bf) " -l" else ""},
                      {if(output_coords_to_tsv) " -t" else ""},
                      {if(encode_cellID_in_pixels) " -m" else ""},
                      {if(fill_interior_pixels) {if(encode_cellID_in_pixels) " -i" else " -m -i"} else ""}
    )

    if(!dry & ignore.stdout) warning("Running CellID through a system call ignoring standard output messages (ignore.stdout = T). This is discouraged!")
    if(!dry) {
      command.output <- system(command = command,
                               wait = T,
                               ignore.stdout = ignore.stdout & !intern,
                               intern = intern)
      if(intern) write(command.output,
                       tempfile(tmpdir = arguments_pos$output[1],
                                fileext = ".txt",
                                pattern = "cellid_log."))
    }
    
    print("---- Done with this position.")
    print(command)
    command
  }
  
  parallel::stopCluster(cl)
  
  cat("\nDone, please examine logs above if anything seems strange :)")
  
  return(invisible(bind_rows(commands = unlist(sent_commands))))
}

#' Cluster test
cluster_test <- function(){
  cl <- parallel::makeCluster(2)
  
  doParallel::registerDoParallel(cl)
  
  # Prueba con base
  parallel::parLapply(cl, list(1,2), function(x) print(x))
  
  # Prueba con foreach
  library(foreach)
  foreach(x=list(1,2)) %dopar% print(x)
  
  parallel::stopCluster(cl)
}

#' Obtener argumentos para CellID
#'
#' @param path directory where images are stored, full path.
#' @param parameters path to the parameters file or a data.frame with "pos" (position number) and "parameter" (path) columns.
#' @param BF.pattern regex pattern to BF images. Defaults to: \code{"^BF"}
#' @param file.pattern regex pattern for all tif files, with one group for each of \code{c("ch", "pos", "t.frame")} in \code{file.pattern.groups.order}. Uses \code{"^(BF|[A-Z]FP)_Position(\\d+)_time(\\d+).tif$"} by default. To omit time, use an empty group for the t.frame in the regex, for example: \code{"^(BF|[A-Z]FP)_Position(\\d+)().tif$"}.
#' @param file.pattern.groups.order a character vector of components \code{c("ch", "pos", "t.frame")} with order corresponding to the order of groups in \code{file.pattern}.
#' @param out.dir name for output directories paths "out"
#' @param tiff.ext regex pattern for the tif file extension
#' @return a data.frame with all the information needed to run CellID
#' @import dplyr tidyr
# @examples
# cell.args <- cellArgs(path = path)
#' @export
cellArgs2 <- function(path,
                      parameters,
                      BF.pattern = "^BF",
                      file.pattern = "^(BF|[A-Z]FP)_Position(\\d+)_time(\\d+).tif$",
                      file.pattern.groups.order = c("ch", "pos", "t.frame"),
                      out.dir = "out",
                      tiff.ext = "tif$"
) {
  
  if(!identical(sort(file.pattern.groups.order), sort(c("ch", "pos", "t.frame")))) 
    stop('file.pattern.groups.order must contain c("ch", "pos", "t.frame")')
  
  path <- normalizePath(path)
  
  pic_files <- dir(path, pattern = file.pattern)
  if(length(pic_files) == 0) stop(paste("cellArgs2 error: no image files retrieved using file pattern:", file.pattern))
  
  pics_df <- data.frame(image = pic_files,
                        path = path) %>% 
    tidyr::extract(col = image, 
                   into = file.pattern.groups.order,
                   regex = file.pattern, 
                   remove = F)
  
  fluor_pics <- pics_df %>% 
    filter(str_detect(string = image,
                               pattern = BF.pattern, 
                               negate = T))
  if(nrow(fluor_pics) == 0) stop("Fluorescence images missing, Check your directories and file.pattern.")
  
  brihtfield_pics <- pics_df %>% 
    dplyr::filter(str_detect(string = image,
                             pattern = BF.pattern)) %>% 
    dplyr::rename(bf = image) %>% 
    dplyr::select(pos, t.frame, bf)
  if(nrow(fluor_pics) == 0) stop("Brightfield images missing, Check your directories and file.pattern.")
  
  
  arguments <- dplyr::left_join(
    fluor_pics,
    brihtfield_pics,
    by = c("pos", "t.frame")
  )
  # Add output column and arrange by position and t.frame
  arguments <- arguments %>% 
    mutate(output = paste0(path, "/Position", pos)) %>% 
    mutate(pos = as.integer(pos),
           t.frame = as.integer(t.frame)) %>% 
    arrange(pos, t.frame)
  
  if(length(parameters) == 1) {arguments$parameters <- parameters} else{
    arguments <- left_join(
      arguments,
      select(parameters, pos, parameters),
      by = "parameters"
    )
  }
  
  arguments <- arguments %>% mutate(parameters = normalizePath(parameters))
  
  if(all(is.na(arguments$t.frame))){
    warning("cellArgs2 warning: No t.frame data extracted, replacing all NAs with '1'. Check your directories and file.pattern if this is unexpected.")
    arguments$t.frame <- 1
  } else if(any(is.na(arguments))){
    print(arguments)
    stop("cellArgs2 error: at least one of the values in the arguments dataframe is missing, check your directories and file.pattern")
  }
  
  return(arguments)
}

#' Cargar el output de cell-id
#'
#' @param path Path to images directory
#' @param pdata Path to metadata CSV file
#' @param position.pattern Regex describing what the position string looks like (default ".*Position(\\d+).*") including a capturing group for the position ID number (coerced to integer).
#' @param ucid.zero.pad amount of decimal digits for the cellID (defaults 4, corresponding to a maximum of 9.999 cellIDs and 9999 positions)
#' @param ... Arguments passed on to \code{cargar.out_all}. Patterns for "out" files, fluorescence channel, and other options may be changed here.
#' @return A list of dataframes: data (CellID data), images (images metadata and paths), image_maping (extra mapping metadata from CellID: BF to FL correspondence, channel flag, bf_as_fl flag, and one-letter channel encoding).
# @examples
# cell.data <- cell.load(path = path, pdata = pdata)
#' @import dplyr stringr tidyr readr
#' @importFrom purrr map
#' @export
cell.load.alt <- function(path,
                          pdata = NULL,
                          position.pattern = ".*Position(\\d+).*",
                          ucid.zero.pad = 4,
                          ...){
  
  path <- normalizePath(path, mustWork = T)

  # Cargar datos out_all y juntar en un solo dataframe usando metadata de "out_bf_fl_mapping"
  cat("\n\nLoading CellID output files...\n")
  d.list <- cargar.out_all(.carpeta = path,
                           position.pattern = position.pattern,
                           ...)  # https://stackoverflow.com/questions/40794780/r-functions-passing-arguments-with-ellipsis/40794874
  
  cat("\rCreating ucid column...                            ")
  d.list$d <- d.list$d %>%
    # By padding to the same number of digits, cells from positions as 90 and 9 could have the same UCID.
    # The next lines fix that bug, which would cells to not be filtered correctly, and be plotted anyways, or other problems.
    # Also, the padding should not be very large, or it will "overflow" R's integer class.
    mutate(cellid.pad = ucid.zero.pad - nchar(as.character(cellID))) %>%
    mutate(
      ucid = as.integer(
        paste0(pos,
               sapply(cellid.pad, FUN = function(pad) paste0(rep(x="0", 
                                                                 times=max(0, pad)), 
                                                             collapse = "")),
               cellID))) %>%
    select(-cellid.pad)
  
  # Check uniqueness of ucid-t.frame combinations
  if(nrow(unique(d.list$d[,c("ucid", "t.frame")])) < nrow(d.list$d)) stop("\nERROR: There are repeated cellID's in the out_all file!")

  # ellipse.perim = perimeter of theoretical ellipse, calculated using each
  # cell's axis values.
  # el.p = ratio of ellipse perim over the perimeter measured by cellID.
  # If this number is small ( < ~0.7) it's probably not a cell.
  cat("\rCreating el.p column...                            ")
  d.list$d <- dplyr::mutate(d.list$d,
                            ellipse.perim = pi *
                              (3 * (maj.axis / 2 + min.axis / 2) -
                                 sqrt((3 * maj.axis / 2 + min.axis / 2) *
                                        (maj.axis / 2 + 3 * min.axis / 2))),

                            el.p = ellipse.perim / perim)

  # Mergear con pdata
  # if(exists("pdata", inherits = F)){
  cat("\rJoining pdata if specified...")
  if(!is.null(pdata)){
    pdata <- read_csv(pdata) %>% mutate(pos = pos %>% as.numeric)
    d.list$d <- d.list$d %>% left_join(pdata, by = "pos")
    cat(" and it was :)                            ")
  } else cat(" but it was not :(                            ")


  # Create paths dataframe and add three-letter code for channel
  cat("\rCreating image paths dataframe...                            ")
  paths <- d.list$d.map %>%
    mutate(channel = paste0(toupper(channel), "FP"))

  paths <- bind_rows( # Bind it with itself, to get entries for BF as well.
    paths %>%  # Get FL image paths
      select(pos, t.frame, channel, fluor) %>%
      rename(file = fluor),
    paths %>%  # Get BF image paths
      select(pos, t.frame, channel, bright) %>%
      rename(file = bright) %>%
      mutate(channel = "BF") %>% 
      unique()  # There may be BF path duplicates in the "bright" column, so keep the unique set
  ) %>%
    mutate(path = dirname(file),  # Add the directory path
           is.out = FALSE)        # Add the is.out column

  paths <- bind_rows(paths,
                     paths %>%  # bind it with part of itself, out files are named exactly the same, but with an extra ".out.tif"
                       mutate(file = paste0(file, ".out.tif"),
                              channel = paste0(channel, ".out"),
                              is.out = TRUE)) %>%
    mutate(image = basename(file))

  d.list$d.paths <- paths
  
  # Make output list
  cell.data <- list(data = d.list$d,
                    images = d.list$d.paths,
                    mapping = d.list$d.map,
                    channels = unique(d.list$flag.channel.mapping))
  
  cat("\rDone loading CellID data!                            \n")
  return(cell.data)
}


#' Una función que lea un .csv y les agregue una columna con un id del archivo (pos)
read_tsv.con.pos <- function(.nombre.archivo, .carpeta, position.pattern, col_types = "c"){
  cat(paste0("\rReading: ", .nombre.archivo), "\033[K")
  .archivo <- normalizePath(paste0(.carpeta, "/", .nombre.archivo))
  .pos <- str_replace(.nombre.archivo, position.pattern, "\\1") %>% as.numeric()
  
  d <-  read_tsv(.archivo, col_types = col_types, trim_ws = T) %>%
    mutate(pos = as.integer(.pos))  # La columna de ID es "pos"
  
  
  if("con.vol_1" %in% names(d)) {
    cat(paste0("\nRemoving 'con.vol_1' column from position: ", .pos, ". Use CellID version > 1.4.6 to stop seeing this message.\n"))
    d <- dplyr::select(d, -con.vol_1)
  }
  
  return(d)
}

#' Una función para leer y combinar todos los archivos "out".
#'
#' @param .nombre.archivos paths de los "out_all"
#' @param .nombre.archivos.map paths de los "bf_fl_mapping"
#' @param .carpeta path a la carpeta donde están os output de CellID
#' @param position.pattern Regex describing what the position string looks like (default 'Position\\d+')
#' @param fluorescence.pattern Regex describing what the fluorescence string looks like (default ".*([GCYRT])FP_Position.*")
#' @import dplyr tidyr readr
#' @importFrom purrr map
#' @return A list of two dataframes: `d` contains the actual output, and `out.map` contains image paths and metadata.
cargar.out_all <- function(#.nombre.archivos, .nombre.archivos.map,
                           .carpeta = "data/images2/",
                           position.pattern = ".*Position(\\d+).*",
                           out_file_pattern = "^out_all$",
                           out_mapping_pattern = "^out_bf_fl_mapping$",
                           fluorescence.pattern = ".*([A-Z])FP_Position.*"){

  # Migrated from cell.load()
  .nombre.archivos <- list.files(path = .carpeta, pattern = out_file_pattern, recursive = T, include.dirs = T)
  .nombre.archivos.map <- list.files(path = .carpeta, pattern = out_mapping_pattern, recursive = T, include.dirs = T)
  # A bit of error handling
  if(length(.nombre.archivos) == 0) stop("Error in cargar.out_all: no CellID output files found, check your path, options and files.")
  if(length(.nombre.archivos.map) == 0) stop("Error in cargar.out_all: no CellID mapping files found, check your path, options and files.")
  if(length(.nombre.archivos) != length(.nombre.archivos.map)) stop("Error in cargar.out_all: different amount of mapping and cellid output files.")
  
  # Cargo y junto los "out_all"
  cat("\rLoading datasets...\033[K")
  d.out <- purrr::map(.x = .nombre.archivos,
                      .f = read_tsv.con.pos, 
                      .carpeta = .carpeta,
                      "iiiiiddddddddddddddddddddddddddddddddddddddddddddddddddd", # types: 5 int columns, 51 double columns
                      position.pattern = position.pattern) %>%
    bind_rows()
  cat("\n Done loading out files!\n")

  # # Cargo y junto los "out_bf_fl_mapping"
  cat("\rLoading mapping...              ")
  d.map <- purrr::map(.x = .nombre.archivos.map,
                      .f = read_tsv.con.pos,  # Una función para leer los archivos "out" y agregarles "pos" segun la carpeta que los contiene
                      .carpeta = .carpeta,
                      col_types = "ciicl",  # types: 2 char columns, 2 int columns, 1 logical column
                      position.pattern = position.pattern) %>%
    bind_rows() %>%
    mutate(channel = str_replace(string = fluor,
                                 pattern = fluorescence.pattern,
                                 replacement = "\\1")) %>%
    mutate(channel = tolower(channel))
  cat("\n Done loading mapping files!\n")
  
  # keep flag-channel mapping
  flag.channel.mapping <- d.map %>% select(flag, channel) %>% unique()

  # Return join (discard flag variable)
  cat("\rJoining data and mapping...\033[K")
  d.out.map <- left_join(d.out,
                         unique(select(d.map, flag, t.frame , pos, channel)),
                         by = c("flag", "t.frame", "pos")) %>% 
    select(-flag)
  
  if(nrow(d.out.map) > nrow(d.out)) stop("Error acargar.out_all: while joining output and mapping, at least one output row matche multiple mappings.")
  
  # Add f.tot columns to data
  d.out.map <- mutate(d.out.map,
                      f = f.tot - (a.tot * f.bg),
                      cf = f / a.tot)

  # Right now the out_all is in a "long" format for the (one-letter) channel variable
  # Spread it to out expectations:
  cat("\rSpreading data from channels...\033[K")
  cdata <- d.out.map %>%
    tidyr::pivot_wider(
      names_from = channel, names_sep = ".",
      id_cols = c(cellID,
                  # flag,  # dropped above
                  t.frame,
                  time,
                  pos,
                  xpos,
                  ypos,
                  a.tot,
                  num.pix,
                  fft.stat,
                  perim,
                  maj.axis,
                  min.axis,
                  rot.vol,
                  con.vol,
                  a.tot.p1,
                  a.tot.m1,
                  a.tot.m2,
                  a.tot.m3,
                  a.local.bg,
                  a.local,
                  a.local2.bg,
                  a.local2,
                  a.surf,
                  # con.vol_1,  # duplicated, removed by read_tsv.con.pos and in recent CellID versions
                  sphere.vol),

      values_from = c(f.tot,
                      f,
                      cf,
                      f.nucl,
                      a.nucl,
                      a.vacuole,
                      f.vacuole,
                      f.bg,
                      f.tot.p1,
                      f.tot.m1,
                      f.tot.m2,
                      f.tot.m3,
                      xpos.nucl,
                      ypos.nucl,
                      f.nucl1,
                      f.nucl.tag1,
                      a.nucl1,
                      f.nucl2,
                      f.nucl.tag2,
                      a.nucl2,
                      f.nucl3,
                      f.nucl.tag3,
                      a.nucl3,
                      f.nucl4,
                      f.nucl.tag4,
                      a.nucl4,
                      f.nucl5,
                      f.nucl.tag5,
                      a.nucl5,
                      f.nucl6,
                      f.nucl.tag6,
                      a.nucl6,
                      f.local.bg,
                      f.local2.bg))
  
  cdata <- mutate(cdata, 
                  cellID = as.integer(cellID),
                  t.frame = as.integer(t.frame))

  return(list(
    "d" = cdata,
    "d.map" = d.map,
    "flag.channel.mapping" = flag.channel.mapping
    ))
}

#' cellArgs2 Summaries
#' 
#' A function to print some summaries, to check cellArgs2 output.
#' 
cellArgs2.summary <- function(arguments){
  arguments %>% group_by(ch) %>% summarise(n_count = n(), .groups = "drop") %>% print()
  arguments %>% select(bf) %>% summarise(unique_BF = "", n_count = length(unique(bf)), .groups = "drop") %>% print()
  arguments %>% group_by(t.frame) %>% summarise(n_count = n(), .groups = "drop") %>% print()
  arguments %>% group_by(pos) %>% summarise(n_count = n(), .groups = "drop") %>% print()
}

#' Pipe
#'
#' Put description here
#'
#' @importFrom purrr %>%
#' @name %>%
#' @rdname pipe
#' @export
#' @param lhs,rhs specify what lhs and rhs are
#' @examples
#' # some examples if you want to highlight the usage in the package
NULL
