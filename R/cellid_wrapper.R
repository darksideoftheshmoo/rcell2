# Previously named cell3.R
#
# library(tidyverse)
# library(foreach)
# library(doParallel)

##### DEFINITIONS ####

#' Correr cellid desde C
#'
#' @param args el comando de cellid entero, tal como se ejecutaria en bash "cell -p ..."
#' @param label_cells Set to 0 to disable labeling cells with their number.
#' @useDynLib rcell2 CellID
#' @export
#' @return Nada, el output está en los directorios.
cellid <- function(args, label_cells=1) {

  # args <- "~/Software/cellID-linux/cell -p /home/nicomic/Projects/Colman/HD/scripts/cellMagick/data/images/parameters.txt -b /tmp/Rtmp7fjlFo/file2b401093d715 -f /tmp/Rtmp7fjlFo/file2b402742f6ef -o /home/nicomic/Projects/Colman/HD/uscope/20200130_Nico_screen_act1_yfp/1/Position001/out"
  argv <- strsplit(args, " ")[[1]]
  argc <- length(argv)

  # print(argv)
  # print(argc)

  exit_code = .C(CellID, 
                 as.integer(argc),         # Argument count
                 as.character(argv),       # Argument character vector
                 integer(1),               # Return variable
                 as.integer(label_cells)   # Option to disable
                 )[[3]]
  
  exit_code
}


#' Buscar BFs
#'
#' @param path directory where images are stored, full path.
#' @param pattern BF pattern "BF_Position\\d+\\.tif$"
#' @return A vector of directories.
# @examples
# cell.images.BF(path, pattern="BF_Position\\d+\\.tif$", out.dir = "out")
cell.images.BF <- function(path, BF.pattern="BF_Position\\d+\\.tif$"){
  f <- sub(x = dir(path = path, pattern = BF.pattern, full.names = T),
           pattern = "//", replacement = "/")
  return(normalizePath(f))
}


#' Buscar ?FPs
#'
#' @param path directory where images are stored, full path.
#' @param pattern BF pattern "FP_Position\\d*\\.tif$"
#' @return A vector of directories.
# @examples
# cell.images.FP(path, pattern="FP_Position\\d*\\.tif$")
cell.images.FP <- function(path, FP.pattern="FP_Position\\d*\\.tif$"){
  # f <- sub(x = paste0(dir(path = path, pattern = pattern, full.names = T)),
  f <- sub(x = dir(path = path, pattern = FP.pattern, full.names = T),
           pattern = "//", replacement = "/")
  return(normalizePath(f))
}


#' Filtrar cdata usando gráficos y dibujando regiones
#'
#' @param path directory where images are stored, full path.
#' @param BF.pattern regex pattern to ?FP images
#' @param O.pattern regex pattern to output directories
#' @param out.dir name for output directories paths "out"
#' @return A vector of directories.
# @examples
# cell.images.out(path, pattern="BF_Position\\d*\\.tif$", out.dir = "out")
cell.images.out <- function(path,
                            BF.pattern="BF_Position\\d+\\.tif$",
                            O.pattern=".*(Position\\d+).*\\.tif",
                            out.dir = "out"){
  .path <- normalizePath(path)

  .pics <- dir(path = .path, pattern = BF.pattern, full.names = T)

  f <- sub(x=.pics,
           pattern = O.pattern,
           replacement = "\\1")
  d <- sub(
    x = paste(.path, f, out.dir, sep = "/"),
    pattern = "//", replacement = "/"
  )

  names(d) <- f

  return(d)
}


#' Filtrar cdata usando gráficos y dibujando regiones
#'
#' @param path directory where images are stored, full path.
#' @param parameters absolute path to the parameters file.
#' @param BF.pattern regex pattern to BF images
#' @param FP.pattern regex pattern to ?FP images
#' @param O.pattern regex pattern to output directories
#' @param out.dir name for output directories paths "out"
# @param channels Default c("b", "f"), don't change, used to create temporal file lists for BF and ?FP.
#' @return Nothing.
# @examples
# cell.args <- cellArgs(path = path)
#' @export
cellArgs <- function(path,
                     parameters,
                     BF.pattern="BF_Position\\d+\\.tif$",
                     FP.pattern="FP_Position\\d*\\.tif$",
                     O.pattern=".*(Position\\d+)\\.tif",
                     out.dir = "out"
                     ) {

  parameters <- normalizePath(parameters)

  b <- cell.images.BF(path, BF.pattern)
  f <- cell.images.FP(path, FP.pattern)
  o <- cell.images.out(path, BF.pattern, O.pattern, out.dir)

  list(p=parameters, b=b, f=f, o=o)
}


#' Print cell.args in a more understandable way.
#' 
#' Note that the output basename ("o" parameter) actually has as many different directories as positions.
#'
#' @param cell.args cell.args list.
#' @param which.args which arguments should be printed
#' @return Nothing.
# @examples
# cellArgs.print(cell.args)
#' @export
cellArgs.print <- function(cell.args, which.args = c("p", "b", "f", "o")) for(i in which.args) {
  print(i, quote = F)
  print("  Dirs:", quote = F)
  print(cell.args[[i]] %>% dirname() %>% unique(), quote = F)
  print("  Files:", quote = F)
  for(arg in cell.args[[i]]) print(basename(arg), quote = F)
  print("----", quote = F)
}


#' Cargar el output de cell-id
#'
#' @param path Path to images directory
#' @param position.pattern Regex describing what the position string looks like (default 'Position\\d+')
#' @param ... pass "pdata" with the path to the pdata file to merge it internally
#' @return a data frame
# @examples
# cell.data <- cell.load(path = path, pdata = pdata)
#' @import dplyr stringr tidyr readr
#' @importFrom purrr map
#' @export
cell.load.alt <- function(path = "data/images2/",
                          pdata = NULL,
                          position.pattern = "Position\\d+",
                          ...){

  if(F){  # TEST
    path <- "/home/nicomic/Projects/Colman/HD/uscope/20200130_Nico_screen_act1_yfp/1/"
    position.pattern = "Position\\d+"
  }

  # Cargar datos out_all y juntar en un solo dataframe usando metadata de "out_bf_fl_mapping"
  d.list <- cargar.out_all(.carpeta = path,
                           position.pattern = position.pattern)

  d.list$d <- d.list$d %>%
    # Previously cells from positions like 10 and 1 could have the same UCID.
    #mutate(ucid = 10 - pos %>% as.character %>% nchar - cellID %>% as.character %>% nchar) %>% #glimpse()
    # Next line fixes that bug, which would cells to not be filtered correctly, and be plotted anyways.
    mutate(ucid = 10 - 1 - cellID %>% as.character %>% nchar) %>%
    mutate(
      ucid = paste0(
        pos,
        sapply(ucid, FUN = function(ucid) rep.int(x = 0, times = ucid) %>% paste0(collapse = "")),
        cellID
      )) %>%
    select(pos, cellID, ucid, everything()) #%>% glimpse()
  # d$ucid %>% nchar() %>% unique()

  # ellipse.perim = perimeter of theoretical ellipse, calculated using each
  # cell's axis values.
  # el.p = ratio of ellipse perim over the perimeter measured by cellID.
  # If this number is small ( < ~0.7) it's probably not a cell.
  # f.x = total fluorescence - background for channel x
  # cf.x = concentration of f.x (divided by cell area)
  d.list$d <- dplyr::mutate(d.list$d,
                            ellipse.perim = pi *
                              (3 * (maj.axis / 2 + min.axis / 2) -
                                 sqrt((3 * maj.axis / 2 + min.axis / 2) *
                                        (maj.axis / 2 + 3 * min.axis / 2))),

                            el.p = ellipse.perim / perim)
  # d.list$d <- d.list$d %>%
  if("f.tot.y" %in% names(d.list$d)) d.list$d <- d.list$d %>%
    dplyr::mutate(f.y = f.tot.y - (a.tot * f.bg.y),
                  cf.y = f.y / a.tot)
  if("f.tot.r" %in% names(d.list$d)) d.list$d <- d.list$d %>%
    dplyr::mutate(f.r = f.tot.r - (a.tot * f.bg.r),
                  cf.r = f.r / a.tot)
  if("f.tot.c" %in% names(d.list$d)) d.list$d <- d.list$d %>%
    dplyr::mutate(f.c = f.tot.c - (a.tot * f.bg.c),
                  cf.c = f.c / a.tot)

  # Mergear con pdata
  # if(exists("pdata", inherits = F)){
  if(!is.null(pdata) & exists("pdata", inherits = F)){
    pdata <- read_csv(pdata) %>% mutate(pos = pos %>% as.numeric)
    d.list$d <- d.list$d %>% left_join(pdata, by = "pos")
  } else print("Either pdata was not supplied or the file was not found in the supplied path.")


  # Creathe paths dataframe
  paths <- d.list$d.map %>%
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

  d.list$d.paths <- paths

  return(d.list)
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
                           position.pattern = "Position\\d+",
                           fluorescence.pattern = ".*([GCYRT])FP_Position.*"){

  if(F){  #TEST
    .carpeta = path
    position.pattern = "Position\\d+"
    fluorescence.pattern = ".*([GCYRT])FP_Position.*"
  }

  # Migrated from cell.load()
  .nombre.archivos <- dir(path = .carpeta, pattern = "^out_all$", recursive = T)
  .nombre.archivos.map <- dir(path = .carpeta, pattern = "^out_bf_fl_mapping$", recursive = T)

  # Una función que lea un .csv y les agregue una columna con un id del archivo (pos)
  read_tsv.con.id <- function(.nombre.archivo,
                              .carpeta = "data/images2/",
                              position.pattern2 = ".*Position(\\d+).*"){

    .archivo <- paste0(.carpeta, "/", .nombre.archivo) %>% sub("//", "/", x = .)
    .pos <- .nombre.archivo %>% str_replace(position.pattern2, "\\1") %>% as.numeric()

    d <- .archivo %>%
      read_tsv(col_types = cols()) %>%
      mutate(pos = .pos)  # La columna de ID es "pos"

    print("Removing 'con.vol_1' duplicated column if it exists.")
    if("con.vol_1" %in% names(d)) d <- d %>% select(-con.vol_1)

    d
  }

  # Cargo y junto los "out_all"
  d.out <- purrr::map(.x = .nombre.archivos,
                      .f = read_tsv.con.id,
                        .carpeta = .carpeta) %>%
    bind_rows() %>%
    mutate(pos = as.integer(pos))
  if(length(d.out$cellID %>% unique()) < length(d.out)){
    stop("There are repeated cellID's in the out_all file!")
  }

  # # Cargo y junto los "out_bf_fl_mapping"
  d.map <- purrr::map(.x = .nombre.archivos.map,
                      .f = read_tsv.con.id,
                        .carpeta = .carpeta) %>%
    bind_rows() %>%
    mutate(pos = as.integer(pos)) %>%
    # select(-bright, -bf.as.fl) %>%
    mutate(channel = str_replace(string = fluor,
                               pattern = fluorescence.pattern,
                               replacement = "\\1")) %>%
    mutate(channel = tolower(channel))

  # Return join
  d.out.map <- d.map %>%
    # select(-bright, -bf.as.fl) %>%
    select(channel, flag, t.frame , pos) %>%
    left_join(d.out, by = c("flag", "t.frame", "pos")) %>%
    select(-flag)
  # d %>% select(channel, flag) %>% unique()  # Check uniqueness of flag-fluor combinations

  d <- d.out.map %>%
    pivot_wider(
      names_from = channel, names_sep = ".",
      id_cols = c(cellID,
                  # flag,
                  t.frame,
                  time,
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
                  # con.vol_1,  # duplicated, removed by read_tsv.con.id
                  sphere.vol,
                  pos),

      values_from = c(f.tot,
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

  return(list(
    "d" = d,
    "d.map" = d.map
    ))
}


#' Filtrar cdata usando gráficos y dibujando regiones
#'
#' @param cell.args An argument list, as built by cellArgs().
#' @param position.pattern a regular expression that recognizes the position in the image name, such as: "Position\\d+" for Position01, Position100, etc.
#' @param position.time.pattern a regular expression that recognizes the time in the image name, such as: "time\\d+" for time01, time007, etc.
#' @param position.time.pattern.sep a regular expression that matches the pos-time separator characters (default is an optional underscore), such as: "_" for Position01_time02
#' @param fluorescence.pattern a regular expression that recognizes the fluorescence tag in the image name, such as "^(.FP)" for TFP, YFP, etc.
#' @param cell.command the CellID command, either "cellBUILTIN" for the builtin binary, a path to the binary executable (get if from https://github.com/naikymen/cellID-linux).
#' @param channels Default c("b", "f"), don't change, used to create temporal file lists for BF and ?FP and pass cell id arguments.
#' @param no_cores Position-wise parallelization,internally capped to number of positions in cell.args.
#' @param dry Do everything without actually running CellID.
#' @return Nothing.
# @examples
# cell(cell.args, path = path)
#' @import purrr dplyr stringr tidyr doParallel readr parallel
#' @rawNamespace import(foreach, except = c("when", "accumulate"))
#' @importFrom purrr map
#' @export
cell <- function(cell.args,
                 position.pattern =  "Position\\d+", 
                 position.time.pattern = NULL, #position.time.pattern = "time\\d+",
                 position.time.pattern.sep = "_?",
                 fluorescence.pattern = "^(.FP)", 
                 cell.command = "cellBUILTIN",
                 # cell.command = "~/Software/cellID-linux/cell",
                 channels = c("b", "f"),
                 old_dirs_path = NULL, old_dirs_pattern = "^Position\\d\\d\\d$",
                 no_cores = NULL, dry = F){
  
  # Optional: remove old dirs
  if(!is.null(old_dirs_path)){
    dir(path = path, 
        pattern = old_dirs_pattern, 
        full.names = T, 
        include.dirs = T) %>% 
      unlink(recursive = T)
  }

  # Setup
  {
    # BF.images <- cell.args$b %>% {.[order(str_extract(., position.pattern))]} # Order FPs by position
    BF.positions <- str_extract(cell.args$b, position.pattern)  # Grab their position identifier

    n_times <- 1
    if(!is.null(position.time.pattern)) n_times <- str_extract(basename(cell.args$b),
                                                               position.time.pattern) %>%
      unique() %>% length() # Count number of times
    u_positions <- unique(BF.positions)
    n_positions <- length(u_positions) # Count number of positions

    f_channels <- basename(cell.args$f) %>% str_extract(fluorescence.pattern) %>% unique()
    n_channels <- length(f_channels)

    # cell.args$f <- cell.args$f %>% {.[order(str_extract(., position.pattern))]}  # Order FPs by position
    # FP.images <- cell.args$f
    # FP.positions <- str_extract(FP.images, position.pattern)  # Grab their position identifier
    # names(cell.args$f) <- FP.positions %>% str_replace(fluorescence.pattern, "\\1")  # name them by pos/time
    if(is.null(position.time.pattern)){
      names(cell.args$f) <- str_extract(cell.args$f, position.pattern)
      names(cell.args$b) <- str_extract(cell.args$b, position.pattern)  # Name the BF vector with their positions
    } else {
      names(cell.args$f) <- str_extract(cell.args$f,
                                        paste0(position.pattern, position.time.pattern.sep, position.time.pattern))
      names(cell.args$b) <- str_extract(cell.args$b,
                                        paste0(position.pattern, position.time.pattern.sep, position.time.pattern))
    }

    # CellID parameters
    parameters <- cell.args$p[1]
    # Lo siguiente no es necesario ahora. Todas las fotos van a compartir el mismo parameters.txt
    # cell.args$p <- rep(cell.args$p[1], length(FP.images)) # Repeat parameters as many times as needed for FP.images

    # Create output directories
    for(k in 1:length(cell.args$b)) dir.create(cell.args$o[k], recursive = T, showWarnings = F)

    # Name output prefix with position
    names(cell.args$o) <- cell.args$o %>% str_extract(position.pattern)  # Name the output vector with their positions
    # cell.args$o <- cell.args$o[FP.positions]  # Repeat output as many times as needed for FP.images

    # Repeat BFs as many times as needed for FP.images
    # cell.args$b <- cell.args$b[names(cell.args$f)]
    # cell.args$b <- rep(cell.args$b, each = n_channels)
    cell.args$b <- rep(cell.args$b, times = n_channels)

    if(!all(names(cell.args$b) == names(cell.args$f))) stop("BF list and ?FP list have a problem...")
  }

  # Run CellID
  if(is.null(no_cores)) no_cores <- min(round(detectCores()/2), 1)  # Problema rarísimo: se repiten rows cada "no_cores" posiciones
  # cl <- makeCluster(no_cores, outfile="tests/foreach.log")
  cl <- makeCluster(
                    min(n_positions,
                        no_cores)
                   )
  # cl <- makeForkCluster()
  # cl <- makeCluster(no_cores, type = "FORK")
  registerDoParallel(cl)

  # registerDoParallel(no_cores)
  commands <- foreach(pos=1:n_positions) %dopar% {
  # for(i in 1:n_positions){

    cell.args.tmp <- c("p" = unname(parameters),  # there is a "feature" here worthy of the R Inferno
                       "o" = unname(cell.args$o[pos*n_times])) # unnaming is necesary or the elements name is joined to the assigned name

    for(channel in channels) {
      tmp <- tempfile(tmpdir = cell.args.tmp["o"],
                      fileext = ".txt",
                      pattern = paste("pos", pos,
                                      "ch", channel,
                                      "param_", sep = "_"))
      # paths <- cell.args[[j]][names(cell.args[[j]]) %in% BF.positions[i]]
      # paths <- cell.args[[channel]][str_detect(names(cell.args[[channel]]), u_positions[pos])]
      paths <- cell.args[[channel]][ grepl(pattern = u_positions[pos],
                                           x = names(cell.args[[channel]])) ]
      write(x = paths, file = tmp)  # readLines(tmp)
      cell.args.tmp[channel] <- normalizePath(tmp)
    }

    command <- paste(cell.command,
                     paste0("-",names(cell.args.tmp), " ", cell.args.tmp,
                            collapse = " ")
    )

    if(cell.command == "cellBUILTIN") {
      if(!dry) cellid(command)
    } else {
      if(!dry) system(command = command, wait = T)
    }

    return(command)
  }

  stopCluster(cl)
  print("Done, please examine logs above if anything seems strange :)")
  return(commands)
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
