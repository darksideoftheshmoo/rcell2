# Previously named cell3.R
#
# library(tidyverse)
# library(foreach)
# library(doParallel)

##### DEFINITIONS ####

#' Correr cellid desde C
#'
#' @param args el comando de cellid entero, tal como se ejecutaria en bash "cell -p ..."
#' @useDynLib rcell2 main_
#' @export
#' @return Nada, el output está en los directorios.
cellid <- function(args) {

  # args <- "~/Software/cellID-linux/cell -p /home/nicomic/Projects/Colman/HD/scripts/cellMagick/data/images/parameters.txt -b /tmp/Rtmp7fjlFo/file2b401093d715 -f /tmp/Rtmp7fjlFo/file2b402742f6ef -o /home/nicomic/Projects/Colman/HD/uscope/20200130_Nico_screen_act1_yfp/1/Position001/out"
  argv <- strsplit(args, " ")[[1]]
  argc <- length(argv)

  # print(argv)
  # print(argc)

  .C(main_, as.integer(argc), as.character(argv), integer(1))[[3]]
}


#' Buscar BFs
#'
#' @param path directory where images are stored, full path.
#' @param pattern BF pattern "BF_Position\\d+\\.tif$"
#' @return A vector of directories.
# @examples
# cell.images.BF(path, pattern="BF_Position\\d+\\.tif$", out.dir = "out")
cell.images.BF <- function(path, pattern="BF_Position\\d+\\.tif$"){
  f <- sub(x = dir(path = path, pattern = pattern, full.names = T),
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
cell.images.FP <- function(path, pattern="FP_Position\\d*\\.tif$"){
  # f <- sub(x = paste0(dir(path = path, pattern = pattern, full.names = T)),
  f <- sub(x = dir(path = path, pattern = pattern, full.names = T),
           pattern = "//", replacement = "/")
  return(normalizePath(f))
}


#' Filtrar cdata usando gráficos y dibujando regiones
#'
#' @param path directory where images are stored, full path.
#' @param pattern BF pattern "BF_Position\\d*\\.tif$"
#' @param out.dir name for output directories paths "out"
#' @return A vector of directories.
# @examples
# cell.images.out(path, pattern="BF_Position\\d*\\.tif$", out.dir = "out")
cell.images.out <- function(path, pattern="BF_Position\\d*\\.tif$", out.dir = "out"){
  .path <- normalizePath(path)
  
  f <- sub(x=dir(path = .path, pattern = pattern, full.names = T),
           pattern = ".*(Position\\d+)\\.tif",
           replacement = "\\1")
  d <- sub(
    x = paste(.path, f, out.dir, sep = "/"),
    pattern = "//", replacement = "/"
  )
  
  return(d)
}


#' Filtrar cdata usando gráficos y dibujando regiones
#'
#' @param path directory where images are stored, full path.
#' @param parameters absolute path to the parameters file.
#' @param b list of arguments: paths to BF images
#' @param f list of arguments: paths to ?FP images
#' @param o list of arguments: paths to output directories
#' @param channels Default c("b", "f"), don't change, used to create temporal file lists for BF and ?FP.
#' @return Nothing.
# @examples
# cell.args <- cellArgs(path = path)
#' @export
cellArgs <- function(path,
                     parameters) {

  parameters <- normalizePath(parameters)

  b <- cell.images.BF(path)
  f <- cell.images.FP(path)
  o <- cell.images.out(path)

  list(p=parameters, b=b, f=f, o=o)
}


#' Print cell.args
#'
#' @param cell.args cell.args list.
#' @param which.args which arguments should be printed
#' @return Nothing.
# @examples
# cellArgs.print(cell.args)
#' @export
cellArgs.print <- function(cell.args, which.args = c("p", "b", "f", "o")) for(i in which.args) {
  print(i, quote = F)
  print(cell.args[[i]])
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
cell.load <- function(path = "data/images2/",
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
#' @param path directory where images are stored, full path.
#' @param position.pattern a regular expression that recognizes the position in the image name
#' @param cell.command the CellID command, either "cellBUILTIN" for the builtin binary, a path to it.
#' @param channels Default c("b", "f"), don't change, used to create temporal file lists for BF and ?FP.
#' @return Nothing.
# @examples
# cell(cell.args, path = path)
#' @import purrr dplyr stringr tidyr foreach doParallel readr parallel
#' @importFrom purrr map
#' @export
cell <- function(cell.args,
                 path = "data/images/",
                 position.pattern =  "Position\\d+",
                 # cell.command = "cellBUILTIN",
                 cell.command = "~/Software/cellID-linux/cell",
                 channels = c("b", "f"),
                 no_cores = NULL){

  # Remove PositionNNN directories
  dir(path = path, pattern = "^Position\\d\\d\\d$", full.names = T) %>% unlink(recursive = T)


  if(F){ # TESTS
    i=1
    path <- paste0("/home/nicomic/Projects/Colman/HD/uscope/20200130_Nico_screen_act1_yfp/", i ,"/")
    path.pdata <- paste0("~/Projects/Colman/HD/uscope/20200130_Nico_screen_act1_yfp/", i, "/pdata.csv")

    position.pattern =  "Position\\d+"
    cell.command = "~/Software/cellID-linux/cell"
    channels = c("b", "f")
    i = 1

    # Delete outputs
    dir(path = path, pattern = "out\\.tif$", full.names = T) %>% unlink()
  }

  # Setup
  {
    parameters <- cell.args$p[1]

    BF.images <- cell.args$b %>% {.[order(str_extract(., position.pattern))]} # Order FPs by position
    BF.positions <- BF.images %>% str_extract(position.pattern)  # Grab their position identifier
    n_positions <- length(BF.images)

    cell.args$f <- cell.args$f %>% {.[order(str_extract(., position.pattern))]}  # Order FPs by position
    FP.images <- cell.args$f
    FP.positions <- FP.images %>% str_extract(position.pattern)  # Grab their position identifier
    names(cell.args$f) <- FP.positions

    cell.args$p <- rep(cell.args$p[1], length(FP.images)) # Repeat parameters as many times as needed for FP.images

    # Create output directories
    for(k in 1:length(cell.args$b)) dir.create(cell.args$o[k], recursive = T, showWarnings = F)

    names(cell.args$o) <- cell.args$o %>% str_extract(position.pattern)  # Name the output vector with their positions
    # cell.args$o <- cell.args$o[FP.positions]  # Repeat output as many times as needed for FP.images

    names(cell.args$b) <- cell.args$b %>% str_extract(position.pattern)  # Name the BF vector with their positions
    cell.args$b <- cell.args$b[FP.positions]  # Repeat BFs as many times as needed for FP.images
  }

  # Run CellID
  if(is.null(no_cores)) no_cores <- min(round(detectCores()/2), 1)  # Problema rarísimo: se repiten rows cada "no_cores" posiciones
  # cl <- makeCluster(no_cores, outfile="tests/foreach.log")
  cl <- makeCluster(no_cores)
  # cl <- makeForkCluster()
  # cl <- makeCluster(no_cores, type = "FORK")
  registerDoParallel(cl)

  # registerDoParallel(no_cores)
  foreach(i=1:n_positions) %dopar% {
  # for(i in 1:n_positions){

    cell.args.tmp <- c(p = parameters)

    for(j in channels) {
      tmp <- tempfile()
      paths <- cell.args[[j]][names(cell.args[[j]]) %in% BF.positions[i]]
      write(x = paths, file = tmp)  # readLines(tmp)
      cell.args.tmp[j] <- tmp
    }

    cell.args.tmp["o"] <- cell.args$o[i]

    command <- paste(cell.command,
                     paste0("-",names(cell.args.tmp), " ", cell.args.tmp,
                            collapse = " ")
    )

    print(command)

    if(cell.command == "cellBUILTIN") {
      cellid(command)
    } else {
      system(command = command, wait = T)
    }
  }

  stopCluster(cl)

  return("Done, please examine logs above if anything seems strange :)")
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
