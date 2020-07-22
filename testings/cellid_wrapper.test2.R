##### USAGE 3 ####
# devtools::document()
devtools::load_all(recompile = T)
# devtools::check()
# devtools::build()

{
  library(dplyr)
  library(stringr)

  path <- "./data/image_samples/"
  path.pdata <- paste0(path, "pdata.csv")
  parameters <-  paste0(path, "parameters.txt")

  # Generate CellID input paths
  cell.args <- cellArgs(path = path, parameters = parameters)
  cellArgs.print(cell.args)  # Para revisar inputs

  position.pattern =  "Position\\d+"
  cell.command = "~/Software/cellID-linux/cell"
  channels = c("b", "f")

  # Delete outputs
  dir(path = path, pattern = "out\\.tif$", full.names = T) %>% unlink()
  dir(path = path, pattern = "^Position\\d\\d\\d$", full.names = T) %>% unlink(recursive = T)

  ####
  BF.images <- cell.args$b %>% {.[order(str_extract(., position.pattern))]} # Order FPs by position
  BF.positions <- BF.images %>% str_extract(position.pattern)  # Grab their position identifier
  n_positions <- length(BF.images)

  cell.args$f <- cell.args$f %>% {.[order(str_extract(., position.pattern))]}  # Order FPs by position
  FP.images <- cell.args$f
  FP.positions <- FP.images %>% str_extract(position.pattern)  # Grab their position identifier
  names(cell.args$f) <- FP.positions

  cell.args$p <- rep(cell.args$p[1], length(FP.images)) # Repeat parameters as many times as needed for FP.images

  # Create output directories
  for(k in 1:length(cell.args$b)) dir.create(cell.args$o[k], recursive = T)

  # names(cell.args$o) <- cell.args$o %>% str_extract(position.pattern)  # Name the output vector with their positions
  # cell.args$o <- cell.args$o[FP.positions]  # Repeat output as many times as needed for FP.images

  names(cell.args$b) <- cell.args$b %>% str_extract(position.pattern)  # Name the BF vector with their positions
  cell.args$b <- cell.args$b[FP.positions]  # Repeat BFs as many times as needed for FP.images

  cell.args.tmp <- c(p = parameters)


  for(j in channels) {
    tmp <- tempfile()
    paths <- cell.args[[j]] # [names(cell.args[[j]]) %in% BF.positions[i]]
    write(x = paths, file = tmp)  # readLines(tmp)
    cell.args.tmp[j] <- tmp
  }

  command <- paste(cell.command,
                   paste0("-",names(cell.args.tmp), " ", cell.args.tmp,
                          collapse = " ")
  )

  command <- strsplit(command, " ")[[1]]
  command[1] <- "cell"
  command <- paste(command, collapse = " ")
}

command

#source("tests/cellid_wrapper.test2.R")
cellid(command)

cell.data <- cell.load(path = path,
                       pdata = path.pdata)
