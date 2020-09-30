#*************************************************************************#
## tidycell: R package for analysis of cellID datasets inside the tidyverse
#*************************************************************************#

# cellID returns table with one row per cell per channel. So, the same
# cell appears "repeated" n times if n channels were used. In general,
# n <= 3 (TFP, RFP, YFP, CFP; in no particular order).
# The variable "flag" (0-3) differentiates the channels.
# The file out_bf_fl_mapping maps a XFP image to its flag number, and
# to its corresponding BF image.

##################### Package Constants #################################
.conflicts.OK = TRUE
.CELLID_ID_VARS = c("pos", "t.frame", "cellID")
.CELLID_ID_VARS_DERIV = c(.CELLID_ID_VARS, "ucid", "time")
.CELLID_DROP_VARS = c("flag", "num.pix", "con.vol.1")

if(getRversion() >= "2.15.1") {
    utils::globalVariables(c("a.tot",
                             "bright",
                             "cellID",
                             "channel",
                             "ellipse.perim",
                             "f.bg.c",
                             "f.bg.r",
                             "f.bg.y",
                             "f.c",
                             "f.r",
                             "f.tot.c",
                             "f.tot.r",
                             "f.tot.y",
                             "f.y",
                             "flag",
                             "fluor",
                             "maj.axis",
                             "min.axis",
                             "perim",
                             "pos"))
}



#*************************************************************************#
## TODO
#*************************************************************************#

# ToDo: documentation on cell.data object

# ToDo: add and check FRET analysis related functions:
#       - .restructure.split.image
#       - .append.identifier

#*************************************************************************#
## DEPENDENCIES
#*************************************************************************#

# do not load plyr after dplyr! namespace conflicts; we prefer dplyr.
#library(plyr)
#library(tibble)
#library(dplyr)
#library(readr)

# g: I want to use skimr as summary for cell.data objects
# library(skimr)

#*************************************************************************#
## FUNCTIONS
#*************************************************************************#


#' Load cellID data
#'
#' \code{load_cell_data} searches a specified directory (the working directory by default)
#' for folders that match a customizable pattern, usually PositionXXX where XXX is the
#' position number. This folders should contain the Cell-ID output files output_all
#' and the output_bf_fl_mapping for each position. The function
#' loads this files and generates a data structure suitable for filtering and
#' plotting. The function
#' returns a cell.data object that contains all the required
#' information for the analysis. All the functions included in the package
#' operate over this object, and its components should not be modified directly,
#' but through the provided functions. Remember to assign the returned value to a
#' variable (e.g. X<-load.cellID.data() )
#'
#' @param path string, path to folder containing PositionXXX folders.
#' @param pattern string, Regex specifying 'PositionXXX' name format.
#' @param basename string, cellID output extension.
#' @param select character vector, defines which variables to include in the cell.data object
#' @param exclude character vector, defines which variables to exclude in the cell.data object.
#' @param load.vars string, character specifying which variables or group of variables of the Cell-ID
#'     out_all file should be loaded.
#' @param split.image boolean, indicates if the images are split and upper cells should be matched
#'     to lower cells. Set to TRUE if analyzing a FRET split image experiment.
#'
#' @details
#' Reads Cell ID output files (basename)_all in folders that match pattern
#' in path and loads them into a cell.data object.
#'
#' It searches for the output_all files in folders of the form specified by
#' pattern (regular expression). If the folder has a numeric value in its name
#' that number is taken as the position index (for example pos01 is given the index 1)
#' If no numeric value is found in the folder name, then a ordinal index is assign.
#'
#' Possible values for load.vars are 'all', 'fl' or 'fluorescence',
#' 'bg' or 'background', 'calc', 'morph' or 'morphological', 'vac' or 'vacuole',
#' 'nucl' or 'nuclear', 'disc'. The group of variables can be specified in either a positive
#' form (i.e. '+fl+bg+morph') or in a negative form (i.e. '-nucl-vac').
#' Combination of positive and negative form is not allowed.
#' A character vector containing the variables names of the out_all file is
#' also allowed. The selection of variables is done before restructuring, so the
#' variable names should correspond to those of the out_all files. Using this argument can be useful
#' if memory issues arise.
#'
#' Alternatively \code{select} and \code{exclude} can be used to subset the dataset.
#' This arguments are applied after the reshaping,
#' so variables names as returned by \code{summary.cell.data} are used. Wildcard patterns (e.g. 'f.*.y')
#' and keywords (e.g. 'all', 'id.vars', 'YFP', etc.) can be used as components of these arguments.
#'
#' @return a cell.data object
#' @export
#'
#' @examples
load_cell_data <-
    function(path = getwd(),
             pattern = "^[Pp]{1}os[:alpha:]*[:digit:]*",
             basename = "out",
             select = NULL,
             exclude = NULL,
             load.vars = "all",
             split.image = FALSE) {

        on.exit(gc())


        # some variable definitions.
        .conflicts.OK = TRUE
        .CELLID_ID_VARS = c("pos", "t.frame", "cellID")
        .CELLID_ID_VARS_DERIV = c(.CELLID_ID_VARS, "ucid", "time")
        .CELLID_DROP_VARS = c("flag", "num.pix", "con.vol.1")

        # transform path to absolute and canonical path in a OS-dependant way
        path <- normalizePath(path)

        # HERE IS THE POSTA
        # Searching for folders that match arg. 'pattern'
        # this makes a char vector with the names of all the position folders
        posdir = dir(pattern = pattern, path = path)
        print(posdir)


        #######################
        ## initialize variables
        #######################

        # vector with the loaded positions, for output
        loaded.pos = c()

        # list with the loaded position directory
        loaded.pos.dir = list()

        # data frame with the position, flag, ch name, number of frames for that flag, and is.bf
        flag.table = data.frame()

        #ToDo: ?
        data <- c()

        # each element corresponds to a single position data
        pos.data <- list()

        # mapping between a bf image to its corresponding fluorescence images.
        # each element corresponds to a position
        bf.fl.mapping <- list()

        # count and array to correct numbers if position folders don't have them
        # count = 0
        # posdir.index = array(-1, dim = c(length(posdir)))

        # variable to assert all output_all have the same columns
        column.names = c()

        ######## ASSERT ########
        #checking if there are Pos folders to be loaded
        if(length(posdir) == 0) {
            stop("No Pos folder found in specified path or working directory.")
        }
        ########

        cat("reading positions...\n")

        for(i in 1:length(posdir)) {

            curr_dir = paste(path, "/", posdir[i], "/", sep = "")


            ######## ASSERT
            # check that current path is a folder, and assign cellID out filenames
            if(!file.info(curr_dir)$isdir) {stop("curr_dir is not a directory")}
            ########


            # paths for out_* files defined
            out_all <- paste(curr_dir, basename, "_", "all", sep = "")
            out_mapping <- paste(curr_dir, basename, "_bf_fl_mapping", sep = "")


            ######## ASSERT
            # check that file exists
            if(!file.exists(out_all)) {
                stop(paste("Missing file:", out_all, "\n Position not loaded\n"))}
            ########


            # rcell had a check for folders with no number, and assigned an ordinal if there
            # wasn't one. I found no use for it, so I removed it


            ####################
            ## reading data
            ####################

            # print position being processed
            cat(gsub("[a-zA-Z_]", "", posdir[i])," ")
            if(i %% 10 == 0) cat("\n")

            pos.data[[i]] <- readr::read_tsv(out_all, col_types = readr::cols())

            # ToDo: check if there's a difference between using Hmisc::import.cleanup or not.


            ######## ASSERT ########
            # asserting that previously loaded positions have the same number and column names
            if(length(column.names) == 0){ #first position
                column.names <- names(pos.data[[i]])

            } else { #not first position
                curr_names <- names(pos.data[[i]])

                if(length(curr_names) != length(column.names)) {
                    stop(out_all," has different number of colums than previous position\n")
                }

                if(sum(column.names == curr_names) != length(column.names)) {
                    stop(out_all," has different column names than previous positions\n")
                }
            }
            ########

            # updates contents of positions loaded.
            loaded.pos <- c(loaded.pos, i)
            loaded.pos.dir[[i]] <- posdir[i]

            #reading output_bf_fl_mapping
            if(file.exists(out_mapping)) {

                # assuming that all mapping files have the same column types
                bf.fl.mapping[[i]] <- readr::read_tsv(out_mapping,
                                                      col_types = "ciici")


                #creating flag table
                pos.flag <- .mk_flag_table(bf.fl.mapping[[i]], pos = i)
                flag.table <- dplyr::bind_rows(pos.flag, flag.table)

            } else warning(out_mapping, "not found")
        }


        ####################
        # creating variables
        ####################

        # creates the variables:
        #    - pos: position number associated to all cells in a position df.
        #    - ucid: unique number associated to each cell
        #    - qc: "quality control", used for filtering


        cat("\ncreating variables...\n")

        for (ipos in loaded.pos) {
            pos.data[[ipos]] <- dplyr::mutate(pos.data[[ipos]],
                                              pos = ipos,
                                              ucid = ipos * 1e6 + cellID,
                                              qc = T)
        }



        ch.names <- .select_channel_names(flag.table)

        ######################
        #Restructuring the data
        ######################


        #ToDo: add selected / removed variables in a tidy way

        rename.non.f = FALSE

        #channels that wont apper in restructured data
        drop.names = .CELLID_DROP_VARS

        #channels that are not renamed, in channel specific manner
        keep.names = .CELLID_ID_VARS_DERIV

        ch.levels = levels(flag.table$channel)
        ch.num = length(ch.levels)

        ######### ASSERT #########
        #Asserting channel names argument
        if(length(ch.names) == 0) {
            ch.names = ch.levels
        } else if (length(ch.names) != ch.num){
            warning("ch.names should have as many elements as channels in the experiment\n",
                    "ch.names=",
                    paste(ch.names),
                    "\n channels=",
                    paste(ch.levels),
                    "\n",
                    "ignoring argument")

            ch.names <- ch.levels
        }
        #########

        ######### ASSERT #########
        #Asserting load.vars
        if(length(load.vars) == 0){
            warning("Loading all variables")
            load.vars <- "all"
        }
        #########

        #Selecting variables to load
        if(length(load.vars) == 1){
            # Rcell::.parse_load_vars
            load.vars <- .parse_load_vars(load.vars,
                                          vars.all = names(pos.data[[loaded.pos[1]]]))

        } else {
            cat("loading variables ", toString(load.vars))
        }

        n.data <- union(union(keep.names,
                              drop.names),
                        load.vars)

        old.ch.header <- c()

        main.header <- c()

        ch.header <- list()

        for (i in ch.levels) {
            ch.header[[i]] = character()
        }

        #generating columns names vector
        for (i in 1:length(n.data)) {

            if(!is.element(n.data[i], keep.names)) {

                if(substr(n.data[i], 1, 2) == "f." |
                   length(grep("[:graph:]*nucl[:graph:]*", n.data[i])) > 0) {

                    # changes the var name
                    old.ch.header <- c(old.ch.header, n.data[i])

                    for (j in 1:ch.num){

                        ch.header[[ch.levels[j]]] <- c(ch.header[[ch.levels[j]]],
                                                       paste(n.data[i], ".", ch.names[j], sep=""))
                    }

                } else if (!is.element(n.data[i], drop.names)){

                    #changes and keeps the name
                    if(rename.non.f) {

                        old.ch.header <- c(old.ch.header, n.data[i])

                        for (j in 1:ch.num)

                            ch.header[[ch.levels[j]]] <- c(ch.header[[ch.levels[j]]],
                                                           paste(n.data[i],".", ch.names[j], sep=""))
                    }

                    main.header <- c(main.header, n.data[i])
                }

            } else
                #keeps the var name unchange
                main.header <- c(main.header, n.data[i])
        }

        output.names <- main.header

        for(i in ch.levels) {
            output.names <- c(output.names, ch.header[[i]])
        }

        data <- c()

        cat("restructuring positions...\n")

        icount <- 0

        for (ipos in loaded.pos) { #loopingin through positions
            posout <- c() #output for this position
            icount <- icount + 1

            cat(formatC(ipos, width = 3), " ")

            if(icount %% 10 == 0) cat("\n")

            #getting flag for each channel in this position
            ch.flag <- subset(flag.table, pos == ipos)$flag

            #using the channel with more t.frames as main channel
            main.flag.index <- which.max(subset(flag.table,pos==ipos)$frame.n)

            curr.pos.data <- subset(pos.data[[ipos]],
                                    flag == ch.flag[main.flag.index],
                                    select = main.header)

            #for(ich in 1:length(ch.flag)){
            for(ich in ch.levels) {
                curr.ch.pos.data <- subset(pos.data[[ipos]],
                                           flag == with(flag.table, flag[channel == ich & pos == ipos]),
                                           select = c(.CELLID_ID_VARS, old.ch.header))

                names(curr.ch.pos.data) <- c(.CELLID_ID_VARS, ch.header[[ich]])

                curr.pos.data <- dplyr::left_join(curr.pos.data,
                                                  curr.ch.pos.data,
                                                  by = c(.CELLID_ID_VARS))
            }

            pos.data[[ipos]] <- curr.pos.data
        }

        cat("\n\n")

        ######### ASSERT #########
        #checking the number of columns after reshaping
        colNum.pos.data <- unlist(lapply(pos.data, function(x) dim(x)[2]))

        if(length(unique(colNum.pos.data)) > 1){
            print(data.frame(variables = colNum.pos.data))
            stop("Positions have different number of variables after reshaping.")
        }
        #########


        # pos.data is a list with a data frame corresponding to each position
        # so, we merge it into a single one by binding rows.
        pos.data <- dplyr::bind_rows(pos.data)


        #################################################################
        # adding variables
        #################################################################

        # geometric variables

        # ellipse.perim = perimeter of theoretical ellipse, calculated using each
        # cell's axis values.

        # el.p = ratio of ellipse perim over the perimeter measured by cellID.
        # If this number is small ( < ~0.7) it's probably not a cell.
        cat('Creating additional variables:\nellipse.perim\nel.p\n')

        pos.data <- dplyr::mutate(pos.data,
                                  ellipse.perim = pi *
                                      (3 * (maj.axis / 2 + min.axis / 2) -
                                           sqrt((3 * maj.axis / 2 + min.axis / 2) *
                                                    (maj.axis / 2 + 3 * min.axis / 2))),

                                  el.p = ellipse.perim / perim)

        # fluorescence variables

        # f.x = total fluorescence - background for channel x
        # cf.x = concentration of f.x (divided by cell area)

        # check if pictures were taken in each channel

        va <- names(pos.data)

        if ("f.tot.y" %in% va) {
            cat("f.y\ncf.y\n")
            pos.data <- dplyr::mutate(pos.data,
                          f.y = f.tot.y - (a.tot * f.bg.y),
                          cf.y = f.y / a.tot)
        }

        if ("f.tot.c" %in% va) {
            cat("f.c\ncf.c\n")
            pos.data <- dplyr::mutate(pos.data,
                          f.c = f.tot.c - (a.tot * f.bg.c),
                          cf.c = f.c / a.tot)
        }

        if ("f.tot.r" %in% va) {
            cat("f.r\ncf.r\n")
            pos.data <- dplyr::mutate(pos.data,
                          f.r = f.tot.r - (a.tot * f.bg.r),
                          cf.r = f.r / a.tot)
        }



        #################################################################
        # Removing duplicates
        #################################################################

        if (identical(pos.data$con.vol, pos.data$con.vol_1)) {
            cat("\nremoving duplicate con.vol\n")
            pos.data <- dplyr::select(pos.data, -con.vol_1)
        }


        #################################################################
        # g: read pdata if it exists
        #################################################################

        pdata_file <- list.files(path = path, pattern = ".*pdata.csv$")

        if (length(pdata_file == 1)) {
            cat("\nJoining pdata!\n\n")
            pdata <- file.path(path, pdata_file)
            pdata <- readr::read_csv(pdata)
            
            if(!"pos" %in% names(pdata)) {
                error_string <- paste0(
                    "Error: aborting because position column 'pos' was not found in pdata file", 
                    " '", normalizePath(paste0(path, "/", pdata_file)), "'. ",
                    "Manually inspect your pdata file, and check it is correctly formatted as a CSV."
                )
                stop(error_string)
            }
            
            pos.data <- dplyr::left_join(pos.data, pdata, by ="pos")
        } else if (length(pdata_file > 1)) {
            cat("\n MULTIPLE PDATA FILES IN EXPERIMENT FOLDER! \n\n\n")
        }


        #################################################################
        # g: hasta aca tengo el DF con los datos crudos: pos.data
        #################################################################


        # g: agrego data de imagenes. paths  y eso.

        for(ipos in loaded.pos){
            bf.fl.mapping[[ipos]] <- transform(bf.fl.mapping[[ipos]], pos = ipos)
        }

        # preparing "hard" image information
        image.info = NULL
        for(i.pos in loaded.pos){


            # in original it was plyr::join
            pii = dplyr::left_join(bf.fl.mapping[[i.pos]],
                                   flag.table[flag.table$pos == i.pos, c("flag", "channel")],
                                   by = "flag")

            pii = transform(pii,
                            fluor = gsub("[\\]", "/", fluor),
                            bright = gsub("[\\]", "/", bright))

            pii = transform(pii,
                            image = basename(fluor),
                            path = dirname(fluor))

            #bf as fluor, aca habria que cambiar algo
            # original was data.frame
            piibf = tibble::tibble(pos = i.pos,
                                   t.frame = pii[pii$flag == 0, "t.frame"],
                                   channel = "BF",
                                   image = basename(pii[pii$flag == 0, "bright"]),
                                   path = dirname(pii[pii$flag == 0, "bright"]))
            #stringsAsFactors = FALSE)

            # in original it was rbind
            pii = dplyr::bind_rows(subset(pii, select = c("pos",
                                                          "t.frame",
                                                          "channel",
                                                          "image",
                                                          "path")),
                                   piibf)

            if(is.null(image.info)) {
                image.info = pii

            } else {
                # original was rbind
                image.info = dplyr::bind_rows(image.info, pii)
            }
        }

        # checking if path in bf_fl_mapping is correct, or should replace with new path
        img.fnames = with(image.info[1:5, ], paste(path, image, sep = "/"))
        img.fnames.exist = file.exists(img.fnames)

        # checking if path argument permorms better
        if(!all(img.fnames.exist)){

            img.fnames2 = paste(path, image.info$image[1:5], sep = "/")
            img.fnames2.exist = file.exists(img.fnames2)

            if(sum(img.fnames2.exist) > sum(img.fnames.exist)) { #replacing path
                image.info$path <- factor(path)

                message("tif files moved since analized with Cell-ID, updating path")

            } else {
                message("tif files moved since analized with Cell-ID, can't find them")
            }
        }

        #adding "out" channels
        img.fnames.out = with(image.info[1:5, ], paste(path, image, ".out.tif", sep="/"))

        if(!all(file.exists(img.fnames.out))){
            # original was rbind
            image.info <- dplyr::bind_rows(transform(image.info,
                                                     is.out = FALSE),
                                           transform(image.info,
                                                     image = paste(image, ".out.tif", sep = ""),
                                                     channel = paste(channel, ".out", sep = ""),
                                                     is.out = TRUE))
        }

        # original was data.frame, with stringsAsFactors=FALSE
        channels = tibble::tibble(posfix = ch.names,
                                  name = levels(flag.table$channel))

        variables = list(id.vars = .CELLID_ID_VARS,
                         id.vars.deriv = .CELLID_ID_VARS_DERIV,
                         morpho = unique(c(setdiff(main.header, c(.CELLID_ID_VARS_DERIV,"qc")),
                                           grep(glob2rx("a.*"), names(pos.data), value = TRUE))),
                         fluor = grep(glob2rx("f.*"), names(pos.data), value = TRUE),
                         qc = "qc",
                         as.factor = c("pos", "cellID", "ucid"),
                         all = names(pos.data))

        for(i in 1:dim(channels)[1])
            variables[[channels[[i,"name"]]]] <- grep(glob2rx(paste("*.",
                                                                    channels[[i, "posfix"]],
                                                                    sep = "")),
                                                      names(pos.data),
                                                      value = TRUE)


        cell.data=
            list(data = pos.data,
                 qc.history = list(),
                 subset.history = list(),
                 transform = list(),
                 channels = channels,
                 variables = variables,
                 images = image.info,
                 software = "Cell-ID",
                 load.date = date(),
                 load.sessionInfo = sessionInfo()
            )

        class(cell.data) <- c("cell.data", "list")

        if(!is.null(select) || !is.null(exclude)) {
            cell.data = subset(cell.data, select = select, exclude = exclude)
        }


        # ToDo: add and check FRET analysis related functions:
        #       - .restructure.split.image
        #       - .append.identifier
        #if(isTRUE(split.image)){
        #    cell.data <- .restructure.split.image(cell.data)
        #}

        # print(summary(cell.data))
        return(cell.data)


    }


#*************************************************************************#
# PUBLIC Functions
#*************************************************************************#

#ToDo: summary.cell.data
#ToDo: transform.cell.data
#ToDo: merge.cell.data






#*************************************************************************#
# PRIVATE Functions
#*************************************************************************#
#private

#' Make Flag Table
#'
#' Generates a table mapping a channel name (3 first characters of the image file)
#' to a flag number for a given bf.fl.mapping data.frame (read from an
#' output_bf_fl_mapping file)
#'
#' @param bf.fl.mapping data.frame of the bf.fl.mapping output file
#' @param pos integer, corresponds to position
#'
#' @return a data.frame containing the bf.fl.mapping of a single position
#'
#' @examples
.mk_flag_table <- function(bf.fl.mapping, pos = NULL){

    flag.name = vector(mode = "character", length = 0)
    flag = c()
    flag.frame.n = c()
    flag.is.bf = c()
    flag.count = 0
    output = data.frame()


    # g: for each fluorescence image (rows in out_bf_fl_mapping)
    for(i in 1:nrow(bf.fl.mapping)){

        # g: separates path to fluorescence image into its components
        part.path <- strsplit(as.character(bf.fl.mapping[i, 1]), "[/\\]")[[1]]

        tmpstr <- substr(part.path[length(part.path)], 1, 3)

        flag.name.index = which(flag.name == tmpstr)


        if(length(flag.name.index) == 0){ #new flag
            tmpflag = as.integer(bf.fl.mapping[i, 2])
            if(length(which(flag == tmpflag)) == 0){#cheking consistency
                flag = c(flag, tmpflag)
                flag.name = c(flag.name, tmpstr)
                flag.frame.n = c(flag.frame.n, 1)
                flag.count = flag.count + 1

                if(tmpstr == substr(bf.fl.mapping[i, 4], 1, 3)){
                    #if(bf.fl.mapping[i,5]==1){
                    flag.is.bf = c(flag.is.bf,TRUE)
                } else {
                    flag.is.bf = c(flag.is.bf,FALSE)
                }
            } else {
                cat(".mk_flag_table: Flag name ambiguity.\n")
            }
        } else if(length(flag.name.index) == 1) { #flag all ready assing
            flag.frame.n[flag.name.index] = flag.frame.n[flag.name.index] + 1
        } else {
            cat(".mk_flag_table: Ambiguous flag name\n")
        }

    }

    output = data.frame(flag = flag,
                        channel = flag.name,
                        frame.n = flag.frame.n,
                        is.bf = flag.is.bf, stringsAsFactors = T)

    if(!is.null(pos)){
        output = data.frame(pos = rep(pos, flag.count), output, stringsAsFactors = T)
    }

    return(output)
}


#*************************************************************************#
#private
#ToDo: improve this, compatibilize with select
#' Parse Names of Variables to Load
#'
#' Parses the input of load.vars and reurns a vector with the elements to be loaded.
#' Possible values for load.vars are 'all', 'fl' or 'fluorescence', 'bg' or 'background', 'calc',
#' 'morph' or 'morphological', 'vac' or 'vacuole', 'nucl' or 'nuclear', 'disc'.
#' The group of variables can be specified in either a positive form (i.e. '+fl+bg+morph')
#' or in a negative form (i.e. '-nucl-vac'). Combination of positive and negative form is not allowed.
#'
#' @param load.vars, pattern of variable names in code
#' @param vars.all, NULL
#'
#' @return character vector containing variable names
#'
#' @examples
.parse_load_vars <- function(load.vars, vars.all = NULL){

    if(length(load.vars) != 1) stop(".parse_load_vars argument should be of length 1\n")

    vars.nucl <- c("f.nucl",
                   "a.nucl",
                   "f.nucl1",
                   "f.nucl.tag1",
                   "a.nucl1",
                   "f.nucl2",
                   "f.nucl.tag2",
                   "a.nucl2",
                   "f.nucl3",
                   "f.nucl.tag3",
                   "a.nucl3",
                   "f.nucl4",
                   "f.nucl.tag4",
                   "a.nucl4",
                   "f.nucl5",
                   "f.nucl.tag5",
                   "a.nucl5",
                   "f.nucl6",
                   "f.nucl.tag6",
                   "a.nucl6",
                   "f.nucl7",
                   "f.nucl.tag7",
                   "a.nucl7",
                   "f.nucl8",
                   "f.nucl.tag8",
                   "a.nucl8")

    vars.nucl2 <- c("f.nucl.tag1",
                    "f.nucl2",
                    "f.nucl.tag2",
                    "a.nucl2",
                    "f.nucl3",
                    "f.nucl.tag3",
                    "a.nucl3",
                    "f.nucl4",
                    "f.nucl.tag4",
                    "a.nucl4",
                    "f.nucl5",
                    "f.nucl.tag5",
                    "a.nucl5",
                    "f.nucl6",
                    "f.nucl.tag6",
                    "a.nucl6",
                    "f.nucl7",
                    "f.nucl.tag7",
                    "a.nucl7",
                    "f.nucl8",
                    "f.nucl.tag8",
                    "a.nucl8")

    vars.vac <- c("a.vacuole", "f.vacuole")

    vars.morph <- c("xpos",
                    "ypos",
                    "a.tot",
                    "num.pix",
                    "fft.stat",
                    "perim",
                    "maj.axis",
                    "min.axis",
                    "rot.vol",
                    "con.vol",
                    "a.surf",
                    "con.vol.1",
                    "sphere.vol")

    vars.fl <- c("f.tot", "a.tot")

    vars.disc <- c("f.tot.p1",
                   "a.tot.p1",
                   "f.tot.m1",
                   "a.tot.m1",
                   "f.tot.m2",
                   "a.tot.m2",
                   "f.tot.m3",
                   "a.tot.m3")

    vars.bg <- c("f.bg",
                 "f.local.bg",
                 "local.bg.num",
                 "local.num",
                 "f.local2.bg",
                 "local2.bg.num",
                 "local2.num")

    if(is.null(vars.all)) {vars.all <- c(vars.bg, vars.fl, vars.morph, vars.vac, vars.nucl)}

    has.plus <- length(grep("[+]", load.vars)) > 0
    has.minus <- length(grep("[-]", load.vars)) > 0

    if(has.plus & has.minus) {
        stop("invalid sintaxis for load.vars")
    }

    if(!has.plus & !has.minus) {
        has.plus <- TRUE
    }

    output = character(0)

    if(has.minus) {
        output = vars.all
    }

    for(i in strsplit(load.vars, split = "[+-]")[[1]]) {
        if(i != "") {
            vars.i <- switch(i,
                             nucl2 = vars.nucl2,
                             nuc2 = vars.nucl2,
                             nucl = vars.nucl,
                             nuc = vars.nucl,
                             nuclear = vars.nucl,
                             vac = vars.vac,
                             vacuole = vars.vac,
                             morph = vars.morph,
                             morphological = vars.morph,
                             fl = vars.fl,
                             fluorescence = vars.fl,
                             bg = vars.bg,
                             background = vars.bg,
                             all = vars.all,
                             disc = vars.disc
            )

            if(is.null(vars.i)) {
                if(is.element(load.vars, vars.all)) {
                    vars.i <- load.vars
                } else {
                    stop("Invalid value for load.vars")
                }
            }

            if(has.plus) {
                output = union(output, vars.i)
            } else {
                output = setdiff(output, vars.i)
            }
        }
    }

    return(output)
}


#*************************************************************************
# private
# selects the channel names from flag.table

.select_channel_names <- function(flag.table) {
    #selecting proper names for the channels
    #atempting to use first letter of channel identifier
    i <- 1
    while(i <= 3){
        ch.names <- substr(levels(flag.table$channel), 1, i)

        #note that ch.names and levels(flag.table$channel) will have the same order

        if (sum(is.na(pmatch(ch.names, levels(flag.table$channel)))) == 0){
            i = 3
            ch.names = tolower(ch.names)

        } else if (i == 3){
            #should never get here
            ch.names = c()
        }
        i = i + 1
    }

    ch.names
}


#*************************************************************************#
#private
#select variables for subsetting
.select <- function(variables,
                    select = NULL,
                    exclude = NULL,
                    warn = TRUE){
    #expanding select
    exp.select = c()
    for(i in select){
        if(substr(i,1,1) == "-"){
            # g: Rcell::.nchar
            exclude <- c(exclude, substr(i, 2, .nchar(i)))

        } else {
            ms <- intersect(i, names(variables))
            if(length(ms) == 1) exp.select <- c(exp.select, variables[[ms]])
            else {
                ms <- grep(glob2rx(i), variables$all, value = TRUE)
                if(length(ms) == 0 && !(select %in% c("", "none")) && warn){
                    warning("unknown selected variable ", i)
                }
                exp.select <- c(exp.select, ms)
            }
        }
    }
    exp.select <- na.omit(unique(exp.select))

    #expanding exclude
    exp.exclude <- c()
    for(i in exclude){
        me <- intersect(i, names(variables))
        if(length(me) == 1) {
            exp.exclude <- c(exp.exclude,variables[[me]])
        } else {
            me = grep(glob2rx(i), variables$all, value = TRUE)
            if(length(me) == 0 && warn) {
                warning("unknown excluded variable ",i)
            }
            exp.exclude = c(exp.exclude, me)
        }
    }
    exp.exclude <- na.omit(unique(exp.exclude))

    if (length(select) == 0 & length(exp.exclude) == 0) {
        return(TRUE)
    } else if (length(select) == 0 & length(exp.exclude) > 0) {
        output <- setdiff(variables$all, exp.exclude)
        return(output[!is.na(output)])
    } else {
        output <- setdiff(exp.select, exp.exclude)
        return(output[!is.na(output)])
    }
}


#*************************************************************************#
#private
#workaround to the change in nchar behavior introduced in R 3.3.0
.nchar<-function(x, type = "chars", allowNA = FALSE){
    if(getRversion() <= "3.2.0"){
        return(nchar(x, type, allowNA))
    } else {
        return(nchar(x, type, allowNA, keepNA = FALSE))
    }
}

#*************************************************************************#
#public
#
#' write tab delimited file
#'
#' Writes a tab delimited file. Wrapper to write.table
#'
#' @param x data frame to write
#' @param file file to output
#' @param quote  a logical value (TRUE or FALSE) or a numeric vector. If TRUE, any character
#'  or factor columns will be surrounded by double quotes. If a numeric vector, its elements
#'  are taken as the indices of columns to quote. In both cases, row and column names are
#'  quoted if they are written. If FALSE, nothing is quoted.
#' @param sep the field separator string. Values within each row of x are separated by this string.
#' @param row.names either a logical value indicating whether the row names of x are to be written
#'  along with x, or a character vector of row names to be written.
#' @param ... arguments to write.table: append, col.names, sep, dec and qmethod cannot be altered.
#'
#' @export
#'
#' @examples
write.delim <- function(x,
                        file = "",
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        ...) {
    write.table(x,
                file = file,
                quote = quote,
                sep = sep,
                row.names = row.names,
                ...)
    }


#ToDo: function within.cell.data
#*************************************************************************#
#private
#' Format sequence of numbers
#'
#' formats sequence of numbers in an short expresion  eg: 1-10, 12-15
#'
#' @param pos last number in sequence to shorten
#'
#' @return string of shortened sequence
#' @export
#'
#' @examples
.format.sequence <- function(pos) {
    if(length(pos) < 2) return(as.character(pos))
    else{
        fs = as.character(pos[1])
        last.pos = pos[1]
        for(i in 2:length(pos)){
            if(last.pos + 1 != pos[i]){
                fs = paste(fs, "-", last.pos, ",", pos[i], sep = "")
                last.pos = pos[i]
            } else if(i == length(pos)) {
                fs = paste(fs, "-", pos[i], sep = "")
            }else{
                last.pos = last.pos + 1
            }
        }
        return(fs)
    }
}

#*************************************************************************#
#public
#
#' is cell data
#'
#' checks if an objects is a cell.data object
#'
#' @param X an object
#'
#' @return boolean indicating if X is a cell.data object
#' @export
#'
#' @examples
is.cell.data <- function(X) {
    inherits(X, "cell.data")
}

#*************************************************************************#
