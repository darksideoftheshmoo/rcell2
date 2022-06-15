#' A function to donwload the latest worflow tempalte in Rmarkdown
#' 
#' Will donwload the .Rmd file to the current working directory.
#' 
#' @param file_name File name for the downloaded wokflow template.
#' @param open.template Try using file.edit to open the file in RStudio after getting it.
#' 
#' @export
get_workflow_template <- function(
    file_name = "rcell2.dev_workflow_template.Rmd",
    open.template = T){
  
  if(file.exists(file_name)) stop("get_workflow_template error: file", file_name, "exists.")
  
  workflow.file <- system.file(
    "rcell2.dev_workflow_template.Rmd",
    package = "rcell2.cellid"
  )
  
  if(file.exists(workflow.file)){
    file.copy(from = workflow.file, to = file_name)
  } else {
    download.file(url = "https://raw.githubusercontent.com/darksideoftheshmoo/rcell2/rcell2.cellid/inst/rcell2.dev_workflow_template.Rmd", 
                  destfile = file_name)
  }
  
  if(open.template){
    file.edit(file_name)
  }
}

