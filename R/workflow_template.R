#' A function to donwload the latest worflow tempalte in Rmarkdown
#' 
#' Will donwload the .Rmd file to the current working directory.
#' 
#' @param file_name File name for the wokflow template.
#' 
#' @export
#' @return 
get_workflow_template <- function(file_name = "rcell2_workflow_template.Rmd"){
  if(file.exists("rcell2_workflow_template.Rmd")) stop("get_workflow_template error: file", file_name, "exists.")
  download.file(url = "https://raw.githubusercontent.com/darksideoftheshmoo/rcell2/master/inst/workflow_template.Rmd",
                destfile = file_name)
}
