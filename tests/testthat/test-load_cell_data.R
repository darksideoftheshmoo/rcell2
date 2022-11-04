mock_data_path <- function(experiment = "") {
  base_path <- ifelse(basename(getwd()) == "rcell2",
                      "./tests/testthat",
                      ".")
  
  file.path(base_path, "mock_data", experiment)
  
}

clean_cdata <- function(cdata) {
  # clean cdata for testing while t.frame bug is not merged
  
  cdata %>%
    dplyr::arrange(cellID) %>%
    dplyr::select(-c(pos, ucid))
  
}

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("load_cell_data loads cell data", {
  experiment <- mock_data_path("sample_time_series")
  res <- load_cell_data(experiment)
  
  clean_data <- clean_cdata(res$data)
  
  tf <- withr::with_tempfile("tf", {
    readr::write_csv(clean_data, tf)
    expect_snapshot_file(tf, name = "clean_cdata.csv")
    })

})