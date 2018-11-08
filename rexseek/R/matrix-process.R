#matrix_path = 'inst/extdata/scirep_sequential_qc.txt'
#classinfo_path = 'scirep_classes.txt'

#' @title read count matrix
#'
#' @param path
#'
#' @return integer matrix
#' @export
#'
#' @examples NULL
read_mat <- function(path) {
	path %>% readr::read_tsv() %>%
		dplyr::rename('transcript' = 1) %>% as.data.frame() %>%
		tibble::column_to_rownames('transcript') %>% as.matrix()
}
