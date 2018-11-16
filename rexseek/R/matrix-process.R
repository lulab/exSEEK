#matrix_path = 'inst/extdata/scirep_sequential_qc.txt'
#classinfo_path = 'scirep_classes.txt'

#' @title read count matrix
#'
#' @param path string.
#' @param ... other arguments passsed on to [readr::read_tsv()]
#'
#' @return integer matrix
#' @export
#'
#' @examples NULL
read_mat <- function(path, ...) {
	path %>% readr::read_tsv(T, readr::cols(.default = 'c'), ...) %>%
		dplyr::mutate_at(-1, readr::parse_integer) %>%
		dplyr::rename('transcript' = 1) %>% as.data.frame() %>%
		tibble::column_to_rownames('transcript') %>% as.matrix()
}

#' Title
#'
#' @param mat integer matrix.
#' @param thres_count integer scalar.
#' @param thres_sample_nums integer scalar.
#'
#' @return scater::calculateQCMetrics()
#' @export
#'
#' @examples NULL
filter_low <- function(mat, thres_count = 2, thres_sample_nums = 5) {
	low_per_row <- rowSums(mat > thres_count)
	keeped_row <- low_per_row > thres_sample_nums
	mat[keeped_row, ]
}

#' @export
plot_highest_exprs <- function(sce) {
	sce %>% {suppressMessages(scater::calculateQCMetrics(.))} %>%
		scater::plotHighestExprs(n = 20)
}

#' @export
plot_PCA <- function(sce) {
	sce %>% scater::plotPCA(
		shape_by = "label", colour_by = "label",
    	run_args = list(exprs_values = "counts")
	)
}

#' @export
plot_TSNE <- function(sce) {
	sce %>% scater::plotTSNE(
	    shape_by = "label", colour_by = "label",
	    run_args = list(exprs_values = "counts")
	)
}





