#' Import data from single EMPOWER2 HPLC files
#'
#' @description
#' reads output files of the EMPOWER 2 SOFTWARE (Waters). Input files must contain data of single samples deposited within the same directory.
#'
#' @param path path to a folder containing input files
#'
#' @param pattern pattern used to select files. By default ".txt"
#'
#' @param skip rows to skip before reading data
#'
#' @param id column containing sample name
#'
#' @inheritParams align_chromatograms
#'
#' @return a list of data frames (each corresponding to a sample)
#'
#' @keywords beta
#'
#' @export
#'
read_empower2 <- function(path = NULL,
                          pattern = ".txt",
                          sep = "\t",
                          skip = 2,
                          id = "SampleName")
    {
    # checks
    if (is.null(path)) stop("Specify path")

    # get files
    files <- list.files(path, pattern = pattern, full.names = T)

    # get variables
    var_names <- readr::read_lines(files[1], n_max = 1)
    var_names <- unlist(stringr::str_split(string = var_names, pattern = sep))
    var_names <- var_names[var_names != ""]
    var_names <- gsub('\"', "", var_names, fixed = TRUE)
    var_names <- stringr::str_trim(var_names)

    # read data
    read_temp <- function(data) {
        dat <- utils::read.table(data, skip = skip + 1, sep = "\t",
                                 stringsAsFactors = F,
                                 col.names = var_names)
        }

    data <- lapply(files, read_temp)

    # get sample names
    sample_names <- unique(do.call("rbind", data)[,id])
    sample_names <- stringr::str_trim(sample_names)

    names(data) <- sample_names

    return(data)
}

