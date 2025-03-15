#' @importFrom data.table rbindlist as.data.table fread
#' @import dplyr
#' @import forcats
#' @import ggplot2
#' @importFrom grDevices dev.off pdf
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom httr status_code GET
#' @importFrom inspectdf inspect_na
#' @importFrom jsonlite fromJSON
#' @import knitr
#' @importFrom lubridate parse_date_time
#' @import naniar
#' @import progress
#' @import purrr
#' @importFrom readr read_lines read_delim
#' @importFrom scales percent
#' @importFrom stats median reorder density
#' @import stringr
#' @import tidyr
#' @importFrom utils URLencode read.csv read.delim write.table
#' @import viridis


utils::globalVariables(
  c("sample_order_num", 
    "total", 
    "missing_pct")
  )

.onLoad <- function(libname, pkgname) {
  invisible()
}