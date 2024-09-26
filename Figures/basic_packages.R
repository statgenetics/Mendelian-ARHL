# R.version
packages_list <- c( "readr", "dplyr", "data.table", "R.utils", "ggplot2", "tidyr")
invisible(lapply(packages_list, library, character.only = TRUE, warn.conflicts = F, quietly = T))
hms_span <- function(start, end) {
    dsec <- as.numeric(difftime(end, start, unit = "secs"))
    hours <- floor(dsec / 3600)
    minutes <- floor((dsec - 3600 * hours) / 60)
    seconds <- dsec - 3600 * hours - 60 * minutes
    time <- paste0(sapply(c(hours, minutes, seconds), function(x) {
        formatC(x, width = 2, format = "d", flag = "0")
    }), collapse = ":")
    print(paste("This process takes", time, "."))
}

