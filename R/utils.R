.Stringtiebin <- function(args = "") {
  if (is.null(args) || args == "") {
    stop("The stringtie binaries nedd to called",
         " with additional arguments")
  }
  args <- gsub("^ *| *$", "", args)
  args <- unlist(strsplit(args, split = " "))
  bin <- file.path(system.file(package = "Rstringtie"), "stringtie")
  output <- system2(command = bin, args = args)
}

.gffcompareBin <- function(args = "") {
  if (is.null(args) || args == "") {
    stop("The stringtie binaries nedd to called",
         " with additional arguments")
  }
  args <- gsub("^ *| *$", "", args)
  args <- unlist(strsplit(args, split = " "))
  bin <- file.path(system.file(package = "Rstringtie"), "gffcompare")
  output <- system2(command = bin, args = args)
}


.parseDots <- function(...) {
  if (...length() == 0) return("")
  dots <- list(...)
  args <- lapply(names(dots), FUN = function(x) {
    paste(x, dots[x], sep = " ")
  })
  args <- paste0(unlist(args), collapse = " ")
  return(args)
}


