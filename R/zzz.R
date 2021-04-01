.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    splitline()
  )
  packageStartupMessage(
    message("Welcome to My package !")
  )
  packageStartupMessage(
    message("Use `load_spat_env Firstly`")
  )
  packageStartupMessage(
    splitline()
  )
}
