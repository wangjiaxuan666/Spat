.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    splitline()
  )
  packageStartupMessage(
    messageline("Welcome to My package !")
  )
  packageStartupMessage(
    messageline("Use `load_spat_env` Firstly")
  )
  packageStartupMessage(
    splitline()
  )
}
