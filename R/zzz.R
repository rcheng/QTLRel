.First.lib <- function(lib, pkg) {
   library.dynam("QTLRel", pkg, lib)
}
.onLoad <- function(lib, pkg) print("R/QTLRel is loaded")
.noGenerics <- TRUE
.onUnload <- function(libpath) library.dynam.unload("QTLRel", libpath)

