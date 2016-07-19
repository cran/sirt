#  zzz.R
#
# This function is simply copied from the mice package.


# on attach sirt
.onAttach <- function(libname,pkgname){
  d <- utils::packageDescription("sirt")
  d1 <- d$Version 
#  nk <- base::paste( base::rep( " " , 20 - base::nchar(d1) ) , collapse="")
  base::packageStartupMessage(
		base::paste("- " , d$Package," " , d1 ," (",d$Date,")",sep="")  )
	}
	
	
version <- function(pkg="sirt"){
  lib <- base::dirname( base::system.file(package = pkg))
  d <- utils::packageDescription(pkg)
  base::return( base::paste(d$Package,d$Version,d$Date,lib))
}

# .First.lib <- function(lib, pkg){
#          library.dynam("sirt", package = pkg, lib.loc = lib)
#          return(invisible(0))
#        } 