
##########################################
# extension of .Call for easy switching between
# R and Cpp versions
CallSwitch <- function( .NAME , ... , PACKAGE = "R" ){
	if ( PACKAGE == "R" ){
		args1 <- base::list(...)
		base::do.call( .NAME , args1 )	
	} else {
		base::.Call( .NAME , ... , PACKAGE = PACKAGE )
	}
}

.CallSwitch <- CallSwitch