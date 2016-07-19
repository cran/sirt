
#########################################################################
xxirt_partable_extract_freeParameters <- function( partable ){
		partable <- partable[ partable$est , ]
		partable <- partable[ base::order(partable$parindex) , ]
		partable <- partable[ ! base::duplicated( partable$parindex) , ]
		x <- partable$value
		base::names(x) <- partable$parlabel
		base::return(x)
}
#########################################################################