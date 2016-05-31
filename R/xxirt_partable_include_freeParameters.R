


xxirt_partable_include_freeParameters <- function( partable , x ){
		vals <- x[ partable$parindex ]
		vals[ is.na(vals) ] <- partable$value[ is.na(vals) ]
		partable$value <- vals
		return(partable)
			}
		