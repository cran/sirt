


xxirt_partable_include_freeParameters <- function( partable , x ){
		vals <- x[ partable$parindex ]
		ind <- base::is.na(vals)
		vals[ ind ] <- partable$value[ ind ]
		partable$value <- vals
		base::return(partable)
}
		