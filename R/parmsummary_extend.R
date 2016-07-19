
####################################################
# extend parameter summary
parmsummary_extend <- function( dfr , level = .95 , 
		est_label = "est" , se_label ="se" ,df_label = "df" )
{
	dfr <- base::as.data.frame(dfr)
	#-- compute t values and p values
	dfr$t <- dfr[,est_label] / dfr[,se_label]
	dfr$p <- 2 * stats::pnorm( - base::abs(dfr$t) )
	#-- compute confidence intervals
	if ( ! base::is.null(level) ){
		if ( "df" %in% base::colnames(dfr) ){
			quant <- - stats::qt( (1-level)/2 , df = dfr$df)
		} else {
			quant <- - stats::qnorm( (1-level)/2 )
		}
		dfr[, paste0("lower" , 100*level ) ] <-
			dfr[,est_label] - quant * dfr[ , se_label ]
		dfr[, paste0("upper" , 100*level ) ] <-
			dfr[,est_label] + quant * dfr[ , se_label ]
	}
	#-- output
	base::return(dfr)
}
#####################################################