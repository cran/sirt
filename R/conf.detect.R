

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Confirmatory DETECT analysis
conf.detect <- function( data , score , itemcluster , bwscale = 1.1 , progress = TRUE ,
                            thetagrid = seq( -3,3,len=200)  ){
    cat("-----------------------------------------------------------\n" )
    cat("Confirmatory DETECT Analysis \n" ) ; flush.console()
    h1 <- is.matrix( score )
    if (h1 ){ PP <- ncol(score) }
    if (! h1  ){  cat("Conditioning on 1 Score\n" )  } else {
            cat(paste("Conditioning on ",PP, " Scores\n" , sep="") ) }
    cat(paste("Bandwidth Scale:" , bwscale , "\n" ) ) 
    utils::flush.console()
    if ( ! h1 ){
        ccovtable <- ccov.np( data , score = score, bwscale = bwscale , 
                                progress= progress , thetagrid = thetagrid )
        res <- detect.index( ccovtable , itemcluster = itemcluster )
                    } else {
            ccovtable.list <- list()
            for (pp in 1:PP){
                cat( paste( "DETECT Calculation Score " , pp , "\n" , sep="") ) ; 
				utils::flush.console()
                ccovtable.list[[pp]] <- ccov.np( data , score = score[,pp], 
                                    bwscale = bwscale , progress= FALSE )
                    }  
        detect.list <- lapply( ccovtable.list , FUN = function( ccovtable ){ 
                    detect.index( ccovtable , itemcluster=itemcluster ) } )
        detect.matrix <- matrix( unlist( lapply( detect.list , FUN = function( ll){ c( ll[1,] , ll[2,] , ll[3,] ) } ) ) , nrow=PP , byrow=T)
        detect.summary <- data.frame( "NScores" = PP , "Mean" = colMeans( detect.matrix ) , 
                    "SD" = apply( detect.matrix , 2 , stats::sd ) , 
                    "Min" = apply( detect.matrix , 2 , min ) , 
                    "Max" = apply( detect.matrix , 2 , max ) 
                    )
        rownames(detect.summary) <- c("DETECT Unweighted" , "DETECT Weighted" , "ASSI Unweighted" , "ASSI Weighted" ,
                        "RATIO Unweighted" , "RATIO Weighted" )
            }
    cat("-----------------------------------------------------------\n" )            
    if ( ! h1){   res <- list(  "detect" = res , "ccovtable" = ccovtable , "detect.summary" = res ) } else
            {     res <- list(  "detect" = detect.list , "ccovtable" = ccovtable.list , "detect.summary" = detect.summary ) }
    print(round(res$detect.summary,3))
    return( res )
    }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

