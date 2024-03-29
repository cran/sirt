## File Name: lsem_fitsem_compute_sufficient_statistics.R
## File Version: 0.109

lsem_fitsem_compute_sufficient_statistics <- function(G, dat, variables_model,
    weights, moderator_variable=NULL, loc_linear_smooth=NULL, moderator.grid=NULL,
    pd=FALSE, residualized_intercepts=NULL,    has_meanstructure=FALSE,
    residualize=TRUE)
{
    wmean <- wcov <- Nobs <- as.list(1:G)
    data_suff <- dat[, variables_model]
    dat_resp <- 1 - is.na(data_suff)
    for (gg in 1:G){
        weights_gg <- weights[,gg]
        # res <- lsem_weighted_mean( x=data_suff, weights=weights_gg, x_resp=dat_resp)
        res <- lsem_weighted_cov( x=data_suff, weights=weights_gg, x_resp=dat_resp,
                    moderator_variable=moderator_variable,
                    loc_linear_smooth=loc_linear_smooth,
                    moderator_value=moderator.grid[gg], pd=pd,
                    residualized_intercepts=residualized_intercepts,
                    has_meanstructure=has_meanstructure, residualize=residualize)
        wmean[[gg]] <- res$mean
        wcov[[gg]] <- res$cov
        Nobs[[gg]] <- round(res$Nobs)
    }

    #** adapt if mean structure is requested
    if ( has_meanstructure & residualize ){
        for (gg in 1:G){
            wmean[[gg]] <- residualized_intercepts[gg,]
        }
    }

    #- output
    res <- list(wmean=wmean, wcov=wcov, Nobs=Nobs)
    return(res)
}
