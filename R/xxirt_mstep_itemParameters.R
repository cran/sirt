## File Name: xxirt_mstep_itemParameters.R
## File Version: 0.35

##############################################################
# M-step item parameters
xxirt_mstep_itemParameters <- function( partable, item_list, items, Theta,
            ncat, partable_index, N.ik, mstep_iter, par0, eps,
            mstep_reltol, mstep_method, item_index, h, use_grad )
{
    #-------------------------------------------------
    #**** define likelihood function
    like_items <- function( x, ... ){
        partable <- xxirt_partable_include_freeParameters( partable=partable, x=x )
        probs1 <- xxirt_compute_itemprobs( item_list=item_list, items=items, Theta=Theta, ncat=ncat,
                            partable=partable, partable_index=partable_index )
        log_probs1 <- log( probs1 + eps )
        # log likelihood value
        ll <- - sum( N.ik * log_probs1 )
        # add prior distributions
        pen <- xxirt_mstep_itemParameters_evalPrior(partable, h=0)
        # add penalty
        ll <- ll + sum(pen)
        return(2*ll)
    }
    #**** end definition likelihood
    #-------------------------------------------------

    #-------------------------------------------------
    #*** define gradient here
    #### x <- partable[ partable$parfree==1, "value"]
    grad_items <- function(x, ...){
        partable0 <- xxirt_partable_include_freeParameters( partable=partable, x=x )
        probs1 <- xxirt_compute_itemprobs( item_list=item_list, items=items, Theta=Theta, ncat=ncat,
                            partable=partable0, partable_index=partable_index )
        NP <- sum( partable$parfree==1)
        grad1 <- rep(0,NP)
        for (pp in 1:NP){
            item_index_pp <- item_index[[pp]]
            xh <- x
            xh[pp] <- xh[pp] + h
            N.ik_pp <- N.ik[item_index_pp,,, drop=FALSE ]
            partable_h <- xxirt_partable_include_freeParameters( partable=partable0, x=xh )
            probs1h <- xxirt_compute_itemprobs( item_list=item_list, items=items, Theta=Theta, ncat=ncat,
                                partable=partable_h, partable_index=partable_index, item_index=item_index_pp )
            ll0 <- - sum( N.ik_pp * log( probs1[ item_index_pp,,,drop=FALSE] + eps ) )
            ll1 <- - sum( N.ik_pp * log( probs1h + eps ) )
            grad1[pp] <- ( ll1 - ll0 )
        }
        #--- prior distributions
        pen <- xxirt_mstep_itemParameters_evalPrior(partable0, h=0)
        pen_h <- xxirt_mstep_itemParameters_evalPrior(partable0, h=h)
        grad1 <- ( grad1 +  pen_h - pen)/h
        return(grad1)
    }
    #**** end definition gradient
    #-------------------------------------------------

    # method <- "BFGS"
    arg_control <- list(maxit=mstep_iter)
    if ( mstep_method=="BFGS"){
        arg_control$reltol <- mstep_reltol
    }
    if ( mstep_method=="L-BFGS-B"){
        arg_control$factr <- mstep_reltol
    }
    # method L-BFGS-B uses 'factr' (and 'pgtol') instead of 'reltol' and 'abstol'
    arg_list <- list( par=par0, fn=like_items, method=mstep_method,
                        control=arg_control)
    if (use_grad){
        arg_list$gr <- grad_items
    }
    if ( mstep_method=="L-BFGS-B"){
        arg_list$lower <- partable[ partable$parfree==1, "lower"]
        arg_list$upper <- partable[ partable$parfree==1, "upper"]
    }
    mod <- do.call( stats::optim, arg_list )
    partable <- xxirt_partable_include_freeParameters( partable=partable, x=mod$par )
    res <- list( ll1=mod$value, partable=partable, par0=mod$par  )
    return(res)
}
###############################################################################
