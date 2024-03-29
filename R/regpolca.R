## File Name: regpolca.R
## File Version: 0.195


#- Regularized polytomous latent class analysis
regpolca <- function(dat, nclasses, weights=NULL, group=NULL,
    regular_type="scad", regular_lam=0, regular_grouped="none",
    sd_noise_init=1, par_item_init=NULL, par_Theta_init=NULL, random_starts=1,
    random_iter=20, random_sd=0, conv=1e-5, h=1e-4, mstep_iter=10, maxit=1000,
    verbose=TRUE, par_item_max=10, set_equal=.01, eps=1e-3)
{
    #*** preliminaries
    CALL <- match.call()
    s1 <- Sys.time()

    #** analyze response patterns
    res <- regpolca_proc_data(dat=dat, group=group)
    ncats <- res$ncats
    lca_dich <- res$lca_dich
    I <- res$I
    N <- res$N
    group <- res$group
    groups <- res$groups
    G <- res$G
    Ni <- res$Ni

    #* initial parameters
    if (! is.null(par_item_init)){
        # random_starts <- 0
    }

    #- define theta class distribution
    K <- nclasses
    par_Theta <- xxirt_classprobs_lca_init_par(K=K, G=G, random_sd=random_sd,
                    par_Theta_init=par_Theta_init)
    customTheta  <- xxirt_createThetaDistribution( par=par_Theta,
                            est=rep(TRUE,G*(K-1)), P=xxirt_classprobs_lca)
    Theta <- diag(K)

    #- define item response functions
    res <- regpolca_define_customItems( ncats=ncats, K=K, dat=dat,
                    par_item_max=par_item_max )
    customItems <- res$customItems
    partable <- res$partable
    itemtype <- res$itemtype

    #-- include penalty function
    penalty_fun_item <- NULL
    if (length(regular_lam)==1){
        regular_lam <- c(regular_lam, 0)
    }
    if (regular_grouped!='none'){
        regular_lam <- c(regular_lam[1], 0)
    }
    if (sum(regular_lam)>0){
        combis <- t( utils::combn(K,2) )
        combis_classes_list <- list()
        combis_categ_list <- list()
        for (ii in 1L:I){
            ncats1_ii <- ncats[ii]-1
            #- classes
            combis_ii <- combis
            if (ncats1_ii>1){
                for (cc in 2:ncats1_ii){
                    combis_ii <- rbind(combis_ii, K*(cc-1) + combis)
                }
            }
            combis_classes_list[[ii]] <- combis_ii
            #- fusing categories
            combis_ii <- NULL
            if (ncats1_ii>1){
                combis0 <- t( utils::combn(ncats1_ii,2) )
                combis2 <- K*( combis0 - 1 )
                for (cc in 1L:K){
                    combis_ii <- rbind( combis_ii, combis2+cc )
                }
                combis_ii <- data.frame(combis_ii)
                colnames(combis_ii) <- c('par1','par2')
                nri <- nrow(combis_ii)
                n1 <- nri/nrow(combis0)
                combis_ii$cat1 <- rep(combis0[,1], n1)
                combis_ii$cat2 <- rep(combis0[,2], n1)
                combis_ii$class <- rep(1:K, each=nri/K)
                cat_pair <- 100*combis_ii$cat1+combis_ii$cat2
                combis_ii$cat_pair <- match(cat_pair, unique(cat_pair))
                #rep(1:(ncats_ii), each=nri
            }
            combis_categ_list[[ii]] <- combis_ii

        }

        if (regular_type=='scad'){ penalty_used <- penalty_D1_scad }
        if (regular_type=='mcp'){ penalty_used <- penalty_D1_mcp }
        if (regular_type=='lasso'){ penalty_used <- penalty_D1_abs }
        fuse_categories <- rep(FALSE,I)
        if (regular_lam[2]>0){
            for (ii in 1L:I){
                if (ncats[ii]>2){
                    fuse_categories[ii] <- TRUE
                }
            }
        }

        #- define penalty function
        penalty_fun_item <- function(x, ...)
        {
            pen <- 0
            pen <- regpolca_penalty_fun( x=x, regular_grouped=regular_grouped, I=I,
                        partable=partable, combis_classes_list=combis_classes_list,
                        regular_lam=regular_lam, eps=eps, penalty_used=penalty_used,
                        Ni=Ni, combis_categ_list=combis_categ_list,
                        fuse_categories=fuse_categories, K=K)
            return(pen)
        }
    }

    #-- create argument list for xxirt
    args <- list( dat=dat, Theta=Theta, partable=partable, customItems=customItems,
                customTheta=customTheta, maxit=random_iter, mstep_iter=mstep_iter,
                penalty_fun_item=penalty_fun_item, h=h, use_grad=TRUE, verbose=2 )

    #-- random starts if required
    args <- regpolca_run_xxirt_random_starts( args=args, random_starts=random_starts,
                    sd_noise_init=sd_noise_init, par_item_init=par_item_init )

    #-- arguments for final xxirt model
    args$verbose <- TRUE
    args$maxit <- maxit

    #-- run xxirt in a final model
    res <- do.call(what=xxirt, args=args)
    res$iter <- res$iter + random_iter*(random_starts>0)

    #- process output
    res$probs_Theta <- regpolca_postproc_prob_Theta(probs_Theta=res$probs_Theta)
    item <- regpolca_postproc_irf(probs_items=res$probs_items, dat=dat,
                    lca_dich=lca_dich)
    res0 <- regpolca_postproc_count_regularized_parameters(item=item,
                    set_equal=set_equal, lca_dich=lca_dich, probs_items=res$probs_items,
                    nclasses=nclasses, ncats=ncats)
    item1_index <- res0$item1_index
    n_reg <- res0$n_reg
    item <- res0$item
    res$probs_items <- res0$probs_items

    #- adapt information criteria
    res$ic <- regpolca_postproc_ic(ic=res$ic, n_reg=n_reg)
    res$ic$opt_val <- res$opt_val
    res$ic$pen_val <- res$pen_val

    #-- arrange output
    res$CALL <- CALL
    res2 <- list(s1=s1, s2=Sys.time(), lca_dich=lca_dich, nclasses=nclasses,
                    item=item, regular_lam=regular_lam, regular_type=regular_type,
                    item1_index=item1_index, n_reg=n_reg,
                    regular_grouped=regular_grouped)
    res <- sirt_add_list_elements(res=res, res2=res2)

    class(res) <- 'regpolca'
    return(res)
}
