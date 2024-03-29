## File Name: lsem_fitsem_joint_estimation_prepare_partable.R
## File Version: 0.318

lsem_fitsem_joint_estimation_prepare_partable <- function(partable, G,
    par_invariant=NULL, par_linear=NULL, par_quadratic=NULL, pw_linear=1,
    pw_quadratic=1)
{
    partable0 <- partable
    partable$id0 <- 1:nrow(partable)
    partable$con <- 0
    label_list <- partable$plabel
    partable1 <- lsem_fitsem_joint_estimation_prepare_partable_include_group_label(
                        partable=partable, gg=1, label_list=label_list)
    partable_mg <- partable1
    for (gg in 2:G){
        partable_gg <- partable
        partable_gg <- lsem_fitsem_joint_estimation_prepare_partable_include_group_label(
                            partable=partable_gg, gg=gg, label_list=label_list)
        partable_gg$group <- partable_gg$block <- gg
        for (vv in c('free','id')){
            partable_gg <- lsem_fitsem_joint_estimation_partable_id(
                                partable_gg=partable_gg, partable_mg=partable_mg, vv=vv)
            partable_gg[,vv][ partable[,vv]==0 ] <- 0
        }
        partable_gg$plabel <- paste0(label_list,'g',gg)
        partable_gg$plabel[ label_list==''] <- ''
        partable_mg <- rbind(partable_mg, partable_gg)
    }

    #- parameter names
    pars <- sirt_lavaan_partable_parnames(partable=partable_mg)
    partable_mg$par <- pars

    # handle constraints
    fixed_invariant <- intersect( paste(partable_mg$par[ partable$free==0]),
                                    par_invariant )
    par_invariant <- setdiff( par_invariant, fixed_invariant)
    par1 <- sirt_define_vector( value='inv', names=par_invariant)
    par2 <- sirt_define_vector( value='lin', names=par_linear)
    par3 <- sirt_define_vector( value='quad', names=par_quadratic)
    par_vec <- c(par1, par2, par3)
    NI <- length(par_vec)
    par_vec_names <- names(par_vec)

    if ( NI > 0 ){
        partable1 <- partable_mg[1,]
        NV <- ncol(partable_mg)
        for (vv in 1:NV){
            if (is.numeric(partable1[1,vv])){
                partable1[1,vv] <- 0
            } else {
                partable1[1,vv] <- ''
            }
        }
        partable1$user <- 2
        partable1$ustart <- NA
        partable1$con <- 0
        partable1$op <- '=='
        partable1c <- partable1

        for (vv in 1:NI){
            par_vec_vv <- par_vec[vv]
            par_vec_names_vv <- par_vec_names[vv]
            ind_vv <- which( paste(partable_mg$par)==par_vec_names[vv] )
            LV2 <- LV <- length(ind_vv)
            N_segments <- 1
            type1 <- 'inv'
            if (par_vec_vv=='lin'){
                LV2 <- LV - 1
                N_segments <- pw_linear
                type1 <- 'lin'
            }
            if (par_vec_vv=='quad'){
                LV2 <- LV - 2
                N_segments <- pw_quadratic
                type1 <- 'quad'
            }
            plabels <- paste(partable_mg$plabel[ind_vv])
            delta <- (LV-1)/N_segments
            par_segments <- data.frame(group=1:LV)
            par_segments$N_segments <- N_segments
            par_segments$type <- type1

            par_segments$segment <- ceiling( ( par_segments$group - 1 ) / delta )
            par_segments$segment <- ifelse( par_segments$segment==0, 1,
                                                par_segments$segment)

            # exclude segments
            par_segments$diff_segment <- c(1,diff(par_segments$segment))
            indices <- par_segments$group[ par_segments$diff_segment==1 ]
            par_segments$con_include <- 1 - par_segments$diff_segment

            if (type1=='lin'){
                indices2 <- setdiff( c(LV,indices-1), seq(-4,1,1) )
                par_segments[indices2, 'con_include'] <- 0
            }
            if (type1=='quad'){
                indices2 <- setdiff( c(LV,indices-1,indices-2), seq(-4,1,1) )
                par_segments[indices2, 'con_include'] <- 0
            }

            for (ll in 2:LV2){
                partable1c$con <- vv
                if (par_vec_vv=='inv'){
                    partable1c$lhs <- plabels[1]
                    partable1c$rhs <- plabels[ll]
                }
                if (par_vec_vv=='lin'){
                    diff1 <- paste0( plabels[ll+1], '-2*', plabels[ll], '+',
                                            plabels[ll-1] )
                    partable1c$lhs <- diff1
                    partable1c$rhs <- 0
                    partable1c$user <- 1
                }
                if (par_vec_vv=='quad'){
                    diff1 <- paste0( plabels[ll+2], '-3*', plabels[ll+1],
                                        '+3*', plabels[ll], '-', plabels[ll-1] )
                    partable1c$lhs <- diff1
                    partable1c$rhs <- 0
                    partable1c$user <- 1
                }
                partable1c$par <- paste0(par_vec_names_vv, '_con', ll-1)
                partable1c$id <- max(partable_mg$id) + 1
                if (par_segments$con_include[ll]==1){
                    partable_mg <- rbind(partable_mg, partable1c)
                }
            }
        }
    }


    #- handle parameter labels
    # same_labels <- setdiff( unique(paste(partable_mg$label)), '')

    #-- output
    return(partable_mg)
}
