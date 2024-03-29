## File Name: lsem_group_moderator.R
## File Version: 0.172


#***** grouping a moderator variable
lsem_group_moderator <- function( data, type, moderator.grid,
            moderator, residualize, h )
{
    moderator.grouped <- NULL
    if (type=='MGM'){
        G1 <- length(moderator.grid)
        moderator.grouped <- data.frame( min=moderator.grid[-G1],
                                            max=moderator.grid[-1] )
        moderator.grouped$mid <- rowMeans( moderator.grouped)
        v1 <- data[, moderator ]
        v2 <- moderator.grouped$mid[1]
        for (gg in 2:G1){
            v2 <- ifelse( v1 > moderator.grouped$max[gg-1],
                        moderator.grouped$mid[gg], v2 )
        }
        data[,moderator] <- v2
        # residualize <- FALSE
        h <- 1E-5
        moderator.grid <- moderator.grouped$mid
    }
    res <- list( data=data, moderator.grouped=moderator.grouped,
                residualize=residualize, h=h,
                moderator.grid=moderator.grid )
    return(res)
}


lsem.group.moderator <- lsem_group_moderator
