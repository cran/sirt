## File Name: mgsem_proc_model_add_specs.R
## File Version: 0.07


mgsem_proc_model_add_specs <- function(model, entry, type, ii, jj, default)
{
    val <- default
    mat <- model[[entry]][[type]]
    if (!is.null(mat)){
        val <- mat[ii,jj]
    }
    return(val)
}
