
# Function to print warnings within a function (even before the top-level function returns)
# This is used in run_labDeaths.R
showMeWarnings <- function(expr) 
{
    .list_of_warnings <- c()
    frame_number <- sys.nframe()
    withCallingHandlers(expr, warning = function(w) 
    {
      .list_of_warnings <<- c(.list_of_warnings, w$message)
      invokeRestart("muffleWarning")
    })
    if (length(.list_of_warnings)>0) cat("\nΠΡΟΣΟΧΗ:\n")
    cat(paste(.list_of_warnings, collapse=""))
    .list_of_warnings
}
 
