#' @export
#'
plot.hhsmm <- function (x, ...) 
{
    tmp = x$model$d
    plot(1:nrow(tmp), tmp[, 1], type = "l", ..., ylab = "d(u)", 
        xlab = "u", ylim = range(tmp))
    for (i in 2:x$J) if(x$model$semi[i]) lines(tmp[, i], type = "l", col = i)
    legend("topright", legend = 1:x$J, col = 1:x$J, lty = 1)
}
