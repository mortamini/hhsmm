.progress <- function (x, max = 100) {
    percent <- x / max * 100
    cat(sprintf('\r[%-10s] %d%%',
                paste(rep('.', percent / 10), collapse = ''),
                floor(percent)))
    if (x == max)
        cat('\n')
}
