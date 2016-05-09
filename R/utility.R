# utility functions for ppmify

# given an object `x` and character vector of classes `expected`,
expectClasses <- function (x,
                             classes,
                             name,
                             some = any,
                             notify = stop) {

  cl <- class(x)
  if (!some(cl %in% classes)) {
    text <- sprintf('%s should have one of these classes: %s, but had class%s: %s',
                    name,
                    paste(classes, collapse = ', '),
                    ifelse(length(cl) > 1, 'es', ''),
                    paste(cl, collapse = ', '))

    notify(text)
  }
}
