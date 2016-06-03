#' @title Linear Equating Analytic Standard Error
#'
#' @description This package is used to calcuate the analytic standard error
#' of linear equating (the delta method) under random groups design.
#'
#' @details According to Kolen and Brennan (2014, p. 250), "the analytic method
#' (referred to as the delta method) that can be used to estimate the standard
#' errors of equating using sample statistics."
#'
#' @param x,y frequency tables of ramdom group equating design
#'
#' @examples
#' #Data prep
#' library("equate")
#' data("ACTmath")
#'
#' #Frequency tables of random group equating
#' x <- as.freqtab(ACTmath[, 1:2])
#' y <- as.freqtab(ACTmath[, c(1, 3)])
#'
#' #Getting std. error of linear equating under random groups design
#' linear_rg(x,y)
#'
#' @export
linear_rg <- function (x, y, ...) {
    if (missing(y))
    y <- margin.table(x, 2)
    if (margins(y) < margins(x))
    x <- margin.table(x, 1)
    xscale <- scales(x)

  # Use stats to get linear equating(revised from equate::linear)
      sigmax <- sd.freqtab(x)
      sigmay <- sd.freqtab(y)
      mux <- mean(x)
      muy <- mean(y)
      slope <- sigmay/sigmax
      intercept <- muy - slope*mux

    yx <- intercept + slope*xscale

  # Get se
    sigmaysq <-var.freqtab(y)
    sux <- summary(x)
    nx <- sux$n
    suy <- summary(y)
    ny <- suy$n
    skx <- skew.freqtab(x)
    sky <- skew.freqtab(y)
    kux <- kurt.freqtab(x, margin = 1)
    kuy <- kurt.freqtab(y, margin = 1)

    sk<- ((skx/nx) + (sky/ny))
    zs <- (xscale - mux)/sigmax
    ku <- ((kux - 1)/(4*nx)) + ((kuy - 1)/(4*ny))

    selrg <- function(nx, ny, sk, ku, zs) {

    sigmaysq * ((1/nx) + (1/ny) + (sk * zs) + (ku * zs^2))

    }
    se <- selrg(nx, ny, sk, ku, zs)

   # if x and y are assumed normally distributed
    selrg1 <- function(nx, ny, zs) {

    (sigmaysq/2) * ((1/nx) + (1/ny)) * (2 + (zs^2))

    }

    se1 <- selrg1(nx, ny, zs)

    # if sample size for x and y are assumed to be equal
    selrg2 <- function(nx, ny, nt, zs) {

    nt <- nx + ny
    ((2*sigmaysq)/(nt)) * (2 + (zs^2))

    }

    se2 <- selrg2(nx, ny, nt, zs)

  out <- cbind(yx, se, se1, se2)
  #out <- list(yx = yx, se = se, se1 = se1, se2 = se2)
  return(out)

}

#' #Checking results
#' round(se, 2)
#' round(se1, 2)
#' round(se2, 2)
#' round(mean(se), 2)
#' round(mean(se1), 2)
#' round(mean(se2), 2)
#'
#' plot(xscale, yx)
#' plot(xscale, se)
#' plot(xscale, se1)
#' plot(xscale, se2)
