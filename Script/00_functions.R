#### sloss function developed followint Lexiguel (https://rdrr.io/github/kamapu/Lexiguel/man/Lexiguel-package.html)
sloss <- function(table, env=data.frame(), area) {
  if(!is.matrix(table)) table <- as.matrix(table)
  area <- substitute(area)
  area <- eval(area, env, parent.frame())
  SLOSS <- list(SL=list(), LS=list())
  # First the calculation from small to large
  SLOSS$SL$area <- c(0, cumsum(area[order(area)]))
  Flor <- apply(table[order(area),], 2, cumsum)
  Flor[Flor > 0] <- 1
  SLOSS$SL$species <- c(0, apply(Flor, 1, sum))
  # Now the calculation from large to small
  SLOSS$LS$area <- c(0, cumsum(area[order(area, decreasing=TRUE)]))
  Flor <- apply(table[order(area, decreasing=TRUE),], 2, cumsum)
  Flor[Flor > 0] <- 1
  SLOSS$LS$species <- c(0, apply(Flor, 1, sum))
  # Calculation of SLOSS index
  SLOSS$Index <- with(SLOSS$SL, curve_area(area,
                                           species))/with(SLOSS$LS, curve_area(area, species))
  # Final object
  class(SLOSS) <- c("SLOSS","list")
  return(SLOSS)
}

curve_area <- function(x, y, bottom=0) {
  D1 <- c(diff(x))
  D2 <- c(diff(y))
  Area <- sum(D1*((y - bottom)[-length(y)]) + D1*D2/2, na.rm=TRUE)
  return(Area)
}

# edited from original function in Lexiguel to include points
plot.SLOSS <- function(x, y=NULL, sl.lty=2, sl.lwd=1, sl.col="black", ls.lty=1,
                       ls.lwd=1, ls.col="black", show.index=TRUE, digits.index=2, cex.index=1,
                       pos.index=c(0.05,0.95), show.legend=FALSE, pos.legend="bottomright",
                       bty.legend="o", main="SLOSS curves",...) {
  with(x$SL, plot(area, species, type="l", lty=sl.lty, lwd=sl.lwd, col=sl.col,
                  main=main, ...))
  with(x$LS, lines(area, species, lty=ls.lty, lwd=ls.lwd, col=ls.col))
  
  with(x$SL, points(area, species, pch = 18, cex = 1))
  with(x$LS, points(area, species, pch = 18, cex = 1))
  
  if(show.legend) {
    legend(pos.legend, lty=c(sl.lty,ls.lty), lwd=c(sl.lwd,ls.lwd),
           legend=c("small to large","large to small"), bty=bty.legend)
  }
  if(show.index) {
    with(x$SL, text(max(area)*pos.index[1], max(species)*pos.index[2],
                    labels=paste("SLOSS-index =",
                                 round(x$Index, digits.index)),
                    cex=cex.index, pos=4))
  }
}


