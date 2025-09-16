
############################
# generic F test functions #
############################

# type = 1 for light red, 2 for light blue
.plot.f.dist <- function(ncp = 1, null.ncp = 0, df1, df2, xlim, type = 1, ticks = TRUE) {


  plot.window.dim <- dev.size("cm")
  cex.axis <- min(plot.window.dim[1] / 15, plot.window.dim[2] / 15)

  ifelse(type == 1,
         color <- adjustcolor(2, alpha.f = 1),
         color <- adjustcolor(4, alpha.f = 1))

  # non-central f function
  funf <- function(x) {
    df(x, df1 = df1, df2 = df2, ncp = ncp)
  }

  # mod of f dist
  x.seq <- seq(xlim[1], xlim[2], length.out = 50)
  y.max.alt <- max(df(x.seq, df1 = df1, df2 = df2, ncp = ncp))
  y.max.null <- max(df(x.seq, df1 = df1, df2 = df2, ncp = null.ncp))
  y.max <- max(c(y.max.alt, y.max.null))

  # plot central f distribution
  ylim <- c(0, y.max + 0.05)
  plot(funf, xlim = xlim, ylim = ylim,
       xaxs = "i", yaxs = "i", bty = "l",
       col = color, lwd = 2, lty = type,
       xlab = "", ylab = "",
       xaxt = "n", yaxt = "n",
       cex.axis = cex.axis)

  if (ticks) {

    axis(side = 1,
         #at = seq(0, size, 20),
         #labels = seq(0, size, 20),
         tick = TRUE,
         col.ticks = "gray30", col.axis = "gray30")

    axis(side = 2,
         #at = seq(0, ymax, 0.02),
         #labels = seq(0, ymax, 0.02),
         tick = TRUE,
         col.ticks = "gray30", col.axis = "gray30")

  } # ticks

}  # end of .plot.f.dist()


# type = 1 for light red shade, 2 for light blue shade, 3 for light black stripes
.paint.f.dist <- function(ncp = 0, df1, df2, xlim, type = 1) {

  color <- switch(type,
                   `1` = adjustcolor(2, alpha.f = 0.3),
                   `2` = adjustcolor(4, alpha.f = 0.3),
                   `3` = adjustcolor(1, alpha.f = 0.3))

  # f denisty function
  funf <- function(x) {
    df(x, df1 = df1, df2 = df2, ncp = ncp)
  }

  x <- seq(min(xlim), max(xlim), by = .001)
  y <- funf(x)
  xs <- c(x, rev(x))
  ys <- c(y, rep(0, length(y)))

  if (type == 1 || type == 2) {
    polygon(x = xs, y = ys, col = color, border = NA)
  } else if (type == 3) {
    polygon(x = xs, y = ys, col = color, density = 25, angle = 45, border = NA)
  }

  prob <- pf(max(xlim), df1 = df1, df2 = df2, ncp = ncp, lower.tail = TRUE) -
    pf(min(xlim), df1 = df1, df2 = df2, ncp = ncp, lower.tail = TRUE)

  return(invisible(prob))

} # end of .paint.f.dist()


.plot.f.t1t2 <- function(ncp, null.ncp = 0, df1, df2, alpha,
                         plot.main = NULL, plot.sub = NULL) {

  if (length(ncp) > 1) stop("Not a valid plotting option.", call. = FALSE)

  f.alpha <- qf(alpha, df1 = df1, df2 = df2, ncp = null.ncp, lower.tail = FALSE)
  yf.alpha <- df(f.alpha, df1 = df1, df2 = df2, ncp = null.ncp)

  # x-axis limits
  ifelse(df1 < 2, prob.lower <- 0.30, prob.lower <- 0.001)
  lower <- qf(prob.lower, df1 = df1, df2 = df2, ncp = null.ncp)
  upper <- qf(0.999, df1 = df1, df2 = df2,  ncp = ncp)
  xlim <- c(lower, upper)

  plot.window.dim <- dev.size("cm")
  cex.legend <- min(plot.window.dim[1] / 18, plot.window.dim[2] / 15)
  cex.title <- min(plot.window.dim[1] / 11, plot.window.dim[2] / 11)
  cex.label <- min(plot.window.dim[1] / 12, plot.window.dim[2] / 12)

  # plots
  .plot.f.dist(df1 = df1, df2 = df2, ncp = ncp, null.ncp = null.ncp,
               xlim = xlim, type = 2, ticks = TRUE)
  par(new = TRUE)
  .plot.f.dist(df1 = df1, df2 = df2, ncp = null.ncp, null.ncp = null.ncp,
               xlim = xlim, type = 1, ticks = FALSE)

  # mod of f dist
  x.seq <- seq(xlim[1], xlim[2], length.out = 50)
  y.seq.alt <- df(x.seq, df1 = df1, df2 = df2, ncp = ncp)
  y.seq.null <- df(x.seq, df1 = df1, df2 = df2, ncp = null.ncp)

  y.max.alt <- max(y.seq.alt)
  y.max.null <- max(y.seq.null)
  x.mod.alt <- x.seq[which(y.max.alt == y.seq.alt)]
  x.mod.null <- x.seq[which(y.max.null == y.seq.null)]

  text(x.mod.alt + 0.50, y.max.alt + 0.01,
       labels = expression(H[1]),
       cex = cex.legend, col = adjustcolor(4, alpha.f = 1))

  text(x.mod.null + 0.50, y.max.null + 0.01,
       labels = expression(H[0]),
       cex = cex.legend, col = adjustcolor(2, alpha.f = 1))

  # draw vertical lines for critical region
  segments(x0 = f.alpha, y0 = rep(0, length(f.alpha)),
           x1 = f.alpha, y1 = yf.alpha, col = 2, lty = 2, lwd = 2)

  # paint regions
  .paint.f.dist(df1 = df1, df2 = df2, ncp = null.ncp,
                xlim = c(f.alpha, max(xlim)), type = 1)
  .paint.f.dist(df1 = df1, df2 = df2, ncp = ncp,
                xlim = c(lower, f.alpha), type = 2)
  power <- .paint.f.dist(df1 = df1, df2 = df2, ncp = ncp,
                         xlim = c(f.alpha, max(xlim)), type = 3)
  # end of paint regions

  # axes labels and subtitle
  title(main = plot.main, line = 2, cex.main = cex.title)
  title(sub = plot.sub, line = 3, cex.sub = cex.title)
  title(ylab = "Probability Density", line = 2.2, cex.lab = cex.label,
        col.lab = adjustcolor(1, alpha.f = 0.8))
  title(xlab = paste0("F Value (df1 = ", round(df1, digits = 2), ", df2 = ", round(df2, digits = 2), ")"),
        line = 2.2, cex.lab = cex.label,
        col.lab = adjustcolor(1, alpha.f = 0.8))

  alpha <- round(alpha, 2)
  beta <- round(1 - power, 2)
  power <- round(power, 2)

  legend("topright", cex = cex.legend,
         c(as.expression(bquote(Power == .(power))),
           as.expression(bquote(alpha == .(alpha))),
           as.expression(bquote(beta == .(beta)))),
         fill = c(adjustcolor(1, alpha.f = 0.3),
                  adjustcolor(2, alpha.f = 0.3),
                  adjustcolor(4, alpha.f = 0.3)),
         border = c(adjustcolor(1, alpha.f = 0.15),
                    adjustcolor(2, alpha.f = 0.15),
                    adjustcolor(4, alpha.f = 0.15)),
         bg = adjustcolor(1, alpha.f = 0.08),
         box.col = adjustcolor(1, alpha.f = 0),
         density = c(30, NA, NA),
         angle = c(45, NA, NA))

} # end of .plot.f.t1t2()
