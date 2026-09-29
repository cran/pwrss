############################
# generic F test functions #
############################

# type = 1 for light red, 2 for light blue
.plot.f.dist <- function(ncp = 1, null.ncp = 0, df1, df2, xlim, type = 1, ticks = TRUE) {


  plot.window.dim <- grDevices::dev.size("cm")
  cex.axis <- min(plot.window.dim[1] / 15, plot.window.dim[2] / 15)

  color <- grDevices::adjustcolor(ifelse(type == 1, 2, 4), alpha.f = 1)

  # non-central f function
  funf <- function(x) stats::df(x, df1 = df1, df2 = df2, ncp = ncp)

  # mod of f dist
  x.seq <- seq(xlim[1], xlim[2], length.out = 50)
  y.max.alt <- max(stats::df(x.seq, df1 = df1, df2 = df2, ncp = ncp))
  y.max.null <- max(stats::df(x.seq, df1 = df1, df2 = df2, ncp = null.ncp))
  y.max <- max(c(y.max.alt, y.max.null))

  # plot central f distribution
  ylim <- c(0, y.max + 0.05)
  graphics::plot(funf, xlim = xlim, ylim = ylim,
                 xaxs = "i", yaxs = "i", bty = "l",
                 col = color, lwd = 2, lty = type,
                 xlab = "", ylab = "",
                 xaxt = "n", yaxt = "n",
                 cex.axis = cex.axis)

  if (ticks) {

    graphics::axis(side = 1,
                   #at = seq(0, size, 20),
                   #labels = seq(0, size, 20),
                   tick = TRUE,
                   col.ticks = "gray30", col.axis = "gray30")

    graphics::axis(side = 2,
                   #at = seq(0, ymax, 0.02),
                   #labels = seq(0, ymax, 0.02),
                   tick = TRUE,
                   col.ticks = "gray30", col.axis = "gray30")

  } # ticks

}  # end of .plot.f.dist()


# type = 1 for light red shade, 2 for light blue shade, 3 for light black stripes
.paint.f.dist <- function(ncp = 0, df1, df2, xlim, type = 1) {

  color <- grDevices::adjustcolor(switch(type, `1` = 2, `2` = 4, `3` = 1), alpha.f = 0.3)

  x <- seq(min(xlim), max(xlim), by = 0.001)
  y <- stats::df(x, df1 = df1, df2 = df2, ncp = ncp)
  xs <- c(x, rev(x))
  ys <- c(y, rep(0, length(y)))

  if (type == 1 || type == 2) {
    graphics::polygon(x = xs, y = ys, col = color, border = NA)
  } else if (type == 3) {
    graphics::polygon(x = xs, y = ys, col = color, density = 25, angle = 45, border = NA)
  }

  invisible(stats::pf(max(xlim), df1 = df1, df2 = df2, ncp = ncp, lower.tail = TRUE) -
            stats::pf(min(xlim), df1 = df1, df2 = df2, ncp = ncp, lower.tail = TRUE))

} # end of .paint.f.dist()


.plot.f.t1t2 <- function(ncp, null.ncp = 0, df1, df2, alpha = 0.05,
                         plot.main = NULL, plot.sub = NULL) {

  check.numeric(ncp, null.ncp)
  if (!is.numeric(df1) || length(df1) != 1 || df1 < 1 || !is.finite(df1))
    stop("`df1` must be numeric, finite, have a value of at least 1 and have a length of 1.", call. = FALSE)
  if (!is.numeric(df2) || length(df2) != 1 || df2 < 1 || !is.finite(df2))
    stop("`df2` must be numeric, finite, have a value of at least 1 and have a length of 1.", call. = FALSE)
  check.proportion(alpha)

  f.alpha <- stats::qf(alpha, df1 = df1, df2 = df2, ncp = null.ncp, lower.tail = FALSE)
  yf.alpha <- stats::df(f.alpha, df1 = df1, df2 = df2, ncp = null.ncp)

  # x-axis limits
  lower <- stats::qf(ifelse(df1 < 2, 0.30, 0.001), df1 = df1, df2 = df2, ncp = null.ncp)
  upper <- stats::qf(0.999, df1 = df1, df2 = df2,  ncp = ncp)
  xlim <- c(lower, upper)

  plot.window.dim <- grDevices::dev.size("cm")
  cex.legend <- min(plot.window.dim[1] / 18, plot.window.dim[2] / 15)
  cex.title <- min(plot.window.dim[1] / 11, plot.window.dim[2] / 11)
  cex.label <- min(plot.window.dim[1] / 12, plot.window.dim[2] / 12)

  # plots
  .plot.f.dist(df1 = df1, df2 = df2, ncp = ncp, null.ncp = null.ncp,
               xlim = xlim, type = 2, ticks = TRUE)
  graphics::par(new = TRUE)
  .plot.f.dist(df1 = df1, df2 = df2, ncp = null.ncp, null.ncp = null.ncp,
               xlim = xlim, type = 1, ticks = FALSE)

  # mod of f dist
  x.seq <- seq(xlim[1], xlim[2], length.out = 50)
  y.seq.alt <- stats::df(x.seq, df1 = df1, df2 = df2, ncp = ncp)
  y.seq.null <- stats::df(x.seq, df1 = df1, df2 = df2, ncp = null.ncp)

  y.max.alt <- max(y.seq.alt)
  y.max.null <- max(y.seq.null)
  x.mod.alt <- x.seq[which(y.max.alt == y.seq.alt)]
  x.mod.null <- x.seq[which(y.max.null == y.seq.null)]

  graphics::text(x.mod.alt + 0.50, y.max.alt + 0.01,
                 labels = expression(H[1]),
                 cex = cex.legend, col = grDevices::adjustcolor(4, alpha.f = 1))

  graphics::text(x.mod.null + 0.50, y.max.null + 0.01,
                 labels = expression(H[0]),
                 cex = cex.legend, col = grDevices::adjustcolor(2, alpha.f = 1))

  # draw vertical lines for critical region
  graphics::segments(x0 = f.alpha, y0 = rep(0, length(f.alpha)),
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
  graphics::title(main = plot.main, line = 2, cex.main = cex.title)
  graphics::title(sub = plot.sub, line = 3, cex.sub = cex.title)
  graphics::title(ylab = "Probability Density", line = 2.2, cex.lab = cex.label,
                  col.lab = grDevices::adjustcolor(1, alpha.f = 0.8))
  graphics::title(xlab = paste0("F Value (df1 = ", round(df1, digits = 2), ", df2 = ", round(df2, digits = 2), ")"),
                  line = 2.2, cex.lab = cex.label,
                  col.lab = grDevices::adjustcolor(1, alpha.f = 0.8))

  alpha <- round(alpha, 2)
  beta <- round(1 - power, 2)
  power <- round(power, 2)

  graphics::legend("topright", cex = cex.legend,
                   c(as.expression(bquote(Power == .(power))),
                     as.expression(bquote(alpha == .(alpha))),
                     as.expression(bquote(beta == .(beta)))),
                   fill = c(grDevices::adjustcolor(1, alpha.f = 0.3),
                            grDevices::adjustcolor(2, alpha.f = 0.3),
                            grDevices::adjustcolor(4, alpha.f = 0.3)),
                   border = c(grDevices::adjustcolor(1, alpha.f = 0.15),
                              grDevices::adjustcolor(2, alpha.f = 0.15),
                              grDevices::adjustcolor(4, alpha.f = 0.15)),
                   bg = grDevices::adjustcolor(1, alpha.f = 0.08),
                   box.col = grDevices::adjustcolor(1, alpha.f = 0),
                   density = c(30, NA, NA),
                   angle = c(45, NA, NA))

} # end of .plot.f.t1t2()
