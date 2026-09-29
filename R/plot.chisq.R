#####################################
# generic Chi-square test functions #
#####################################

# type = 1 for light red, 2 for light blue
.plot.chisq.dist <- function(ncp = 0, null.ncp = 0, df, xlim, type = 1, ticks = TRUE) {


  plot.window.dim <- grDevices::dev.size("cm")
  cex.axis <- min(plot.window.dim[1] / 15, plot.window.dim[2] / 15)

  color <- grDevices::adjustcolor(ifelse(type == 1, 2, 4), alpha.f = 1)

  # non-central f function
  funchisq <- function(x) stats::dchisq(x, df = df, ncp = ncp)

  # mod of f dist
  x.seq <- seq(xlim[1], xlim[2], length.out = 50)
  y.max.alt <- max(stats::dchisq(x.seq, df = df, ncp = ncp))
  y.max.null <- max(stats::dchisq(x.seq, df = df, ncp = null.ncp))
  y.max <- max(c(y.max.alt, y.max.null))

  # plot central f distribution
  ylim <- c(0, y.max + y.max * 0.05)
  graphics::plot(funchisq, xlim = xlim, ylim = ylim,
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

} # end of .plot.chisq.dist()


# type = 1 for light red shade, 2 for light blue shade, 3 for light black stripes
.paint.chisq.dist <- function(ncp = 0, df, xlim, type = 1) {

  color <- grDevices::adjustcolor(switch(type, `1` = 2, `2` = 4, `3` = 1), alpha.f = 0.3)

  x <- seq(min(xlim), max(xlim), by = 0.001)
  y <- stats::dchisq(x, df = df, ncp = ncp)
  xs <- c(x, rev(x))
  ys <- c(y, rep(0, length(y)))

  if (type == 1 || type == 2) {
    graphics::polygon(x = xs, y = ys, col = color, border = NA)
  } else if (type == 3) {
    graphics::polygon(x = xs, y = ys, col = color, density = 25, angle = 45, border = NA)
  }

  invisible(stats::pchisq(max(xlim), df = df, ncp = ncp, lower.tail = TRUE) -
            stats::pchisq(min(xlim), df = df, ncp = ncp, lower.tail = TRUE))

} # end of .paint.chisq.dist()



.plot.chisq.t1t2 <- function(ncp, null.ncp = 0, df, alpha = 0.05,
                             plot.main = NULL, plot.sub = NULL) {

  check.numeric(ncp)
  check.proportion(alpha)
  if (!is.numeric(df) || length(df) != 1 || df < 1 || !is.finite(df))
    stop("`df` must be numeric, finite, have a value of at least 1 and have a length of 1.", call. = FALSE)

  chisq.alpha <- stats::qchisq(alpha, df = df, ncp = null.ncp, lower.tail = FALSE)
  y.chisq.alpha <- stats::dchisq(chisq.alpha, df = df, ncp = null.ncp)

  # x-axis limits
  lower <- stats::qchisq(ifelse(df < 2, 0.30, 0.001), df = df, ncp = 0)
  upper <- stats::qchisq(0.999, df = df, ncp = ncp)
  xlim <- c(lower, upper)

  plot.window.dim <- grDevices::dev.size("cm")
  cex.legend <- min(plot.window.dim[1] / 18, plot.window.dim[2] / 15)
  cex.title <- min(plot.window.dim[1] / 11, plot.window.dim[2] / 11)
  cex.label <- min(plot.window.dim[1] / 12, plot.window.dim[2] / 12)

  # plots
  .plot.chisq.dist(ncp = ncp, df = df, xlim = xlim, type = 2, ticks = TRUE)
  graphics::par(new = TRUE)
  .plot.chisq.dist(ncp = null.ncp, df = df, xlim = xlim, type = 1, ticks = FALSE)

  # mod of f dist
  x.seq <- seq(xlim[1], xlim[2], length.out = 50)
  y.seq.alt <- stats::dchisq(x.seq, df = df, ncp = ncp)
  y.seq.null <- stats::dchisq(x.seq, df = df, ncp = null.ncp)

  y.max.alt <- max(y.seq.alt)
  y.max.null <- max(y.seq.null)
  x.mod.alt <- x.seq[which(y.max.alt == y.seq.alt)]
  x.mod.null <- x.seq[which(y.max.null == y.seq.null)]

  graphics::text(x.mod.alt + 0.50, y.max.alt + y.max.alt * 0.06,
                 labels = expression(H[1]),
                 cex = cex.legend, col = grDevices::adjustcolor(4, alpha.f = 1))

  graphics::text(x.mod.null + 0.50, y.max.null + y.max.null * 0.04,
                 labels = expression(H[0]),
                 cex = cex.legend, col = grDevices::adjustcolor(2, alpha.f = 1))

  # draw vertical lines for critical region
  graphics::segments(x0 = chisq.alpha, y0 = rep(0, length(chisq.alpha)),
                     x1 = chisq.alpha, y1 = y.chisq.alpha, col = 2, lty = 2, lwd = 2)

  # paint regions
  .paint.chisq.dist(ncp = null.ncp, df = df,  xlim = c(chisq.alpha, max(xlim)), type = 1)
  .paint.chisq.dist(ncp = ncp, df = df, xlim = c(lower, chisq.alpha), type = 2)
  power <- .paint.chisq.dist(ncp = ncp, df = df, xlim = c(chisq.alpha, max(xlim)), type = 3)
  # end of paint regions

  # axes labels and subtitle
  graphics::title(main = plot.main, line = 2, cex.main = cex.title)
  graphics::title(sub = plot.sub, line = 3, cex.sub = cex.title)
  graphics::title(ylab = "Probability Density", line = 2.2, cex.lab = cex.label,
                  col.lab = grDevices::adjustcolor(1, alpha.f = 0.8))
  graphics::title(xlab = as.expression(bquote(chi[2] ~ "Value (" * df == .(df) * ")")),
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
