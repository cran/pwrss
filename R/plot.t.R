##################################
# generic t and z test functions #
##################################

# type = 1 for light red, 2 for light blue
.plot.t.dist <- function(ncp = 0, df = Inf,
                         xlim, ylim = c(0, 0.50),
                         type = 1, ticks = TRUE) {


  plot.window.dim <- grDevices::dev.size("cm")
  cex.axis <- min(plot.window.dim[1] / 15, plot.window.dim[2] / 15)

  color <- grDevices::adjustcolor(ifelse(type == 1, 2, 4), alpha.f = 1)

  # non-central t function
  funt <- function(x) stats::dt(x, df = df, ncp = ncp)

  # plot central t distribution
  graphics::plot(funt, xlim = xlim, ylim = ylim,
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

}

# type = 1 for light red shade, 2 for light blue shade, 3 for light black stripes
.paint.t.dist <- function(ncp = 0, df = Inf, xlim, type = 1) {

  color <- grDevices::adjustcolor(switch(type, `1` = 2, `2` = 4, `3` = 1), alpha.f = ifelse(type == "s", 0.6, 0.3))

  x <- seq(min(xlim), max(xlim), by = 0.001)
  y <- stats::dt(x, df = df, ncp = ncp)
  xs <- c(x, rev(x))
  ys <- c(y, rep(0, length(y)))

  if (type == 1 || type == 2) {
    graphics::polygon(x = xs, y = ys, col = color, border = NA)
  } else if (type == 3) {
    graphics::polygon(x = xs, y = ys, col = color, density = 25, angle = 45, border = NA)
  }

  invisible(stats::pt(max(xlim), df = df, ncp = ncp, lower.tail = TRUE) -
            stats::pt(min(xlim), df = df, ncp = ncp, lower.tail = TRUE))

}

.plot.t.t1t2 <- function(ncp, null.ncp = 0, df = Inf, alpha = 0.05,
                         alternative = c("one.sided", "two.sided", "two.one.sided"),
                         plot.main = NULL, plot.sub = NULL, digits = 2) {

  alternative <- match.arg(alternative)

  check.numeric(ncp)
  null.ncp <- check.margins(null.ncp, check.numeric, alternative)
  if (!is.numeric(df) || length(df) != 1 || df < 1)
    stop("`df` must be numeric, have a value of at least 1 and have a length of 1.", call. = FALSE)
  check.proportion(alpha)

  # critical t line segment coordinates
  if (alternative == "two.one.sided") {

    # equivalence test
    if (ncp > min(null.ncp) && ncp < max(null.ncp)) {

      t.alpha.upper <- stats::qt(alpha, df = df, ncp = min(null.ncp), lower.tail = FALSE)
      t.alpha.lower <- stats::qt(alpha, df = df, ncp = max(null.ncp), lower.tail = TRUE)
      t.alpha <- c(t.alpha.lower, t.alpha.upper)

      yt.alpha <- stats::dt(rev(t.alpha), df = df, ncp = null.ncp)

    }

    # minimum effect test
    if (ncp < min(null.ncp) || ncp > max(null.ncp)) {

      t.alpha.lower <- stats::qt(alpha / 2, df = df, ncp = min(null.ncp), lower.tail = TRUE)
      t.alpha.upper <- stats::qt(alpha / 2, df = df, ncp = max(null.ncp), lower.tail = FALSE)
      t.alpha <- c(t.alpha.lower, t.alpha.upper)

      yt.alpha <- stats::dt(t.alpha, df = df, ncp = null.ncp)

    }

  } else if (alternative == "two.sided") {

    t.alpha.upper <- stats::qt(alpha / 2, df = df, ncp = null.ncp, lower.tail = FALSE)
    t.alpha.lower <- stats::qt(alpha / 2, df = df, ncp = null.ncp, lower.tail = TRUE)
    t.alpha <- c(t.alpha.lower, t.alpha.upper)

    yt.alpha <- stats::dt(t.alpha, df = df, ncp = null.ncp)

  } else if (alternative == "one.sided") {

    t.alpha <- stats::qt(alpha, df = df, ncp = null.ncp, lower.tail = ncp <= null.ncp)

    yt.alpha <- stats::dt(t.alpha, df = df, ncp = null.ncp)

  } # alternative

  # x-axis limits
  prob.extreme <- ifelse(df < 20, 0.001, 0.0001)
  lower <- min(min(stats::qt(prob.extreme,     df = df, ncp = ncp,      lower.tail = TRUE)),
                   stats::qt(prob.extreme,     df = df, ncp = null.ncp, lower.tail = TRUE))
  upper <- max(max(stats::qt(1 - prob.extreme, df = df, ncp = ncp,      lower.tail = TRUE)),
                   stats::qt(1 - prob.extreme, df = df, ncp = null.ncp, lower.tail = TRUE))
  xlim <- c(lower, upper)

  plot.window.dim <- grDevices::dev.size("cm")
  cex.legend <- min(plot.window.dim[1] / 18, plot.window.dim[2] / 15)
  cex.title <- min(plot.window.dim[1] / 11, plot.window.dim[2] / 11)
  cex.label <- min(plot.window.dim[1] / 12, plot.window.dim[2] / 12)

  # plots
  if (alternative == "two.one.sided") {

    .plot.t.dist(ncp = null.ncp[1], df = df, xlim = xlim, type = 1, ticks = TRUE)
    graphics::par(new = TRUE)
    .plot.t.dist(ncp = null.ncp[2], df = df, xlim = xlim, type = 1, ticks = FALSE)
    graphics::par(new = TRUE)
    .plot.t.dist(ncp = ncp, df = df, xlim = xlim, type = 2, ticks = FALSE)

    graphics::text(ncp, stats::dt(ncp, df = df, ncp = ncp) + 0.05,
                   labels = expression(H[1]),
                   cex = cex.legend, col = grDevices::adjustcolor(4, alpha.f = 1))

    graphics::text(null.ncp[1], stats::dt(null.ncp[1], df = df, ncp = null.ncp[1]) + 0.05,
                   labels = expression(H[0]),
                   cex = cex.legend, col = grDevices::adjustcolor(2, alpha.f = 1))

    graphics::text(null.ncp[2], stats::dt(null.ncp[2], df = df, ncp = null.ncp[2]) + 0.05,
                   labels = expression(H[0]),
                   cex = cex.legend, col = grDevices::adjustcolor(2, alpha.f = 1))

  } else {

    .plot.t.dist(ncp = ncp, df = df, xlim = xlim, type = 2, ticks = TRUE)
    graphics::par(new = TRUE)
    .plot.t.dist(ncp = null.ncp, df = df, xlim = xlim, type = 1, ticks = FALSE)

    graphics::text(ncp, stats::dt(ncp, df = df, ncp = ncp) + 0.05,
                   labels = expression(H[1]),
                   cex = cex.legend, col = grDevices::adjustcolor(4, alpha.f = 1))

    graphics::text(null.ncp, stats::dt(null.ncp, df = df, ncp = null.ncp) + 0.05,
                   labels = expression(H[0]),
                   cex = cex.legend, col = grDevices::adjustcolor(2, alpha.f = 1))

  } # end of plots

  # draw vertical lines for critical region
  graphics::segments(x0 = t.alpha, y0 = rep(0, length(t.alpha)),
                     x1 = t.alpha, y1 = yt.alpha, col = 2, lty = 2, lwd = 2)

  # paint regions
  if (alternative == "two.one.sided") {

    # equivalence test
    if (ncp > min(null.ncp) && ncp < max(null.ncp)) {

      .paint.t.dist(ncp = null.ncp[1], df = df, xlim = c(t.alpha[2], max(xlim)), type = 1)
      .paint.t.dist(ncp = null.ncp[2], df = df, xlim = c(min(xlim), t.alpha[1]), type = 1)

      .paint.t.dist(ncp = ncp, df = df, xlim = c(t.alpha[1], max(xlim)), type = 2)
      .paint.t.dist(ncp = ncp, df = df, xlim = c(min(xlim), t.alpha[2]), type = 2)

      power <- ifelse(t.alpha[1] > t.alpha[2],
                      .paint.t.dist(ncp = ncp, df = df, xlim = t.alpha, type = 3),
                      0)

    }

    # minimum effect test
    if (ncp < min(null.ncp) || ncp > max(null.ncp)) {

      .paint.t.dist(ncp = null.ncp[1], df = df, xlim = c(t.alpha[1], min(xlim)), type = 1)
      .paint.t.dist(ncp = null.ncp[2], df = df, xlim = c(max(xlim), t.alpha[2]), type = 1)

      if (ncp > max(null.ncp)) {
        .paint.t.dist(ncp = ncp, df = df, xlim = c(t.alpha[2], min(xlim)), type = 2)
      } else {
        .paint.t.dist(ncp = ncp, df = df, xlim = c(t.alpha[1], max(xlim)), type = 2)
      }

      power.upper <- .paint.t.dist(ncp = ncp, df = df, xlim = c(t.alpha[1], min(xlim)), type = 3)
      power.lower <- .paint.t.dist(ncp = ncp, df = df, xlim = c(t.alpha[2], max(xlim)), type = 3)
      power <- power.upper + power.lower

    }

  } else if (alternative == "two.sided") {

    .paint.t.dist(ncp = null.ncp, df = df, xlim = c(t.alpha[1], min(xlim)), type = 1)
    .paint.t.dist(ncp = null.ncp, df = df, xlim = c(t.alpha[2], max(xlim)), type = 1)

    .paint.t.dist(ncp = ncp, df = df, xlim = c(t.alpha[1], t.alpha[2]), type = 2)

    # type S region
    if (ncp > null.ncp) {
      .paint.t.dist(ncp = ncp, df = df, xlim = c(min(t.alpha), min(xlim)), type = 1)
    } else {
      .paint.t.dist(ncp = ncp, df = df, xlim = c(max(t.alpha), max(xlim)), type = 1)
    }

    # type S
    Phi.p <- stats::pt(q = max(t.alpha), df = df, ncp = ncp)
    Phi.m <- stats::pt(q = min(t.alpha), df = df, ncp = ncp)
    type.s <- min(Phi.m, 1 - Phi.p) / (Phi.m + 1 - Phi.p)
    type.s <- round(type.s, digits)

    # type M
    type.m <- suppressWarnings({
      bounds <- stats::qt(c(1e-10, 1 - 1e-10), df = df, ncp = ncp)
      integrand <- function(t) abs(t) * stats::dt(t, df = df, ncp = ncp)
      numerator <- stats::integrate(integrand, min(bounds), min(t.alpha))$value +
                   stats::integrate(integrand, max(t.alpha), max(bounds))$value
      denominator  <- abs(ncp) * (stats::pt(min(t.alpha), df = df, ncp = ncp, lower.tail = TRUE) +
                                  stats::pt(max(t.alpha), df = df, ncp = ncp, lower.tail = FALSE))
      numerator / denominator
    })
    type.m <- round(type.m, digits)

    power.left <- .paint.t.dist(ncp = ncp, df = df, xlim = c(t.alpha[1], min(xlim)), type = 3)
    power.right <- .paint.t.dist(ncp = ncp, df = df, xlim = c(t.alpha[2], max(xlim)), type = 3)
    power <- power.left + power.right

  } else {

    if (ncp < null.ncp) {

               .paint.t.dist(ncp = null.ncp, df = df, xlim = c(t.alpha, min(xlim)), type = 1)
               .paint.t.dist(ncp = ncp,      df = df, xlim = c(t.alpha, max(xlim)), type = 2)
      power <- .paint.t.dist(ncp = ncp,      df = df, xlim = c(t.alpha, min(xlim)), type = 3)

    } else {

               .paint.t.dist(ncp = null.ncp, df = df, xlim = c(t.alpha, max(xlim)), type = 1)
               .paint.t.dist(ncp = ncp,      df = df, xlim = c(t.alpha, min(xlim)), type = 2)
      power <- .paint.t.dist(ncp = ncp,      df = df, xlim = c(t.alpha, max(xlim)), type = 3)

    }

  } # end of paint regions

  # axes labels and subtitle
  graphics::title(main = plot.main, line = 2, cex.main = cex.title)
  graphics::title(sub = plot.sub, line = 3, cex.sub = cex.title)
  graphics::title(ylab = "Probability Density", line = 2.2, cex.lab = cex.label,
                  col.lab = grDevices::adjustcolor(1, alpha.f = 0.8))
  if (is.finite(df)) {
    graphics::title(xlab = paste0("T Value (df = ", round(df, digits = 2), ")"),
                    line = 2.2, cex.lab = cex.label,
                    col.lab = grDevices::adjustcolor(1, alpha.f = 0.8))
  } else {
    graphics::title(xlab = "Z Value", line = 2.2, cex.lab = cex.label,
                    col.lab = grDevices::adjustcolor(1, alpha.f = 0.8))
  }

  power[power < 0] <- 0
  alpha <- round(alpha, digits)
  beta <- round(1 - power, digits)
  power <- round(power, digits)

  if (alternative == "two.sided") {
    graphics::legend("topright", cex = cex.legend,
                     c(as.expression(bquote(Power == .(power))),
                       as.expression(bquote(alpha == .(alpha))),
                       as.expression(bquote(beta == .(beta))),
                       as.expression(bquote(S == .(type.s))),
                       as.expression(bquote(M == .(type.m)))),
                     fill = c(grDevices::adjustcolor(1, alpha.f = 0.3),
                              grDevices::adjustcolor(2, alpha.f = 0.3),
                              grDevices::adjustcolor(4, alpha.f = 0.3),
                              grDevices::adjustcolor(2, alpha.f = 0.6),
                              grDevices::adjustcolor(1, alpha.f = 0)),
                     border = c(grDevices::adjustcolor(1, alpha.f = 0.15),
                                grDevices::adjustcolor(2, alpha.f = 0.15),
                                grDevices::adjustcolor(4, alpha.f = 0.15),
                                grDevices::adjustcolor(2, alpha.f = 0.30),
                                grDevices::adjustcolor(1, alpha.f = 0)),
                     bg = grDevices::adjustcolor(1, alpha.f = 0.08),
                     box.col = grDevices::adjustcolor(1, alpha.f = 0),
                     density = c(30, NA, NA, NA, NA),
                     angle = c(45, NA, NA, NA, NA))
  } else {
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
  }

} # end of .plot.t.t1t2()
