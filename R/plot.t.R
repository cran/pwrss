##################################
# generic t and z test functions #
##################################

# type = 1 for light red, 2 for light blue
.plot.t.dist <- function(ncp = 0, df = Inf,
                         xlim, ylim = c(0, 0.50),
                         type = 1, ticks = TRUE) {


  plot.window.dim <- dev.size("cm")
  cex.axis <- min(plot.window.dim[1] / 15, plot.window.dim[2] / 15)

  ifelse(type == 1,
         color <- adjustcolor(2, alpha.f = 1),
         color <- adjustcolor(4, alpha.f = 1))

  # non-central t function
  funt <- function(x) {
    dt(x, df = df, ncp = ncp)
  }

  # plot central t distribution
  plot(funt, xlim = xlim, ylim = ylim,
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

}

# type = 1 for light red shade, 2 for light blue shade, 3 for light black stripes
.paint.t.dist <- function(ncp = 0, df = Inf, xlim, type = 1) {

  color <- switch(type,
                   `1` = adjustcolor(2, alpha.f = 0.3),
                   `2` = adjustcolor(4, alpha.f = 0.3),
                   `3` = adjustcolor(1, alpha.f = 0.3))

  # non-central t function
  funt <- function(x) {
    dt(x, df = df, ncp = ncp)
  }

  x <- seq(min(xlim), max(xlim), by = .001)
  y <- funt(x)
  xs <- c(x, rev(x))
  ys <- c(y, rep(0, length(y)))

  if (type == 1 || type == 2) {
    polygon(x = xs, y = ys, col = color, border = NA)
  } else if (type == 3) {
    polygon(x = xs, y = ys, col = color, density = 25, angle = 45, border = NA)
  }

  prob <- pt(max(xlim), df = df, ncp = ncp, lower.tail = TRUE) -
    pt(min(xlim), df = df, ncp = ncp, lower.tail = TRUE)

  return(invisible(prob))

}

.plot.t.t1t2 <- function(ncp, null.ncp = 0, df = Inf, alpha = 0.05,
                         alternative = c("one.sided", "two.sided", "two.one.sided"),
                         plot.main = NULL, plot.sub = NULL) {

  alternative <- tolower(match.arg(alternative))

  # critical t line segment coordinates
  if (alternative == "two.one.sided") {

    ifelse(is.numeric(ncp) && length(ncp) == 1,
           valid.ncp <- TRUE,
           valid.ncp <- FALSE)

    ifelse(is.numeric(null.ncp) && length(null.ncp) %in% c(1, 2),
           valid.null.ncp <- TRUE,
           valid.null.ncp <- FALSE)

    if (isFALSE(valid.ncp)) stop("'ncp' must be numeric and of length one for equivalence tests.", call. = FALSE)
    if (isFALSE(valid.null.ncp)) stop("'null.ncp' must be numeric and of length one (absolute value) or length two (with lower and upper bounds) for the equivalence test.", call. = FALSE)

    if (length(null.ncp) == 1) null.ncp <- c(min(c(-null.ncp, null.ncp)), max(-null.ncp, null.ncp))

    # equivalence test
    if (ncp > min(null.ncp) && ncp < max(null.ncp)) {

      t.alpha.upper <- qt(alpha, df = df, ncp = min(null.ncp), lower.tail = FALSE)
      t.alpha.lower <- qt(alpha, df = df, ncp = max(null.ncp), lower.tail = TRUE)
      t.alpha <- c(t.alpha.lower, t.alpha.upper)

      yt.alpha <- dt(rev(t.alpha), df = df, ncp = null.ncp)

    }

    # minimum effect test
    if (ncp < min(null.ncp) || ncp > max(null.ncp)) {

      t.alpha.lower <- qt(alpha / 2, df = df, ncp = min(null.ncp), lower.tail = TRUE)
      t.alpha.upper <- qt(alpha / 2, df = df, ncp = max(null.ncp), lower.tail = FALSE)
      t.alpha <- c(t.alpha.lower, t.alpha.upper)

      yt.alpha <- dt(t.alpha, df = df, ncp = null.ncp)

    }

  } else if (alternative == "two.sided") {

    ifelse(is.numeric(ncp) && length(ncp) == 1,
           valid.ncp <- TRUE,
           valid.ncp <- FALSE)

    ifelse(is.numeric(null.ncp) && length(null.ncp) == 1,
           valid.null.ncp <- TRUE,
           valid.null.ncp <- FALSE)

    if (isFALSE(valid.ncp) || isFALSE(valid.null.ncp)) stop("'ncp' or 'null.ncp' must be numeric, positive (absolute value), and of length one for the two-sided test.", call. = FALSE)
    # if (ncp < null.ncp) stop("'ncp' must be equal or greater than 'null.ncp' for the two-sided test.", .call = FALSE)

    t.alpha.upper <- qt(alpha / 2, df = df, ncp = null.ncp, lower.tail = FALSE)
    t.alpha.lower <- qt(alpha / 2, df = df, ncp = null.ncp, lower.tail = TRUE)
    t.alpha <- rbind(t.alpha.lower, t.alpha.upper)

    yt.alpha <- dt(t.alpha, df = df, ncp = null.ncp)

  } else {

    ifelse(is.numeric(ncp) && length(ncp) == 1,
           valid.ncp <- TRUE,
           valid.ncp <- FALSE)

    ifelse(is.numeric(null.ncp) && length(null.ncp) == 1,
           valid.null.ncp <- TRUE,
           valid.null.ncp <- FALSE)

    if (isFALSE(valid.ncp) || isFALSE(valid.null.ncp)) stop("'ncp' or 'null.ncp' must be numeric and of length one for the one-sided test.", call. = FALSE)

    ifelse(ncp > null.ncp,
           lower.tail <- FALSE,
           lower.tail <- TRUE)
    t.alpha <- qt(alpha, df = df, ncp = null.ncp, lower.tail = lower.tail) # if ncp > null.ncp

    yt.alpha <- dt(t.alpha, df = df, ncp = null.ncp)

  } # alternative

  # x-axis limits
  ifelse(df < 20, prob.extreme <- 0.001, prob.extreme <- 0.0001)
  lower <- min(min(qt(prob.extreme, df = df, ncp = ncp, lower.tail = TRUE)),
               qt(prob.extreme, df = df, ncp = null.ncp, lower.tail = TRUE))
  upper <- max(max(qt(1 - prob.extreme, df = df, ncp = ncp, lower.tail = TRUE)),
               qt(1 - prob.extreme, df = df, ncp = null.ncp, lower.tail = TRUE))
  xlim <- c(lower, upper)

  plot.window.dim <- dev.size("cm")
  cex.legend <- min(plot.window.dim[1] / 18, plot.window.dim[2] / 15)
  cex.title <- min(plot.window.dim[1] / 11, plot.window.dim[2] / 11)
  cex.label <- min(plot.window.dim[1] / 12, plot.window.dim[2] / 12)

  # plots
  if (alternative == "two.one.sided") {

    .plot.t.dist(ncp = null.ncp[1], df = df, xlim = xlim, type = 1, ticks = TRUE)
    par(new = TRUE)
    .plot.t.dist(ncp = null.ncp[2], df = df, xlim = xlim, type = 1, ticks = FALSE)
    par(new = TRUE)
    .plot.t.dist(ncp = ncp, df = df, xlim = xlim, type = 2, ticks = FALSE)

    text(ncp, dt(ncp, df = df, ncp = ncp) + 0.05,
         labels = expression(H[1]),
         cex = cex.legend, col = adjustcolor(4, alpha.f = 1))

    text(null.ncp[1], dt(null.ncp[1], df = df, ncp = null.ncp[1]) + 0.05,
         labels = expression(H[0]),
         cex = cex.legend, col = adjustcolor(2, alpha.f = 1))

    text(null.ncp[2], dt(null.ncp[2], df = df, ncp = null.ncp[2]) + 0.05,
         labels = expression(H[0]),
         cex = cex.legend, col = adjustcolor(2, alpha.f = 1))

  } else {

    .plot.t.dist(ncp = ncp, df = df, xlim = xlim, type = 2, ticks = TRUE)
    par(new = TRUE)
    .plot.t.dist(ncp = null.ncp, df = df, xlim = xlim, type = 1, ticks = FALSE)

    text(ncp, dt(ncp, df = df, ncp = ncp) + 0.05,
         labels = expression(H[1]),
         cex = cex.legend, col = adjustcolor(4, alpha.f = 1))

    text(null.ncp, dt(null.ncp, df = df, ncp = null.ncp) + 0.05,
         labels = expression(H[0]),
         cex = cex.legend, col = adjustcolor(2, alpha.f = 1))

  } # end of plots

  # draw vertical lines for critical region
  segments(x0 = t.alpha, y0 = rep(0, length(t.alpha)),
           x1 = t.alpha, y1 = yt.alpha, col = 2, lty = 2, lwd = 2)

  # paint regions
  if (alternative == "two.one.sided") {

    # equivalence test
    if (ncp > min(null.ncp) && ncp < max(null.ncp)) {

      .paint.t.dist(ncp = null.ncp[1], df = df, xlim = c(t.alpha[2], max(xlim)), type = 1)
      .paint.t.dist(ncp = null.ncp[2], df = df, xlim = c(min(xlim), t.alpha[1]), type = 1)

      .paint.t.dist(ncp = ncp, df = df, xlim = c(t.alpha[1], max(xlim)), type = 2)
      .paint.t.dist(ncp = ncp, df = df, xlim = c(min(xlim), t.alpha[2]), type = 2)

      ifelse(t.alpha[1] > t.alpha[2],
             power <- .paint.t.dist(ncp = ncp, df = df, xlim = t.alpha, type = 3),
             power <- 0)

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

    power.left <- .paint.t.dist(ncp = ncp, df = df, xlim = c(t.alpha[1], min(xlim)), type = 3)
    power.right <- .paint.t.dist(ncp = ncp, df = df, xlim = c(t.alpha[2], max(xlim)), type = 3)
    power <- power.left + power.right

  } else {

    ifelse(ncp < null.ncp,
           .paint.t.dist(ncp = null.ncp, df = df, xlim = c(t.alpha, min(xlim)), type = 1),
           .paint.t.dist(ncp = null.ncp, df = df, xlim = c(t.alpha, max(xlim)), type = 1))

    ifelse(ncp < null.ncp,
           .paint.t.dist(ncp = ncp, df = df, xlim = c(t.alpha, max(xlim)), type = 2),
           .paint.t.dist(ncp = ncp, df = df, xlim = c(t.alpha, min(xlim)), type = 2))

    ifelse(ncp < null.ncp,
           power <- .paint.t.dist(ncp = ncp, df = df, xlim = c(t.alpha, min(xlim)), type = 3),
           power <- .paint.t.dist(ncp = ncp, df = df, xlim = c(t.alpha, max(xlim)), type = 3))

  } # end of paint regions

  # axes labels and subtitle
  title(main = plot.main, line = 2, cex.main = cex.title)
  title(sub = plot.sub, line = 3, cex.sub = cex.title)
  title(ylab = "Probability Density", line = 2.2, cex.lab = cex.label,
        col.lab = adjustcolor(1, alpha.f = 0.8))
  if (is.finite(df)) {
    title(xlab = paste0("T Value (df = ", round(df, digits = 2), ")"),
          line = 2.2, cex.lab = cex.label,
          col.lab = adjustcolor(1, alpha.f = 0.8))
  } else {
    title(xlab = "Z Value", line = 2.2, cex.lab = cex.label,
          col.lab = adjustcolor(1, alpha.f = 0.8))
  }

  if (power < 0) power <- 0
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

} # end of .plot.t.t1t2()
