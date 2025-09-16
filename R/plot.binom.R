# type = 1 for light red, 2 for light blue
.plot.binom.dist <- function(size, prob = 0.50,
                             xlim = c(0, size),
                             ylim = c(0, dbinom(x = size * prob, size = size, prob = prob) * 1.10),
                             type = 1, ticks = FALSE) {


  plot.window.dim <- dev.size("cm")
  cex.axis <- min(plot.window.dim[1] / 15, plot.window.dim[2] / 15)

  ifelse(type == 1,
         color <- adjustcolor(2, alpha.f = 1),
         color <- adjustcolor(4, alpha.f = 1))

  # binom density function
  fun.binom <- function(x) {
    dbinom(x, size = size, prob = prob)
  }

  xseq <- seq(xlim[1], xlim[2])
  yseq <- fun.binom(xseq)

  plot(xseq - 0.50, yseq, type = "s", # x-values are nudged by -0.50 unit to match the bars later
       xlim = xlim, ylim = ylim,
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

} # .plot.binom.dist()


# type = 1 for light red shade, 2 for light blue shade, 3 for light black stripes
.paint.binom.dist <- function(size, prob = 0.50, xlim = c(0, size), type = 1) {

  color <- switch(type,
                   `1` = adjustcolor(2, alpha.f = 0.3),
                   `2` = adjustcolor(4, alpha.f = 0.3),
                   `3` = adjustcolor(1, alpha.f = 0.3))

  # binom density function
  fun.binom <- function(x) {
    dbinom(x, size = size, prob = prob)
  }

  x <- seq(min(xlim), max(xlim), by = 1)
  y <- fun.binom(x)
  xs <- c(x, rev(x))
  ys <- c(y, rep(0, length(y)))


  if (type == 1 || type == 2) {

    # polygon(x = xs, y = ys, col = color, border = NA)

    bar.width <- 0.8
    for (i in seq_along(xs)) {
      x.left <- x[i] - bar.width / 2
      x.right <- x[i] + bar.width / 2
      polygon(x = c(x.left, x.left, x.right, x.right),
              y = c(0, ys[i], ys[i], 0),
              col = color, border = NA)
    }

  } else if (type == 3) {

    # polygon(x = xs, y = ys, col = color, density = 20, angle = 45, border = NA)

    bar.width <- 0.8
    for (i in seq_along(xs)) {
      ifelse(x[i] %% 2 == 0, angle <- 45, angle <- 135)
      x.left <- x[i] - bar.width / 2
      x.right <- x[i] + bar.width / 2
      polygon(x = c(x.left, x.left, x.right, x.right),
              y = c(0, ys[i], ys[i], 0),
              density = 30, angle = angle,
              col = color, border = NA)
    }

  }

  if (min(xlim) > 0) {
    prob.shaded <- pbinom(max(xlim), size = size, prob = prob, lower.tail = TRUE) -
      pbinom(min(xlim), size = size, prob = prob, lower.tail = TRUE)
  } else {
    prob.shaded <- pbinom(max(xlim), size = size, prob = prob, lower.tail = TRUE)
  }

  return(invisible(prob.shaded))

}



.plot.binom.t1t2 <- function(size, prob, null.prob = 0.50, alpha, alternative,
                             plot.main = NULL, plot.sub = NULL) {

  if (size < 10) stop("Number of trials should be greater than 10 for plotting", call. = FALSE)

  # critical t line segment coordinates
  if (alternative == "two.one.sided") {

    if (length(null.prob) == 1) null.prob <- rbind(prob - abs(null.prob - prob), prob + abs(null.prob - prob))
    if (length(null.prob) == 2 && is.vector(null.prob)) null.prob <- rbind(min(null.prob), max(null.prob))
    if (length(null.prob) > 2) stop("Not a valid plotting option.", call. = FALSE)

    if (prob > min(null.prob) && prob < max(null.prob)) {
      # equivalence
      q.binom.upper <- qbinom(alpha, size = size, prob = null.prob[1], lower.tail = FALSE)
      p.binom.upper <- pbinom(q.binom.upper, size = size, prob = null.prob[1], lower.tail = FALSE)
      if (p.binom.upper > alpha) q.binom.upper <- q.binom.upper + 1

      q.binom.lower <- qbinom(alpha, size = size, prob = null.prob[2], lower.tail = TRUE)
      p.binom.lower <- pbinom(q.binom.lower, size = size, prob = null.prob[2], lower.tail = TRUE)
      if (p.binom.lower > alpha) q.binom.lower <- q.binom.lower - 1

      q.binom.alpha <- c(q.binom.upper, q.binom.lower)

      y.binom.alpha <- dbinom(q.binom.alpha, size = size, prob = null.prob)

    } else {
      # minimal effect

      q.binom.lower <- qbinom(alpha / 2, size = size, prob = null.prob[1], lower.tail = TRUE)
      p.binom.lower <- pbinom(q.binom.lower, size = size, prob = null.prob[1], lower.tail = TRUE)
      if (p.binom.lower > alpha / 2) q.binom.lower <- q.binom.lower - 1

      q.binom.upper <- qbinom(alpha / 2, size = size, prob = null.prob[2], lower.tail = FALSE)
      p.binom.upper <- pbinom(q.binom.upper, size = size, prob = null.prob[2], lower.tail = FALSE)
      if (p.binom.upper > alpha / 2) q.binom.upper <- q.binom.upper + 1

      q.binom.alpha <- c(q.binom.lower, q.binom.upper)

      y.binom.alpha <- dbinom(q.binom.alpha, size = size, prob = null.prob)

    }

  } else if (alternative == "two.sided") {

    if (length(prob) > 1) stop("Not a valid plotting option.", call. = FALSE)

    q.binom.lower <- qbinom(alpha / 2, size = size, prob = null.prob, lower.tail = TRUE)
    p.binom.lower <- pbinom(q.binom.lower, size = size, prob = null.prob, lower.tail = TRUE)
    if (p.binom.lower > alpha / 2) q.binom.lower <- q.binom.lower - 1

    q.binom.upper <- qbinom(alpha / 2, size = size, prob = null.prob, lower.tail = FALSE)
    p.binom.upper <- pbinom(q.binom.upper, size = size, prob = null.prob, lower.tail = FALSE)
    if (p.binom.upper > alpha / 2) q.binom.upper <- q.binom.upper + 1

    q.binom.alpha <- c(q.binom.lower, q.binom.upper)

    y.binom.alpha <- dbinom(q.binom.alpha, size = size, prob = null.prob)

  } else if (alternative == "one.sided") {

    if (length(prob) > 1) stop("Not a valid plotting option.", call. = FALSE)

    if (prob < null.prob) {
      # less
      q.binom.alpha <- qbinom(alpha, size = size, prob = null.prob, lower.tail = TRUE)
      p.binom.alpha <- pbinom(q.binom.alpha, size = size, prob = null.prob, lower.tail = TRUE)
      if (p.binom.alpha > alpha) q.binom.alpha <- q.binom.alpha - 1

      y.binom.alpha <- dbinom(q.binom.alpha, size = size, prob = null.prob)

    } else {
      # greater
      q.binom.alpha <- qbinom(alpha, size = size, prob = null.prob, lower.tail = FALSE)
      p.binom.alpha <- pbinom(q.binom.alpha, size = size, prob = null.prob, lower.tail = FALSE)
      if (p.binom.alpha > alpha) q.binom.alpha <- q.binom.alpha + 1

      y.binom.alpha <- dbinom(q.binom.alpha, size = size, prob = null.prob)

    }

  } # alternative

  # x-axis limits
  prob.extreme <- 0.00001
  lower <- min(min(qbinom(prob.extreme, size = size, prob = prob, lower.tail = TRUE)),
               qbinom(prob.extreme, size = size, prob = null.prob, lower.tail = TRUE))
  upper <- max(max(qbinom(1 - prob.extreme, size = size, prob = prob, lower.tail = TRUE)),
               qbinom(1 - prob.extreme, size = size, prob = null.prob, lower.tail = TRUE))
  xlim <- c(max(0, lower - 1), upper + 1)


  # y-axis limits
  ymax.HA <- dbinom(round(prob * size), size = size, prob = prob)
  ymax.H0 <- dbinom(round(null.prob * size), size = size, prob = null.prob)

  ymax <- max(c(ymax.HA, ymax.H0))
  ylim <- c(0, ymax * 1.20)

  plot.window.dim <- dev.size("cm")
  cex.legend <- min(plot.window.dim[1] / 18, plot.window.dim[2] / 15)
  cex.title <- min(plot.window.dim[1] / 11, plot.window.dim[2] / 11)
  cex.label <- min(plot.window.dim[1] / 12, plot.window.dim[2] / 12)

  # plots
  if (alternative == "two.one.sided") {

    .plot.binom.dist(prob = null.prob[1], size = size, xlim = xlim, ylim = ylim, type = 1, ticks = TRUE)
    par(new = TRUE)
    .plot.binom.dist(prob = null.prob[2], size = size, xlim = xlim, ylim = ylim, type = 1)
    par(new = TRUE)
    .plot.binom.dist(prob = prob, size = size, xlim = xlim, ylim = ylim, type = 2)


    text(round(prob * size), ymax.HA * 1.05,
         labels = expression(H[1]),
         cex = cex.legend, col = adjustcolor(4, alpha.f = 1))

    text(round(null.prob * size), ymax.H0 * 1.05,
         labels = expression(H[0]),
         cex = cex.legend, col = adjustcolor(2, alpha.f = 1))

  } else {

    .plot.binom.dist(prob = prob, size = size, xlim = xlim, ylim = ylim, type = 2, ticks = TRUE)
    par(new = TRUE)
    .plot.binom.dist(prob = null.prob, size = size, xlim = xlim, ylim = ylim, type = 1)

    text(round(prob * size), ymax.HA * 1.05,
         labels = expression(H[1]),
         cex = cex.legend, col = adjustcolor(4, alpha.f = 1))

    text(round(null.prob * size), ymax.H0 * 1.05,
         labels = expression(H[0]),
         cex = cex.legend, col = adjustcolor(2, alpha.f = 1))

  } # end of plots

  # draw vertical lines for critical region
  segments(x0 = q.binom.alpha, y0 = rep(0, length(q.binom.alpha)),
           x1 = q.binom.alpha, y1 = y.binom.alpha, col = 2, lty = 2, lwd = 2)

  # paint regions
  if (alternative == "two.one.sided") {

    if (prob > min(null.prob) && prob < max(null.prob)) {
      # equivalence
      .paint.binom.dist(prob = null.prob[1], size = size, xlim = c(q.binom.alpha[1], max(xlim)), type = 1)
      .paint.binom.dist(prob = null.prob[2], size = size, xlim = c(q.binom.alpha[2], min(xlim)), type = 1)

      .paint.binom.dist(prob = prob, size = size, xlim = c(q.binom.alpha[1] - 1, min(xlim)), type = 2)
      .paint.binom.dist(prob = prob, size = size, xlim = c(q.binom.alpha[2] + 1, max(xlim)), type = 2)

      ifelse(q.binom.alpha[2] > q.binom.alpha[1],
             power <- .paint.binom.dist(prob = prob, size = size, xlim = q.binom.alpha, type = 3),
             power <- 0)

    } else {
      # minimal effect
      .paint.binom.dist(prob = null.prob[1], size = size, xlim = c(min(xlim), q.binom.alpha[1]), type = 1)
      .paint.binom.dist(prob = null.prob[2], size = size, xlim = c(q.binom.alpha[2], max(xlim)), type = 1)

      if (prob > max(null.prob)) {
        .paint.binom.dist(prob = prob, size = size, xlim = c(q.binom.alpha[2] - 1, min(xlim)), type = 2)
      } else {
        .paint.binom.dist(prob = prob, size = size, xlim = c(max(xlim), q.binom.alpha[1] + 1), type = 2)
      }

      power.upper <- .paint.binom.dist(prob = prob, size = size, xlim = c(q.binom.alpha[2], max(xlim)), type = 3)
      power.lower <- .paint.binom.dist(prob = prob, size = size, xlim = c(q.binom.alpha[1], min(xlim)), type = 3)
      power <- power.upper + power.lower

    }

  } else if (alternative == "two.sided") {

    .paint.binom.dist(prob = null.prob, size = size, xlim = c(q.binom.alpha[1], min(xlim)), type = 1)
    .paint.binom.dist(prob = null.prob, size = size, xlim = c(q.binom.alpha[2], max(xlim)), type = 1)

    .paint.binom.dist(prob = prob, size = size, xlim = c(q.binom.alpha[1] + 1, q.binom.alpha[2] - 1), type = 2)

    power.left <- .paint.binom.dist(prob = prob, size = size, xlim = c(min(xlim), q.binom.alpha[1]), type = 3)
    power.right <- .paint.binom.dist(prob = prob, size = size, xlim = c(q.binom.alpha[2], max(xlim)), type = 3)
    power <- power.left + power.right

  } else {

    ifelse(prob < null.prob,
           .paint.binom.dist(prob = null.prob, size = size, xlim = c(q.binom.alpha, min(xlim)), type = 1),
           .paint.binom.dist(prob = null.prob, size = size, xlim = c(q.binom.alpha, max(xlim)), type = 1))

    ifelse(prob < null.prob,
           .paint.binom.dist(prob = prob, size = size, xlim = c(q.binom.alpha + 1, max(xlim)), type = 2),
           .paint.binom.dist(prob = prob, size = size, xlim = c(q.binom.alpha - 1, min(xlim)), type = 2))

    ifelse(prob < null.prob,
           power <- .paint.binom.dist(prob = prob, size = size, xlim = c(min(xlim), q.binom.alpha), type = 3),
           power <- .paint.binom.dist(prob = prob, size = size, xlim = c(q.binom.alpha, max(xlim)), type = 3))

  } # end of paint regions

  # axes labels and subtitle
  title(main = plot.main, line = 2, cex.main = cex.title)
  title(sub = plot.sub, line = 3, cex.sub = cex.title)
  title(ylab = "Probability Density", line = 2.2, cex.lab = cex.label,
        col.lab = "grey30")
  title(xlab = paste0("X ~ Binomial(n =", size, ")"),
        line = 2.2, cex.lab = cex.label,
        col.lab = "grey30")

  if (power < 0) power <- 0
  alpha <- round(alpha, 2)
  beta <- round(1 - power, 2)
  power <- round(power, 2)

  legend("topright", cex = cex.legend,
         c(as.expression(bquote(Power == .(power))),
           as.expression(bquote(alpha == .(alpha))), # %~~%
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

} # end of .plot.binom.t1t2()
