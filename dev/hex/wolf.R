library(ggplot2)
# original howling-wolf silhouette, drawn by hand (sitting, muzzle up-left)
wolf <- data.frame(matrix(c(
  0.90, 9.70,   # nose tip
  1.60, 9.75,
  2.40, 9.55,   # top of muzzle
  3.05, 9.30,   # stop
  3.35, 9.55,
  3.75, 10.45,  # ear tip
  4.10, 9.55,
  4.45, 9.20,   # back of head
  5.05, 8.70,
  5.65, 7.90,   # nape
  6.30, 6.85,   # back
  7.05, 5.60,
  7.70, 4.30,
  8.10, 3.05,   # rump
  8.30, 2.20,
  8.80, 2.35,   # bushy tail
  9.60, 2.35,
  10.35, 1.85,
  10.05, 1.55,
  10.40, 1.20,
  9.70, 0.95,
  8.90, 0.85,
  8.05, 0.55,
  7.70, 0.25,   # hind paw
  5.60, 0.25,
  5.55, 0.55,
  6.05, 0.75,
  6.20, 1.25,   # haunch
  5.55, 1.45,
  5.15, 0.70,
  5.05, 0.25,   # front paws
  3.85, 0.25,
  3.90, 0.55,
  4.30, 0.75,
  4.35, 1.70,   # front leg
  4.20, 3.00,
  3.95, 4.10,   # chest ruff
  3.78, 4.45,
  3.70, 4.90,
  3.38, 5.30,
  3.35, 5.75,
  3.02, 6.15,
  2.95, 6.60,   # throat
  2.40, 7.40,
  1.95, 8.05,
  1.35, 8.55,   # lower jaw tip
  1.55, 8.85,
  2.15, 8.90,   # mouth corner
  1.45, 9.30,   # upper jaw underside
  0.95, 9.45
), ncol = 2, byrow = TRUE))
names(wolf) <- c("x", "y")
# Chaikin corner cutting, keeping the ear tip, nose and paws crisp
chaikin <- function(d, keep, n = 2) {
  for (k in seq_len(n)) {
    m <- nrow(d); out <- list(); kk <- integer(0)
    for (i in seq_len(m)) {
      p <- d[i, ]; q <- d[ifelse(i == m, 1, i + 1), ]
      if (i %in% keep) { out[[length(out) + 1]] <- p; kk <- c(kk, length(out)) }
      out[[length(out) + 1]] <- 0.75 * p + 0.25 * q
      out[[length(out) + 1]] <- 0.25 * p + 0.75 * q
    }
    d <- do.call(rbind, out); keep <- kk
  }
  d
}
wolf <- chaikin(wolf, keep = c(1, 6, 18, 20))

arc <- function(r, cx = 0.75, cy = 9.85, from = 95, to = 205) {
  th <- seq(from, to, length.out = 60) * pi / 180
  data.frame(x = cx + r * cos(th), y = cy + r * sin(th), r = r)
}
arcs <- do.call(rbind, lapply(c(0.9, 1.6, 2.3), arc))

circ <- function(cx, cy, r, n = 200) {
  th <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + r * cos(th), y = cy + r * sin(th))
}
moon_df <- circ(5.6, 8.4, 2.5)

wolf_plot <- function(fill = "#F4EBD9", wave = "#F2A541", moon = "#3C5A80",
                      lw = 1.2) {
  ggplot() +
    geom_polygon(data = moon_df, aes(x, y), fill = moon, colour = NA) +
    geom_path(data = arcs, aes(x, y, group = r), colour = wave,
              linewidth = lw, lineend = "round") +
    geom_polygon(data = wolf, aes(x, y), fill = fill, colour = fill,
                 linewidth = 0.3, linejoin = "round") +
    coord_fixed(xlim = c(-1.9, 10.6), ylim = c(-0.2, 12.4), expand = FALSE) +
    theme_void() + theme(legend.position = "none")
}
