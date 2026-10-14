# Run from the package root: source("dev/hex/make_hex_80s_v2.R")
# crywolf hex sticker, 80s version with the PhyloPic wolf
# Wolf silhouette: Canis lupus monstrabilis by Tracy A. Heath, PhyloPic, CC0 1.0
# https://www.phylopic.org/images/e4e306cd-73b6-4ca3-a08c-753a856f7f12
library(hexSticker); library(ggplot2); library(showtext)
sysfonts::font_add_google("Righteous", "righteous")
showtext_auto()

bg <- "#241042"; magenta <- "#FF3EA5"; cyan <- "#2DE2E6"; sun1 <- "#FFB627"
cream <- "#F4EBD9"

wolf_png <- rsvg::rsvg("dev/hex/pp/wolf_cream.svg", width = 1536)   # RGBA array, cream fill
wolf_box <- c(xmin = 0.4, xmax = 10.2, ymin = 0, ymax = 9.0)  # keeps 1536:1409

circ <- function(cx, cy, r, n = 200, from = 0, to = 360) {
  th <- seq(from, to, length.out = n) * pi / 180
  data.frame(x = cx + r * cos(th), y = cy + r * sin(th))
}
sun <- circ(6.3, 7.6, 3.1)
bands <- data.frame(ymin = c(4.9, 5.65, 6.3, 6.85), ymax = c(5.3, 5.95, 6.5, 6.95))
arcs <- do.call(rbind, lapply(c(0.8, 1.45, 2.1), function(r)
  transform(circ(0.55, 8.75, r, 60, 95, 200), r = r)))
set.seed(7)
eq <- data.frame(x = seq(0.4, 10.2, by = 0.5))
eq$h <- 0.25 + 0.9 * abs(sin(eq$x * 1.3)) * runif(nrow(eq), 0.5, 1)

p <- ggplot() +
  geom_polygon(data = sun, aes(x, y), fill = sun1, colour = NA) +
  geom_rect(data = bands, aes(xmin = 2.9, xmax = 9.7, ymin = ymin, ymax = ymax),
            fill = bg, colour = NA) +
  geom_rect(data = eq, aes(xmin = x - 0.17, xmax = x + 0.17,
                           ymin = -1.0, ymax = -1.0 + h), fill = cyan) +
  geom_path(data = arcs, aes(x, y, group = r), colour = magenta,
            linewidth = 1.3, lineend = "round") +
  annotation_raster(wolf_png, wolf_box["xmin"], wolf_box["xmax"],
                    wolf_box["ymin"], wolf_box["ymax"], interpolate = TRUE) +
  coord_fixed(xlim = c(-2.0, 11.0), ylim = c(-1.6, 11.4), expand = FALSE) +
  theme_void() + theme_transparent()

make <- function(file, size, dpi = NULL) {
  if (!is.null(dpi)) showtext_opts(dpi = dpi)
  args <- list(p, package = "crywolf",
               p_family = "righteous", p_size = size, p_color = cream,
               p_x = 1, p_y = 0.42,
               s_x = 1.0, s_y = 1.15, s_width = 1.15, s_height = 1.15,
               h_fill = bg, h_color = magenta, h_size = 1.6,
               filename = file)
  if (!is.null(dpi)) args$dpi <- dpi
  do.call(sticker, args)
}
make("man/figures/logo.png", 7, dpi = 600)
make("man/figures/logo.svg", 7 * 96 / 600 * 5.4, dpi = 96)
