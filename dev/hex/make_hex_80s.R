library(hexSticker); library(ggplot2); library(showtext)
source("wolf.R")
sysfonts::font_add_google("Righteous", "righteous")
showtext_auto(); showtext_opts(dpi = 600)

bg <- "#241042"; magenta <- "#FF3EA5"; cyan <- "#2DE2E6"; sun1 <- "#FFB627"
cream <- "#F4EBD9"

# 80s sun: solid disc cut by horizontal bands of the background colour
sun <- circ(5.7, 7.6, 3.0)
bands <- data.frame(ymin = c(5.0, 5.75, 6.4, 6.95),
                    ymax = c(5.4, 6.05, 6.6, 7.05))
# equalizer, a nod to the acoustic data
set.seed(7)
eq <- data.frame(x = seq(0.2, 10.2, by = 0.5))
eq$h <- 0.25 + 0.9 * abs(sin(eq$x * 1.3)) * runif(nrow(eq), 0.5, 1)

p <- ggplot() +
  geom_polygon(data = sun, aes(x, y), fill = sun1, colour = NA) +
  geom_rect(data = bands, aes(xmin = 2.4, xmax = 9.0, ymin = ymin, ymax = ymax),
            fill = bg, colour = NA) +
  geom_rect(data = eq, aes(xmin = x - 0.17, xmax = x + 0.17,
                           ymin = -1.6, ymax = -1.6 + h), fill = cyan) +
  geom_path(data = arcs, aes(x, y, group = r), colour = magenta,
            linewidth = 1.3, lineend = "round") +
  geom_polygon(data = wolf, aes(x, y), fill = cream, colour = cream,
               linewidth = 0.3, linejoin = "round") +
  coord_fixed(xlim = c(-1.9, 10.6), ylim = c(-1.8, 12.4), expand = FALSE) +
  theme_void() + theme_transparent()

sticker(p, package = "crywolf",
        p_family = "righteous", p_size = 7, p_color = cream,
        p_x = 1, p_y = 0.42,
        s_x = 1.0, s_y = 1.17, s_width = 1.12, s_height = 1.12,
        h_fill = bg, h_color = magenta, h_size = 1.6,
        filename = "crywolf_hex_80s.png", dpi = 600)
showtext_opts(dpi = 96)
sticker(p, package = "crywolf",
        p_family = "righteous", p_size = 7 * 96 / 600 * 5.4, p_color = cream,
        p_x = 1, p_y = 0.42,
        s_x = 1.0, s_y = 1.17, s_width = 1.12, s_height = 1.12,
        h_fill = bg, h_color = magenta, h_size = 1.6,
        filename = "crywolf_hex_80s.svg")
