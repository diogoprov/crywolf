library(hexSticker); library(ggplot2); library(showtext)
source("wolf.R")
sysfonts::font_add_google("Fredoka", "fredoka")
showtext_auto(); showtext_opts(dpi = 600)

bg <- "#1D2D44"; border <- "#F2A541"
p <- wolf_plot(fill = "#F4EBD9", wave = "#F2A541", moon = "#3E5C76", lw = 1.3) + theme_transparent()

sticker(p, package = "crywolf",
        p_family = "fredoka", p_fontface = "bold", p_size = 6.8, p_color = "#F4EBD9",
        p_x = 1, p_y = 0.43,
        s_x = 1.0, s_y = 1.14, s_width = 1.0, s_height = 1.0,
        h_fill = bg, h_color = border, h_size = 1.6,
        filename = "crywolf_hex.png", dpi = 600)
