# Font theme
theme_nature_fonts <- function(base_size = 5) {
  theme(
    axis.text = element_text(size = base_size),
    axis.title = element_text(size = base_size),
    plot.title = element_text(size = base_size, hjust = 0.5),
    plot.subtitle = element_text(size = base_size, hjust = 0.5),
    legend.title = element_text(size = base_size),
    legend.text = element_text(size = base_size),
    strip.text = element_text(size = base_size)
  )
}

# Small legend theme
theme_small_legend <- function(...) {
  theme(
    legend.key.height = unit(0.3, "cm"),  
    legend.key.width = unit(0.3, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    ...
  )
}

# White background
theme_white_background <- function(...) {
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),   
    panel.grid.major = element_blank(),                           
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(0.5, "mm"),                    
    ...
  )
}

# White strip background
theme_white_strip <- function(...) {
  theme(
    strip.background = element_rect(fill = "white", color = NA),  
    strip.text = element_text(color = "black"),   
    strip.placement = "inside",                                  
    ...
  )
}
