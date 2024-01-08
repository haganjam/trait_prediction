#' @title theme_meta
#' @description A customised plotting theme to equalise the formatting of
#'  all plots plotted with ggplot2
#' @author James G. Hagan (james_hagan(at)outlook.com)
theme_meta <- function(base_size = 12, base_family = "") {
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(
      fill = "NA",
      color = "black",
      linewidth = 0.75,
      linetype = "solid"
    ),
    axis.line.x = element_line(color = "black", linewidth = 0.2),
    axis.line.y = element_line(color = "black", linewidth = 0.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(-0.16, "cm"),
    axis.title.x = element_text(
      colour = "black",
      size = 12,
      face = "plain",
      margin = margin(5, 0, 0, 0, "pt")
    ),
    axis.title.y = element_text(
      colour = "black",
      size = 12,
      face = "plain",
      margin = margin(0, 5, 0, 0, "pt")
    ),
    axis.text.x = element_text(
      colour = "black",
      size = 11,
      face = "plain",
      margin = margin(10, 0, 0, 0, "pt")
    ),
    axis.text.y = element_text(
      colour = "black",
      size = 11,
      face = "plain",
      margin = margin(0, 10, 0, 0, "pt")
    ),
    axis.ticks.x = element_line(colour = "black", size = 0.4),
    axis.ticks.y = element_line(colour = "black", size = 0.4),
    legend.text = element_text(colour = "black", size = 10, face = "plain"),
    legend.title = element_text(colour = "black", size = 10, face = "plain")
  )
}

#' @title theme_transparent
#' @description A customised plotting theme to make the plots transparent
#' for use in presentations.
#' @author James G. Hagan (james_hagan(at)outlook.com)

theme_transparent <- function() {
  
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )
  
  }

# model fit statistics
lm_eqn <- function(m){
  eq <- substitute(italic(r)^2~"="~r2,
                   list(r2 = format(summary(m)$r.squared, digits = 2)))
  as.character(as.expression(eq))
}







