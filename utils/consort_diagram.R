# =============================================================================
# CONSORT-style inclusion flow diagram (ggplot2 only)
# =============================================================================
# Renders the cohort attrition funnel from the structured attrition tibble
# produced by utils/attrition_log.R, using only ggplot2 primitives so no
# DiagrammeR / consort package dependency is introduced.
#
# The main vertical chain shows n_remaining at each step; right-hand-side branch
# boxes show the number excluded and the reason at each transition.

library(tidyverse)

render_consort <- function(attrition_tbl, title = NULL) {
  required_columns <- c("step_order", "step_label", "n_remaining", "n_excluded", "exclusion_reason")
  missing_columns <- setdiff(required_columns, names(attrition_tbl))
  
  if (length(missing_columns) > 0L) {
    stop("render_consort: missing required column(s): ", paste(missing_columns, collapse = ", "))
  }
  
  at <- attrition_tbl %>% arrange(step_order)
  n_steps <- nrow(at)
  
  if (n_steps == 0L) {
    stop("render_consort: empty attrition table")
  }
  
  if (!is.na(at$n_excluded[[1]]) && at$n_excluded[[1]] > 0) {
    stop("render_consort: the first step cannot have a preceding exclusion branch.")
  }
  
  ys <- seq(1, 0, length.out = n_steps)
  main_x <- 0.27
  excl_x <- 0.78
  half_gap <- if (n_steps > 1L) (ys[1] - ys[2]) / 2 * 0.55 else 0.1
  
  main <- tibble(
    x = main_x, y = ys,
    label = paste0(
      stringr::str_wrap(at$step_label, width = 38),
      "\nn = ", format(at$n_remaining, big.mark = ",")
    )
  )
  
  vertical <- tibble(
    x = main_x, xend = main_x,
    y = head(ys, -1L) - half_gap,
    yend = tail(ys, -1L) + half_gap
  )
  
  exclusion_indices <- which(!is.na(at$n_excluded) & at$n_excluded > 0)
  
  p <- ggplot() +
    geom_segment(
      data = vertical, aes(x = x, y = y, xend = xend, yend = yend),
      arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed"),
      linewidth = 0.4, colour = "grey30"
    )
  
  if (length(exclusion_indices) > 0L) {
    excluded <- tibble(
      x = excl_x,
      y = (ys[exclusion_indices] + ys[exclusion_indices - 1L]) / 2,
      label = paste0(
        stringr::str_wrap(coalesce(at$exclusion_reason[exclusion_indices], "Excluded"), width = 38),
        "\nn = ", format(at$n_excluded[exclusion_indices], big.mark = ",")
      )
    )
    
    horizontal <- tibble(
      x = main_x + 0.02, xend = excl_x - 0.18,
      y = excluded$y, yend = excluded$y
    )
    
    p <- p +
      geom_segment(
        data = horizontal, aes(x = x, y = y, xend = xend, yend = yend),
        arrow = grid::arrow(length = grid::unit(0.12, "cm"), type = "closed"),
        linewidth = 0.4, colour = "#D55E00"
      ) +
      geom_label(
        data = excluded, aes(x = x, y = y, label = label),
        fill = "#FBE7DD", colour = "#D55E00", linewidth = 0.3,
        size = 2.8, lineheight = 0.95,
        label.padding = grid::unit(0.35, "lines")
      )
  }
  
  p <- p +
    geom_label(
      data = main, aes(x = x, y = y, label = label),
      fill = "grey95", colour = "black", linewidth = 0.3,
      size = 3.1, lineheight = 0.95,
      label.padding = grid::unit(0.45, "lines")
    ) +
    coord_cartesian(xlim = c(0, 1.08), ylim = c(-0.10, 1.10), clip = "off") +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.margin = margin(t = 20, r = 35, b = 20, l = 35)
    )
  
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  
  p
}
