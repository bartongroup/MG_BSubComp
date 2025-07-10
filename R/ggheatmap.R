ggheatmap <- function(tab, order.col = TRUE, order.row = TRUE, dendro.line.size = 0.5,
                      text.size = 12, legend.text.size = 12, legend.name = "value",
                      dist.method = "euclidean", clust.method = "complete", title = NULL,
                      with.y.text = FALSE, with.x.text = TRUE, palette = "distiller",
                      max.fc = NULL) {
  
  d <- tab |>
    as_tibble(rownames = "rowname") |>
    mutate(rowname = factor(rowname, levels = rownames(tab))) |>
    pivot_longer(-rowname, names_to = "variable", values_to = "value") |>
    mutate(variable = factor(variable, levels = colnames(tab))) |>
    mutate(value = as.numeric(value))
  
  dd <- function(X) {
    X |> dist(method = dist.method) |> hclust(method = clust.method) |> as.dendrogram()
  }
  
  # Cluster rows
  if (order.row) {
    dd.row <- dd(tab)
    row.ord <- order.dendrogram(dd.row)
    ordered_row_names <- row.names(tab[row.ord, ])
    d$rowname <- factor(d$rowname, levels = ordered_row_names)
  }
  
  # Cluster columns
  if (order.col) {
    dd.col <- dd(t(tab))
    col.ord <- order.dendrogram(dd.col)
    ordered_col_names <- colnames(tab[, col.ord])
    d$variable <- factor(d$variable, levels = ordered_col_names)
  }
  
  if (is.null(max.fc)) {
    mx <- max(abs(d$value), na.rm = TRUE)
  } else {
    mx <- max.fc
  }
  heat_plot <- ggplot(d, aes(x = variable, y = rowname, fill = value)) +
    geom_tile() +
    theme_classic() +
    theme(
      axis.line = element_blank(),
      text = element_text(size = text.size)
    ) +
    scale_y_discrete(position = "right", expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, title = title)
  
  if (palette == "viridis") {
    heat_plot <- heat_plot +
      scale_fill_viridis_c(option = "cividis", name = legend.name)
  } else {
    heat_plot <- heat_plot +
      scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-mx, mx), name = legend.name)
  }
  
  if (with.x.text) {
    heat_plot <- heat_plot +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  } else {
    heat_plot <- heat_plot +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  if (with.y.text) {
    heat_plot <- heat_plot +
      theme(axis.text.y = element_text(size = 8))
  } else {
    heat_plot <- heat_plot +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  
  final_plot <- heat_plot
  if (order.row) {
    dendro_data_row <- ggdendro::dendro_data(dd.row, type = "rectangle")
    dendro_row <- axis_canvas(heat_plot, axis = "y", coord_flip = TRUE) +
      geom_segment(data = ggdendro::segment(dendro_data_row), aes(y = -y, x = x, xend = xend, yend = -yend), linewidth = dendro.line.size) +
      coord_flip()
    final_plot <- final_plot |>
      cowplot::insert_yaxis_grob(dendro_row, grid::unit(0.2, "null"), position = "left")
  }
  
  if (order.col) {
    dendro_data_col <- ggdendro::dendro_data(dd.col, type = "rectangle")
    dendro_col <- cowplot::axis_canvas(heat_plot, axis = "x") +
      geom_segment(data = ggdendro::segment(dendro_data_col), aes(x = x, y = y, xend = xend, yend = yend), linewidth = dendro.line.size)
    final_plot <- final_plot |>
      cowplot::insert_xaxis_grob(dendro_col, grid::unit(0.2, "null"), position = "top")
  }
  
  ggdraw(final_plot)
}
