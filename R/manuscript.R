mn_plot_volcano_toxins <- function(de_toxan, toxan, fdr_limit = 0.01, logfc_limit = 1) {
  d <- de_toxan |> 
    mutate(
      toxan = id %in% toxan$id,
      sig = FDR < fdr_limit & logFC < -logfc_limit
    )
  
  d_ntox <- d |> filter(!toxan)
  d_tox <- d |> filter(toxan)
  d_toxsig <- d |> filter(toxan & sig)
  
  ggplot(d, aes(x = logFC, y = -log10(PValue), label = gene_symbol)) +
    th +
    theme(legend.position = "none") +
    geom_point(data = d_ntox, colour = "grey80", size = 0.6, alpha = 0.3) +
    geom_point(data = d_tox, aes(colour = sig), size = 1.3) +
    scale_colour_manual(values = c("deepskyblue2", "black")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    geom_text_repel(data = d_toxsig, size = 4.0, min.segment.length = 0.03, force = 20, segment.colour = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
    #facet_wrap(~ contrast) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P))
}


mn_plot_browser <- function(de, genes, operons, range, fdr_limit = 0.01, arrow_limit = 150000, text_limit = 300000,
                         operon_limit = 100000) {
  dat <- de |> 
    left_join(genes, by = join_by(id, gene_symbol)) |> 
    mutate(gene_symbol = if_else(is.na(gene_symbol), id, gene_symbol))
  
  width <- range$end - range$start
  mid <- range$start + width / 2

  g <- dat |> 
    filter(chr == range$chr & end > range$start & start < range$end) |> 
    mutate(
      sig = FDR < fdr_limit,
      sig = factor(sig, levels = c(FALSE, TRUE))
    )
  mx <- max(abs(g$logFC)) * 1.3
  dm <- tibble(
    dx = c(range$start, range$end),
    dy = c(-mx, mx)
  )

  ar <- g |> 
    mutate(
      xstart = if_else(strand == "+", start, end),
      xend = if_else(strand == "+", end, start),
      yarr = if_else(strand == "+", mx * 0.01, -mx * 0.01)
    )
  
  op <- g |> 
    mutate(pos = (start + end) / 2) |> 
    distinct(id, pos) |> 
    left_join(operons, by = join_by(id)) |> 
    drop_na() 

  wdt <- case_when(
    width < 10000 ~ sprintf("%d bp", width),
    width >= 10000 & width < 1000000 ~ sprintf("%.1f kbp", width / 1000),
    width >= 1000000 ~ sprintf("%.2f Mbp", width / 1000000)
  )
  
  pl <- g |> 
    ggplot() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(size = 12),
      legend.position = "none"
    ) +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = logFC, fill = sig), alpha = 0.6, colour = "black", linewidth = 0.3) +
    geom_blank(data = dm, aes(x = dx, y = dy)) +
    geom_hline(yintercept = 0, linewidth = 0.6, alpha = 0.6) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.2) +
    geom_vline(xintercept = mid, linetype = "dashed", alpha = 0.2) +
    #scale_colour_manual(values = c("royalblue3", "springgreen4", "grey"), drop = FALSE) +
    scale_fill_manual(values = c("grey80", "brown"), drop = FALSE) +
    scale_x_continuous(expand = c(0, 0), labels = scales::label_comma()) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(range$start, range$end)) +
    labs(x = paste("Position along", range$chr), y = expression(log[2]~FC))

  if(width < arrow_limit) {
    pl <- pl + geom_segment(data = ar, aes(x = xstart, xend = xend, y = logFC, yend = logFC), 
                           arrow = arrow(length = unit(0.15, "inches"), type = "closed"), linewidth = 1.3)
  }

  if(width < text_limit) {
    pl <- pl + geom_text(aes(label = gene_symbol, x = (start + end) / 2, y = logFC, vjust = -sign(logFC) + 0.5))
  }

  if(width < operon_limit) {
     pl <- pl + 
      geom_point(data = op, aes(x = pos, y = 0), size = 3, colour = "darkgreen") +
      geom_line(data = op, aes(x = pos, y = 0, group = operon), linewidth = 2, colour = "darkgreen")
  }

  pl
}



mn_plot_gse_nes <- function(gse, fdr_limit = 0.05) {
  list_rbind(gse, names_to = "ontology") |> 
    filter(fdr < fdr_limit) |> 
    mutate(
      term_name = term_name |>
        str_remove(" - .*") |> 
        str_wrap(width = 30) |> 
        fct_reorder(nes)
    ) |> 
    ggplot(aes(x = nes, y = term_name)) +
    th +
    theme(
      legend.position = "right"
    ) +
    geom_segment(aes(xend = 0, yend = term_name), colour = "grey70") +
    geom_point(aes(size = size, fill = fdr), shape = 21, colour = "grey70") +
    geom_vline(xintercept = 0, alpha = 0.5) +
    scale_fill_viridis_c(option = "cividis") +
    scale_size_area() +
    scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    labs(x = "Normalized Enrichment Score (NES)", y = NULL, fill = "FDR", size = "Count")
}

mn_plot_gse_heatmap <- function(gse, de, genes, fdr_limit = 0.05, log_fc_limit = 1) {
  tab <- list_rbind(gse, names_to = "ontology") |> 
    mutate(
      term_name = str_remove(term_name, " - .*"),
      term_name = fct_reorder(term_name, nes)
    ) |> 
    filter(fdr < fdr_limit) |>
    unnest(leading_edge) |> 
    left_join(de |> select(id, log_fc = logFC, de_fdr = FDR), by = join_by(leading_edge == id)) |>
    left_join(genes |> select(id, gene_symbol), by = join_by(leading_edge == id)) |> 
    filter(de_fdr < fdr_limit & abs(log_fc) > log_fc_limit) |> 
    mutate(term_name = fct_reorder(term_name, nes)) |>
    pivot_wider(id_cols = gene_symbol, names_from = term_name, values_from = log_fc, values_fill = 0) |> 
    column_to_rownames("gene_symbol")

  ggheatmap(tab, with.y.text = TRUE, order.col = FALSE, with.x.dendro = FALSE, with.y.dendro = FALSE,
            legend.name = expression(log[2]~FC), legend.text.size = 8)
}