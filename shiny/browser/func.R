expand_range <- function(range, margin = 0.5) {
  width <- range[2] - range[1]
  delta <- width * margin
  range[1] <- range[1] - delta
  range[2] <- range[2] + delta
  range
}

parse_str_range <- function(srange) {
  s1 <- str_split_1(srange, ":")
  s2 <- str_split_1(s1[2], "-")
  list(
    chr = s1[1],
    start = as.integer(s2[1]),
    end = as.integer(s2[2])
  )
}

make_str_range <- function(range) {
  start <- round(range$start)
  end <- round(range$end)
  str_glue("{range$chr}:{start}-{end}")
}

make_de_table <- function(de, fdr_limit) {
  de |> 
    select(id, gene_symbol, description, log_fc, fdr) |> 
    filter(fdr < fdr_limit) |> 
    arrange(log_fc)
}

plot_browser <- function(genes, de, range, fdr_limit = 0.01) {
  g <- genes |> 
    filter(chr == range$chr & end > range$start & start < range$end) |> 
    inner_join(select(de, -gene_symbol), by = join_by(id)) |> 
    mutate(
      sig = fdr < fdr_limit,
      sig = factor(sig, levels = c(FALSE, TRUE))
    )
  mx <- max(abs(g$log_fc)) * 1.3
  dm <- tibble(
    dx = c(range$start, range$end),
    dy = c(-mx, mx)
  )
  print(g)
  print(dm)

  g |> 
    ggplot() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.ticks.y = element_blank(),
      ##axis.text.y = element_blank(),
      legend.position = "none"
    ) +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = log_fc, fill = sig), alpha = 0.6, colour = "black") +
    geom_text(aes(label = gene_symbol, x = (start + end) / 2, y = log_fc, vjust = -sign(log_fc) + 0.5)) +
    geom_blank(data = dm, aes(x = dx, y = dy)) +
    geom_hline(yintercept = 0, linewidth = 1.2, alpha = 0.6) +
    scale_fill_manual(values = c("grey80", "brown"), drop = FALSE) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(range$start, range$end)) +
    labs(x = NULL, y = expression(log[2]~FC))

}