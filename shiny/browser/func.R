expand_range <- function(range, chrlens, margin = 0.5) {
  width <- range$end - range$start
  delta <- width * margin
  chlen <- chrlens[[range$chr]]
  range$start <- max(range$start - delta, 0)
  range$end <- min(range$end + delta, chlen)
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


plot_browser <- function(dat, ops, range, fdr_limit = 0.01, arrow_limit = 150000, text_limit = 300000,
                         operon_limit = 100000) {
  width <- range$end - range$start
  mid <- range$start + width / 2

  g <- dat |> 
    filter(chr == range$chr & end > range$start & start < range$end) |> 
    mutate(
      sig = fdr < fdr_limit,
      sig = factor(sig, levels = c(FALSE, TRUE))
    )
  mx <- max(abs(g$log_fc)) * 1.3
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
    left_join(ops, by = join_by(id)) |> 
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
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = log_fc, fill = sig), alpha = 0.6, colour = "black", linewidth = 0.3) +
    geom_blank(data = dm, aes(x = dx, y = dy)) +
    geom_hline(yintercept = 0, linewidth = 0.6, alpha = 0.6) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", alpha = 0.2) +
    geom_vline(xintercept = mid, linetype = "dashed", alpha = 0.2) +
    #scale_colour_manual(values = c("royalblue3", "springgreen4", "grey"), drop = FALSE) +
    scale_fill_manual(values = c("grey80", "brown"), drop = FALSE) +
    scale_x_continuous(expand = c(0, 0), labels = label_comma()) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(xlim = c(range$start, range$end)) +
    labs(x = NULL, y = expression(log[2]~FC), title = wdt)

  if(width < arrow_limit) {
    pl <- pl + geom_segment(data = ar, aes(x = xstart, xend = xend, y = log_fc, yend = log_fc), 
                           arrow = arrow(length = unit(0.15, "inches"), type = "closed"), linewidth = 1.3)
  }

  if(width < text_limit) {
    pl <- pl + geom_text(aes(label = gene_symbol, x = (start + end) / 2, y = log_fc, vjust = -sign(log_fc) + 0.5))
  }

  if(width < operon_limit) {
     pl <- pl + 
      geom_point(data = op, aes(x = pos, y = 0), size = 3, colour = "darkgreen") +
      geom_line(data = op, aes(x = pos, y = 0, group = operon), linewidth = 2, colour = "darkgreen")
  }

  pl
}