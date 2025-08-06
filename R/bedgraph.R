read_bedgraphs <- function(config, meta) {
  path <- file.path(config$data_path, config$dir$bedgraph)
  map2(meta$raw_sample, meta$sample, function(rsam, sam) {
    bfile <- file.path(path, glue::glue("{rsam}.{config$window_size}.bedgraph"))
    read_tsv(bfile, col_names = c("chr", "start", "end", "score"), col_types = "ciid") |> 
      add_column(sample = sam)
  }) |> 
    list_rbind()
}

plot_bedgraphs_region <- function(bg, .chr, .start, .end, log_scale = FALSE) {
  b <- bg |>
    mutate(pos = start / 1e6) |> 
    filter(chr == .chr & pos >= .start & pos <= .end) |> 
    select(-c(chr, start, end))
  ylab <- "Count"
  if(log_scale) {
    b <- b |>
      mutate(score = log10(score + 1))
    ylab = expression(log[10]~(count+1))
  }
  b |> 
    ggplot() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
    ) +
    geom_step(aes(x = pos, y = score)) +
    facet_grid(sample ~ .) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = "Position (Mb)", y = ylab)
}


plot_bedgraph_chromosomes <- function(bg, samp, chroms=CHROMOSOMES, ymax=NULL, brks=waiver(), log_scale=FALSE, highlight_peaks=FALSE, scales="fixed") {
  bin_size <- (bg$start[2] - bg$start[1]) 
  b <- bg |> 
    filter(chr %in% chroms & sample == samp) |> 
    mutate(
      chr =factor(chr, levels=chroms),
      pos = start / 1e6,
      y = score / bin_size
    )
  if(log_scale) b <- b |> filter(y > 0)
  
  g <- ggplot(b, aes(x = pos, y = y)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank(),
      panel.border = element_blank(),
      #axis.line.x.bottom = element_line(size=0.1),
      axis.line.y.left = element_line(colour = "grey70")
    ) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_step() +
    facet_grid(chr ~ ., scales=scales) +
    labs(x = "Position (Mb)", y = "Count density") 
  if(!is.null(ymax)) g <- g + coord_cartesian(ylim = c(0, ymax))
  if(log_scale) {
    g <- g + scale_y_log10(breaks = brks)
  } else {
    g <-  g + scale_y_continuous(breaks = brks, trans="sqrt")
  }
  if(highlight_peaks) g <- g + geom_point(data=b |> filter(peak), colour="red")
  g
}


plot_bedgraphs_depth_distribution <- function(bg, what="control") {
  bg |> 
    mutate(val = !!sym(what)) |> 
    filter(val > 0) |>
    group_by(val) |>
    tally() |>
  ggplot(aes(x=val, y=n)) +
    theme_bw() +
    geom_step() +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    labs(x=what, y="Frequency")
}

