okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#AAAAAA", "#000000")

th <- ggplot2::theme_bw() +
  ggplot2::theme(panel.grid = ggplot2::element_blank())

scientific_10 <- function(v, limit = 3.01) {
  sc <- map_chr(v, function(x) {
    if (is.na(x)) {
      as.character(x)
    } else if (x == 0 | abs(log10(x)) < limit) {
      format(x, scientific = FALSE)
    } else {
      format(x, scientific = TRUE) |> 
        str_replace("e", "%*%10^") |> 
        str_replace("\\+", "") |> 
        str_replace("\\^0+", "^")
    }
  })
  parse(text = sc)
}


gs <- function(gg, name, width, height) {
  ggplot2::ggsave(filename = file.path("fig", paste0(name, ".png")), plot = gg, device = "png",
         width = width, height = height, dpi = 300, limitsize = FALSE)
}

gp <- function(gg, name, width, height) {
  ggplot2::ggsave(filename = file.path("fig", paste0(name, ".pdf")), plot = gg, device = "pdf",
                  width = width, height = height)
}




plot_volma <- function(res, point_size, point_alpha) {
  r <- res |>
    select(x, y, sel, group)
  r_sel <- r |> filter(sel)
  r_nsel <- r |> filter(!sel)
  
  rm(res, r)  # Minimise environment for serialisation
  
  g <- ggplot() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_point(data = r_nsel, aes(x = x, y = y), colour = "grey70",
               size = point_size, alpha = point_alpha) +
    facet_wrap( ~ group)
  
  if (nrow(r_sel) > 0) {
    g <- g + geom_point(data = r_sel, aes(x = x, y = y), colour = "black",
                        size = point_size, alpha = point_alpha)
  }
  g
}

plot_ma <- function(res, a = "AveExpr", fc = "logFC", p = "PValue", fdr = "FDR", group = "contrast",
                    fdr_limit = 0.05, logfc_limit = 0, point_size = 0.5, point_alpha = 0.5) {
  res |>
    mutate(
      x = get(a),
      y = get(fc),
      group = get(group),
      sel = get(fdr) < fdr_limit & abs(get(fc)) >= logfc_limit
    ) |>
    plot_volma(point_size, point_alpha) +
    geom_hline(yintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(Mean~log[2]~CPM), y = expression(log[2]~FC))
}

plot_volcano <- function(res, fc = "logFC", p = "PValue", fdr = "FDR", group = "contrast",
                         fdr_limit = 0.05, logfc_limit = 0, point_size = 0.5, point_alpha = 0.5) {
  res |>
    mutate(
      x = get(fc),
      y = -log10(get(p)),
      group = get(group),
      sel = get(fdr) < fdr_limit & abs(get(fc)) >= logfc_limit
    ) |>
    plot_volma(point_size, point_alpha) +
    geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
}

plot_pdist <- function(res, p = "PValue", group = "contrast", n_bins = 50) {
  brks <- seq(0, 1, length.out = n_bins)
  r <- res |>
    mutate(p = get(p), grp = get(group)) |>
    select(p, grp)
  rm(res)
  ggplot(r, aes(x = p, y = after_stat(density))) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_histogram(breaks = brks) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    facet_wrap(~grp)
}

plot_vmp <- function(res, a = "logCPM", fc = "logFC", p = "PValue", fdr = "FDR", group = "contrast",
                     fdr_limit = 0.05, logfc_limit = 0, point_size = 0.5, point_alpha = 0.5, n_bins = 50){
  list(
    volcano = plot_volcano(res, fc, p, fdr, group, fdr_limit, logfc_limit, point_size, point_alpha),
    ma = plot_ma(res, a, fc, p, fdr, group, fdr_limit, logfc_limit, point_size, point_alpha),
    pdist = plot_pdist(res, p, group, n_bins)
  )
}


mn_plot_volcano <- function(de, group_ids = NULL, id_var = "id", group_var = "group", 
                            fdr_limit = 0.01, logfc_limit = 0, with_names = TRUE,
                            palette = okabe_ito_palette) {
  
  d <- de |> 
    mutate(id = get(id_var))
  
  if(is.null(group_ids)) {
    d <- d |> 
      mutate(
        sig = FDR < fdr_limit & abs(logFC) > logfc_limit,
        sel = sig,
        group = sig
      )
  } else {
    gi <- group_ids |> 
      transmute(id = get(id_var), group = get(group_var))
    d <- d |> 
      left_join(gi, by = "id") |> 
      mutate(sel = !is.na(group))
  }
  
  d_nsel <- d |> filter(!sel)
  d_sel <- d |> filter(sel)
  
  g <- ggplot(d, aes(x = logFC, y = -log10(PValue), label = gene_symbol)) +
    th +
    geom_point(data = d_nsel, colour = "grey80", size = 0.6, alpha = 0.3) +
    geom_point(data = d_sel, aes(colour = group), size = 1.3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
    facet_wrap(~ contrast) +
    #geom_vline(xintercept = c(-logfc_limit, logfc_limit), linetype = "dashed", linewidth = 0.2, colour = "darkgreen") +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P))
  
  if(is.null(group_ids)) {
    g <- g +
      theme(legend.position = "none") +
      scale_colour_manual(values = "black")
  } else {
    g <- g +
      scale_colour_manual(values = palette)
  }
  
  if(with_names) {
    g <- g + geom_text_repel(data = d_sel, size = 2.5, max.overlaps = 100) 
  }
  
  sigfilt_01 <- d |> filter(FDR < 0.01)
  sigfilt_05 <- d |> filter(FDR < 0.05)
  if(nrow(sigfilt_01) > 0) {
    siglim <- -log10(sigfilt_01$PValue) |> min()
    g <- g + geom_hline(yintercept = siglim, linetype = "dashed", linewidth = 0.2, colour = "darkgreen")
  }
  if(nrow(sigfilt_05) > 0) {
    siglim <- -log10(sigfilt_05$PValue) |> min()
    g <- g + geom_hline(yintercept = siglim, linetype = "dashed", linewidth = 0.2, colour = "darkgreen")
  }
  g
}

make_fc_fc <- function(de, ctr1, ctr2, fdr_limit = 0.01, sig_sel = c("any", "both", "x", "y")) {
  sig_sel = match.arg(sig_sel)

  fcfc <- de |> 
    filter(contrast %in% c(ctr1, ctr2)) |> 
    mutate(ctr = case_match(
      contrast,
      ctr1 ~ "x",
      ctr2 ~ "y"
    )) |> 
    pivot_wider(id_cols = id, names_from = ctr, values_from = c(logFC, FDR)) 

  if(sig_sel == "any") {
    fcfc <- fcfc |> mutate(sel = FDR_x < fdr_limit | FDR_y < fdr_limit)
  } else if(sig_sel == "both") {
    fcfc <- fcfc |> mutate(sel = FDR_x < fdr_limit & FDR_y < fdr_limit)
  } else if(sig_sel == "x") {
    fcfc <- fcfc |> mutate(sel = FDR_x < fdr_limit)
  } else if(sig_sel == "y") {
    fcfc <- fcfc |> mutate(sel = FDR_y < fdr_limit)
  } 
}

plot_fc_fc <- function(fcfc, x_lab, y_lab) {
  fcfc |> 
    mutate(
      group = "Fold change vs fold change"
    ) |> 
    select(x = logFC_x, y = logFC_y, group, sel) |> 
    plot_volma(point_size = 1, point_alpha = 1) +
    labs(x = x_lab, y = y_lab) +
    geom_hline(yintercept = 0, colour = "brown", alpha = 0.3) +
    geom_vline(xintercept = 0, colour = "brown", alpha = 0.3)
}


make_hist <- function(x, bins) {
  h <- hist(x, breaks = bins, plot = FALSE)
  tibble(
    x = h$mids,
    y = h$density
  )
}


plot_sample_quasirandom <- function(set, what = "count", colour_var = "treatment", log_scale = TRUE) {
  df <- set$dat |> 
    left_join(set$metadata, by = "sample") |> 
    mutate(val = get(what), colvar = get(colour_var)) 
  if(log_scale)
    df <- df |> filter(val > 0)
  dm <- df |> 
    group_by(sample) |> 
    summarise(M = median(val)) |> 
    mutate(x = as.integer(sample))
  w <- 0.3
  
  g <- df |> 
    ggplot() +
    th +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "bottom"
    ) +
    geom_quasirandom(aes(x = sample, y = val, colour = colvar), size = 0.4, varwidth = TRUE) +
    geom_segment(data = dm, aes(x = x - w, xend = x + w, y = M, yend = M)) +
    scale_color_manual(values = okabe_ito_palette) +
    labs(x = NULL, y = what, colour = colour_var)
  if(log_scale)
    g <- g + scale_y_log10(labels = scientific_10)
  g
}

plot_sample_distributions <- function(set, bins = 100, ncol = 15, what = "rlog", colour_var = "group",
                                      x_breaks = c(-5, 0, 5), text_size = 10, x_lim = c(-5, 5), log_scale = FALSE,
                                      scales = "free_y") {
  d <- set$dat |> 
    left_join(set$metadata, by = "sample") |> 
    mutate(val = get(what), colvar = get(colour_var)) |> 
    select(sample, val, colvar)
  if(log_scale) {
    d <- d |> 
      filter(val > 0) |> 
      mutate(val = log10(val))
    what <- str_glue("log {what}")    
  }
  d <- d |> 
    nest(data = val) |> 
    mutate(hist = map(data, ~make_hist(.x$val, bins = bins))) |> 
    select(-data) |> 
    unnest(hist)
  w <- d$x[2] - d$x[1]
  
  rm(set)
  
  ggplot(d, aes(x = x, y = y, fill = colvar)) +
    theme_bw() +
    theme(
      text = element_text(size = text_size),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.background = element_blank()
    ) +
    geom_col(width = w) +
    geom_hline(yintercept = 0, colour = "black") +
    #geom_vline(xintercept = 0, colour = "brown") +
    facet_wrap(~sample, scales = scales, ncol = ncol) +
    scale_fill_manual(values = okabe_ito_palette) +
    scale_x_continuous(breaks = x_breaks, limits = x_lim) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = what, y = "Density", fill = colour_var)
}

plot_kernels <- function(set, what = "abu_med", xlab = "Normalised abundance", log_scale = FALSE) {
  env <- new.env(parent = globalenv())
  env$what <- what
  env$xlab <- xlab

  d <- set$dat
  if(log_scale)
    d <- d |> mutate(!!what := log10(get(what)))

  env$d <- d |>
    mutate(val = get(what)) |>
    select(sample, val) |>
    nest(data = val) |>
    mutate(
      krn = map(data, ~density(.x$val, bw = "SJ") |> tidy())
    ) |>
    select(-c(data)) |>
    unnest(krn)

  with(env, {
    ggplot(d, aes(x = x, y = y, group = sample)) +
      theme_bw() +
      theme(
        panel.grid = element_blank()
      ) +
      geom_line(alpha = 0.1) +
      labs(x = xlab, y = "Density") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
  })
}


plot_clustering <- function(set, text_size = 10, what = "rlog", id_sel = NULL, names = "sample", dist.method = "euclidean",
                            clust.method = "complete", colour_var = "group", skip_selection = FALSE) {

  env <- new.env(parent = globalenv())
  env$text_size <- text_size

  #dat <- select_good_data(set, skip_selection)
  dat <- set$dat
  if(!is.null(id_sel)) {
    dat <- dat |> filter(id %in% id_sel)
  }

  tab <- dat2mat(dat, value_col = what, name_col = names)

  dendr <- t(tab) |>
    dist(method = dist.method) |>
    hclust(method = clust.method) |>
    dendsort::dendsort() |>
    ggdendro::dendro_data()

  env$seg <- ggdendro::segment(dendr)
  meta <- set$metadata |>
    mutate(
      colvar = get(colour_var),
      namevar = get(names) |> as.character()
    )
  env$labs <- left_join(dendr$labels |> mutate(label = as.character(label)), meta, by = c("label" = "namevar")) |>
    mutate(colour = okabe_ito_palette[as_factor(colvar)])
  env$theme.d <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = text_size, colour = env$labs$colour),
    axis.line.y = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(linewidth = 0.5),
    axis.ticks.y = ggplot2::element_blank()
  )

  with(env, {
    ggplot() +
      theme.d +
      coord_flip() +
      geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend)) +
      scale_x_continuous(breaks = seq_along(labs$label), labels = labs$label, expand = expansion(add = 1)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, max(seg$y) * 1.03)) +
      scale_colour_manual(values = okabe_ito_palette) +
      labs(x = NULL, y = "Distance")
  })
}


plot_clustering_circular <- function(set, what = "rlog", dist.method = "euclidean",
                                     clust.method = "complete", colour_var = "group", shape_var = "group") {
  env <- new.env(parent = globalenv())
  env$shape_var <- shape_var
  env$colour_var <- colour_var

  dat <- select_good_data(set) |>
    left_join(set$metadata, by = "sample")

  tab <- dat2mat(dat, value_col = what)
  gd <- which(rowSums(is.na(tab)) == 0)
  tab <- tab[gd, ]

  ph <- t(tab) |>
    dist(method = dist.method) |>
    hclust(method = clust.method) |>
    ape::as.phylo() |>
    ape::unroot()

  m <- set$metadata |>
    filter(sample %in% ph$tip.label)
  colours <- split(m$sample, m[[colour_var]])
  shapes <- split(m$sample, m[[shape_var]])

  phg <- ph |>
    ggtree::groupOTU(colours, group_name = "colour") |>
    ggtree::groupOTU(shapes, group_name = "shape")
  
  env$m <- m
  env$phg <- phg
  with(env, {
    ggtree::ggtree(phg, aes(colour = colour, shape = factor(shape, levels = levels(m[[shape_var]]))), layout = "daylight", branch.length = "branch.length") +
      theme(legend.position = "bottom") +
      ggtree::geom_tippoint() +
      scale_colour_manual(values = okabe_ito_palette) +
      scale_shape_manual(values = c(15:18, 0, 1, 2, 5, 6), na.value = 4, name = shape_var) +
      labs(colour = colour_var, shape = shape_var) +
      guides(colour = guide_legend(nrow = 4), shape = guide_legend(nrow = 4))
  })

}



plot_distance_matrix <- function(set, what = "rlog", text_size = 10, skip_selection = FALSE) {
  env <- new.env(parent = globalenv())
  env$text_size <- text_size

  #dat <- select_good_data(set, skip_selection)
  dat <- set$dat
  tab <- dat2mat(dat, value_col = what)

  env$d <- cor(tab, use = "complete.obs") |>
    as_tibble(rownames = "sample") |>
    pivot_longer(-sample) |>
    mutate(sample = factor(sample, levels = set$metadata$sample)) |>
    mutate(name = factor(name, levels = set$metadata$sample))

  with(env, {
    ggplot(d, aes(x = sample, y = name)) +
      geom_tile(aes(fill = value)) +
      scale_fill_viridis_c(option = "cividis") +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = text_size),
        axis.text.y = element_text(size = text_size),
        legend.position = "top"
      ) +
      labs(x = NULL, y = NULL, fill = "Correlation")
  })
}

plot_median_correlation <- function(set, what = "rlog", colour_var = "group", shape_var = "group",
                                  text_size = 10, skip_selection = FALSE) {
  env <- new.env(parent = globalenv())
  env$text_size <- text_size
  env$colour_var <- colour_var
  env$shape_var <- shape_var

  #dat <- select_good_data(set, skip_selection)
  dat <- set$dat
  tab <- dat2mat(dat, value_col = what)

  env$d <- cor(tab, use = "complete.obs") |>
    as_tibble(rownames = "sample") |>
    pivot_longer(-sample) |>
    filter(sample > name) |>
    group_by(sample) |>
    summarise(
      M = median(value),
      Q1 = quantile(value, 0.25),
      Q2 = quantile(value, 0.75)
    ) |>
    left_join(set$metadata, by = "sample") |> 
    mutate(
      colour = get(colour_var),
      shape = get(shape_var)
    )

  with(env, {
    d |> 
      ggplot(aes(x = fct_reorder(sample, M), y = M, colour = colour, shape = shape)) +
      th +
      geom_errorbar(aes(ymin = Q1, ymax = Q2), width = 0, colour = "grey70") +
      geom_point() +
      coord_flip() +
      scale_color_manual(values = okabe_ito_palette, name = colour_var) +
      scale_shape_manual(values = c(15:18, 0, 1, 2, 5, 6), na.value = 4, name = shape_var) +
      labs(x = NULL, y = "Mean correlation")
  })
}





plot_mean_var <- function(set, point_size = 0.5, what = "count_norm", group_var = "group") {
  env <- new.env(parent = globalenv())
  env$point_size <- point_size

  env$d <- set$dat |> 
    filter(!bad) |> 
    left_join(set$metadata, by = "sample") |> 
    mutate(val = get(what), group = get(group_var)) |> 
    group_by(id, group) |> 
    summarise(M = mean(val), V = var(val))
  with(env, {
    d |> 
      ggplot(aes(x = log10(M), y = log10(V))) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      ggpointdensity::geom_pointdensity(size = point_size) +
      scale_colour_viridis_c(option = "cividis") +
      geom_abline(slope = 1, intercept = 0, colour = "red", alpha = 0.3) +
      facet_wrap(~group) +
      labs(x = expression(log[10]~mean), y = expression(log[10]~variance))
  })
}


plot_xy <- function(dat, colour_var, shape_var, point_size = 1) {
  dat |> 
    mutate(colvar = get(colour_var), shapevar = get(shape_var)) |> 
    ggplot(aes(x = x, y = y, colour = colvar, shape = shapevar)) +
    theme_bw() +
    geom_point(size = point_size) +
    scale_color_manual(values = okabe_ito_palette, name = colour_var) +
    scale_shape_manual(values = c(15:18, 0, 1, 2, 5, 6), na.value = 4, name = shape_var) +
    labs(x = NULL, y = NULL) +
    theme(
      panel.grid = element_blank()
    ) +
    ggrepel::geom_text_repel(aes(label = sample), size = 3, colour = "grey40")
}

pca2xy <- function(pc, meta) {
  pc$x |> 
    as_tibble(rownames = "sample") |> 
    select(x = PC1, y = PC2, sample) |> 
    left_join(meta, by = "sample")
}

umap2xy <- function(um, meta) {
  colnames(um) <- c("x", "y")
  um |> 
    as_tibble() |> 
    mutate(sample = meta$sample) |> 
    left_join(meta, by = "sample")
}

plot_pca <- function(set, point_size = 2, what = "rlog", colour_var = "group", shape_var = "day") {
  tab <- dat2mat(set$dat, value_col = what)
  
  # remove rows with zero variance
  tab <- tab[apply(tab, 1, function(v) sum(is.na(v)) == 0), ]
  tab <- tab[apply(tab, 1, sd) > 0, ]
  
  pca <- prcomp(t(tab), scale. = TRUE, center = TRUE)
  var.perc <- 100 * (pca$sdev)^2 / sum((pca$sdev)^2)
  pca1 <- sprintf("PCA1 (%5.1f%%)", var.perc[1])
  pca2 <- sprintf("PCA2 (%5.1f%%)", var.perc[2])
  pca2xy(pca, set$metadata) |> 
    plot_xy(colour_var, shape_var, point_size) +
    labs(x = pca1, y = pca2)
}


plot_umap <- function(set, what = "rlog", point_size = 2, seed = 1,
                      n_neighbours = 15, min_dist = 0.01,
                      colour_var = "group", shape_var = "group") {
  tab <- dat2mat(set$dat, value_col = what)
  
  set.seed(seed)
  tab <- tab[apply(tab, 1, function(v) sum(is.na(v)) == 0), ]
  uwot::umap(t(tab), n_neighbors = n_neighbours, min_dist = min_dist) |> 
    umap2xy(set$metadata) |> 
    plot_xy(colour_var, shape_var, point_size)
}


plot_up_down <- function(res, groupvar = "contrast") {
  res |>
    filter(sig) |>
    mutate(direction = if_else(logFC > 0, "up", "down") |> factor(levels = c("up", "down"))) |>
    mutate(gr = get(groupvar)) |> 
    group_by(gr, direction) |>
    tally() |>
    mutate(n = if_else(direction == "down", -n, n)) |>
  ggplot(aes(x = gr, y = n, fill = direction)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    geom_col() +
    scale_fill_manual(values = okabe_ito_palette) +
    labs(x = NULL, y = "Num significant genes")
}

plot_pair_outliers <- function(set, pair, genes, amp = 2, tau = 3, label.size = 2.5, value = "count_norm", sel = NULL) {
  env <- new.env(parent = globalenv())
  env$pair <- pair
  env$label.size <- label.size
  env$sel <- sel
  
  efun <- function(x) {amp * exp(-x / tau)}
  genes <- genes |> 
    select(id, gene_symbol) |> 
    distinct()
  env$d <- set$dat |>
    filter(!bad) |> 
    filter(sample %in% pair) |>
    mutate(val = !!sym(value)) |> 
    mutate(val = log10(val + 0.5)) |> 
    pivot_wider(id_cols = id, names_from = sample, values_from = val) |> 
    rename(v1 = pair[1], v2 = pair[2]) |> 
    mutate(x = v1 + v2, y = v2 - v1) |> 
    mutate(f = efun(x)) |> 
    select(id, x, y, f)
  env$dsel <- env$d |> 
    filter(abs(y) > f) |> 
    left_join(genes, by = "id")

  with(env, {
    pl <- ggplot(d, aes(x = x, y = y)) +
      theme_bw() +
      theme(legend.position = "none") +
      geom_point(alpha = 0.3, colour = "black", size = 0.5) +
      geom_line(aes(y = f), colour = "grey90") +
      geom_line(aes(y = -f), colour = "grey90") +
      scale_fill_viridis_c(option = "cividis", trans = "sqrt") +
      labs(x = paste(pair[1], "+", y = pair[2]), y = paste(pair[2], "-", y = pair[1])) +
      geom_text_repel(data = dsel, aes(label = gene_symbol), size = label.size, max.overlaps = 30)
    if(!is.null(sel)) {
      dcol <- d |> 
        filter(id %in% sel)
      pl <- pl +
        geom_point(data = dcol, colour = "red", size = 1)
    }
    pl
  })  
}

plot_pair_outliers_rlog <- function(set, pair, gene_info, limit = 2, label.size = 2.5) {
  gene_info <- gene_info |> 
    select(id, gene_symbol) |> 
    distinct()
  d <- set$dat |>
    filter(sample %in% pair) |>
    pivot_wider(id_cols = id, names_from = sample, values_from = rlog) |> 
    rename(v1 = pair[1], v2 = pair[2]) |> 
    mutate(x = v1 + v2, y = v2 - v1) |> 
    select(id, x, y)
  dsel <- d |> 
    filter(abs(y) > limit) |> 
    left_join(gene_info, by = "id") 
  ggplot(d, aes(x = x, y = y)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_point(alpha = 0.3, colour = "black", size = 0.5) +
    geom_hline(yintercept = c(-limit, limit), colour = "grey90") +
    scale_fill_viridis_c(option = "cividis", trans = "sqrt") +
    labs(x = paste(pair[1], "+", y = pair[2]), y = paste(pair[2], "-", y = pair[1])) +
    geom_text_repel(data = dsel, aes(label = gene_symbol), size = label.size, max.overlaps = 30)
}


plot_biotype_de <- function(res, groupvar = "contrast") {
  baseline <- res |>
    group_by(gene_biotype) |>
    tally() |>
    mutate(base_prop = n / sum(n)) |>
    rename(base_n = n)
  
  d <- res |>
    filter(FDR < 0.05) |>
    mutate(gr = !!sym(groupvar)) |> 
    mutate(change = if_else(logFC>0, "up", "down") |> as_factor()) |>
    group_by(gr, change, gene_biotype) |>
    summarise(n = n()) |>
    left_join(baseline, by = "gene_biotype") |>
    mutate(rat = n / base_n, rat_se = sqrt(rat * (1 - rat) / base_n)) |>
    mutate(rat = if_else(change == "down", -rat, rat)) |>
    mutate(tp = str_remove(gr, "DMSO-ML")) |>
    arrange(gr) |>
    mutate(tp = as_factor(tp)) |>
    ungroup() |>
    group_by(gr, change) |>
    mutate(prop = n / sum(n), prop_se = sqrt(prop * (1 - prop) / sum(n))) |>
    mutate(prop = if_else(change == "down", -prop, prop))
  
  g0 <- baseline |>
    ggplot(aes(x = base_n, y = gene_biotype)) +
    theme_bw() +
    geom_col() +
    labs(x = "Count", y = NULL) +
    scale_x_continuous(expand = c(0,0), limits = c(0,max(baseline$base_n)*1.03))
  
  g1 <- d |>
    filter(base_n > 30) |>
    ggplot(aes(x = tp, y = rat)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), legend.position = "none") +
    geom_col(aes(fill = change)) +
    geom_errorbar(aes(ymin = rat-rat_se, ymax = rat+rat_se), width = 0.1) +
    scale_fill_manual(values = okabe_ito_palette) +
    facet_wrap(~ gene_biotype, scales = "fixed", ncol = 5) +
    labs(y = "Proportion of total", x = NULL)
  
  g2 <- d |>
    filter(gene_biotype != 'protein_coding') |>
    ggplot(aes(x = gene_biotype, y = prop)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), legend.position = "none") +
    geom_col(aes(fill = change)) +
    geom_errorbar(aes(ymin = prop-prop_se, ymax = prop+prop_se), width = 0.1) +
    scale_fill_manual(values = okabe_ito_palette) +
    facet_wrap(~ tp, ncol = 4) +
    labs(y = "Composition", x = NULL) +
    coord_flip()
  #scale_y_continuous(trans = "log10")
  
  list(baseline = g0, proportion = g1, composition = g2)
}


plot_gene_count <- function(set, gid, val = "count_norm", log.scale = FALSE, ncol = 2, symbol.size = 3, fill = "group", shape = "zero") {
  d <- set$dat |> 
    filter(id %in% gid) |> 
    mutate(value = get(val)) |> 
    mutate(zero = value == 0) |> 
    left_join(set$metadata, by = "sample") |> 
    left_join(set$genes, by = "id") |> 
    mutate(gene_symbol = if_else(is.na(gene_symbol), id, gene_symbol)) |> 
    mutate(sample = factor(sample, levels = set$metadata$sample))
  
  if(log.scale) {
    d$value <- log10(d$value)
    ylab <- expression(log[10]~count)
  } else {
    ylab <- "Count"
  }
  
  g <- ggplot(d, aes_string(x = "sample", y = "value", fill = fill, shape = shape)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      #legend.position = "none"
    ) +
    geom_point(size = symbol.size) +
    labs(x = NULL, y = ylab) +
    facet_wrap(~gene_symbol, scales = "free_y", ncol = ncol) +
    scale_shape_manual(values = c(21,23)) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    scale_fill_manual(values = okabe_ito_palette)
  if(!log.scale) g <-  g + scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA)) 
  g
}

plot_gene_groups <- function(set, gid, val = "count_norm", log_scale = FALSE, zero_scale = TRUE,
                             ncol = 4, symbol_size = 3, cex = 2, group_var = "group", w_mean = 0.3,
                             scales = "free_y") {
  env <- new.env(parent = globalenv())
  env$ncol <- ncol
  env$symbol_size <- symbol_size
  env$cex <- cex
  env$w_mean <- w_mean
  env$scales <- scales
  env$log_scale <- log_scale
  env$zero_scale <- zero_scale
  
  d <- set$dat |> 
    filter(id %in% gid & !bad) |>
    select(-bad) |> 
    left_join(set$metadata, by = "sample") |> 
    left_join(set$genes, by = "id") |> 
    filter(!bad) |> 
    mutate(value = get(val), group = get(group_var)) |> 
    mutate(zero = value == 0) |> 
    mutate(gene_symbol = if_else(is.na(gene_symbol), id, gene_symbol)) |> 
    mutate(sample = factor(sample, levels = set$metadata$sample))

  if (log_scale) {
    d$value <- log10(d$value + 0.25)
    ylab <- paste("log_10", val)
  } else {
    ylab <- val
  }
  
  dm <- d |> 
    group_by(gene_symbol, group) |> 
    summarise(M = mean(value))  |> 
    mutate(ix = group |> factor() |> as.integer())

  env$d <- d
  env$dm <- dm
  env$ylab <- ylab
  with(env, {
    g <- ggplot(d) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
        #legend.position = "none"
      ) +
      geom_beeswarm(aes(x = group, y = value, fill = as.factor(replicate), shape = zero), size = symbol_size, cex = cex) +
      geom_segment(data = dm, aes(x = ix - w_mean, xend = ix + w_mean, y = M, yend = M)) +
      labs(x = NULL, y = ylab, fill = "Replicate") +
      facet_wrap(~gene_symbol, labeller = label_wrap_gen(), scales = scales, ncol = ncol) +
      scale_shape_manual(values = c(21,23)) +
      guides(fill = guide_legend(override.aes = list(shape = 21))) +
      scale_fill_manual(values = okabe_ito_palette)
    if (!log_scale & zero_scale) g <-  g + scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA)) 
    g
  })
}


plot_term_genes <- function(set, terms, trm, de, ctr, n_top = 24, n_col = 4, sign = 1, val = "count_norm") {
  all_gids <- terms$term2feature[[trm]]
  gids <- de |> 
    filter(contrast == ctr & id %in% all_gids) |> 
    arrange(-sign * logFC) |> 
    pull(id) |> 
    unique() |> 
    head(n_top)
  plot_gene_groups(set, gids, val, ncol = n_col)
}

plot_fc_comparison <- function(res, groupvar = "contrast") {
  mx <- max(abs(res$logFC))
  d <- res |> 
    pivot_wider(id_cols = id, names_from = sym(groupvar), values_from = logFC)
    
  groups <- levels(res[[groupvar]])
  pairs <-  expand_grid(x = as_factor(groups), y = as_factor(groups)) |>
    filter(as.integer(x) < as.integer(y)) |> 
    mutate(across(everything(), as.character))
  map(1:nrow(pairs), function(i) {
    r <- pairs[i, ]
    d |> 
      rename(x = sym(r$x), y = sym(r$y)) |> 
    ggplot(aes(x = x, y = y)) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      geom_point(size = 0.3, alpha = 0.3) +
      geom_hline(yintercept = 0, alpha = 0.4) +
      geom_vline(xintercept = 0, alpha = 0.4) +
      labs(x = r$x, y = r$y) +
      xlim(-mx, mx) +
      ylim(-mx, mx)
  }) |> 
    plot_grid(plotlist = ., nrow = 1)
    
}

plot_fc_heatmap <- function(set, what = "rlog", min_n = 3, max_fc = 2,
                            id_sel = NULL, sample_sel = NULL, order_col = TRUE,
                            with_x_text = FALSE, with_y_text = FALSE) {
  d <- set$dat |> 
    mutate(val = get(what))
  if (!is.null(sample_sel))
    d <- d |> filter(sample %in% sample_sel)
  if (!is.null(id_sel))
    d <- d |> filter(id %in% id_sel)
  
  d <- d |> 
    group_by(id) |> 
    mutate(
      n = n(),
      fc = val - mean(val)
    ) |> 
    filter(n > min_n)
  tab <- dat2mat(d, value_col = "fc")
  
  smpls <- set$metadata |> 
    filter(sample %in% colnames(tab)) |> 
    arrange(group) |> 
    pull(sample)
  tab <- tab[, smpls]
  
  ggheatmap(tab, order.col = order_col, with.x.text = with_x_text, with.y.text = with_y_text, dendro.line.size = 0.2, max.fc = max_fc, legend.name = "logFC")
}


plot_de_dist <- function(set, de_genes) {
  set$dat |> 
    filter(id %in% de_genes) |> 
    select(id, sample, count) |> 
    left_join(set$metadata, by = c("sample")) |> 
    group_by(id, group) |> 
    summarise(M = mean(count)) |> 
  ggplot(aes(x = M)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_histogram(bins = 100) +
    scale_x_log10() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    labs(x = "Mean read count across replicates", y = "Frequency")
}


plot_group_profiles <- function(set, groups, ncol = 3) {
  set$dat |> 
    mutate(gr = groups[id]) |> 
    group_by(gr, id) |>
    mutate(M = mean(rlog), rm = rlog - M) |>
    ungroup() |> 
    group_by(gr, sample) |>
    summarise(m = mean(rm)) |>
    filter(gr > 0) |> 
    left_join(set$metadata, by = "sample") |> 
  ggplot(aes(x = sample, y = m, colour = group)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_segment(aes(xend = sample, yend = 0), colour = "grey80") +
    geom_point() +
    facet_wrap(~gr, ncol = ncol) +
    scale_colour_manual(values = okabe_ito_palette)
}



plot_clustergram <- function(enr, fdr_limit = 0.05) {
  # different format from fgsea
  if (!("ids") %in% names(enr)) {
    enr <- enr |> 
      rowwise() |> 
      mutate(ids = str_c(leading_edge, collapse = ", "))
  }
  
  d <- enr |> 
    select(term_name, ids, p_adjust) |> 
    filter(p_adjust < fdr_limit) |> 
    separate_rows(ids, sep = ",\\s+")
  
  did <- d |>
    group_by(ids) |> 
    tally() |> 
    arrange(desc(n))
  dtm <- d |> 
    group_by(term_name) |> 
    tally() |> 
    arrange(desc(n))
  
  divs <- seq(0.5, max(length(unique(d$term_name)), length(unique(d$ids))))
  
  d |> 
    mutate(
      x = factor(ids, levels = did$ids),
      y = factor(term_name, levels = dtm$term_name),
      z = factor(1)
    ) |> 
    ggplot(aes(x = x, y = y, fill = z)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      legend.position = "none"
    ) +
    geom_tile() +
    geom_vline(xintercept = divs, colour = "grey") +
    geom_hline(yintercept = divs, colour = "grey") +
    scale_fill_manual(values = "darkred") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL)
}



plot_network <- function(edges, clr, min_weight = 0.8) {
  edg <- edges |> 
    filter(colour == clr & weight > min_weight) 
  
  genes <- tibble(
    gene_symbol = c(edg$fromAltName, edg$toAltName) |> 
      unique()
  ) |> 
    mutate(id = seq_along(gene_symbol) - 1)
  gene2id <- set_names(genes$id, genes$gene_symbol)
  
  links <- edg |> 
    mutate(
      source = gene2id[fromAltName],
      target = gene2id[toAltName],
      value =  (weight - min(weight)) / (max(weight) - min(weight))
    ) |> 
    select(source, target, value)
  
  nodes <- genes |> 
    arrange(id) |> 
    mutate(group = 1, size = 1)
  
  networkD3::forceNetwork(
    links,
    nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "gene_symbol",
    Nodesize = "size",
    Group = "group",
    opacity = 0.9,
    fontSize = 10,
    opacityNoHover = 1,
    zoom = TRUE
  )
}

plot_volcano_term <- function(res, gso, ctr, termids, ncol = 3) {
  d <- res |> 
    filter(contrast == ctr) |> 
    mutate(
      x = logFC,
      y = -log10(PValue)
    ) |> 
    select(id, gene_symbol, x, y)
  
  termids <- intersect(termids, unique(gso$term_id))
  
  dat <- map_dfr(termids, function(termid) {
    gns <- gso |> 
      filter(contrast == ctr & term_id %in% termid) |> 
      pull(leading_edge) |> 
      unlist()
    termname <- gso |> 
      filter(contrast == ctr & term_id == termid) |> 
      pull(term_name)
    term <- str_glue("{termid}: {termname}")
    
    d |> 
      mutate(sel = id %in% gns) |> 
      add_column(term = term)
  })
  
  sel <- dat |> 
    filter(sel)
  
  
  ggplot(dat, aes(x = x, y = y)) +
    theme_bw() +
    theme(
      panel.grid = element_blank()
    ) +
    geom_hline(yintercept = 0, colour = "darkgreen") +
    geom_vline(xintercept = 0, colour = "darkgreen") +
    geom_point(colour = "grey80", size = 0.5) +
    geom_point(data = sel, colour = "black", size = 1) +
    geom_text_repel(data = sel, aes(label = gene_symbol)) +
    facet_wrap(~ term, labeller = label_wrap_gen(), ncol = ncol) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P))
}


plot_volcano_house <- function(res, house) {
  classes <- unique(house$class)
  map_dfr(classes, function(cl) {
    h <- filter(house, class == cl)
    res |> 
      mutate(
        contrast = cl,
        sig = gene_symbol %in% h$gene_symbol
      )
  }) |> 
    plot_volcano()
}


plot_condition_normalization <- function(set) {
  set$dat |> 
    left_join(set$metadata, by = "sample") |> 
    pivot_longer(starts_with("logratio_")) |> 
    mutate(name = str_remove(name, "logratio_")) |> 
    group_by(sample, condition, day, name) |> 
    summarise(m = mean(value)) |> 
  ggplot(aes(x = condition, y = m, colour = day)) +
    th +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_point() +
    geom_line(aes(group = day)) +
    facet_wrap(~name) +
    scale_colour_manual(values = okabe_ito_palette) +
    labs(x = NULL, y = "Mean log-ratio")
}

plot_de_heatmap <- function(set, sites, what = "value_med", max.scale = NULL) {
  if (nrow(sites) == 0) return(NULL)
  mat <- set$dat %>%
    mutate(val = get(what)) %>%
    right_join(sites, by = c("id", "multi")) %>%
    mutate(sample = factor(sample, levels = set$metadata$sample)) %>%
    group_by(id, multi) %>%
    mutate(M = mean(val), logfc = log2(val / M)) %>%
    pivot_wider(id_cols = c(id, multi), names_from = sample, values_from = logfc) %>%
    left_join(select(set$info, id, gene_name), by = "id") %>%
    unite("gid", gene_name, id, multi, sep = "-") %>%
    column_to_rownames("gid")
  ggheatmap(mat, legend.name = expression(log[2]~FC), with.y.text = TRUE, order.col = FALSE)
}