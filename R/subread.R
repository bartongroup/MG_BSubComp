

read_and_process_subread <- function(config, meta, gene_info, min_count = 10) {
  subread_path <- file.path(config$data_path, config$dir$subread)

  parse_subread_counts(subread_path, meta) |> 
    add_subread_log(subread_path, meta) |> 
    add_gene_names(gene_info) |> 
    normalize_subread_counts(gene_info) |> 
    filter_star_min_count(min_count)
}


parse_one_subread_count <- function(file, smpl) {
  if (!file.exists(file)) {
    warning(str_glue("File {file} not found."))
    return(NULL)
  }
  read_tsv(file, comment = "#", show_col_types = FALSE) |> 
    select(id = Geneid, count = last_col()) |> 
    add_column(raw_sample = smpl, .after = "id")
}


parse_subread_counts <- function(path, meta) {
  s2n <- set_names(meta$sample, meta$raw_sample)
  dat <- meta$raw_sample |> 
    map(~parse_one_subread_count(file.path(path, paste0(.x, ".txt")), .x)) |> 
    list_rbind() |> 
    mutate(sample = as.character(s2n[raw_sample])) |> 
    group_by(id) |>
    mutate(gene_count = sum(count)) |>
    ungroup() |>
    filter(gene_count > 0) |>
    select(id, sample, count) |> 
    mutate(sample = factor(sample, levels = meta$sample))
  tab <- dat |> 
    dat2mat(value_col = "count", id_col = "id", name_col = "sample")
  list(dat = dat, tab = tab, metadata = meta, sel = rownames(tab))
}


parse_one_subread_log <- function(file, smpl) {
  if (!file.exists(file)) {
    warning(paste("File", file, "not found"))
    return(NULL)
  }
  tb <- tibble(x = readLines(file)) |> 
    separate(x, c("key", "value"), sep = "\\t", fill = "right") |> 
    drop_na() |>
    filter(key != "Status") |> 
    mutate(
      value = as.integer(value),
      raw_sample = smpl
    )
  tot <- tibble(
    key = "Total",
    value = sum(tb$value),
    raw_sample = smpl
  )
  bind_rows(tb, tot)
}


parse_subread_logs <- function(path, meta) {
  s2n <- set_names(meta$sample, meta$raw_sample)
  meta$raw_sample |> 
    map(~parse_one_subread_log(file.path(path, paste0(.x, ".txt.summary")), .x)) |> 
    list_rbind() |> 
    mutate(sample = as.character(s2n[raw_sample]))
}


add_subread_log <- function(set, path, meta) {
  set$subread_log <- parse_subread_logs(path, meta)
  return(set)
}



normalize_subread_counts <- function(set, gene_info) {
  libsize <- set$subread_log |> 
    filter(key == "Total") |> 
    mutate(size = as.numeric(value)) |> 
    select(sample, size)
  
  set |> 
    normalise_to_library(gene_info, libsize) |> 
    #normalise_edger() |> 
    regularised_log()
}



plot_subread_count <- function(set, nrow = 1) {
  env <- new.env(parent = globalenv())
  env$nrow <- nrow

  counts <- set$subread_log |> 
    filter(key %in% c("Total", "Assigned")) |> 
    mutate(key = recode(key, "Total" = "Input reads", "Assigned" = "Uniquely mapped reads")) |> 
    select(sample, value, key) |> 
    mutate(value = as.numeric(value)) |> 
    mutate(value = value / 1e6) |> 
    mutate(key = fct_relevel(key, c("Input reads", "Uniquely mapped reads")))
  
  sample_input <- counts |> 
    filter(key == "Input reads") |>
    rename(input = value) |> 
    arrange(input) |> 
    mutate(sample = sample |> as_factor()) |> 
    select(-key)
  
  perc <- counts |> 
    left_join(sample_input, by = "sample") |> 
    mutate(value = 100 * value / input) |> 
    select(-input)
  
  rank_by <- function(w, sel) {
    r <- w |> 
      filter(key == sel) |> 
      mutate(rank = rank(value)) |> 
      select(sample, rank)
    w |> 
      left_join(r, by = "sample")
  }
  
  dat <- bind_rows(
      counts |>
        rank_by("Input reads") |> 
        add_column(what = "Count"),
      perc |>
        rank_by("Uniquely mapped reads") |> 
        add_column(what = "Percentage") |> 
        mutate(rank = rank + 1000)
    )
  
  labs <- dat |> 
    select(sample, rank) |> 
    distinct() |> 
    mutate(rank = as.character(rank))
  
  env$dat <- dat
  env$labs <- labs
  with(env, {
    get_labels <- function(rnk) {
      tibble(rank = rnk) |> 
        left_join(labs, by = "rank") |> 
        pull(sample)
    }
    dat |> 
      ggplot(aes(x = as_factor(rank), y = value, colour = key)) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      ) +
      geom_segment(aes(xend = as_factor(rank), yend = 0), colour = "grey90") +
      geom_point() +
      scale_colour_manual(values = okabe_ito_palette) +
      scale_x_discrete(labels = get_labels) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
      facet_wrap(~ what, scales = "free", nrow = nrow) +
      labs(x = NULL, y = "Read count (millions)", colour = "Legend")
  })
}
