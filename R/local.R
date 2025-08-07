read_toxins_and_antibiotics <- function(fname, genes, sw_genes, max_diff = 50,
                                        manual_match = c("yfjE" = "B4U62_RS04555")) {
  g <- readxl::read_excel(fname, col_names = FALSE, skip = 1) |> 
    rename(x = 1) |> 
    filter(!(x %in% c("Chromosome", "Plasmid"))) |> 
    mutate(row_id = row_number()) |> 
    pivot_longer(-row_id) |> 
    drop_na() |> 
    pull(value) |> 
    unique()

  tb <- tibble(
    gene_symbol = g
  ) |> 
    distinct()

  # Add manual matches to genes
  genes[genes$id %in% manual_match, ]$gene_symbol <- names(manual_match)

  # match directly to NCBI genes
  tb_ncbi <- tb |> 
    left_join(select(genes, id, gene_symbol, start, end), by = join_by(gene_symbol))

  # missing from NCBI
  mis_ncbi <- tb_ncbi |> 
    filter(is.na(id)) |> 
    select(gene_symbol)

  # match SubtiWiki genes
  sw <- sw_genes |> 
    select(id, name, synonyms, genomic_annotations.start, genomic_annotations.end)
  tb_sw <- mis_ncbi |> 
    left_join(sw, by = c("gene_symbol" = "name"))

  # Missing from SubtiWiki first attemps
  mis_sw <- tb_sw |> 
    filter(is.na(id)) |> 
    select(gene_symbol)

  # Try matching synonyms
  tb_sw_syn <- mis_sw |> 
    left_join(sw, by = c("gene_symbol" = "synonyms"))
  # Will have to rename to synonyms
  name2synomym <- tb_sw_syn |>
    filter(!is.na(id)) |> 
    select(synonym = gene_symbol, name)

  # Match all together
  tb1 <- tb_ncbi |> 
    drop_na() |> 
    add_column(ncbi = TRUE)

  tb_sw_all <- bind_rows(tb_sw, tb_sw_syn) |> 
    filter(!is.na(id))
  sw_sel <- sw_genes |> 
    filter(id %in% tb_sw_all$id)
  # Rename name to synonym
  sw_sel[sw_sel$name %in% name2synomym$name, ]$name <- name2synomym$synonym

  # Match SubtiWiki with NCBI by coordinates
  sw_matched <- match_genes_by_coordinates(genes, sw_sel, max_diff) |> 
    filter(!is.na(sw_id)) |> 
    select(name, id, start = start.x, end = end.x)
  tb2 <- mis_ncbi |>
    left_join(sw_matched, by = c("gene_symbol" = "name")) |> 
    add_column(ncbi = FALSE)

  bind_rows(tb1, tb2) |> 
    left_join(select(genes, id, ncbi_gene = gene_symbol, chr), by = "id") |> 
    relocate(ncbi_gene, .after = gene_symbol) |> 
    mutate(gene_symbol = make.unique(gene_symbol))
}




match_genes_by_coordinates <- function(genes, sw, max_diff = 50) {
  g1 <- genes |> 
    filter(chr == "NZ_CP020102.1") |> 
    select(chr, start, end, id, gene_symbol)
  g2 <- sw |> 
    select(sw_id = id, name, synonyms, start = genomic_annotations.start, end = genomic_annotations.end) |> 
    mutate(
      s_min = start - max_diff,
      s_max = start + max_diff,
      e_min = end - max_diff,
      e_max = end + max_diff
    ) 
  full_join(g1, g2,
    by = join_by(
      between(start, s_min, s_max), 
      between(end, e_min, e_max)
    )
  ) |> 
    select(-c(s_min, s_max, e_min, e_max))
}


compare_ncbi_sw_genes <- function(genes, sw_genes, limit_line = 50) {
  d <- match_genes_by_coordinates(genes, sw_genes, max_diff = 10000) |>
    filter(gene_symbol == name) |> 
    mutate(
      diff_start = abs(start.y - start.x),
      diff_end = abs(end.y - end.x)
    )
  pl <- d |> 
    ggplot(aes(x = start.x / 1e6, y = diff_end)) +
    th +
    geom_point() +
    scale_y_log10() +
    geom_hline(yintercept = limit_line, colour = "brown") +
    labs(x = "Position in chromosome (Mbp)", y = "Difference in start posision (bp)")
  list(
    df = d,
    plot = pl
  )
}


match_genes <- function(gene_list, genes, sw_genes, max_diff = 50) {
  # match directly to NCBI genes
  tb <- tibble(
    gene_symbol = gene_list
  ) |> 
    left_join(select(genes, id, gene_symbol, start, end))
}


rename_toxan_genes <- function(df, toxan) {
  df |> 
    left_join(select(toxan, id, gene_symbol), by = "id") |>
    mutate(gene_symbol = dplyr::if_else(is.na(gene_symbol.y), gene_symbol.x, gene_symbol.y)) |> 
    select(-c(gene_symbol.x, gene_symbol.y))
}

rename_toxan_genes_in_set <- function(set, toxan) {
  set$genes <- rename_toxan_genes(set$genes, toxan)
  set
}


get_forward_reverse <- function(star, star_opposite) {
  full_join(star$dat, star_opposite$dat, by = join_by(id, sample)) |> 
    select(id, sample, fwd_count = count.x, rev_count = count.y) |> 
    mutate(across(c(fwd_count, rev_count), ~replace_na(.x, 0))) |> 
    left_join(select(star$metadata, sample, group), by = "sample") |> 
    left_join(select(star$genes, id, gene_symbol), by = "id")
}