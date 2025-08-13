write_shiny_file <- function(obj, method = c("qs", "rds"), shiny_dir = "shiny") {
  method <- match.arg(method)
  path <- file.path(shiny_dir, "data")
  if (!dir.exists(path)) dir.create(path)
  obj_name <- deparse(substitute(obj))
  cat(paste("  Writing", obj_name))
  if(method == "qs") {
    file_name <- file.path(path, str_glue("{obj_name}.qs"))
    qs2::qs_save(obj, file_name)
  } else if(method == "rds") {
    file_name <- file.path(path, str_glue("{obj_name}.rds"))
    write_rds(obj, file_name, compress = "xz")
  }
  cat("\n")
}

save_data_for_shiny <- function(set, de, fterms, gse, what = "count_norm", shiny_dir = "shiny") {

  de <- de |>
    select(id, log_fc = logFC, expr = logCPM, p_value = PValue, fdr = FDR, contrast)
  data <- set$dat |>
    #filter(!bad) |>
    mutate(value = get(what)) |>
    select(id, sample, value)
  metadata <- set$metadata |>
    filter(!bad) |>
    select(-bad) |> 
    select(sample, replicate, group)

  dexset <- dexdash::dexdash_set(de, data, metadata, "RNA-seq")

  features <- set$genes |>
    select(id, description, name = gene_symbol) |> 
    distinct()

  write_shiny_file(dexset, shiny_dir = shiny_dir)
  write_shiny_file(features, shiny_dir = shiny_dir)
  write_shiny_file(fterms, shiny_dir = shiny_dir)
  write_shiny_file(gse, shiny_dir = shiny_dir)
}

save_data_for_shiny_browser <- function(de, genes, operons, shiny_dir = "shiny") {
  genes <- genes |> 
    select(id, gene_symbol, description, chr, start, end, strand) |> 
    mutate(chr = if_else(chr == "NZ_CP020102.1", "chromosome", "plasmid"))
  operons <- operons |> 
    select(id, operon)
  de <- de |>
    select(id, log_fc = logFC, fdr = FDR)

  write_shiny_file(de, shiny_dir = shiny_dir)
  write_shiny_file(genes, shiny_dir = shiny_dir)
  write_shiny_file(operons, shiny_dir = shiny_dir)
}