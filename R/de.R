add_genes <- function(res, info) {
  g <- info |>
    select(id, gene_symbol) |> 
    distinct()
  res |>
    left_join(g, by = "id")
}

#' Differential Expression Analysis with Formula
#'
#' This function performs differential expression analysis using the edgeR
#' package with a specified formula.
#'
#' @param set A list containing the dataset, with `dat` and `metadata` elements.
#' @param formula A formula specifying the model for the analysis.
#' @param what A character string specifying the data type to analyze. Default
#'   is `"count_norm"`.
#' @param filt A character string representing a logical expression to filter
#'   the metadata. Default is `"TRUE"`.
#' @param id_col The name of the column with identifiers, default is `"id"`.
#' @param sample_var A character string specifying the sample names variable. Default
#'   is `"sample"`.
#'
#' @return A tibble with the differential expression results.
edger_de_formula <- function(set, formula, what = "count", filt = "TRUE", id_col = "id", sample_var = "sample") {
  meta <- set$metadata |>
    filter(!bad & !!rlang::parse_expr(filt)) |>
    droplevels()
  dat <- set$dat |> 
    filter(good)
  
  tab <- dat2mat(dat, value_col = what, id_col = id_col, name_col = sample_var)[, as.character(meta[[sample_var]])]
  design_mat <- model.matrix(as.formula(formula), data = meta)
  coefs <- colnames(design_mat)[-1]
  
  fit <- tab |>
    DGEList() |>
    normLibSizes() |>
    glmQLFit(design = design_mat)
  
  coefs |> 
    map(function(coef) {
      fit |> 
        glmQLFTest(coef = coef) |> 
        topTags(n = 1e16, sort.by = "none") |> 
        pluck("table") |> 
        as_tibble(rownames = "id") |>
        add_column(contrast = coef)
   }) |> 
    list_rbind() |> 
    select(-F) |>
    mutate(contrast = factor(contrast, levels = coefs)) |> 
    add_genes(set$genes)
}


#' Differential Expression Analysis for Selected Contrasts
#'
#' This function performs differential expression analysis using the limma
#' package for selected contrasts.
#'
#' @param set A list containing the dataset, with `dat` and `metadata` elements.
#' @param contrasts A character vector of contrasts to analyze. If `NULL`, all
#'   pairs of contrasts are analyzed. Default is `NULL`.
#' @param group_var A character string specifying the grouping variable in the
#'   metadata. Default is `"treatment"`.
#' @param what A character string specifying the data type to analyze. Default
#'   is `"abu_norm"`.
#' @param filt A character string representing a logical expression to filter
#'   the metadata. Default is `"TRUE"`.
#' @param names A character string specifying the sample names variable. Default
#'   is `"sample"`.
#'
#' @return A tibble with the differential expression results.
edger_de_contrasts <- function(set, contrasts = NULL, group_var = "group", what = "count", filt = "TRUE",
                     id_col = "id", sample_var = "sample") {
  meta <- set$metadata |>
    filter(!bad & !!rlang::parse_expr(filt)) |>
    mutate(group = get(group_var)) |>
    filter(!is.na(group)) |>
    droplevels()
  dat <- set$dat |> 
    filter(!bad)
  groups <- unique(as.character(meta$group)) |>
    janitor::make_clean_names()
  design_mat <- model.matrix(~ 0 + group, data = meta)
  colnames(design_mat) <- groups
  
  tab <- dat2mat(dat, value_col = what, id_col = id_col, name_col = sample_var)[, as.character(meta[[sample_var]])]

  if (is.null(contrasts)) {
    contrasts <- expand_grid(x = as_factor(groups), y = as_factor(groups)) |>
      filter(as.integer(x) < as.integer(y)) |>
      unite(contrast, c(y, x), sep = "-") |>
      pull(contrast)
  }
  contrast_mat <- limma::makeContrasts(contrasts = contrasts, levels = design_mat)
  
  fit <- tab |>
    DGEList(group = meta$group) |>
    normLibSizes() |>
    glmQLFit(design = design_mat)
  
  map(contrasts, function(ctr) {
    fit |> 
      glmQLFTest(contrast = contrast_mat[, ctr]) |>
      topTags(n = 1e16, adjust.method = "BH", sort.by = "none") |>
      pluck("table") |>
      as_tibble(rownames = "id") |>
      mutate(contrast = ctr)
  }) |> 
    list_rbind() |> 
    mutate(contrast = contrast |> str_replace(" - ", "-") |> str_remove_all("group") |> as_factor()) |> 
    select(-F) |>
    add_genes(set$genes)
}

select_top_genes <- function(de, ctr, direction = "up", n_genes = 6,
                             fdr_limit = 0.01, log_cpm_limit = 4) {
  d <- de |> 
    filter(contrast == ctr & FDR < fdr_limit & logCPM > log_cpm_limit)
  if(direction == "up") {
    d <- d |> arrange(-logFC)
  } else {
    d <- d |> arrange(logFC)
  }
  d |> 
    head(n_genes) |> 
    pull(id)
}


write_de <- function(de, genes, file) {
  genes <- genes |> 
    select(id, description)
  de |> 
    left_join(genes, by = "id") |> 
    mutate(across(where(is.numeric), ~signif(.x, 4))) |> 
    write_csv(file)
}


