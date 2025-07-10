one_count_file <- function(config, meta, n) {
  path <- file.path(config$data_path, config$dir$read_count)
  rs <- meta |> 
    slice(n) |> 
    pull(raw_sample) 
  file.path(path, str_glue("{rs}.txt"))
}

#' Read and Process STAR Output
#'
#' Reads STAR output files, processes them by adding gene names, normalizing counts,
#' filtering based on minimum count, and applying additional processing such as
#' log transformations.
#'
#' @param config A list containing configuration options, including paths.
#' @param meta Metadata data frame for the samples.
#' @param gene_info Data frame containing gene information.
#' @param suffix (Optional) File suffix for read count files, defaults to ".txt".
#' @param min_count (Optional) Minimum count threshold for filtering, defaults to 10.
#' @param fix_names_fun (Optional) A function to fix gene names, defaults to NULL.
#' @return A processed data set after applying all transformations and filters.
read_and_process_star <- function(config, meta, gene_info, suffix = ".txt", min_count = 10,
                                  fix_names_fun = NULL) {
  readcount_path <- file.path(config$data_path, config$dir$read_count)
  starmap_path <- file.path(config$data_path, config$dir$mapping)
  
  star_column <- star_strand_column(readcount_path, meta, suffix)
  strand <- names(star_column)
  msg <- str_glue("Detected strand - {strand}. Reading STAR column {star_column}.")
  message(msg)
  
  parse_star_counts(readcount_path, meta, star_column, suffix, fix_names_fun) |> 
    add_star_log(starmap_path, meta) |> 
    add_gene_names(gene_info) |> 
    normalize_star_counts(gene_info) |> 
    filter_star_min_count(min_count)
}

#' Add STAR Log Information
#'
#' Adds log information from STAR alignment to the dataset.
#'
#' @param set A dataset to which the STAR log information will be added.
#' @param path Path to the directory containing STAR log files.
#' @param meta Metadata dataframe for the samples.
#' @return The dataset with added STAR log information.
add_star_log <- function(set, path, meta) {
  set$star_log <- parse_star_logs(path, meta)
  return(set)
}

#' Parse a Single STAR Log File
#'
#' Reads and parses a single STAR log file, extracting key-value pairs of log
#' information.
#'
#' @param file Path to a STAR log file.
#' @param smpl Sample name associated with the log file.
#' @return A tibble with the parsed log information, or NULL if the file does
#'   not exist.
parse_one_star_log <- function(file, smpl) {
  if (!file.exists(file)) {
    warning(paste("File", file, "not found"))
    return(NULL)
  }
  tibble(x = readLines(file)) |> 
    separate(x, c("key", "value"), sep = "\\t", fill = "right") |> 
    drop_na() |> 
    mutate(key = str_remove(key, " \\|")) |> 
    mutate(key = str_remove(key, "^\\s+")) |> 
    mutate(raw_sample = smpl)
}

#' Parse All STAR Log Files
#'
#' Reads and parses all STAR log files specified in the metadata, combining them
#' into a single dataframe.
#'
#' @param path Path to the directory containing all STAR log files.
#' @param meta Metadata data frame for the samples.
#' @return A data frame with all parsed STAR log information.
parse_star_logs <- function(path, meta) {
  s2n <- set_names(meta$sample, meta$raw_sample)
  meta$raw_sample |> 
    map(~parse_one_star_log(file.path(path, paste0(.x, "_Log.final.out")), .x)) |> 
    list_rbind() |> 
    mutate(sample = as.character(s2n[raw_sample]))
}

# tabulate star log for printing
tabulate_star_log <- function(slog) {
  slog |> 
    select(-raw_sample) |> 
    pivot_wider(names_from = "sample", values_from = "value") |> 
    rename(Description = key) 
}

#' Parse a Single STAR Count File
#'
#' Reads a single STAR count file, extracting gene counts for a specific sample.
#'
#' @param file Path to a STAR count file.
#' @param smpl Sample name associated with the count file.
#' @param column (Optional) The column number to extract counts from, defaults to 2.
#' @param fix_names_fun (Optional) A function to fix gene names, defaults to NULL.
#' @return A data frame with gene counts for the specified sample.
parse_one_star_count <- function(file, smpl, column = 2, fix_names_fun = NULL) {
  if (!file.exists(file)) {
    warning(str_glue("File {file} not found."))
    return(NULL)
  }
  d <- read_tsv(file, col_names = FALSE, skip = 4, col_types = "ciii") |> 
    select(1, all_of(column)) |> 
    set_names("id", "count") |> 
    add_column(raw_sample = smpl, .after = "id")
  if (!is.null(fix_names_fun)) {
    d <- d |> mutate(id = fix_names_fun(id))
  }
  d
}

#' Parse All STAR Count Files
#'
#' Parses all STAR count files specified in the metadata, combining them into a
#' comprehensive dataset that includes wide matrix and long tibble formats,
#' along with sample metadata and selected gene identifiers.
#'
#' @param path Path to the directory containing STAR count files.
#' @param meta Metadata data frame for the samples, including sample names and
#'   raw sample identifiers.
#' @param column (Optional) Column number to extract counts from, defaults to 2.
#' @param suffix (Optional) Suffix for the count files, defaults to ".txt".
#' @param fix_names_fun (Optional) Function to correct gene names, defaults to
#'   NULL.
#' @return A list containing the data in both tab (wide matrix) and dat (long
#'   tibble) formats, along with metadata and selected row identifiers.
parse_star_counts <- function(path, meta, column = 2, suffix = ".txt", fix_names_fun) {
  s2n <- set_names(meta$sample, meta$raw_sample)
  dat <- meta$raw_sample |> 
    map(~parse_one_star_count(file.path(path, paste0(.x, suffix)), .x, column, fix_names_fun)) |> 
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

#' Add Gene Names to Dataset
#'
#' Adds human-readable gene names to the dataset based on gene identifiers.
#'
#' @param set A dataset to which gene names will be added.
#' @param gene_info Data frame containing gene identifiers and their
#'   corresponding symbols.
#' @return The dataset with added gene symbols.
add_gene_names <- function(set, gene_info) {
  set$genes <- tibble(id = rownames(set$tab)) |> 
    left_join(gene_info |> select(id, gene_symbol, description) |> distinct(), by = "id") |> 
    mutate(gene_symbol = if_else(is.na(gene_symbol), id, gene_symbol))
  set
}

#' Normalize Counts to Library Size
#'
#' Calculates normalization factors based on library sizes and applies them to
#' count data.
#'
#' @param set A dataset containing count data.
#' @param libsize (Optional) A pre-computed library size vector, defaults to
#'   NULL, in which case the library size is calculated from the data.
#' @return A list containing normalization factors and size factors for further
#'   adjustments.
#' @examples
normalise_to_size <- function(set, libsize = NULL) {
  if (is.null(libsize)) {
    libsize <- set$dat |> 
      group_by(sample) |> 
      summarise(size = sum(count, na.rm = TRUE)) 
  }
  
  libsize <- libsize |> mutate(normfac = size / mean(size))
  
  # look-up tables are much faster than left_join
  list(
    norm_fac = set_names(libsize$normfac, libsize$sample),
    size_fac = set_names(libsize$size / 1e6, libsize$sample),
    normfac_tab = libsize
  )
  
}

#' Normalize Dataset to Library Size and RPKM
#'
#' Applies library size normalization and calculates RPKM (Reads Per Kilobase
#' Million) for the dataset.
#'
#' @param set A dataset containing count data and gene information.
#' @param gene_info Data frame containing gene length information.
#' @param input_size A data frame or vector specifying the input library sizes
#'   for normalization.
#' @return The dataset with normalized count data and RPKM values.
normalise_to_library <- function(set, gene_info, input_size) {
  
  # len_fac = set_names(gene_info$length / 1e3, gene_info$id)
  
  # normalised to total mapped and counted reads
  mapped <- normalise_to_size(set)
  # normalised to total input reads
  input <- normalise_to_size(set, input_size)
  
  
  set$dat <- set$dat |> 
    mutate(
      count_norm = count / mapped$norm_fac[sample] |> unname()
      #rpkm = (count + 1) / (mapped$size_fac[sample] * len_fac[id] |> unname())
      #count_inputnorm = count / input$norm_fac[sample] |> unname,
      #rpkm_inputnorm = (count + 1) / (input$size_fac[sample] * len_fac[id] |> unname),
    )
  
  set$mapped_normfac <- mapped$normfac_tab
  set$input_normfac <- input$normfac_tab
  
  return(set)
}

#' Normalize Counts Using edgeR
#'
#' Applies TMM (Trimmed Mean of M-values) normalization from the edgeR package
#' to the count data.
#'
#' @param set A dataset containing count data.
#' @return The dataset with TMM-normalized count data.
normalise_edger <- function(set) {
  ed <- edgeR::DGEList(set$tab) |> 
    edgeR::calcNormFactors() |> 
    pluck("samples") |> 
    rownames_to_column("sample") |> 
    rename(normfac = norm.factors) |> 
    select(sample, normfac) |> 
    as_tibble()
  
  # look-up tables are much faster than left_join
  norm_fac <- set_names(ed$normfac, ed$sample)
  
  set$dat <- set$dat |> 
    mutate(count_tmm = count / norm_fac[sample] |> unname())
  set$edger_normfac <- ed
  
  return(set)
}

#' Apply Regularized Logarithm Transformation
#'
#' Transforms count data using the regularized log transformation from the
#' DESeq2 package, converting to log10 scale.
#'
#' @param set A dataset containing count data.
#' @return The dataset with count data transformed using the regularized log
#'   transformation.
regularised_log <- function(set) {
  rdat <- DESeq2::rlog(set$tab) |> 
    as_tibble(rownames = "id") |>
    pivot_longer(-id, names_to = "sample", values_to = "rlog") |> 
    mutate(
      rlog = rlog / log2(10),
      sample = factor(sample, levels = set$metadata$sample)
    )
  set$dat <- set$dat |> 
    left_join(rdat, by = c("id", "sample"))
 
  return(set) 
}


#' Normalize STAR Counts and Apply Transformations
#'
#' Normalizes STAR counts, applies library size and RPKM normalizations, and
#' adds regularized logarithm transformations. This function wraps several
#' normalization steps into one for convenience.
#'
#' @param set A dataset containing count data and STAR log information.
#' @param gene_info Data frame containing gene information, including lengths for
#'   RPKM calculation.
#' @return The dataset after applying all normalization and transformation
#'   steps.
normalize_star_counts <- function(set, gene_info) {
  libsize <- set$star_log |> 
    filter(key == "Number of input reads") |> 
    mutate(size = as.numeric(value)) |> 
    select(sample, size)
  
  set |> 
    normalise_to_library(gene_info, libsize) |> 
    #normalise_edger() |> 
    regularised_log()
}

# Normalise to one condition assuming replicates are matched
normalize_to_condition <- function(set, ref_cond, min_count = 10) {
  # Select genes where each replicate in reference is > min_count
  sel_genes <- set$dat |> 
    left_join(set$metadata, by = "sample") |> 
    filter(condition == ref_cond) |> 
    group_by(id) |> 
    summarise(m = min(count)) |> 
    filter(m >= min_count) |> 
    pull(id)
  
  d <- set$dat |> 
    filter(id %in% sel_genes)
  
  nrm <- set$metadata |> 
    filter(condition == ref_cond) |> 
    select(sample, replicate) |> 
    left_join(d, by = "sample") |> 
    select(id, replicate, ref = count_norm)
  
  lr <- d |> 
    left_join(set$metadata |> select(sample, replicate, condition), by = "sample") |> 
    filter(condition != ref_cond) |> 
    droplevels() |> 
    left_join(nrm, by = c("id", "replicate")) |> 
    mutate(logratio = log2((count_norm + 0.5) / ref)) |> 
    select(id, sample, logratio) |> 
    set_names("id", "sample", str_c("logratio_", ref_cond))
  
  set$metadata <- set$metadata |> 
    filter(condition != ref_cond) |> 
    droplevels()
  set$dat <- d |> 
    right_join(lr, by = c("id", "sample")) |> 
    droplevels()
  set$sel <- set$dat$id |> unique()
  
  set
}

#' Plot STAR Log Statistics
#'
#' Generates plots for selected statistics from STAR log files, such as the
#' number of input reads and percentages of uniquely mapped reads, across
#' different samples and conditions or groups.
#'
#' @param set A dataset containing count data and STAR log information.
#' @param group_var (Optional) The name of the variable in `meta` that defines
#'   the grouping of samples, defaults to "condition".
#' @param descs (Optional) A character vector specifying which statistics from
#'   the STAR log to plot. Commonly includes "Number of input reads", "Uniquely
#'   mapped reads %", and other mapping statistics. Defaults to a pre-defined
#'   list of important metrics.
#' @return A ggplot object visualizing the selected statistics from the STAR log
#'   files, separated by the specified group variable and for each sample.
plot_star_log <- function(set, group_var = "condition",
                        descs = c(
                          "Number of input reads",
                          "Uniquely mapped reads %",
                          "% of reads mapped to multiple loci",
                          "% of reads mapped to too many loci",
                          "% of reads unmapped: too short")
                        ) {
  env <- new.env(parent = globalenv())
  env$meta <- set$metadata
  
  env$s <- set$star_log |> 
    filter(key %in% descs) |>
    mutate(value = str_remove(value, "%") |> as.numeric()) |> 
    mutate(key = factor(key, levels = descs)) |> 
    select(-raw_sample) |> 
    left_join(set$metadata, by = "sample") |> 
    mutate(group = get(group_var)) |> 
    mutate(value = if_else(key ==  "Number of input reads", value / 1e6, value)) |> 
    mutate(key = recode(key,  "Number of input reads" =  "Number of input reads (millions)",))
  env$su <- env$s |> 
    filter(key == "Uniquely mapped reads %") |> 
    arrange(value)
  
  with(env, {
    s |> 
      mutate(sample = factor(sample, levels = meta$sample)) |> 
      ggplot(aes(x = sample, y = value, fill = group)) +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position = "bottom"
      ) +
      #geom_segment(aes(xend = sample, yend = 0), colour = "grey70") +
      geom_col(colour = "grey30") +
      facet_wrap(~key, scales = "free_y", ncol = 1) +
      scale_y_continuous(expand = c(0, 0), labels = scientific_10) +
      scale_fill_manual(values = okabe_ito_palette) +
      geom_text(aes(y = value * 1.05, label = "")) + # blank geom to expand axis
      labs(x = NULL, y = "Value")
  })
}

plot_star_log_map <- function(set) {
  env <- new.env(parent = globalenv())
  
  env$d <- set$star_log |> 
    filter(str_detect(key, "%")) |> 
    mutate(value = str_remove(value, "%") |> as.numeric()) |> 
    mutate(value = na_if(value, 0)) |> 
    mutate(sample = factor(sample, levels = set$metadata$sample)) |> 
    mutate(Description = factor(key))
  with(env, {
    ggplot(d, aes(x = sample, y = Description, fill = value)) + 
      theme_bw() +
      geom_tile() +
      scale_fill_viridis_c(option = "cividis", trans = "log10", limits = c(0.1, 100), breaks = c(0.1,1,10,100), labels = c(0.1,1,10,100)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
      labs(x = NULL, y = NULL, fill = "Percentage")
  })
}


plot_star_sense <- function(file) {
  read_tsv(file, skip = 4, col_names = c("gene", "Unstranded", "First", "Second"), col_types = cols()) |> 
    pivot_longer(-gene, names_to = "column", values_to = "count") |>
    filter(count > 0) |>
  ggplot(aes(x = log10(count), fill = column, group = column)) +
    theme_bw() +
    geom_density(alpha = 0.3)
}

#' Summarize Column Counts from a STAR Output File
#'
#' Reads a STAR output file and summarizes counts for 'unstranded', 'first', and
#' 'second' columns, filtering out genes with fewer than 10 'unstranded' counts.
#'
#' @param file Path to the STAR output file.
#' @param smpl Sample name associated with the file.
#' @return A dataframe summarizing the total counts for 'unstranded', 'first',
#'   and 'second' columns, along with the sample name.
star_col_count <- function(file, smpl) {
  read_tsv(file, skip = 4, col_names = c("gene", "unstranded", "first", "second"), col_types = cols()) |> 
    filter(unstranded > 10) |> 
    summarise(unstranded = sum(unstranded), first = sum(first), second = sum(second)) |> 
    mutate(raw_sample = smpl)
}

#' Determine Strand Specificity from STAR Outputs
#'
#' Analyzes STAR output files to determine the strand specificity of sequencing
#' based on count proportions and a predefined limit.
#'
#' @param path Path to the directory containing STAR output files.
#' @param meta Metadata data frame for the samples, including sample names and
#'   raw sample identifiers.
#' @param suffix (Optional) Suffix for the output files, defaults to ".txt".
#' @param prop.limit (Optional) Proportion limit to determine strand
#'   specificity, defaults to 0.8.
#' @return The determined strand specificity as "first", "second", or
#'   "unstranded". Will stop if multiple strand specificities are detected
#'   across samples.
star_strand <- function(path, meta, suffix = ".txt", prop.limit = 0.8) {
  s2n <- set_names(meta$sample, meta$raw_sample)
  m <- meta$raw_sample |> 
    map_dfr(~star_col_count(file.path(path, paste0(.x, suffix)), .x)) |> 
    mutate(sample = as.character(s2n[raw_sample])) |>
    mutate(r1 = first / unstranded, r2 = second / unstranded) |> 
    mutate(strand = if_else(r1 > prop.limit, "first", if_else(r2 > prop.limit, "second", "unstranded")))
  ms <- m |> distinct(strand)
  if (nrow(ms) == 1) {
    return(ms$strand)
  } else {
    print("Different stranding detected")
    print(m)
    stop()
  }
}

#' Wrapper Around `star_strand` to Return STAR Column Number
#'
#' Determines the strand specificity using `star_strand` and maps it to the
#' corresponding STAR output column number.
#'
#' @param ... Arguments passed on to the `star_strand` function.
#' @return Column number corresponding to the determined strand specificity.
star_strand_column <- function(...) {
  str2col <- set_names(
    c(3, 4, 2),
    c("first","second", "unstranded")
  )
  
  strand <- star_strand(...)
  str2col[strand]
}

#' Filter Genes Based on Minimum Count Threshold
#'
#' Filters genes based on a minimum count threshold in at least one sample,
#' within the specified count column.
#'
#' @param set A dataset containing count data.
#' @param min_count (Optional) Minimum count threshold, defaults to 10.
#' @param count_column (Optional) The count column to apply the threshold to,
#'   defaults to "count_norm".
#' @return The dataset filtered based on the minimum count threshold.
filter_star_min_count <- function(set, min_count = 10, count_column = "count_norm"){
  set$dat <- set$dat |> 
    group_by(id) |> 
    mutate(bad = max(get(count_column)) < min_count) |> 
    ungroup()
    
  set$sel <- set$dat |> 
    filter(!bad) |> 
    pull(id) |> 
    unique()

  set  
}

#' Filter Samples Based on Expression Criteria
#'
#' Filters samples from the dataset based on specified expression criteria.
#'
#' @param set A dataset containing count data and sample metadata.
#' @param expr An expression to evaluate for filtering samples.
#' @return The dataset filtered based on the expression criteria.
#' @examples
#' # Assuming `set` contains your dataset and `expr` is your criteria.
#' filtered_set <- filter_star_samples(set, "condition == 'treated'")
filter_star_samples <- function(set, expr) {
  meta <- set$metadata |> 
    filter(eval(rlang::parse_expr(expr))) |> 
    droplevels()
  smpl_sel <- as.character(meta$sample)
  
  set$dat <- set$dat |> filter(sample %in% smpl_sel)
  set$tab <- set$tab[, smpl_sel]
  set$metadata <- meta  
  
  set
}

#' Find Genes with Zero Counts in Any Group
#'
#' Identifies genes with zero counts in at least one group and summarizes counts
#' per group for each gene.
#'
#' @param set A dataset containing count data and metadata.
#' @param group_var (Optional) The variable name in metadata to define groups,
#'   defaults to "group".
#' @return A wide-format dataframe with gene IDs, and their sum of counts per
#'   group, highlighting genes with zero counts in any group.
find_zeroes <- function(set, group_var = "group") {
  set$dat |>
    left_join(set$metadata, by = c("sample", "raw_sample")) |>
    mutate(group = get(group_var)) |> 
    group_by(id, group) |>
    summarise(S = sum(count)) |>
    ungroup() |>
    group_by(id) |>
    mutate(m = min(S)) |>
    ungroup() |>
    filter(m == 0) |> 
    pivot_wider(id_cols = id, names_from = group, values_from = S)
}


merge_star_dat <- function(dat, meta, genes, columns = c("group", "replicate")) {
  dat |> 
    left_join(meta, by = "sample") |> 
    left_join(genes, by = "id") |> 
    mutate(gene_symbol = na_if(gene_symbol, "")) |> 
    select(id, gene_symbol, description, sample, all_of(columns), count, count_norm, rpkm)
}


plot_qualities <- function(qc) {
  qc |> 
    ggplot(aes(x = Base, y = Mean, group = root, colour = pair)) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_line(alpha = 0.3, linewidth = 0.6) +
    scale_colour_manual(values = okabe_ito_palette) +
    labs(x = "Base", y = "Mean quality score")
}

plot_cluster_qualities <- function(qc, text.size = 10) {
  tab <- qc |> 
    pivot_wider(id_cols = c(pair, Base), names_from = sample, values_from = Mean) |> 
    select(-c(pair, Base)) |> 
    as.matrix()
  
  hc <- t(tab) |> 
    dist() |> 
    hclust()
  
  dendr <- ggdendro::dendro_data(hc)
  seg <- ggdendro::segment(dendr)
  theme.d <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = text.size),
    axis.line.y = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(linewidth = 0.5),
    axis.ticks.y = ggplot2::element_blank()
  )
  ggplot() +
    theme.d +
    coord_flip() +
    geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_x_continuous(breaks = seq_along(dendr$labels$label), labels = dendr$labels$label) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(seg$y) * 1.03)) +
    scale_colour_manual(values = okabe_ito_palette) +
    labs(x = NULL, y = "Distance")
  
}


plot_map_qual <- function(qc, slog, base = 20) {
  qcf <- qc |> 
    filter(Base == base) |> 
    group_by(raw_sample) |> 
    summarise(qual = mean(Mean))
  slogf <- slog |> 
    filter(key == "Uniquely mapped reads %") |> 
    mutate(mapped = str_remove(value, "%") |> as.numeric())
  qcf |> 
    left_join(slogf, by = "raw_sample") |> 
  ggplot(aes(x = mapped, y = qual)) +
    theme_bw() +
    geom_point() +
    labs(y = paste("Mean read quality at", base), x = "Uniquely mapped reads")
}







despike <- function(set) {
  set$sel <- set$sel |> 
    str_subset("SPIKE", negate = TRUE)
  
  set
}

#' Plot Mapped Count and Percentages for Samples
#'
#' Visualizes the number of input reads, uniquely mapped reads, and reads
#' counted in genes for each sample in a dataset. It presents data in both
#' absolute counts (millions) and percentages to give a comprehensive view of
#' sequencing and mapping efficiency.
#'
#' @param set A dataset containing sequencing and mapping information. The
#'   dataset must have `$tab` for count data and `$star_log` for log
#'   information, including the number of input reads and uniquely mapped reads.
#' @return A ggplot object that displays the mapped count and percentages for
#'   samples. The visualization includes separate facets for absolute counts and
#'   percentages, facilitating easy comparison and interpretation.
plot_mapped_count <- function(set, nrow = 1) {
  env <- new.env(parent = globalenv())
  env$nrow <- nrow

  sm <- colSums(set$tab)
  cnt <- tibble(
    sample = names(sm),
    value = sm,
    key = "Reads counted in genes"
  )
  counts <- set$star_log |> 
    filter(key %in% c("Number of input reads", "Uniquely mapped reads number")) |> 
    mutate(key = recode(key, "Number of input reads" = "Input reads", "Uniquely mapped reads number" = "Uniquely mapped reads")) |> 
    select(sample, value, key) |> 
    mutate(value = as.numeric(value)) |> 
    bind_rows(cnt) |> 
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


select_strong_genes <- function(set, limit = 100) {
  set$dat |> 
    select(id, sample, count) |> 
    left_join(set$metadata, by = "sample") |> 
    group_by(id, group) |> 
    summarise(min_group_count = min(count)) |> 
    ungroup() |> 
    group_by(id) |> 
    summarise(best_count = max(min_group_count)) |> 
    filter(best_count > limit) |> 
    pull(id) |> 
    unique()
}


save_count_data <- function(set, clock_genes, file, what = "count_norm") {
  set$dat |> 
    mutate(val = get(what)) |> 
    left_join(set$genes, by = "id") |> 
    pivot_wider(id_cols = c(id, gene_symbol), names_from = sample, values_from = val) |> 
    mutate(across(where(is.numeric), ~signif(.x, 4))) |> 
    left_join(clock_genes |> select(-gene_symbol), by = "id") |> 
    relocate(class, .after = "gene_symbol") |> 
    mutate(class = as.character(class)) |> 
    mutate(class = replace_na(class, "-")) |> 
    rename(clock_class = class) |> 
    write_tsv(file)
}

dat2mat <- function(dat, id_col = "id", value_col = "abu_norm", name_col = "sample") {
  dat |>
    pivot_wider(id_cols = !!id_col, names_from = !!name_col, values_from = !!value_col) |>
    column_to_rownames(id_col) |>
    as.matrix()
}


star_mean_expr <- function(set, what = "count", id_var = "id", group_var = "condition") {
  set$dat |> 
    left_join(set$metadata, by = "sample") |> 
    mutate(val = get(what), group = get(group_var), id = get(id_var)) |> 
    group_by(id, group) |> 
    summarise(M = mean(val)) |> 
    ungroup() |> 
    rename(!!id_var := id, !!group_var := group) 
}


read_fastq_length <- function(path, meta) {
  map2_dfr(meta$raw_sample, meta$sample, function(rsam, sam) {
    cfile <- file.path(path, str_glue("{rsam}_count.txt"))
    if(file.exists(cfile)) {
      read_tsv(cfile, col_names = "read_count", show_col_types = FALSE) |> 
        mutate(sample = sam, raw_sample = rsam, .before = 1)
    }
  })
}

read_fastq_lengths <- function(config, meta) {
  trimmed_path <- file.path(config$data_path, config$dir$fastq_trimmed)
  f_trim <- read_fastq_length(trimmed_path, meta) |> 
    add_column(type = "Trimmed")
  clean_path <- file.path(config$data_path, config$dir$fastq_clean)
  f_clean <- read_fastq_length(clean_path, meta) |> 
    add_column(type = "Clean")
  bind_rows(f_trim, f_clean)
}



mat2dat <- function(tab, to_what = "abu_norm", names = "sample", id_var = "id") {
  tab |> 
    as.data.frame() |> 
    rownames_to_column(id_var) |> 
    pivot_longer(-{{ id_var }}, names_to = names, values_to = to_what)
}



#' Remove Batch Effects from Dataset
#'
#' Applies batch effect correction to the expression data in the dataset usingx
#' the limma method from the HarmonizR package. It's designed to adjust for
#' systematic differences between batches that might skew the analysis results.
#'
#' @param set A list containing the dataset to be adjusted, including `$dat` for
#'   expression data, `$metadata` for sample metadata, and optionally, other
#'   elements used in downstream analyses.
#' @param what The name of the column in `$dat` containing the expression values
#'   to be adjusted.
#' @param names The name of the column in `$metadata` representing the sample
#'   identifiers.
#' @param id_var The name of the column in `$dat` representing gene identifiers.
#' @param batch_var The name of the column in `$metadata` representing batch
#'   identifiers.
#' @param formula A formula specifying the model to adjust for batch effects,
#'   typically including the condition of interest and the batch.
#' @param filt A logical expression as a character string to filter the metadata
#'   for samples to be included in the adjustment process.
#' @return The modified dataset with batch effects removed from the specified
#'   expression column. The function adds the corrected expression values as a
#'   new column in `$dat` and updates `$metadata` to include only the samples
#'   that were included in the adjustment.
remove_batch_effects <- function(set, what = "count_norm", names = "sample", id_var = "id",
                                 batch_var = "day", formula = "~ condition + day",
                                 filt = "TRUE") {
  
  file_dat <- tempfile("data")
  file_desc <- tempfile("description")
  
  meta <- set$metadata |> 
    filter(!!rlang::parse_expr(filt)) |> 
    droplevels() |> 
    mutate(batch = get(batch_var)) |> 
    mutate(ibatch = as.integer(batch)) |> 
    arrange(ibatch)
  design_mat <- model.matrix(as.formula(formula), data = meta)
  
  dat <- set$dat |> 
    pivot_wider(id_cols = {{ id_var }}, names_from = !!names, values_from = !!what)
  dat <- dat[, c(id_var, as.character(meta$sample))]
  tab <- dat |> 
    column_to_rownames(id_var) |> 
    as.matrix()
  
  desc <- meta |> 
    rename(ID = sample) |> 
    mutate(sample = row_number(), batch = ibatch) |> 
    select(ID, sample, batch)
  
  write_tsv(dat, file_dat, na = "NaN")
  write_csv(desc, file_desc)
  
  # hr_combat <- HarmonizR::harmonizR(file_dat, file_desc, algorithm = "ComBat")
  hr_limma <- HarmonizR::harmonizR(file_dat, file_desc, algorithm = "limma")
  # bat <- limma::removeBatchEffect(tab, batch = meta[[batch_var]], design = design_mat)
  
  # dat_bat <- mat2dat(bat, "count_batch", names)
  # dat_combat <- mat2dat(hr_combat, "count_combat", names)
  dat_limma <- mat2dat(hr_limma, "count_limma", names)
  set$dat <- set$dat |>   
    # inner_join(dat_combat, by = c(id_var, names)) |> 
    inner_join(dat_limma, by = c(id_var, names)) 
    # inner_join(dat_bat, by = c(id_var, names))
  set$metadata <- meta
  
  set
}


write_counts <- function(set, file, what = "count_norm") {
  i2n <- set$genes |> 
    select(id, gene_symbol) |> 
    distinct()
  set$dat |>
    pivot_wider(id_cols = id, names_from = sample, values_from = !!what) |> 
    left_join(i2n, by = "id") |> 
    relocate(gene_symbol, .after = 1) |> 
    mutate(across(where(is.numeric), ~signif(.x, 4))) |> 
    write_csv(file)
}


group_counts_operons <- function(set, operons, min_count = 10) {
  dat <- set$genes |> 
    left_join(operons, by = "gene_symbol", relationship = "many-to-many") |> 
    filter(!is.na(operon_id)) |> 
    select(id, operon) |> 
    left_join(set$dat, by = "id", relationship = "many-to-many") |> 
    group_by(operon, sample) |>
    summarise(count = sum(count))  |> 
    ungroup() |> 
    rename(id = operon)
  genes <- dat |> 
    select(id) |> 
    distinct() |> 
    mutate(gene_symbol = id)

  list(
    metadata = set$metadata,
    genes = genes,
    dat = dat,
    tab = dat2mat(dat, id_col = "id", value_col = "count"),
    star_log = set$star_log
  ) |> 
    normalize_star_counts(genes) |> 
    filter_star_min_count(min_count)
}