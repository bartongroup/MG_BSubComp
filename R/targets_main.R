targets_main <- function() {

  setup_experiment = tar_plan(
    tar_target(config_file, "config.json", format = "file"),
    config = jsonlite::read_json(config_file),
    tar_target(yaml_file, config$yaml_file, format = "file"),
    yaml = yaml::read_yaml(yaml_file),
    tar_target(meta_file, config$metadata_file, format = "file"),
    metadata = make_metadata(meta_file)
  )

  get_annotations = tar_plan(
    gtf = rtracklayer::import(yaml$gtf_url),
    operons = get_operons(),
    genes = get_gene_annotations_from_gtf(gtf),
    chromosomes = genes$chr |> unique() |> as.character(),
    terms = get_functional_terms(gtf),
    fterms = prepare_terms_fenr(terms, genes$id)
  )

  qc <- tar_plan(
    fastq_len = read_fastq_lengths(config, metadata),
    fig_ribo_prop = plot_ribo_proportion(fastq_len),
    
    fscrn = parse_fscreens(config, metadata),
    qcs = parse_qcs(config, metadata, paired = TRUE),
    idxstats = parse_idxstats(config, metadata),
    
    fig_fscreen = plot_fscreen_map(fscrn),
    fig_read_qual = plot_qualities(qcs),
    fig_read_qual_clust = plot_cluster_qualities(qcs),
    
    png_sample_dist = plot_sample_quasirandom(star, colour_var = "group") |> gs("sample_quasirandom", 8, 4),
    fig_mean_var = plot_mean_var(star, group_var = "group"),
    fig_distance_mat = plot_distance_matrix(star),
    fig_clustering = plot_clustering(star, colour_var = "group"),
    fig_pca = plot_pca(star, colour_var = "group", shape_var = "group")
  )

   star <- tar_plan(
    star = read_and_process_star(config, metadata, genes, min_count = 10),
    bg = read_bedgraphs(config, metadata),
    
    tab_star_log = star$star_log,
    fig_star_log = plot_star_log(star, group_var = "group"),
    fig_star_log_map = plot_star_log_map(star),
    fig_map_count = plot_mapped_count(star),
    
    example_count_file = one_count_file(config, metadata, 1),
    fig_star_sense = plot_star_sense(example_count_file),

    fig_bedgraphs = plot_bedgraphs_region(bg, "NZ_CP020102.1", 0, Inf)
  )
  
  diff_expr <- tar_plan(
    my_contrasts = c("mut-ctrl"),
    de = edger_de_contrasts(star, contrasts = my_contrasts),  
    figs_de = plot_vmp(de, fdr_limit = config$fdr_limit, logfc_limit = config$logfc_limit),
    gse = fgsea_all_terms(de, fterms, rank_expr = "-logFC * log10(PValue)"),

    der = edger_de_contrasts(star, contrasts = my_contrasts, filt = 'sample != "mut_5"'),  
    figs_der = plot_vmp(der, fdr_limit = config$fdr_limit, logfc_limit = config$logfc_limit),
    gser = fgsea_all_terms(der, fterms, rank_expr = "-logFC * log10(PValue)")
  )

  group_operons <- tar_plan(
    ops = group_counts_operons(star, operons),

    deo = edger_de_contrasts(ops, contrasts = my_contrasts),  
    figs_deo = plot_vmp(deo, fdr_limit = config$fdr_limit, logfc_limit = config$logfc_limit),
  )

  for_report <- tar_plan(
    rrna_genes = genes |> dplyr::filter(biotype == "rRNA") |> dplyr::pull(id),
    fig_mut_1_5 = plot_pair_outliers(star, c("mut_1", "mut_5"), genes, amp = 2.5, tau = 4, sel = rrna_genes),
    fig_mut_1_2 = plot_pair_outliers(star, c("mut_1", "mut_2"), genes, amp = 2.5, tau = 4, sel = rrna_genes),

    de_genes = de |> dplyr::filter(FDR < config$fdr_limit & abs(logFC) >= config$logfc_limit) |> select(id, gene_symbol)
  )

  shiny <- tar_plan(
    sav_shiny = save_data_for_shiny(star, de, fterms, gse)
  )

  get_environment <- tar_plan(
    tar_target(version_file, file.path(config$data_path, "environment.txt"), format = "file"),
    software_versions = get_software_versions(
      version_file,
      selection = c("fastqc", "fastq-screen", "fastp", "ribodetector", "multiqc", "star", "samtools", "bedtools", "snakemake", "drmaa", "deeptools")
    )
  )
  
  c(
    get_annotations,
    setup_experiment,
    get_environment,
    qc,
    star,
    diff_expr,
    group_operons,
    for_report,
    shiny
  )
}