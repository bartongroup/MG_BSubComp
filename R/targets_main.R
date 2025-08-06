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
    genome = Biostrings::readDNAStringSet(config$genome_file),
    gtf = rtracklayer::import(config$gtf_file),
    # Following https://pmc.ncbi.nlm.nih.gov/articles/PMC3754741,
    # plasmid GFF was downloaded from https://ncbi.nlm.nih.gov/nucleotide/KF365913
    gff_KF365913 = rtracklayer::import(config$KF365913_file),
    genes_gtf = get_gene_annotations_from_gtf(gtf),
    chromosomes = genes_gtf$chr |> unique() |> as.character(),
    plasmid_annot = ncbi_plasmid_annotations(genes_gtf, gff_KF365913),
    genes = annotate_plasmid_genes(genes_gtf, plasmid_annot),
    # Downloaded from https://subtiwiki.uni-goettingen.de/v5/data-export, using full export
    sw_genes = readr::read_csv(config$subti_wiki_file),
    sw_gene_match = match_genes_by_coordinates(genes, sw_genes, max_diff = 50),
    operons = get_operons(genes),
    terms = get_functional_terms(gtf, kg_spec = "bsu"),
    fterms = prepare_terms_fenr(terms, genes$id)
  )

  qc <- tar_plan(
    fastq_len = read_fastq_lengths(config, metadata),
    fig_ribo_prop = plot_ribo_proportion(fastq_len),

    tbl_strandedness = parse_strandedness(config, metadata),
    
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
    star_opposite = parse_star_counts(file.path(config$data_path, config$dir$read_count), metadata, column = 3),
    bg = read_bedgraphs(config, metadata),
    
    tab_star_log = star$star_log,
    fig_star_log = plot_star_log(star, group_var = "group"),
    fig_star_log_map = plot_star_log_map(star),
    fig_map_count = plot_mapped_count(star),
    
    example_count_file = one_count_file(config, metadata, 1),
    fig_star_sense = plot_star_sense(example_count_file),

    fig_bedgraphs = plot_bedgraphs_region(bg, "NZ_CP020102.1", 0, Inf, log_scale = TRUE)
  )
  
  diff_expr <- tar_plan(
    my_contrasts = c("mut-ctrl"),
    de = edger_de_contrasts(star, contrasts = my_contrasts),
    de_genes = de |> dplyr::filter(FDR < config$fdr_limit & abs(logFC) >= config$logfc_limit) |> select(id, gene_symbol),

    figs_de = plot_vmp(de, fdr_limit = config$fdr_limit, logfc_limit = config$logfc_limit),
    gse = fgsea_all_terms(de, fterms, rank_expr = "-logFC * log10(PValue)"),

    der = edger_de_contrasts(star, contrasts = my_contrasts, filt = 'sample != "mut_5"'),  
    figs_der = plot_vmp(der, fdr_limit = config$fdr_limit, logfc_limit = config$logfc_limit),
    gser = fgsea_all_terms(der, fterms, rank_expr = "-logFC * log10(PValue)"),

    top_sig = de |> dplyr::filter(id %in% de_genes$id) |> dplyr::arrange(PValue) |> head(50) |> dplyr::pull(id),
    fig_top_sig_heatmap = plot_fc_heatmap(star, id_sel = top_sig, max_fc = NA, with_x_text = TRUE, with_y_text = TRUE)
  )

  group_operons <- tar_plan(
    star_ops = group_counts_operons(star, operons),
    terms_ops = group_terms_operons(terms, operons),
    fterms_ops = prepare_terms_fenr(terms_ops, star_ops$genes$id),

    de_ops = edger_de_contrasts(star_ops, contrasts = my_contrasts),  
    de_operons = de_ops |> dplyr::filter(FDR < config$fdr_limit & abs(logFC) >= config$logfc_limit) |> select(id, gene_symbol),

    figs_de_ops = plot_vmp(de_ops, fdr_limit = config$fdr_limit, logfc_limit = config$logfc_limit),

    cmp_genes_ops = compare_de_genes_operons(de, de_ops, "mut-ctrl", fdr_limit = config$fdr_limit, logfc_limit = config$logfc_limit),

    gse_ops = fgsea_all_terms(de_ops, fterms_ops, rank_expr = "-logFC * log10(PValue)"),

    top_sig_ops = de_ops |> dplyr::filter(id %in% de_operons$id) |> dplyr::arrange(PValue) |> dplyr::pull(id),
    fig_top_sig_ops_heatmap = plot_fc_heatmap(star_ops, id_sel = top_sig_ops, max_fc = NA, with_x_text = TRUE, with_y_text = TRUE, text_size = 10)
  )

  for_report <- tar_plan(
    forward_reverse = get_forward_reverse(star, star_opposite),
    fig_fwd_rev_example = plot_forward_reverse(forward_reverse, c("S1449", "S1450", "B4U62_RS20160", "B4U62_RS20165", "B4U62_RS20155")),

    rrna_genes = genes |> dplyr::filter(biotype == "rRNA") |> dplyr::pull(id),
    fig_mut_1_5 = plot_pair_outliers(star, c("mut_1", "mut_5"), genes, amp = 2.5, tau = 4, sel = rrna_genes),
    fig_mut_1_2 = plot_pair_outliers(star, c("mut_1", "mut_2"), genes, amp = 2.5, tau = 4, sel = rrna_genes),

    tbl_gse = print_gse(gse, genes, fdr_limit = 0.05),
    gse_example = list(ontology = "kg", term_id = "bsu02010", contrast = "mut-ctrl"),
    gse_random = list(ontology = "kg", term_id = "bsu03010", contrast = "mut-ctrl"),
    gse_example_stats = make_gse_example(de, terms, gse, gse_example),
    gse_random_stats = make_gse_example(de, terms, gse, gse_random),
    
    fig_fg_example = plot_fgsea_enrichment(de, fterms, gse_example, rank_expr = "-logFC * log10(PValue)"),
    fig_fg_example_random = plot_fgsea_enrichment(de, fterms, gse_random, rank_expr = "-logFC * log10(PValue)"),

    dea = print_de(de, genes),
    dea_ops = print_de(de_ops, genes),
    sav_de = write_table(dea),
    sav_de_ops = write_table(dea_ops),

    gsea = print_gse(gse, star$genes, fdr_limit = 1),
    gsea_ops = print_gse(gse_ops, star_ops$genes, fdr_limit = 1),
    sav_gse = write_table(gsea),
    sav_gse_ops = write_table(gsea_ops),

    cmp_ncbi_sw = compare_ncbi_sw_genes(genes, sw_genes, limit_line = 50)
  )

  for_manuscript <- tar_plan(
    pdf_volcano_toxins = mn_plot_volcano_toxins(de_toxan, toxan) |> gp("volcano_toxins", 4, 4)
  )

  toxins_and_antibiotics <- tar_plan(
    tar_target(toxan_file, "info/toxins_and_antibiotics.xlsx", format = "file"),
    toxan = read_toxins_and_antibiotics(toxan_file, genes, sw_genes, max_diff = 50),
    de_toxan = rename_toxan_genes(de, toxan),
    de_toxan_sel = de_toxan |> filter(id %in% toxan$id),
    star_toxan = rename_toxan_genes_in_set(star, toxan),
    fig_volcano_toxan = mn_plot_volcano(de_toxan, group_ids = toxan |> add_column(group = "Toxins and antibiotics")),
    fig_heatmap_toxan = plot_fc_heatmap(star_toxan, id_sel = toxan$id, with_x_text = TRUE, with_y_text = TRUE, max_fc = NA)
  )

  shiny <- tar_plan(
    sav_shiny_genes = save_data_for_shiny(star, de, fterms, gse, shiny_dir = "shiny/de_genes"),
    sav_shiny_operons = save_data_for_shiny(star_ops, de_ops, fterms_ops, gse_ops, shiny_dir = "shiny/de_operons")
  )

  get_environment <- tar_plan(
    tar_target(version_file, file.path(config$data_path, "environment.txt"), format = "file"),
    software_versions = get_software_versions(
      version_file,
      selection = c("fastqc", "fastq-screen", "fastp", "ribodetector", "multiqc", "star", "subread", "rseqc", "samtools", "bedtools", "snakemake", "drmaa", "deeptools")
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
    toxins_and_antibiotics,
    for_report,
    for_manuscript,
    shiny
  )
}