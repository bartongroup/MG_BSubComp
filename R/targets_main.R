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
    gff = rtracklayer::import(yaml$gff_url),
    genes = get_gene_annotations_from_gff(gff),
    terms = get_functional_terms(gff, kg_spec = "bsu"),
    fterms = prepare_terms_fenr(terms, genes$gene_symbol)
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
    
    png_sample_dist = plot_sample_quasirandom(star) |> gs("sample_quasirandom", 8, 4),
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
    fig_star_sense = plot_star_sense(example_count_file)
  )
  
  de <- tar_plan(
    my_contrasts = c("mut-ctlr"),
    de = edger_de_contrasts(star, contrasts = my_contrasts),
    
    figs_de = plot_vmp(de, fdr_limit = config$fdr_limit),
    
    gse = fgsea_all_terms(de, fterms)
  )

  get_environment <- tar_plan(
    tar_target(version_file, file.path(config$data_path, "environment.txt"), format = "file"),
    software_versions = get_software_versions(
      version_file,
      selection = c("fastqc", "fastq-screen", "fastp", "ribodetector", "multiqc", "star", "samtools", "bedtools", "snakemake", "drmaa")
    )
  )
  
  c(
    get_annotations,
    setup_experiment,
    get_environment
  )
}