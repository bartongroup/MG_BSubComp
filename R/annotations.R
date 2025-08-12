
get_gene_annotations_from_gtf <- function(gtf) {
  coords <- gtf |> 
    as_tibble() |> 
    select(
      chr = seqnames,
      start,
      end,
      strand
    )
  mc <- mcols(gtf) |> 
    as_tibble()

 # We use "gene" selection in STAR alignment, so need to follow here
 genes <- bind_cols(coords, mc) |> 
  filter(type == "gene")  |> 
  select(
    chr,
    start,
    end,
    strand,
    id = gene_id,
    gene_symbol = gene,
    biotype = gene_biotype
  )

  # only rows marked as CDS contain gene desscription
  CDS <-  bind_cols(coords, mc) |> 
  filter(type == "CDS")  |> 
  select(
    id = gene_id,
    description = product
  )

  left_join(genes, CDS, by = "id")
}

get_go_mapping_from_gtf <- function(gtf) {
  mapping <- mcols(gtf) |> 
    as_tibble() |> 
    filter(type == "CDS") |> 
    select(
      id = gene_id,
      term_id = Ontology_term
    ) |> 
    unnest(term_id) |> 
    distinct() |> 
    drop_na()
  terms <- fenr:::fetch_go_terms(use_cache = TRUE, on_error = "stop") |> 
    filter(term_id %in% mapping$term_id)

  list(
    mapping = mapping,
    terms = terms
  )
}


get_functional_terms <- function(gtf, kg_spec) {
  cat("Loading GO term data\n")
  go <- get_go_mapping_from_gtf(gtf)
  cat("Loading KEGG data\n")
  kg <- fenr::fetch_kegg(species = kg_spec)

  gns <- get_gene_annotations_from_gtf(gtf) |> 
    select(id, gene_symbol) |>
    distinct()

  kg$mapping <- kg$mapping |> 
    inner_join(gns, by = "gene_symbol", relationship = "many-to-many") |> 
    select(id, term_id) |> 
    distinct()

  terms = list(
    go = go,
    kg = kg
  )
}


prepare_terms_fenr <- function(terms, all_features) {
  ontologies <- names(terms)
  
  map(ontologies, function(ont) {
    trm <- terms[[ont]]
    fenr::prepare_for_enrichment(trm$terms, trm$mapping, all_features, feature_name = "id")
  }) |> 
    set_names(ontologies)
}


# Download gene annotations from SubtiWiki

get_subtiwiki_genes <- function() {
  # First, get all gene IDs
  url <- "https://subtiwiki.uni-goettingen.de/v5/api/gene/"
  req <- httr2::request(url)
  resp <- httr2::req_perform(req)
  js <- httr2::resp_body_json(resp)

  data <- js$data
  ids <- map(data, function(d) {
    tibble(id = d$id, name = d$name)
  }) |> 
    list_rbind()

  gns <- map(ids$id, function(id) {
    qry <- str_glue("{url}{id}?representation=default")
    req <- httr2::request(qry)
    resp <- httr2::req_perform(req)
    js <- httr2::resp_body_json(resp)
    
    js$data |> 
      unlist() |> 
      enframe() |>
      filter(
        name %in% c("synonyms", "id", "legacy_id", "name", "description", "function", "product", "essential") |
        str_detect(name, "^genomic_annotations")
      ) |> 
      pivot_wider()
  }, .progress = TRUE) |> 
    list_rbind() |> 
    mutate(across(c(genomic_annotations.start, genomic_annotations.end), as.integer))
}



# Download operon clustering from SubtiWiki
get_operons <- function(genes) {
  url <- "https://subtiwiki.uni-goettingen.de/v5/api/operon/"
  req <- httr2::request(url)
  resp <- httr2::req_perform(req)
  js <- httr2::resp_body_json(resp)

  data <- js$data
  ops <- map(data, function(d) {
    id <- d$id
    map(d$genes, function(gene) {
      tibble(gene_symbol = gene$name)
    }) |> 
      list_rbind() |> 
      mutate(operon_id = id)
  }) |> 
    list_rbind() 

  ops_names <- ops |> 
    group_by(operon_id) |> 
    summarise(operon = str_c(sort(gene_symbol), collapse = "-"))

  id2name <- genes |> 
    select(id, gene_symbol) |> 
    distinct() |> 
    drop_na()
  ops |> 
    left_join(ops_names, by = "operon_id") |> 
    left_join(id2name, by = "gene_symbol") |> 
    select(id, gene_symbol, operon)
}


group_terms_operons <- function(terms, operons) {
  ontologies <- names(terms)
  map(ontologies, function(ont) {
    trm <- terms[[ont]]

    mapping <- trm$mapping |> 
      left_join(operons, by = "id") |> 
      select(-id) |> 
      select(id = operon, term_id)

    list(
      terms = trm$terms,
      mapping = mapping
    )
  }) |> 
    set_names(ontologies)
}

#' Annotate plasmid genes
#'
#' Annotate plasmid genes in the NCBI GTF with gene symbols from
#'  https://pmc.ncbi.nlm.nih.gov/articles/PMC3754741.
#' 
#' Here we assume that the relative gene coordinates are the same in both annotations, only 
#' the origin is different. Hence, the "anchor" gene has to be provided, a gene for which locus ID and
#' gene symbol are know. It will be used to link the two annotations.
#' 
#' 
#' @param genes Gene information from NCBI, extracted from a GTF file, lacking gene symbols for the plasmid.
#' @param gff GRanges object representing the GFF file from https://pmc.ncbi.nlm.nih.gov/articles/PMC3754741.
#' @param anchor_gene A named string with a known gene, linking its id to the gene symbol.
#' @param len Length of the plasmid.
#'
#' @returns
#'
#' @export
#' @examples
ncbi_plasmid_annotations <- function(genes, gff, anchor_gene = c("B4U62_RS22265" = "rapP"), len = 84215) {
  g <- genes |> 
    filter(chr == "NZ_CP020103.1")
  f <- gff |> 
    as_tibble() |> 
    filter(type == "gene") |> 
    select(start, Name)
  # Anchor gene start in NCBI
  s1 <- g |> 
    filter(id == names(anchor_gene)) |> 
    pull(start)
  # Anchor gene start in GFF
  s2 <- f |> 
    filter(Name == anchor_gene) |> 
    pull(start)
  # Offset between origins in the two annotations
  dif <- s2 - s1
  # Add gene symbols to NCBI, merge by gene start with offset
  g |> 
    mutate(s = start + dif) |> 
    mutate(s = if_else(s > len, s - len, s))  |> 
    left_join(f, by = c("s" = "start")) |> 
    select(-c(gene_symbol, s)) |> 
    rename(gene_symbol = Name) |> 
    relocate(gene_symbol, .after = id)
}

annotate_plasmid_genes <- function(genes, plasmid_annot) {
  chrom <- genes |> 
    filter(chr == "NZ_CP020102.1")
  bind_rows(chrom, plasmid_annot)
}






match_sequences_to_genome <- function(genome, sequences, max_mismatch = 0, with_indels = FALSE) {
  sequences |> 
    dplyr::mutate(
      hits = purrr::map2(sequence, strand, ~{
        s <- Biostrings::DNAString(.x)
        if(.y == "-")
          s <- Biostrings::reverseComplement(s)
        mp <- Biostrings::matchPattern(s, genome, max.mismatch = max_mismatch, with.indels = with_indels)
        if(length(mp) > 0) {
          tibble::tibble(
            start = BiocGenerics::start(mp),
            end = BiocGenerics::end(mp),
            width = BiocGenerics::width(mp)
          )
        } else {
          tibble::tibble(start = NA, end = NA, width = NA)
        }
        
      }, .progress = TRUE)
    ) |> 
    tidyr::unnest(hits) |> 
    dplyr::select(-sequence) |>
    dplyr::arrange(start)
}

match_subtiwiki_names <- function(genes, sw_genes, genome, max_mismatch = 0, with_indels = FALSE) {
  sequences <- sw_genes |> 
    select(
      id,
      name,
      sequence = genomic_annotations.dna_sequence,
      strand = genomic_annotations.orientation
    ) |> 
    mutate(strand = if_else(strand == "Positive", "+", "-")) |> 
    drop_na()

  # All matches/non-matches
  seq_match <- match_sequences_to_genome(genome[[1]], sequences, max_mismatch, with_indels)
  
  # Unique matches
  seq_good_match <- seq_match |> 
    drop_na() |> 
    distinct(start, end, strand, .keep_all = TRUE) |> # Remove multiple sequences matching the same locus
    distinct(id, .keep_all = TRUE) |> # Remove the same sequence mapping multiple loci
    select(sw_id = id, sw_gene_symbol = name, name, start, end, strand)

  # Gene table with new column
  new_genes <- genes |> 
    rename(ncbi_gene_symbol = gene_symbol) |> 
    left_join(seq_good_match, by = join_by(start, end, strand)) |> 
    mutate(gene_symbol = if_else(is.na(sw_gene_symbol), ncbi_gene_symbol, sw_gene_symbol))

  list(
    sw_genes = sequences  |> select(-c(sequence, strand)),
    match_results = seq_match,
    good_matches = seq_good_match,
    genes = new_genes
  )
}

match_subtiwiki_stats <- function(sw_match) {
  chgenes <- sw_match$genes |>
       filter(chr == "NZ_CP020102.1") 
  tribble(
    ~Genes, ~Count, ~tot_base,
    "SubtiWiki sequences", sw_match$sw_genes |> nrow(), "sw",
    "Sequences matches", sw_match$match_results |> drop_na() |> distinct(id) |> nrow(), "sw",
    "Unique matches", sw_match$good_matches |> nrow(), "sw",
    "Mapped sequence loci match NCBI coordinates", sw_match$genes |> filter(!is.na(sw_id)) |> nrow(), "sw",
    "NCBI annotated genes", chgenes |> nrow(), "ncbi",
    "NCBI genes that changed name", chgenes |> filter(ncbi_gene_symbol != sw_gene_symbol) |> nrow(), "ncbi",
    "Non-ribosomal NCBI genes that gained name", chgenes |> filter(!(biotype %in% c("tRNA", "rRNA")) & is.na(ncbi_gene_symbol) & !is.na(sw_gene_symbol) & !str_detect(sw_gene_symbol, "BSU")) |> nrow(), "ncbi",
    "Protein coding genes not matched between NCBI and SubtiWiki", chgenes |> filter(is.na(sw_id) & biotype == "protein_coding") |> nrow(), "ncbi"
  ) |> 
    mutate(total = case_match(
      tot_base,
      "sw" ~ sw_match$sw_genes |> nrow(),
      "ncbi" ~ chgenes |> nrow()
    )) |> 
    mutate(Percentage = 100 * Count / total) |> 
    select(-c(tot_base, total))
}


manually_curate_genes <- function(genes, curated_genes) {
  genes |> 
    left_join(select(curated_genes, id, curated_gene_symbol = gene_symbol), by = join_by(id)) |> 
    mutate(gene_symbol = if_else(is.na(curated_gene_symbol), gene_symbol, curated_gene_symbol))
}