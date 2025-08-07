
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




