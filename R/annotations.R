
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
    gene_symbol = gene
  )

  # only rows marked as CDS contain gene desscription
  CDS <-  bind_cols(coords, mc) |> 
  filter(type == "CDS")  |> 
  select(
    id = gene_id,
    description = product
  )

  left_join(genes, CDS)
}

get_go_mapping_from_gtf <- function(gtf) {
  mapping <- mcols(gtf) |> 
    as_tibble() |> 
    filter(type == "CDS") |> 
    select(
      gene_symbol = gene,
      term_id = Ontology_term
    ) |> 
      unnest(term_id) |> 
    distinct()
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
  terms = list(
    go = go,
    kg = kg
  )
}


prepare_terms_fenr <- function(terms, all_features) {
  ontologies <- names(terms)
  
  map(ontologies, function(ont) {
    trm <- terms[[ont]]
    fenr::prepare_for_enrichment(trm$terms, trm$mapping, all_features, feature_name = "gene_symbol")
  }) |> 
    set_names(ontologies)
}
