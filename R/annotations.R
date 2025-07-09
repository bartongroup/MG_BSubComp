
get_gene_annotations_from_gff <- function(gff) {
  coords <- gff |> 
    as_tibble() |> 
    select(
      chr = seqnames,
      start,
      end,
      strand
    )
  mc <- mcols(gff) |> 
    as_tibble()

 bind_cols(coords, mc) |> 
  filter(type == "CDS")  |> 
  select(
    chr,
    start,
    end,
    strand,
    id = Name,
    gene_symbol = gene,
    description = product
  )
}

get_go_mapping_from_gff <- function(gff) {
  mapping <- mcols(gff) |> 
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


get_functional_terms <- function(gff, kg_spec) {
  cat("Loading GO term data\n")
  go <- get_go_mapping_from_gff(gff)
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
