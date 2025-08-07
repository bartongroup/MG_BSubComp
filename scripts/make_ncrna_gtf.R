suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(readr)
  library(Biostrings)
})

get_ncrna_sequences <- function(sw_genes) {
  s <- sw_genes |> 
    dplyr::filter(str_detect(name, "^S\\d+$")) |> 
    dplyr::select(
      name,
      locus_tag = genomic_annotations.locus_tag,
      product = product,
      sw_start = genomic_annotations.start,
      sw_end = genomic_annotations.end,
      strand = genomic_annotations.orientation,
      sequence = genomic_annotations.dna_sequence
    ) |> 
      dplyr::mutate(strand = dplyr::if_else(strand == "Positive", "+", "-"))
}

match_ncrna_genome <- function(sw_genes, genome, max_mismatch = 0, with_indels = FALSE) {
  ncrna <- get_ncrna_sequences(sw_genes)

  ncrna |> 
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


collapse_attributes <- function(df) {
  df |>
    dplyr::mutate(value = paste0('"', value, '"')) |>
    tidyr::unite(attr, c(key, value), sep = " ") |>
    dplyr::pull(attr) |>
    stringr::str_c(collapse = "; ")
}

make_ncnra_gtf <- function(ncrna, gtf_file, seqname = "NZ_CP020102.1") {
  ncrna |> 
    tidyr::drop_na() |> 
    dplyr::rowwise() |> 
    dplyr::group_split() |> 
    purrr::map(function(g) {
      gene_attr <- tibble::tribble(
        ~key, ~value,
        "gene_id", g$name,
        "transcript_id", "",
        "gbkey", "Gene",
        "gene", g$name,
        "gene_biotype", "ncRNA",
        "locus_tag", g$locus_tag
      ) |> collapse_attributes()
      cds_attribute <- tibble::tribble(
        ~key, ~value,
        "gene_id", g$name,
        "transcript_id", stringr::str_glue("transcript_{g$name}"),
        "gbkey", "CDS",
        "gene", g$name,
        "product", g$product,
        "transl_table", "1",
        "exon_number", "1"
      ) |> collapse_attributes()

      tibble(
        seqname = c(seqname, seqname),
        source = c("SubtiWiki", "SubtiWiki"),
        feature = c("gene", "CDS"),
        start = c(g$start, g$start),
        end = c(g$end, g$end),
        score = c(".", "."),
        strand = c(g$strand, g$strand),
        frame = c(0, 0),
        attribute = c(gene_attr, cds_attribute)
      )
    }) |> 
      purrr::list_rbind() |> 
      readr::write_tsv(gtf_file, escape = "none", col_names = FALSE)
}


#######################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Three arguments must be supplied: genome_file, subtiwiki_file and output_file")
}

#args <- c(
#  "rna_seq/genome/GCF_002055965.1_ASM205596v1_genomic.fna",
#  "info/2025-07-21_09-26-20_subti_wiki_export_genes.csv",
#  "info/ncrna.gtf"
#)

genome_file <- args[1]
subtiwiki_file <- args[2]
output_file <- args[3]

genome <- Biostrings::readDNAStringSet(genome_file)
sw_genes <- readr::read_csv(subtiwiki_file, guess_max = 6000, show_col_types = FALSE)

ncrna <- match_ncrna_genome(sw_genes, genome[[1]], max_mismatch = 2, with_indels = TRUE)
sav_ncrna <- make_ncnra_gtf(ncrna, output_file)
