make_metadata <- function(mfile) {
  readxl::read_excel(mfile) |>
    janitor::clean_names() |>
    mutate(group = if_else(genotype == "NRS6942_Control", "ctrl", "mut")) |> 
    separate_wider_delim(sample_name, delim = "_", names = c("group_id", "replicate"), cols_remove = FALSE) |> 
    unite(sample, c(group, replicate), remove = FALSE) |> 
    mutate(
      sample = as_factor(sample),
      group = as_factor(group)
    ) |> 
    rename(raw_sample = sample_id) |> 
    mutate(bad = FALSE)
}

