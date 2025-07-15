write_table <- function(obj, dir = "tab") {
  path <- file.path(dir)
  if (!dir.exists(path)) dir.create(path)
  obj_name <- deparse(substitute(obj))
  cat(paste("  Writing", obj_name))
  file_name <- file.path(path, str_glue("{obj_name}.csv"))
  obj |> 
    mutate(across(where(is.numeric), ~signif(.x, 4))) |> 
    write_csv(file_name)
  cat("\n")
}
