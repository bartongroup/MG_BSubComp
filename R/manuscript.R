mn_plot_volcano_toxins <- function(de_toxan, toxan, fdr_limit = 0.01, logfc_limit = 1) {
  d <- de_toxan |> 
    mutate(
      toxan = id %in% toxan$id,
      sig = FDR < fdr_limit & logFC < -logfc_limit
    )
  
  d_ntox <- d |> filter(!toxan)
  d_tox <- d |> filter(toxan)
  d_toxsig <- d |> filter(toxan & sig)
  
  ggplot(d, aes(x = logFC, y = -log10(PValue), label = gene_symbol)) +
    th +
    theme(legend.position = "none") +
    geom_point(data = d_ntox, colour = "grey80", size = 0.6, alpha = 0.3) +
    geom_point(data = d_tox, aes(colour = sig), size = 1.3) +
    scale_colour_manual(values = c("deepskyblue2", "black")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    geom_text_repel(data = d_toxsig, size = 4.0, min.segment.length = 0.03, force = 20, segment.colour = "grey50") +
    geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
    #facet_wrap(~ contrast) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P))
}