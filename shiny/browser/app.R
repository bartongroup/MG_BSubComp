library(shiny)
library(bslib)
library(qs2)
library(dplyr)
library(tibble)
library(stringr)
library(purrr)
library(ggplot2)

source("func.R")

de <- qs_read("./data/de.qs")
genes <- qs_read("./data/genes.qs")

de <- de |> 
  left_join(select(genes, id, description), by = "id")

chroms <- unique(genes$chr)
chrlens <- map(chroms, function(chrom) {
  genes |> 
    filter(chr == chrom) |> 
    pull(end) |> 
    max()
}) |> 
  set_names(chroms)
  

card_de <- card(
  card_header(
    "DE results",
    class = "d-flex justify-content-between"
  ),
  card_body(
    min_height = 600,
    DT::dataTableOutput(
      outputId = "de_table"
    )
  )
)

range_str <- textInput("range_str", "Range", value = "")
ab_left <- actionButton("left", HTML("<span class='small'><i class='glyphicon glyphicon-arrow-left'></i></span>"))
ab_in <- actionButton("zoomin", HTML("<span class='small'><i class='glyphicon glyphicon-plus'></i></span>"))
ab_out <- actionButton("zoomout", HTML("<span class='small'><i class='glyphicon glyphicon-minus'></i></span>"))
ab_right <- actionButton("right", HTML("<span class='small'><i class='glyphicon glyphicon-arrow-right'></i></span>"))

card_browser <- card(
  card_header(
    "Browser",
    shiny::span(range_str, ab_left, ab_in, ab_out, ab_right),
    class = "d-flex justify-content-between"
  ),
  card_body(
    min_height = 500,
    max_height = 500,
    shiny::plotOutput(
      outputId = "browser_plot",
      width = "100%",
      height = "500px",
      dblclick = "dblclick",
      brush = "brush"
    )
  )
)


 ui <- page_fillable(
    theme = bs_theme(bootswatch = "spacelab", spacer = "0.8rem") |>
      bs_add_rules(".info-pop { max-width: 500px; background-color: #f0f9e8; }"),
    title = "BROWSER",

    layout_columns(
      col_widths = c(5, 7),
      style = "height: 100%;", 
      card_de,
      card_browser
    )  
)


server <- function(input, output, session) {
  # Prevents RStudio from crashing when Shiny window closed manually
  session$onSessionEnded(function() {
    shiny::stopApp()
  })
 
  make_table <- shiny::reactive({
    fdr_limit <- 0.01
    req(fdr_limit)
    make_de_table(de, fdr_limit)
  })
  
  output$de_table <- DT::renderDataTable({
    dt <- make_table() |> 
      select(-id)
    DT::datatable(
      dt,
      options = list(paging = FALSE, dom = "ft"),
      style = "bootstrap",
      selection = "single",
      rownames = FALSE
    ) |>
      DT::formatStyle(columns = colnames(dt), fontSize = "80%") |> 
      DT::formatSignif(columns = c("fdr", "log_fc"), digits = 3)
  })

  observeEvent(input$de_table_rows_selected, ignoreNULL = FALSE, {
    rows_sel <- input$de_table_rows_selected
    if(is.null(rows_sel)) {
      sr <- ""
    } else {
      dt <- make_table()
      ids <- dt[rows_sel, ]$id
      r <- genes |> 
        filter(id == ids)
      range <- list(chr = r$chr, start = r$start, end = r$end) |> 
        expand_range(chrlens, margin = 0.5)
      sr <- make_str_range(range)
    }
    updateTextInput(session, "range_str", value = sr)
  })

  output$browser_plot <- renderPlot({
    srange <- input$range_str
    req(srange)
    if(srange != "") {
      range <- parse_str_range(srange)
      plot_browser(genes, de, range, fdr_limit = 0.01)
    }
  })
    

  observeEvent(input$dblclick, {
    centre <- input$dblclick$x
    crange <- input$range_str
    range <- parse_str_range(crange)
    delta <- centre - 0.5 * (range$start + range$end)
    range$start <- range$start + delta
    range$end <- range$end + delta
    sr <- make_str_range(range)
    updateTextInput(session, "range_str", value = sr)
  })
  
  observeEvent(input$brush, {
    crange <- input$range_str
    range <- parse_str_range(crange)
    range$start <- input$brush$xmin
    range$end <- input$brush$xmax
    if(range$end < range$start + 10) range$end <- range$start + 10
    sr <- make_str_range(range)
    updateTextInput(session, "range_str", value = sr)
    session$resetBrush("brush")
  })

  observeEvent(input$left, {
    crange <- input$range_str
    range <- parse_str_range(crange)
    width <- range$end - range$start
    range$start <- max(range$start - 0.2 * width, 0)
    range$end <- max(range$end - 0.2 * width, 0)
    sr <- make_str_range(range)
    updateTextInput(session, "range_str", value = sr)
  })

  observeEvent(input$right, {
    crange <- input$range_str
    range <- parse_str_range(crange)
    width <- range$end - range$start
    chlen <- chrlens[[range$chr]]
    range$start <- min(range$start + 0.2 * width, chlen)
    range$end <- min(range$end + 0.2 * width, chlen)
    sr <- make_str_range(range)
    updateTextInput(session, "range_str", value = sr)
  })

  observeEvent(input$zoomin, {
    crange <- input$range_str
    range <- parse_str_range(crange)
    width <- range$end - range$start
    if(width > 10) {
      range$start <- range$start + width / 4
      range$end <- range$end - width / 4
      sr <- make_str_range(range)
      updateTextInput(session, "range_str", value = sr)
    }
  })
  
  observeEvent(input$zoomout, {
    crange <- input$range_str
    range <- parse_str_range(crange)
    width <- range$end - range$start
    chlen <- chrlens[[range$chr]]
    range$start <- max(range$start - width / 4, 0)
    range$end <- min(range$end + width / 4, chlen)
    sr <- make_str_range(range)
    updateTextInput(session, "range_str", value = sr)
  })



}


# Run the application
shinyApp(ui = ui, server = server)
