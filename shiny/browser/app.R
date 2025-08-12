library(shiny)
library(bslib)
library(qs2)
library(dplyr)
library(tibble)
library(stringr)
library(purrr)
library(ggplot2)
library(DT)

source("func.R")

# Prepare UI elements

fdr_limit_input <- numericInput(
  inputId = "fdr_limit",
  label = "FDR limit",
  value = 0.01,
  min = 0,
  max = 1
)

card_de <- card(
  card_header(
    "DE results",
    span(fdr_limit_input),
    class = "d-flex justify-content-between"
  ),
  card_body(
    min_height = 600,
    dataTableOutput(
      outputId = "de_table"
    )
  )
)

range_str <- textInput("range_str", "Range", value = "")
ab_left <- actionButton("left", icon("arrow-left"), class = "small")
ab_in <- actionButton("zoomin", icon("plus"), class = "small")
ab_out <- actionButton("zoomout", icon("minus"), class = "small")
ab_right <- actionButton("right", icon("arrow-right"), class = "small")

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


### UI

ui <- page_fillable(
  theme = bs_theme(bootswatch = "spacelab", spacer = "0.8rem") |>
          bs_add_rules(".info-pop { max-width: 500px; background-color: #f0f9e8; }"),
  title = "DE BROWSER",

  tags$style(HTML("
    .small {
      padding: 8px 8px 9px !important; 
      font-size: 10px; !important;
      height: 28px !important;
      line-height: 1.0 !important;
    }
  ")),

  # Wheel → Shiny input
  tags$script(HTML("
    document.addEventListener('DOMContentLoaded', function() {
      const el = document.getElementById('browser_plot');        // plotOutput id
      if (!el) return;
      el.addEventListener('wheel', function(e){
        e.preventDefault();                            // don't scroll the page
        const rect = el.getBoundingClientRect();
        const fracX = (e.clientX - rect.left) / rect.width;  // 0..1 across the plot
        Shiny.setInputValue('p_wheel', {
          deltaY: e.deltaY,                            // >0 = down, <0 = up
          fracX: Math.min(Math.max(fracX, 0), 1),      // clamp
          ctrl:  e.ctrlKey,
          shift: e.shiftKey,
          time:  Date.now()
        }, {priority:'event'});
      }, {passive:false});
    });
  ")),

  layout_columns(
    col_widths = c(5, 7),
    style = "height: 100%;", 
    card_de,
    card_browser
  )  
)

### Server

server <- function(input, output, session) {
  # Prevents RStudio from crashing when Shiny window closed manually
  session$onSessionEnded(function() {
    shiny::stopApp()
  })
 
  genes <- reactive({
    qs_read("data/genes.qs")
  })

  de <- reactive({
    qs_read("data/de.qs")
  })

  gene_data <- reactive({
    de() |> 
      select(-gene_symbol) |> 
      left_join(genes(), by = join_by(id)) |> 
      mutate(gene_symbol = if_else(is.na(gene_symbol), id, gene_symbol))
  })

  chrlens <- reactive({
    gns <- genes()
    chroms <- unique(gns$chr)
    map(chroms, function(chrom) {
      gns |> 
        filter(chr == chrom) |> 
        pull(end) |> 
        max()
    }) |> 
    set_names(chroms)
  })
 
  make_de_table <- shiny::reactive({
    fdr_limit <- input$fdr_limit
    dat <- gene_data()
    req(fdr_limit, dat)
    dat |> 
      select(id, gene_symbol, description, log_fc, fdr) |> 
      filter(fdr < fdr_limit) |> 
      arrange(log_fc)
  })

  output$de_table <- DT::renderDataTable({
    dt <- make_de_table() |> 
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
      dt <- make_de_table()
      ids <- dt[rows_sel, ]$id
      r <- genes() |> 
        filter(id == ids)
      range <- list(chr = r$chr, start = r$start, end = r$end) |> 
        expand_range(chrlens(), margin = 0.5)
      sr <- make_str_range(range)
    }
    updateTextInput(session, "range_str", value = sr)
  })

  output$browser_plot <- renderPlot({
    srange <- input$range_str
    dat <- gene_data()
    req(srange, dat)
    if(srange != "") {
      range <- parse_str_range(srange)
      plot_browser(dat, range, fdr_limit = 0.01)
    }
  })
    
  observeEvent(input$dblclick, {
    centre <- input$dblclick$x
    range <- parse_str_range(input$range_str)
    delta <- centre - 0.5 * (range$start + range$end)
    range$start <- range$start + delta
    range$end <- range$end + delta
    sr <- make_str_range(range)
    updateTextInput(session, "range_str", value = sr)
  })
  
  observeEvent(input$brush, {
    range <- parse_str_range(input$range_str)
    range$start <- input$brush$xmin
    range$end <- input$brush$xmax
    if(range$end < range$start + 10) range$end <- range$start + 10
    sr <- make_str_range(range)
    updateTextInput(session, "range_str", value = sr)
    session$resetBrush("brush")
  })

  observeEvent(input$left, {
    range <- parse_str_range(input$range_str)
    width <- range$end - range$start
    range$start <- max(range$start - 0.2 * width, 0)
    range$end <- max(range$end - 0.2 * width, 0)
    sr <- make_str_range(range)
    updateTextInput(session, "range_str", value = sr)
  })

  observeEvent(input$right, {
    range <- parse_str_range(input$range_str)
    width <- range$end - range$start
    chlen <- chrlens()[[range$chr]]
    range$start <- min(range$start + 0.2 * width, chlen)
    range$end <- min(range$end + 0.2 * width, chlen)
    sr <- make_str_range(range)
    updateTextInput(session, "range_str", value = sr)
  })

  observeEvent(input$zoomin, {
    range <- parse_str_range(input$range_str)
    width <- range$end - range$start
    if(width > 10) {
      range$start <- range$start + width / 4
      range$end <- range$end - width / 4
      sr <- make_str_range(range)
      updateTextInput(session, "range_str", value = sr)
    }
  })
  
  observeEvent(input$zoomout, {
    range <- parse_str_range(input$range_str)
    width <- range$end - range$start
    chlen <- chrlens()[[range$chr]]
    range$start <- max(range$start - width / 4, 0)
    range$end <- min(range$end + width / 4, chlen)
    sr <- make_str_range(range)
    updateTextInput(session, "range_str", value = sr)
  })


  observeEvent(input$p_wheel, ignoreInit = TRUE, {
    w <- input$p_wheel
    range <- parse_str_range(input$range_str)
    width <- range$end - range$start
    chlen <- chrlens()[[range$chr]]
    if (width <= 0) return()

    # zoom step (tweak to taste)
    base <- if (isTRUE(w$shift)) 1.35 else if (isTRUE(w$ctrl)) 1.08 else 1.20
    factor <- if (w$deltaY > 0) base else 1/base   # down=zoom out, up=zoom in

    min_width <- 10L
    new_width <- max(round(width * factor), min_width)

    # center at mouse position within the current window
    center <- range$start + as.integer(round(w$fracX * width))
    range$start <- max(center - new_width %/% 2, 0)
    range$end <- min(range$start + new_width, chlen)
    sr <- make_str_range(range)

    updateTextInput(session, "range_str", value = sr)
    
  })



}


# Run the application
shinyApp(ui = ui, server = server)
