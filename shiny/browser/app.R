library(shiny)
library(bslib)
library(qs2)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(purrr)
library(ggplot2)
library(scales)
library(DT)

source("func.R")

# Prepare UI elements

style_small <- tags$style("
  .small {
    padding: 8px 8px 9px !important; 
    font-size: 10px; !important;
    height: 28px !important;
    line-height: 1.0 !important;
  }
")

# Wheel → Shiny input (written by ChatGPT 5.0)
script_wheel <- tags$script("
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
")

# Differential expression gene table card
card_de <- card(
  card_header(
    "DE results",
    class = "d-flex justify-content-between",
    div(
      class = "d-flex align-items-center",
      tags$label("FDR limit", `for` = "fdr_limit", class = "me-2"),
      div(
        style = "width: 100px;",
        numericInput(
          inputId = "fdr_limit",
          label = NULL,
          value = 0.01,
          min = 0,
          max = 1
        )
      )
    )
  ),
  card_body(
    min_height = 600,
    dataTableOutput(
      outputId = "de_table"
    )
  )
)

# Genomic browser card
card_browser <- card(
  card_header(
    "Browser",
    class = "d-flex justify-content-between",
    div(
      class = "d-flex align-items-center",
      tags$label("Range", `for` = "range_str", class = "me-2"),
      textInput("range_str", NULL, value = "")
    )
  ),
  card_body(
    min_height = 500,
    max_height = 500,
    plotOutput(
      outputId = "browser_plot",
      width = "100%",
      height = "500px",
      dblclick = "dblclick",
      brush = "brush"
    ),
    uiOutput("nav_buttons")
  )
)



### UI

ui <- page_fillable(
  theme = bs_theme(
    bootswatch = "spacelab",
    spacer = "0.8rem"
  ) |>
    bs_add_rules(".info-pop { max-width: 500px; background-color: #f0f9e8; }"),
  title = "DE BROWSER",
  style_small,
  script_wheel,
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

  operons <- reactive({
    qs_read("data/operons.qs")
  })

  de <- reactive({
    qs_read("data/de.qs")
  })

  gene_data <- reactive({
    de() |> 
      left_join(genes(), by = join_by(id)) |> 
      mutate(gene_symbol = if_else(is.na(gene_symbol), id, gene_symbol))
  })

  chrlens <- reactive({
    gns <- genes()
    req(gns)
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
        expand_range(chrlens(), margin = 0.8)
      sr <- make_str_range(range)
    }
    updateTextInput(session, "range_str", value = sr)
  })

  plot_ready <- reactive({
    nzchar(input$range_str)
  })

  output$browser_plot <- renderPlot({
    req(plot_ready())
    dat <- gene_data()
    ops <- operons()
    req(dat, ops)
    srange <- input$range_str
    range <- parse_str_range(srange)
    plot_browser(dat, ops, range, fdr_limit = input$fdr_limit)
  })

  output$nav_buttons <- renderUI({
    req(plot_ready())
    div(
      class = "d-flex justify-content-center mt-2",
      actionButton("left",  icon("arrow-left"), class = "small"),
      actionButton("zoomin", icon("plus"),       class = "small"),
      actionButton("zoomout",icon("minus"),      class = "small"),
      actionButton("right", icon("arrow-right"), class = "small")
    )
  })
    
  observeEvent(input$dblclick, {
    range <- parse_str_range(input$range_str)
    chlen <- chrlens()[[range$chr]]
    delta <- input$dblclick$x - 0.5 * (range$start + range$end)
    range$start <- max(range$start + delta, 0)
    range$end <- min(range$end + delta, chlen)
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

  observeEvent(input$left, ignoreInit = TRUE, {
    range <- parse_str_range(input$range_str)
    width <- range$end - range$start
    range$start <- max(range$start - 0.2 * width, 0)
    range$end <- max(range$end - 0.2 * width, 0)
    sr <- make_str_range(range)
    updateTextInput(session, "range_str", value = sr)
  })

  observeEvent(input$right, ignoreInit = TRUE, {
    range <- parse_str_range(input$range_str)
    width <- range$end - range$start
    chlen <- chrlens()[[range$chr]]
    range$start <- min(range$start + 0.2 * width, chlen)
    range$end <- min(range$end + 0.2 * width, chlen)
    sr <- make_str_range(range)
    updateTextInput(session, "range_str", value = sr)
  })

  observeEvent(input$zoomin, ignoreInit = TRUE, {
    range <- parse_str_range(input$range_str)
    width <- range$end - range$start
    if(width > 10) {
      range$start <- range$start + width / 4
      range$end <- range$end - width / 4
      sr <- make_str_range(range)
      updateTextInput(session, "range_str", value = sr)
    }
  })
  
  observeEvent(input$zoomout, ignoreInit = TRUE, {
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
    base <- 1.2
    factor <- ifelse(w$deltaY > 0, 1 / base, base)   # down=zoom out, up=zoom in

    frac <- w$fracX

    # fracX is the fraction in a window that includes ticks and labels on the left;
    # needs correcting
    wm <- 1.06
    frac <- 1 - (1 - frac) * wm

    # Zoom in/out around the mouse position
    zoom_point <- range$start + frac * width
    range$start <- max(zoom_point - (zoom_point - range$start) / factor, 0)
    range$end <- min(zoom_point + (range$end - zoom_point) / factor, chlen)

    sr <- make_str_range(range)

    updateTextInput(session, "range_str", value = sr)
  })

}


# Run the application
shinyApp(ui = ui, server = server)
