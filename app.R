

# =====================================================================
# North Carolina SIDS - Spatial Statistics Dashboard
# Methods: (1) Global Moran's I (+ Geary's C)  (2) Local Moran's I (LISA)
#          (3) Spatial regression (OLS vs spatial lag vs spatial error)
# Data:    nc.shp shipped with the sf package (100 counties)
#
# Install once:
# install.packages(c("shiny","bslib","sf","spdep","spatialreg",
#                    "leaflet","ggplot2","plotly"))
# Run:  shiny::runApp("nc_sids_dashboard.R")   (or press "Run App" in RStudio)
# =====================================================================

library(shiny)
library(bslib)
library(sf)
library(spdep)
library(spatialreg)
library(leaflet)
library(ggplot2)
library(plotly)

# spdep also exports a function called card(); make sure bslib's version is used
card <- bslib::card

# ---------------------------- Data ---------------------------------
nc_raw  <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
nc_proj <- st_transform(nc_raw, 32119)          # NC state plane (metres) for neighbours
nc_ll   <- st_transform(nc_raw, 4326)           # lon/lat for leaflet
cents   <- st_coordinates(st_centroid(st_geometry(nc_proj)))

var_choices <- c("SIDS rate (per 1,000 births)" = "rate",
                 "Non-white birth proportion"   = "nwprop",
                 "Log of births"                = "logbirths")

make_dat <- function(p) {
  d   <- nc_ll
  bir <- d[[paste0("BIR", p)]]
  d$rate      <- d[[paste0("SID", p)]] / bir * 1000
  d$nwprop    <- d[[paste0("NWBIR", p)]] / bir
  d$logbirths <- log(bir)
  d
}

cluster_cols <- c("High-High" = "#d7191c", "Low-Low" = "#2c7bb6",
                  "High-Low"  = "#fdae61", "Low-High" = "#abd9e9",
                  "Not significant" = "#e0e0e0")

# ------------------------ Helper functions -------------------------
poly_map <- function(d, vals, pal, title, lab_name = title) {
  labs <- lapply(sprintf("<strong>%s County</strong><br/>%s: %s",
                         d$NAME, lab_name,
                         ifelse(is.numeric(vals), formatC(vals, digits = 3, format = "g"), as.character(vals))),
                 HTML)
  leaflet(d, options = leafletOptions(zoomControl = TRUE)) |>
    addTiles() |>
    addPolygons(fillColor = pal(vals), fillOpacity = 0.85, color = "white", weight = 1,
                label = labs,
                highlightOptions = highlightOptions(weight = 3, color = "#333", bringToFront = TRUE)) |>
    addLegend(pal = pal, values = vals, title = title, position = "bottomright")
}

interpret <- function(stat, expected, p, alpha, what = "Moran's I") {
  if (p < alpha) {
    if (stat > expected) sprintf("Significant positive spatial autocorrelation (p = %.4f < %.2f): similar values cluster in space.", p, alpha)
    else sprintf("Significant negative spatial autocorrelation (p = %.4f < %.2f): neighbouring counties tend to be dissimilar.", p, alpha)
  } else {
    sprintf("No significant spatial autocorrelation at the %.0f%% level (p = %.4f): the pattern is consistent with spatial randomness.", alpha * 100, p)
  }
}

# ------------------------------ UI ---------------------------------
ui <- page_sidebar(
  title = "North Carolina SIDS Dashboard",
  theme = bs_theme(version = 5, bootswatch = "flatly", primary = "#2C5F7C"),
  fillable = FALSE,
  sidebar = sidebar(
    width = 320,
    h5("1. Data"),
    selectInput("period", "Time period", c("1974-78" = "74", "1979-84" = "79")),
    selectInput("var", "Variable (Moran's I and LISA tabs)", var_choices),
    hr(),
    h5("2. Spatial weights"),
    selectInput("nbtype", "Neighbour definition",
                c("Queen contiguity" = "queen", "Rook contiguity" = "rook",
                  "k-nearest neighbours" = "knn")),
    conditionalPanel("input.nbtype == 'knn'",
                     sliderInput("k", "Number of neighbours (k)", 1, 10, 4, step = 1)),
    selectInput("style", "Weights style",
                c("Row-standardised (W)" = "W", "Binary (B)" = "B")),
    hr(),
    h5("3. Inference"),
    sliderInput("alpha", "Significance level (\u03b1)", 0.01, 0.10, 0.05, step = 0.01),
    helpText("Tip: row-standardised weights are recommended, especially for the regression models.")
  ),
  navset_card_tab(
    id = "method",
    # ---------- Method 1 ----------
    nav_panel("Global Moran's I",
              layout_columns(col_widths = c(12), uiOutput("moran_boxes")),
              layout_columns(col_widths = c(6, 6),
                             card(card_header("Choropleth map"), leafletOutput("map_var", height = 420)),
                             card(card_header("Moran scatterplot"), plotlyOutput("moran_plot", height = 420))),
              card(card_header("Interpretation"),
                   uiOutput("moran_interp"),
                   checkboxInput("geary", "Also compute Geary's C", FALSE),
                   uiOutput("geary_out"))
    ),
    # ---------- Method 2 ----------
    nav_panel("Local Moran's I (LISA)",
              layout_columns(col_widths = c(4, 8),
                             card(card_header("Options"),
                                  selectInput("padj", "Multiple-testing adjustment",
                                              c("None" = "none", "Bonferroni" = "bonferroni",
                                                "False discovery rate (BH)" = "fdr")),
                                  helpText("LISA runs one test per county, so adjusting p-values reduces false positives."),
                                  tableOutput("lisa_tab")),
                             card(card_header("LISA cluster map"), leafletOutput("map_lisa", height = 480)))
    ),
    # ---------- Method 3 ----------
    nav_panel("Spatial regression",
              layout_columns(col_widths = c(4, 8),
                             card(card_header("Model setup"),
                                  p(strong("Response:"), "SIDS rate per 1,000 births"),
                                  checkboxGroupInput("preds", "Predictors",
                                                     c("Non-white birth proportion" = "nwprop",
                                                       "Log of births" = "logbirths"),
                                                     selected = "nwprop"),
                                  selectInput("model", "Model to inspect",
                                              c("OLS (non-spatial)" = "ols", "Spatial lag (SAR)" = "lag",
                                                "Spatial error (SEM)" = "err")),
                                  helpText("OLS assumes independence. If its residuals are spatially autocorrelated, a lag or error model is more appropriate.")),
                             card(card_header("Model comparison"), tableOutput("cmp_tab"))),
              layout_columns(col_widths = c(6, 6),
                             card(card_header("Model summary"), verbatimTextOutput("model_txt")),
                             card(card_header("Residual map"), leafletOutput("map_resid", height = 420)))
    ),
    # ---------- About ----------
    nav_panel("About the data",
              card(
                h4("Dataset"),
                p("County-level data for the 100 counties of North Carolina (USA) on Sudden Infant Death Syndrome (SIDS), from Cressie (1993), distributed with the", code("sf"), "package."),
                tags$ul(
                  tags$li(strong("SID74 / SID79:"), " SIDS deaths in 1974-78 / 1979-84"),
                  tags$li(strong("BIR74 / BIR79:"), " live births in each period"),
                  tags$li(strong("NWBIR74 / NWBIR79:"), " non-white births in each period")),
                p("The SIDS rate is computed per 1,000 births so that counties of different sizes are comparable. Associations shown here are not causal; the non-white birth proportion is a crude demographic indicator."),
                h4("Methods"),
                tags$ul(
                  tags$li(strong("Global Moran's I / Geary's C:"), " test for overall spatial autocorrelation."),
                  tags$li(strong("Local Moran's I (LISA):"), " locates clusters (High-High, Low-Low) and spatial outliers."),
                  tags$li(strong("Spatial regression:"), " compares OLS with spatial lag and spatial error models, and checks residuals for remaining autocorrelation."))
              )
    )
  )
)

# ---------------------------- Server -------------------------------
server <- function(input, output, session) {
  
  dat <- reactive(make_dat(input$period))
  x   <- reactive(dat()[[input$var]])
  var_label <- reactive(names(var_choices)[match(input$var, var_choices)])
  
  lw <- reactive({
    nb <- switch(input$nbtype,
                 queen = poly2nb(nc_proj, queen = TRUE),
                 rook  = poly2nb(nc_proj, queen = FALSE),
                 knn   = knn2nb(knearneigh(cents, k = input$k)))
    nb2listw(nb, style = input$style, zero.policy = TRUE)
  })
  
  # ---------- Method 1: Global Moran's I ----------
  mor <- reactive(moran.test(x(), lw(), zero.policy = TRUE))
  
  output$moran_boxes <- renderUI({
    m <- mor()
    layout_columns(
      value_box("Moran's I", sprintf("%.3f", m$estimate[1]), theme = "primary"),
      value_box("Expected I", sprintf("%.3f", m$estimate[2]), theme = "secondary"),
      value_box("z-score", sprintf("%.2f", m$statistic), theme = "info"),
      value_box("p-value", format.pval(m$p.value, digits = 3, eps = 1e-4),
                theme = if (m$p.value < input$alpha) "success" else "warning"))
  })
  
  output$moran_interp <- renderUI({
    m <- mor()
    p(interpret(m$estimate[1], m$estimate[2], m$p.value, input$alpha))
  })
  
  output$geary_out <- renderUI({
    req(input$geary)
    g <- geary.test(x(), lw(), zero.policy = TRUE)
    p(sprintf("Geary's C = %.3f (expected %.3f), p = %.4f. Values below 1 indicate positive autocorrelation.",
              g$estimate[1], g$estimate[2], g$p.value))
  })
  
  output$map_var <- renderLeaflet({
    v <- x()
    poly_map(dat(), v, colorNumeric("YlOrRd", v), var_label())
  })
  
  output$moran_plot <- renderPlotly({
    v  <- x()
    df <- data.frame(county = dat()$NAME, value = v,
                     lagged = lag.listw(lw(), v, zero.policy = TRUE))
    g <- ggplot(df, aes(value, lagged)) +
      geom_hline(yintercept = mean(df$lagged), linetype = "dashed", colour = "grey60") +
      geom_vline(xintercept = mean(df$value),  linetype = "dashed", colour = "grey60") +
      geom_point(aes(text = county), colour = "#2C5F7C", alpha = 0.8, size = 2) +
      geom_smooth(method = "lm", se = FALSE, colour = "#d7191c", formula = y ~ x) +
      labs(x = var_label(), y = "Spatially lagged value") +
      theme_minimal(base_size = 12)
    ggplotly(g, tooltip = c("text", "x", "y"))
  })
  
  # ---------- Method 2: LISA ----------
  lisa <- reactive({
    v  <- x()
    z  <- as.numeric(scale(v))
    lz <- lag.listw(lw(), z, zero.policy = TRUE)
    li <- localmoran(v, lw(), zero.policy = TRUE)
    p  <- p.adjust(li[, 5], method = input$padj)
    q  <- ifelse(z >= 0 & lz >= 0, "High-High",
                 ifelse(z <  0 & lz <  0, "Low-Low",
                        ifelse(z >= 0 & lz <  0, "High-Low", "Low-High")))
    q[p >= input$alpha] <- "Not significant"
    list(cluster = factor(q, levels = names(cluster_cols)), Ii = li[, 1], p = p)
  })
  
  output$map_lisa <- renderLeaflet({
    l   <- lisa()
    d   <- dat()
    pal <- colorFactor(unname(cluster_cols), levels = names(cluster_cols))
    labs <- lapply(sprintf("<strong>%s County</strong><br/>Cluster: %s<br/>Local I: %.3f<br/>p: %.4f",
                           d$NAME, l$cluster, l$Ii, l$p), HTML)
    leaflet(d) |>
      addTiles() |>
      addPolygons(fillColor = pal(l$cluster), fillOpacity = 0.85, color = "white", weight = 1,
                  label = labs,
                  highlightOptions = highlightOptions(weight = 3, color = "#333", bringToFront = TRUE)) |>
      addLegend(colors = unname(cluster_cols), labels = names(cluster_cols),
                title = "LISA cluster", position = "bottomright")
  })
  
  output$lisa_tab <- renderTable({
    t <- table(lisa()$cluster)
    data.frame(Cluster = names(t), Counties = as.integer(t))
  })
  
  # ---------- Method 3: Spatial regression ----------
  fits <- reactive({
    validate(need(length(input$preds) > 0, "Select at least one predictor."))
    f <- as.formula(paste("rate ~", paste(input$preds, collapse = " + ")))
    d <- st_drop_geometry(dat())
    tryCatch(
      list(ols = lm(f, data = d),
           lag = lagsarlm(f, data = d, listw = lw(), zero.policy = TRUE),
           err = errorsarlm(f, data = d, listw = lw(), zero.policy = TRUE)),
      error = function(e) NULL)
  })
  
  output$cmp_tab <- renderTable({
    f <- fits()
    validate(need(!is.null(f), "The spatial models could not be fitted with these weights. Try row-standardised weights or a different neighbour definition."))
    nm <- c(ols = "OLS", lag = "Spatial lag (SAR)", err = "Spatial error (SEM)")
    do.call(rbind, lapply(names(f), function(k) {
      m <- moran.test(residuals(f[[k]]), lw(), zero.policy = TRUE)
      data.frame(Model = nm[[k]], AIC = AIC(f[[k]]),
                 `Residual Moran's I` = unname(m$estimate[1]),
                 `Residual p-value` = m$p.value, check.names = FALSE)
    }))
  }, digits = 4)
  
  output$model_txt <- renderPrint({
    f <- fits()
    validate(need(!is.null(f), "Model could not be fitted."))
    summary(f[[input$model]])
  })
  
  output$map_resid <- renderLeaflet({
    f <- fits()
    validate(need(!is.null(f), "Model could not be fitted."))
    r   <- as.numeric(residuals(f[[input$model]]))
    lim <- max(abs(r))
    pal <- colorNumeric("RdBu", domain = c(-lim, lim), reverse = TRUE)
    poly_map(dat(), r, pal, "Residual", "Residual")
  })
}

shinyApp(ui, server)

