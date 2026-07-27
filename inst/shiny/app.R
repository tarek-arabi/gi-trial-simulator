# Explorer for GI trial designs. Launch with gitrialsim::run_explorer().
# ponytail: single file, no modules. Split only if this grows past a few hundred lines.

library(shiny)
library(gitrialsim)

packs <- list_packs()

ui <- fluidPage(
  title = "GI trial design explorer",
  tags$head(tags$style(HTML("
    body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif; }
    .result-n { font-size: 2.4rem; font-weight: 600; line-height: 1.1; }
    .result-label { color: #666; font-size: 0.85rem; text-transform: uppercase;
                    letter-spacing: 0.06em; }
    .note { color: #666; font-size: 0.85rem; }
    .panel-box { border: 1px solid #e0e0e0; padding: 1rem 1.25rem; margin-bottom: 1rem; }
  "))),
  h3("GI trial design explorer"),
  p(class = "note",
    "Sample sizes and boundaries come from rpact. Bayesian operating characteristics are simulated.
     Every scenario below is parameterised from published aggregate results; see the citation shown
     for the selected endpoint."),
  sidebarLayout(
    sidebarPanel(
      width = 4,
      selectInput("pack", "Scenario", choices = stats::setNames(packs$id, packs$title)),
      uiOutput("endpoint_ui"),
      hr(),
      numericInput("control_rate", "Control event rate", value = 0.0658,
                   min = 0.001, max = 0.999, step = 0.005),
      numericInput("treatment_rate", "Treatment event rate", value = 0.0395,
                   min = 0.001, max = 0.999, step = 0.005),
      actionButton("reset_rates", "Reset to published rates", class = "btn-sm"),
      hr(),
      selectInput("design", "Design",
                  c("Fixed" = "fixed",
                    "Group sequential" = "gs",
                    "Bayesian adaptive" = "bayes")),
      conditionalPanel(
        "input.design == 'gs'",
        numericInput("k", "Number of analyses", value = 3, min = 2, max = 6, step = 1),
        selectInput("spending", "Boundary type",
                    c("O'Brien-Fleming spending" = "asOF",
                      "Pocock spending" = "asP",
                      "O'Brien-Fleming classical" = "OF",
                      "Pocock classical" = "P")),
        selectInput("futility", "Futility",
                    c("None" = "none", "Non-binding O'Brien-Fleming" = "nonbinding_obf"))
      ),
      conditionalPanel(
        "input.design == 'bayes'",
        numericInput("bayes_n", "Maximum total sample size", value = 1000, min = 50, step = 50),
        numericInput("bayes_looks", "Number of analyses", value = 3, min = 1, max = 6, step = 1),
        numericInput("bayes_nsim", "Simulation replications", value = 2000,
                     min = 200, max = 20000, step = 500),
        p(class = "note", "Bayesian operating characteristics require simulation, so this is slower.")
      ),
      hr(),
      sliderInput("alpha", "One-sided type I error", value = 0.025,
                  min = 0.005, max = 0.10, step = 0.005),
      sliderInput("power", "Target power", value = 0.9, min = 0.5, max = 0.99, step = 0.01)
    ),
    mainPanel(
      width = 8,
      div(class = "panel-box",
          div(class = "result-label", "Maximum total sample size"),
          div(class = "result-n", textOutput("n_total", inline = TRUE)),
          uiOutput("n_detail")),
      div(class = "panel-box", verbatimTextOutput("design_print")),
      conditionalPanel("input.design == 'gs'",
                       div(class = "panel-box",
                           h5("Stopping boundaries"),
                           tableOutput("boundaries"))),
      div(class = "panel-box",
          h5("Sample size against treatment effect"),
          plotOutput("power_curve", height = "320px"),
          p(class = "note",
            "Total sample size required at the selected power, across treatment event rates.
             The vertical line marks the currently selected treatment rate.")),
      div(class = "panel-box",
          h5("Source"),
          textOutput("citation"))
    )
  )
)

server <- function(input, output, session) {
  pack <- reactive(load_pack(input$pack))

  output$endpoint_ui <- renderUI({
    eps <- pack()$endpoints
    labels <- vapply(eps, function(e) e$label %||% "", character(1))
    roles <- vapply(eps, function(e) e$role %||% "", character(1))
    selectInput("endpoint", "Endpoint",
                choices = stats::setNames(names(eps), labels),
                selected = names(eps)[roles == "primary"][1])
  })

  published <- reactive({
    req(input$endpoint)
    scenario(pack(), endpoint = input$endpoint)
  })

  observeEvent(published(), {
    updateNumericInput(session, "control_rate", value = published()$control_rate)
    updateNumericInput(session, "treatment_rate", value = published()$treatment_rate)
  })

  observeEvent(input$reset_rates, {
    updateNumericInput(session, "control_rate", value = published()$control_rate)
    updateNumericInput(session, "treatment_rate", value = published()$treatment_rate)
  })

  current <- reactive({
    req(input$control_rate, input$treatment_rate)
    validate(
      need(input$control_rate > 0 && input$control_rate < 1,
           "Control rate must be strictly between 0 and 1."),
      need(input$treatment_rate > 0 && input$treatment_rate < 1,
           "Treatment rate must be strictly between 0 and 1."),
      need(abs(input$treatment_rate - input$control_rate) > 1e-6,
           "The two event rates are equal, so there is no effect to power against.")
    )
    scenario(pack(), endpoint = input$endpoint,
             control_rate = input$control_rate,
             treatment_rate = input$treatment_rate)
  })

  design <- reactive({
    s <- current()
    switch(input$design,
      fixed = design_fixed(s, alpha = input$alpha, power = input$power),
      gs = design_group_sequential(
        s, alpha = input$alpha, power = input$power,
        k = input$k, type_of_design = input$spending, futility = input$futility
      ),
      bayes = {
        withProgress(message = "Simulating operating characteristics", value = 0.4, {
          design_bayesian(s, n_max = input$bayes_n, looks = input$bayes_looks,
                          nsim = input$bayes_nsim, seed = 1)
        })
      }
    )
  })

  output$n_total <- renderText(format(design()$n_total, big.mark = ","))

  output$n_detail <- renderUI({
    d <- design()
    parts <- paste0(format(d$n_per_arm, big.mark = ","), " per arm")
    en <- d$detail$expected_n_h1 %||% d$detail$expected_n_alt
    if (!is.null(en) && is.finite(en)) {
      parts <- paste0(parts, "; expected ", format(round(en), big.mark = ","),
                      " under the alternative")
    }
    p(class = "note", parts)
  })

  output$design_print <- renderPrint(print(design()))

  output$boundaries <- renderTable(gs_boundaries(design()), digits = 4)

  output$power_curve <- renderPlot({
    s <- current()
    control <- s$control_rate
    lower <- if (s$direction == "lower_is_better") control * 0.30 else control * 1.02
    upper <- if (s$direction == "lower_is_better") control * 0.98 else min(control * 2.5, 0.98)
    rates <- seq(lower, upper, length.out = 25)

    n <- vapply(rates, function(r) {
      out <- try(
        design_fixed(
          scenario(pack(), endpoint = input$endpoint,
                   control_rate = control, treatment_rate = r),
          alpha = input$alpha, power = input$power
        )$n_total,
        silent = TRUE
      )
      if (inherits(out, "try-error")) NA_real_ else out
    }, numeric(1))

    old <- graphics::par(mar = c(4.5, 5, 1, 1), bty = "n", las = 1)
    on.exit(graphics::par(old), add = TRUE)
    plot(rates, n, type = "l", lwd = 2, col = "#1f4e79", log = "y",
         xlab = "Treatment arm event rate", ylab = "Total sample size required")
    graphics::abline(v = s$treatment_rate, col = "#b03a2e", lty = 2)
    graphics::points(s$treatment_rate, design_fixed(
      s, alpha = input$alpha, power = input$power
    )$n_total, pch = 19, col = "#b03a2e")
  })

  output$citation <- renderText({
    published()$source %||% "No citation recorded for this endpoint."
  })
}

`%||%` <- function(x, y) if (is.null(x)) y else x

shinyApp(ui, server)
