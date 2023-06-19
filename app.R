library(shiny)
library(tidyr)
library(ggplot2)
library(hash)
library(mcr)

make_ba <- function(dF, title, y_ax, x_ax, clia_pt) {
  
  mean_bias <- mean(dF$Bias)
  
  if (clia_pt > 0) {
    ggplot(data = dF, mapping = aes(x = Predicate, y = Bias)) +
      geom_point(aes(alpha = 0.3), pch = 21, fill = 'pink', color = 'black', show.legend = FALSE) +
      geom_hline(aes(yintercept = mean_bias, color = 'Mean Bias'), linewidth = 0.4) +
      geom_hline(aes(yintercept = clia_pt, color = 'CLIA Limit'), linetype = 'dashed', linewidth = 0.3) +
      geom_hline(aes(yintercept = -clia_pt, color = 'CLIA Limit'), linetype = 'dashed', linewidth = 0.3) +
      labs(title = title) +
      labs(y = y_ax, x = x_ax) +
      scale_color_manual(name = 'Legend',
                         breaks = c('Mean Bias', 'CLIA Limit'),
                         values = c('green', 'red')
      ) +
      xlim(0, 1.2*max(dF$Predicate)) +
      theme_light(base_size = 12)
  } else {
    ggplot(data = dF, mapping = aes(x = Predicate, y = Bias)) +
      geom_point(aes(alpha = 0.3), pch = 21, fill = 'pink', color = 'black', show.legend = FALSE) +
      geom_hline(aes(yintercept = mean_bias, color = 'Mean Bias'), linewidth = 0.4) +
      labs(title = title) +
      labs(y = y_ax, x = x_ax) +
      scale_color_manual(name = 'Legend',
                         breaks = 'Mean Bias',
                         values = 'green'
                         ) +
      xlim(0, 1.2*max(dF$Predicate)) +
      theme_light(base_size = 12)
  }
}
make_pb <- function(dF, title, y_ax, x_ax, clia_pt) {
  pb_reg <- mcreg(x = dF$Predicate, y = dF$Test, method.reg = 'PaBa')
  pearson_r <- round(cor(x = dF$Predicate, y = dF$Test, method = 'pearson'),2)
  slope <- round(pb_reg@para[2,1], 2)
  y_int <- round(pb_reg@para[1,1], 2)
  eq <- paste('Y = ', slope, 'X +', y_int)
  r_sq <- paste('R^2 = ', pearson_r)

  if (clia_pt > 0) {
    ggplot(data <- dF, mapping = aes(x = Predicate, y = Test)) +
      geom_point(aes(alpha = 0.3), pch =21, fill = 'pink', color = 'black', show.legend = FALSE) +
      geom_abline(aes(slope = slope, intercept = y_int, color = 'Regression Line'), linewidth = 0.4) +
      geom_abline(aes(slope = (1 + clia_pt/100), intercept = 0, color = 'CLIA Limit'), linetype = 'dashed', linewidth = 0.1) +
      geom_abline(aes(slope = (1 - clia_pt/100), intercept = 0, color = 'CLIA Limit'), linetype = 'dashed', linewidth = 0.1) +
      geom_abline(aes(slope = 1, intercept = 0, color = 'Identity Line'), linetype = 'dashed', linewidth = 0.1) +
      labs(title = title) +
      labs(y = y_ax, x = x_ax) +
      ylim(0, 1.2*max(dF$Predicate)) +
      xlim(0, 1.2*max(dF$Predicate)) +
      scale_color_manual(name = 'Legend',
                         breaks = c('Regression Line', 'Identity Line', 'CLIA Limit'),
                         values = c('green', 'black', 'red')
      ) +
      annotate(geom = 'text', x = as.integer(max(dF$Predicate)/1.11), y = as.integer(min(dF$Predicate)*2), label = eq) +
      annotate(geom = 'text', x = as.integer(max(dF$Predicate)/1.11), y = 0, label = r_sq) +
      theme_light(base_size = 12)
    
  } else {
    ggplot(data <- dF, mapping = aes(x = Predicate, y = Test)) +
      geom_point(aes(alpha = 0.3), pch =21, fill = 'pink', color = 'black', show.legend = FALSE) +
      geom_abline(aes(slope = slope, intercept = y_int, color = 'Regression Line'), linewidth = 0.4) +
      geom_abline(aes(slope = 1, intercept = 0, color = 'Identity Line'), linetype = 'dashed', linewidth = 0.1) +
      labs(title = title) +
      labs(y = y_ax, x = x_ax) +
      ylim(0, 1.2*max(dF$Predicate)) +
      xlim(0, 1.2*max(dF$Predicate)) +
      scale_color_manual(name = 'Legend',
                         breaks = c('Regression Line', 'Identity Line'),
                         values = c('green', 'black')
      ) +
      annotate(geom = 'text', x = as.integer(max(dF$Predicate)/1.11), y = as.integer(min(dF$Predicate)*2), label = eq) +
      annotate(geom = 'text', x = as.integer(max(dF$Predicate)/1.11), y = 0, label = r_sq) +
      theme_light(base_size = 12)
  }
}
hover_function <- function(data, hover){
  point <- nearPoints(data, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
  if (nrow(point) == 0) return(NULL)
  
  left_px <- hover$coords_css$x
  top_px <- hover$coords_css$y
  
  # create style property for tooltip
  # background color is set so tooltip is a bit transparent
  # z-index is set so we are sure tooltip will be on top
  style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                  "left:", left_px + 2, "px; top:", top_px + 2, "px;")
  
  # actual tooltip created as wellPanel
  wellPanel(
    style = style,
    p(HTML(paste0("<b> Point: </b>", rownames(point), "<br/>",
                  "<b> Test: </b>", point$Test, "<br/>",
                  "<b> Predicate: </b>", point$Predicate, "<br/>"
                  )))
  )
}
extract_nums <- function(text) {
  text <- gsub(" ", "", text)
  split <- strsplit(text, ",", fixed = FALSE)[[1]]
  as.numeric(split)
}
{analytes <- c(
  "Don't show",
  'Fibrinogen',
  'Hemoglobin',
  'Hematocrit',
  'Partial thromboplastin time',
  'Platelet count',
  'Prothrombin time',
  'Red blood cell count',
  'White blood cell count'
)
} #analytes
{clia_map <- hash()
clia_map[["Don't show"]] <- 0
clia_map[['Hemoglobin']] <- 7
clia_map[['Hematocrit']] <- 6
clia_map[['Red blood cell count']] <- 6
clia_map[['Platelet count']] <- 25
clia_map[['White blood cell count']] <- 15
clia_map[['Fibrinogen']] <- 20
clia_map[['Partial thromboplastin time']] <- 15
clia_map[['Prothrombin time']] <- 15
} #clia hash

ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = 'journal'),
  fluidRow(
    column(10,
           h1("Method comparison dashboard")
           )
  ),
  fluidRow(
    column(12,
           p('The purpose of this dashboard is to help people better visualize their data when running
             method comparison studies. Today, several regulatory bodies (ISO, CLIA, and CLSI) outline
             standardized methods to properly test medical applications built by industry. Because of 
             this standardization, industry players are held to the same expectations when developing 
             their particular product.'),
           p("Many today likely collect their data, then quickly run their own scripts to analyze that data.
              Since we're all essentially doing the same thing each with some small variation, this tool was
             built to be a low-touch, comprehensive and standardized way to assess our data in a flash while
             producing publication-ready figures."),
           p("Enjoy!")
    )
  ),
  fluidRow(
    column(10,
           fileInput("upload", "Please upload your data", buttonLabel = 'Upload*...', accept = '.csv'),
           tags$h6('*your data must be a .csv file with only 2 columns. One column labeled "Predicate" and the other "Test".'),
           tags$h6('Please wait a few seconds for the regression to complete. The passing bablok is computationally expensive.')
           )
    ),
  fluidRow(
    column(5,
           div(
             style = "position:relative",
             plotOutput("ba_tooltip",
                        hover = hoverOpts("plot_hover", delay = 100, delayType = "debounce")),
             uiOutput("hover_info_ba", style = "pointer-events: none")
           ),
           uiOutput('show_ba_title'),
           uiOutput('show_ba_x'),
           uiOutput('show_ba_y'),
           uiOutput('show_clia_pt_ba'),
           uiOutput('show_ba_analysis'),
           textOutput('ba_mean'),
           textOutput('ba_sd')
           ),
    column(5,
           div(
             style = "position:relative",
             plotOutput("pb_tooltip",
                        hover = hoverOpts("plot_hover", delay = 100, delayType = "debounce")),
             uiOutput("hover_info_pb", style = "pointer-events: none")
           ),
           uiOutput('show_pb_title'),
           uiOutput('show_pb_x'),
           uiOutput('show_pb_y'),
           uiOutput('show_clia_pt_pb'),
           uiOutput('show_pb_analysis'),
           uiOutput('enter_md_levels'),
           uiOutput('show_pb_param_analysis'),
           verbatimTextOutput('pb_param_analysis'),
           uiOutput('show_pb_md_analysis'),
           verbatimTextOutput('pb_md_analysis')
           )

  )
)

server <- function(input, output, session) {
  
  data <- reactive({
    loader <- req(input$upload)
    if (is.null(loader)) {
      return(NULL)
    }
    read.csv(loader$datapath)
    ext <- tools::file_ext(input$upload$name)
    switch(ext,
           csv = vroom::vroom(input$upload$datapath, delim = ","),
           validate("Invalid file... please upload a .csv file")
    )
  })
  
  data_trans <- reactive({
    
    if (sum(is.na(data())) > 0) {
      data_new <- na.omit(data())
    } else {
      data_new <- data()
    }
    Bias <- ((data_new$Test - data_new$Predicate) / data_new$Predicate) * 100
    data_new <- cbind(data_new, Bias)
    
  })
  
  output$show_ba_title <- renderUI({
    if (is.null(data())) return()
    textInput('ba_title', "Title", value = 'Bland-Altman Analysis')
  })
  ba_title <- reactive({
    if (is.null(input$ba_title)) return()
    input$ba_title
  })
  
  output$show_ba_x <- renderUI({
    if (is.null(data())) return()
    textInput('ba_x', "X-axis label", value = 'Predicate XS-1000i (K/µL)')
  })
  ba_x <- reactive({
    if (is.null(input$ba_x)) return()
    input$ba_x
  })
  
  output$show_ba_y <- renderUI({
    if (is.null(data())) return()
    textInput('ba_y', 'Y-axis label', value = 'Proportional Bias (%)')
  })
  ba_y <- reactive({
    if (is.null(input$ba_y)) return()
    input$ba_y
  })
  
  output$show_ba_analysis <- renderUI({
    if (is.null(data())) return()
    h5('Bland Altman Analysis')
  })
  
  output$show_clia_pt_ba <- renderUI({
    if (is.null(data())) return()
    selectInput('clia_pt_ba', 'Show CLIA proficiency testing limits (as of 06/10/2023) for', analytes)
  })
  
  pt_limit_ba <- eventReactive(input$clia_pt_ba, {
    clia_map[[input$clia_pt_ba]]
  })
  
  output$ba_tooltip <- renderPlot({
    make_ba(data_trans(), title = ba_title(), x_ax = ba_x(), y_ax = ba_y(), clia_pt = pt_limit_ba())
  })
  
  output$hover_info_ba <- renderUI({
    hover <- input$plot_hover
    hover_function(data_trans(), hover)
  })
  
  output$ba_mean <- renderText({
    mean_bias <- mean(data_trans()$Bias)
    paste0('Mean proportional bias: ', round(mean_bias,2))
  })
  
  output$ba_sd <- renderText({
    sd_bias <- sd(data_trans()$Bias)
    paste0('2 standard deviations of proportional bias: ', round(sd_bias,2))
  })
  
  output$show_pb_title <- renderUI({
    if (is.null(data())) return()
    textInput('pb_title', "Title", value = 'Passing Bablok Regression')
  })
  pb_title <- reactive({
    if (is.null(input$pb_title)) return()
    input$pb_title
  })
  
  output$show_pb_x <- renderUI({
    if (is.null(data())) return()
    textInput('pb_x', "X-axis label", value = 'Predicate XS-1000i (u)')
  })
  pb_x <- reactive({
    if (is.null(input$pb_x)) return()
    input$pb_x
  })
  
  output$show_pb_y <- renderUI({
    if (is.null(data())) return()
    textInput('pb_y', 'Y-axis label', value = 'Athelas Home (u)')
  })
  pb_y <- reactive({
    if (is.null(input$pb_y)) return()
    input$pb_y
  })
  
  output$show_clia_pt_pb <- renderUI({
    if (is.null(data())) return()
    selectInput('clia_pt_pb', 'Show CLIA proficiency testing limits (as of 06/10/2023) for', analytes)
  })
  
  pt_limit_pb <- eventReactive(input$clia_pt_pb, {
    clia_map[[input$clia_pt_pb]]
  })
  
  output$show_pb_analysis <- renderUI({
    if (is.null(data())) return()
    h5('Passing Bablok Analysis')
  })
  
  output$pb_tooltip <- renderPlot({
    make_pb(data_trans(), title = pb_title(), x_ax = pb_x(), y_ax = pb_y(), clia_pt = pt_limit_pb())
  })
  
  output$hover_info_pb <- renderUI({
    hover <- input$plot_hover
    hover_function(data_trans(), hover)
  })
  
  output$enter_md_levels <- renderUI({
    if (is.null(data())) return()
    textInput('md_levels', 'Enter medical decision levels (if entering multiple, separate with a comma)', value = '0,1,2')
  })
  
  output$show_pb_param_analysis <- renderUI({
    if (is.null(data())) return()
    h5('Parameter Estimates')
  })
  
  output$pb_param_analysis <- renderPrint({
    pb_reg <- mcreg(x = data_trans()$Predicate, y = data_trans()$Test, method.reg = 'PaBa')
    pb_reg@para
  })
  
  output$show_pb_md_analysis <- renderUI({
    if (is.null(data())) return()
    h5('Systematic Bias at Specified Medical Decision Levels')
  })
  
  output$pb_md_analysis <- renderPrint({
    pb_reg <- mcreg(x = data_trans()$Predicate, y = data_trans()$Test, method.reg = 'PaBa')
    calcBias(pb_reg, x.levels = extract_nums(input$md_levels))
  })
  
}

shinyApp(ui, server)