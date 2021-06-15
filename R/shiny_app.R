#' Plot spatial maps from FLAMe system
#'
#' Run this function (with no arguments) to launch a shiny application that allows for visualizing spatial measurements for user-chosen dates and variables.
#'
#' @return opens a shiny app
#' @import shiny
#' @import ggplot2
#' @export
#'
#' @examples
#' ## Only run this example in interactive R sessions
#' if (interactive()) {
#'   plot_FLAMe_maps()
#' }
plot_FLAMe_maps <- function(){

# Define UI for application that draws a histogram
  ui <- fluidPage(

    # Application title
    titlePanel("Spatial Data"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
      sidebarPanel(
        checkboxGroupInput("lakes", label=h3("Select Lake(s)"), choices=c(Peter="R", Paul="L")),
        uiOutput("getVariable"),
        uiOutput("getDates")
      ),

      # Show a plot of the generated distribution
      mainPanel(
        plotOutput("flamePlot", height="700px")
      )
    )
  )

  # Define server logic required to draw a histogram
  server <- function(input, output) {
    # read in data, setup variables for UI
    data = flame_data
    unique_lakesDates = unique(data[, c("Lake", "Date")])
    dates =  unique(unique_lakesDates[, "Date"] %>% dplyr::pull())
    tau_variables = c(`BGA (ug/L)` = "BGApc_ugL_tau", `D.O. sat (%)` = "ODO_percent_tau", `pH` = "pH_tau")
    # create UI elements for selecting data/variables
    output$getDates = renderUI({
      checkboxGroupInput("selectedDates", label = h3("Choose Date(s)"), choices = dates, inline=TRUE)
    })
    output$getVariable = renderUI({
      selectInput("selectedVariable", label=h5("Select Variable"), choices = names(tau_variables), selected=1, multiple=FALSE)
    })
    # create the plot
    plotData = reactive({
      data %>%
        dplyr::filter(Lake %in% input$lakes & Date %in% as.Date(input$selectedDates))
    })

    output$flamePlot <- renderPlot({
      req(nrow(plotData()) > 0)
      ggplot(plotData(), aes(x=longitude, y=latitude, color=get(tau_variables[input$selectedVariable]))) +
        facet_wrap(~Date) +
        geom_point(size = 2 + (length(input$selectedDates))^-1) +
        theme_bw() +
        viridis::scale_color_viridis(option="viridis") +
        labs(color = input$selectedVariable)
    })
  }

  # Run the application
  shinyApp(ui = ui, server = server)
}

