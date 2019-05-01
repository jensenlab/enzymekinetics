#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- splitLayout(cellWidths = c("50%", "50%"),
         fluidPage(
           fluidRow(h3("  Enzyme Kinetics Fitting Tool")),
           fluidRow(column(width=4, textInput("t1", "Time [min]",0)),
                   column(width=4, textInput("x1", "Substrate 1 [nM]",54.8)),
                   column(width=4, textInput("y1", "Substrate 2 [nM]",54.8))),
           
           fluidRow(column(width=4, textInput("t2", NULL,1)),
                    column(width=4, textInput("x2", NULL,50.3)),
                    column(width=4, textInput("y2", NULL,53.1))),
           
           fluidRow(column(width=4, textInput("t3", NULL,5)),
                    column(width=4, textInput("x3", NULL,39.9)),
                    column(width=4, textInput("y3", NULL,47.1))),
           
           fluidRow(column(width=4, textInput("t4", NULL,15)),
                    column(width=4, textInput("x4", NULL,27.9)),
                    column(width=4, textInput("y4", NULL,38.1))),
           
           fluidRow(column(width=4, textInput("t5", NULL,30)),
                    column(width=4, textInput("x5", NULL,14.2)),
                    column(width=4, textInput("y5", NULL,26.0))),
           
           fluidRow(column(width=4, textInput("t6", NULL,45)),
                    column(width=4, textInput("x6", NULL,0)),
                    column(width=4, textInput("y6", NULL,22.1))),
           
           fluidRow(column(width=4, textInput("t7", NULL,60)),
                    column(width=4, textInput("x7", NULL,0)),
                    column(width=4, textInput("y7", NULL,6.4))),
           
           fluidRow(column(width=12, actionButton("goButton", "Fit Data"))),
           
           fluidRow(column(width=12, verbatimTextOutput(outputId = "result")))
           
         ),
         fluidPage(
           fluidRow(plotOutput("distPlot")),
           fluidRow(plotOutput("velocityPlot"))
         )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  data <- eventReactive(input$goButton, {
    n <- 7
    d <- data.frame(t=numeric(n), S1=numeric(n), S2 = numeric(n))
    for (i in 1:n) {
      d$t[i] <- as.numeric(input[[paste0("t",i)]])
      d$S1[i] <- as.numeric(input[[paste0("x",i)]])
      d$S2[i] <- as.numeric(input[[paste0("y",i)]])
    }
    d
  })
  
  reportData <- function(data) {
    fit1 <- fit_MM(data()$t, data()$S1)
    fit2 <- fit_MM(data()$t, data()$S2)
    fmt <- function(x) sprintf("%.3f", x)
    paste0(
      "Enzyme 1\n     Km = ", fmt(fit1$Km), " nM\n   Vmax = ", fmt(fit1$Vmax), " nM/min\n\n",
      "Enzyme 2\n     Km = ", fmt(fit2$Km), " nM\n   Vmax = ", fmt(fit2$Vmax), " nM/min\n"
    )
  }
  
  output$result <- shiny::renderText({reportData(data())})
  
   output$distPlot <- renderPlot({
      plot_MM_comparison(fit_MM(data()$t, data()$S1), fit_MM(data()$t, data()$S2))
   })
   
   output$velocityPlot <- renderPlot({
     plot_velocity_comparison(fit_MM(data()$t, data()$S1), fit_MM(data()$t, data()$S2))
   })
}

solve_S <- function(Km, Vmax, t, S0, max_iter=100) {
  # for initial guess, solve for \phi = 1 - (rxn extent)
  # i.e. S(t) = \phi*S0
  # assume that \phi << 1 near the endpoint
  phi <- exp((S0 - Vmax*t)/Km)
  S <- phi*S0
  for (i in seq_len(max_iter)) {
    S <- S - (Km*log(S0/S) + S0 - S - Vmax*t)/(-Km/S-1)
  }
  return(S)
}

fit_MM <- function(t, S, phi=1e-3) {
  MM_int_formula <- t ~ (Km*log(S0/S)+S0-S)/Vmax
  
  # we cannot have any S = 0; these are replaced with \phi*S0
  S[S < phi*max(S)] <- phi*max(S)
  
  data <- data.frame(
    t = t,
    S = S,
    S0 = max(S)
  )
  
  fit <- nls(MM_int_formula, data=data, start=list(Km=1, Vmax=1))
  Km <- unname(coef(fit)["Km"])
  Vmax <- unname(coef(fit)["Vmax"])
  
  return(list(
    Km = Km,
    Vmax = Vmax,
    S0 = data$S0[1],
    St = solve_S(Km, Vmax, max(t), max(S)),
    t = t,
    S = S
  ))
}

plot_MM_comparison <- function(fit1, fit2) {
  plot(fit1$t, fit1$S, xlab="time [min]", ylab="substrate [nM]")
  S <- seq(from=fit1$S0, to=fit1$St, length.out=100)
  tplot <- (fit1$Km*log(fit1$S0/S)+fit1$S0-S)/fit1$Vmax
  lines(tplot, S, lwd=2)
  
  points(fit2$t, fit2$S, pch=2, col="red")
  S <- seq(from=fit2$S0, to=fit2$St, length.out=100)
  tplot <- (fit2$Km*log(fit2$S0/S)+fit2$S0-S)/fit2$Vmax
  lines(tplot, S, lwd=2, col="red")
}

plot_velocity_comparison <- function(fit1, fit2) {
  S <- seq(from=0, to=fit1$S0, length.out=100)
  v <- fit1$Vmax*S/(fit1$Km+S)
  plot(S, v, type="l", lwd=3, xlab="substrate [nM]", ylab="velocity [nM/min]", ylim=c(0,max(fit1$Vmax, fit2$Vmax)))
  abline(h=fit1$Vmax, lty="dashed", lwd=2)
  lines(c(fit1$Km, fit1$Km), c(0, fit1$Vmax*fit1$Km/(fit1$Km+fit1$Km)), lty="dotted", lwd=2)
  
  v <- fit2$Vmax*S/(fit2$Km+S)
  lines(S, v, lwd=3, col="red")
  abline(h=fit2$Vmax, lty="dashed", lwd=2, col="red")
  lines(c(fit2$Km, fit2$Km), c(0, fit2$Vmax*fit2$Km/(fit2$Km+fit2$Km)), lty="dotted", col="red", lwd=2)
}

# Run the application 
shinyApp(ui = ui, server = server)

