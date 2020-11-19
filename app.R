#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

t <- c(0, 1, 5, 15, 30, 45, 60)
S1 <- c(54.8, 50.3, 39.9, 27.9, 14.2, 0, 0)
S2 <- c(54.8, 53.1, 47.1, 38.1, 26, 22.1, 6.4)

t <- c(0, 15, 30, 45, 55, 59, 60)
S1 <- c(399, 256.2, 225.4, 222, 166.5, 66.7, 56.4)
S2 <- c(251.9, 208.6, 205.4, 204.9, 123.1, 86.4, 56.4)

# Define UI for application that draws a histogram
ui <- splitLayout(cellWidths = c("50%", "50%"),
         fluidPage(
           fluidRow(h3("  Enzyme Kinetics Fitting Tool")),
           fluidRow(column(width=4, textInput("t1", "Time [min]", t[1])),
                   column(width=4, textInput("x1", "Substrate 1 [nM]", S1[1])),
                   column(width=4, textInput("y1", "Substrate 2 [nM]", S2[1]))),
           
           fluidRow(column(width=4, textInput("t2", NULL, t[2])),
                    column(width=4, textInput("x2", NULL, S1[2])),
                    column(width=4, textInput("y2", NULL, S2[2]))),
           
           fluidRow(column(width=4, textInput("t3", NULL, t[3])),
                    column(width=4, textInput("x3", NULL, S1[3])),
                    column(width=4, textInput("y3", NULL, S2[3]))),
           
           fluidRow(column(width=4, textInput("t4", NULL, t[4])),
                    column(width=4, textInput("x4", NULL, S1[4])),
                    column(width=4, textInput("y4", NULL, S2[4]))),
           
           fluidRow(column(width=4, textInput("t5", NULL, t[5])),
                    column(width=4, textInput("x5", NULL, S1[5])),
                    column(width=4, textInput("y5", NULL, S2[5]))),
           
           fluidRow(column(width=4, textInput("t6", NULL, t[6])),
                    column(width=4, textInput("x6", NULL, S1[6])),
                    column(width=4, textInput("y6", NULL, S2[6]))),
           
           fluidRow(column(width=4, textInput("t7", NULL, t[7])),
                    column(width=4, textInput("x7", NULL, S1[7])),
                    column(width=4, textInput("y7", NULL, S2[7]))),
           
           fluidRow(column(width=4, actionButton("goButton", "Fit Data")),
                    column(width=8, radioButtons("solver", "Solver:", c("Nonlinear" = "nonlinear", "Absolute" = "abs")))),
           
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
  
  solvef <- reactive({
    switch(input$solver,
           nonlinear = fit_MM,
           linear = fit_MM_linear,
           abs = fit_MM_abs)
  })
  
  reportData <- function(data) {
    fit1 <- solvef()(data()$t, data()$S1)
    fit2 <- solvef()(data()$t, data()$S2)
    fmt <- function(x) sprintf("%.3f", x)
    paste0(
      "Enzyme 1\n     Km = ", fmt(fit1$Km), " nM\n   Vmax = ", fmt(fit1$Vmax), " nM/min\n\n",
      "Enzyme 2\n     Km = ", fmt(fit2$Km), " nM\n   Vmax = ", fmt(fit2$Vmax), " nM/min\n"
    )
  }
  
  output$result <- shiny::renderText({reportData(data())})
  
   output$distPlot <- renderPlot({
      plot_MM_comparison(solvef()(data()$t, data()$S1), solvef()(data()$t, data()$S2))
   })
   
   output$velocityPlot <- renderPlot({
     plot_velocity_comparison(solvef()(data()$t, data()$S1), solvef()(data()$t, data()$S2))
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

find_start <- function(t, S) {
  Km <- (max(S) - min(S))/2
  n <- length(t)
  X = matrix(nrow=n, ncol=3, data=1)
  X[ ,2] <- t
  X[ ,3] <- t^2
  beta = MASS::ginv(X) %*% S
  t_Km <- max(
    (-beta[2] + sqrt(beta[2]^2 - 4*beta[2]*beta[1]))/(2*beta[3]),
    (-beta[2] - sqrt(beta[2]^2 - 4*beta[2]*beta[1]))/(2*beta[3])
  )
  Vmax <- beta[2] + 2*beta[3]*t_Km
  if (is.nan(Vmax) || Vmax < 0) {
    Vmax <- max(abs(S[2:n] - S[1:(n-1)]) / (t[2:n] - t[1:(n-1)]))
  }
  if (is.na(Vmax) || Vmax < 0) {
    Vmax <- 1
  }
  list(Km=Km, Vmax=Vmax)
}

fit_MM <- function(t, S, phi=1e-3, exclude_zeros=FALSE) {
  MM_int_formula <- t ~ (Km*log(S0/S)+S0-S)/Vmax
  
  # we cannot have any S = 0; these are replaced with \phi*S0
  if (exclude_zeros) {
    ex <- S < phi*max(S)
    S <- S[!ex]
    t <- t[!ex]
  } else {
    S[S < phi*max(S)] <- phi*max(S)
  }
  
  data <- data.frame(
    t = t,
    S = S,
    S0 = max(S)
  )
  
  #fit <- nls(MM_int_formula, data=data, start=list(Km=1, Vmax=1))
  fit <- nls(MM_int_formula, data=data, start=find_start(t,S))
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

fit_MM_abs <- function(t, S, phi=1e-3, exclude_zeros=FALSE) {
  if (exclude_zeros) {
    ex <- S < phi*max(S)
    S <- S[!ex]
    t <- t[!ex]
  } else {
    S[S < phi*max(S)] <- phi*max(S)
  }
  
  S0 <- max(S)
  
  f <- function(beta) sum(abs(beta[1]*log(S0/S)+S0-S - beta[2]*t))
  st <- find_start(t,S)
  beta_0 <- c(st$Km, st$Vmax)
  sol <- optim(beta_0, f, lower=c(0,0), method="L-BFGS-B")
  beta <- sol$par
  
  return(list(
    Km = beta[1],
    Vmax = beta[2],
    S0 = S0,
    St = solve_S(beta[1], beta[2], max(t), max(S)),
    t = t,
    S = S
  ))
}

fit_MM_linear <- function(t, S, phi=1e-3) {
  # we cannot have any S = 0; these are replaced with \phi*S0
  S[S < phi*max(S)] <- 1
  
  S0 = max(S)
  n <- length(t)
  X = matrix(nrow=n, ncol=2)
  X[ ,1] <- log(S0/S)
  X[ ,2] <- -t
  y <- S - S0
  beta <- MASS::ginv(X) %*% y
  
  return(list(
    Km = beta[1],
    Vmax = beta[2],
    S0 = S0,
    St = solve_S(beta[1], beta[2], max(t), max(S)),
    t = t,
    S = S
  ))
}

plot_MM_comparison <- function(fit1, fit2) {
  ylim <- c(
    min(min(fit1$S), min(fit2$S)), 
    max(fit1$S0, fit2$S0)
  )
  plot(fit1$t, fit1$S, xlab="time [min]", ylab="substrate [nM]", ylim=ylim)
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

# t <- c(0,1,5,15,30,45,60)
# S_b1 <- c(57.24,
#           61.49,
#           56.58,
#           53.89,
#           38.47,
#           41.64,
#           27.74)
# S_b2 <- c(57.24,
#           41.06,
#           20.25,
#           24.13,
#           5.45,
#           0,
#           0)
# S_c1 <- c(54.9,
#           53.1,
#           46.4,
#           43.1,
#           33.5,
#           36.2,
#           23.9)
# S_c2 <- c(54.9,
#           37.8,
#           22.4,
#           23.2,
#           7.17,
#           0,
#           0)
# 
# print(fit_MM(t, S_b1))
# print(fit_MM(t[1:6], S_b2[1:6]))
# print(fit_MM(t, S_c1))
# print(fit_MM(t[1:6], S_c2[1:6]))
# 
# print(fit_MM_abs(t, S_b1))
# print(fit_MM_abs(t[1:6], S_b2[1:6]))
# print(fit_MM_abs(t, S_c1))
# print(fit_MM_abs(t[1:6], S_c2[1:6]))
