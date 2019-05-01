
# times are in minutes
# concentrations are in nM

times <- c(0,1,5,15,30,45,60)
S_HF <- c(54847.11367, 50367.61578, 39989.17491, 27865.21746, 14233.57999, 6, 6)/1000
S_reg <- c(54847.11367, 53188.19125, 47128.91215, 38164.35909, 25979.01489, 22187.11001, 6441.779328)/1000

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

fit_MM_linear <- function(t, S, phi=1e-3) {
  # we cannot have any S = 0; these are replaced with \phi*S0
  S[S < phi*max(S)] <- phi*max(S)
  
  # Km*log(S0/S) + S0 - S = Vmax*t
  # log(S0/S)*Km - t*Vmax = S - S0
  n <- length(t)
  A <- matrix(0, nrow=n, ncol=2)
  A[ ,1] <- log(max(S)/S)
  A[ ,2] <- -t
  print(A)
  print(MASS::ginv(A))
  x <- MASS::ginv(A) %*% (S - max(S))
  print(x)
  
  return(list(
    Km = x[1],
    Vmax = x[2],
    S0 = max(S),
    St = solve_S(x[1], x[2], max(t), max(S)),
    t = t,
    S = S
  ))
}

fit_MM_analytical <- function(t, S, phi=1e-3, max_iter=100) {
  # we cannot have any S = 0; these are replaced with \phi*S0
  S[S < phi*max(S)] <- phi*max(S)
  
  # \beta_1 = Vmax, \beta_2 = Km
  n <- length(t)
  J <- matrix(0, nrow=n, ncol=2)
  J[ ,1] <- -t
  J[ ,2] <- log(max(S)/S)
  pinvJ <- MASS::ginv(J)
  beta = c(1,1)
  for (i in seq_len(max_iter)) {
    r <- max(S) - S + beta[2]*log(max(S)/S) - beta[1]*t
    beta <- beta - pinvJ %*% r
  }
  
  return(list(
    Km = beta[2],
    Vmax = beta[1],
    S0 = max(S),
    St = solve_S(beta[2], beta[1], max(t), max(S)),
    t = t,
    S = S
  ))
}

fit_MM <- function(t, S, phi=1e-3, nonnegative=FALSE) {
  MM_int_formula <- t ~ (Km*log(S0/S)+S0-S)/Vmax
  
  # we cannot have any S = 0; these are replaced with \phi*S0
  S[S < phi*max(S)] <- phi*max(S)
  
  data <- data.frame(
    t = t,
    S = S,
    S0 = max(S)
  )
  
  if (nonnegative) {
    fit <- nls(MM_int_formula, data=data, start=list(Km=1, Vmax=1))
  } else {
    fit <- nls(MM_int_formula, data=data, start=list(Km=1, Vmax=1), lower=c(0,0), algorithm="port")
  }
  
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

plot_MM_fit <- function(fit) {
  plot(fit$t, fit$S)
  S <- seq(from=fit$S0, to=fit$St, length.out=100)
  tplot <- (fit$Km*log(fit$S0/S)+fit$S0-S)/fit$Vmax
  lines(tplot, S)
}

plot_MM_comparison <- function(fit1, fit2) {
  plot(fit1$t, fit1$S, xlab="time [min]", ylab="substrate [nM]")
  S <- seq(from=fit1$S0, to=fit1$St, length.out=100)
  tplot <- (fit1$Km*log(fit1$S0/S)+fit1$S0-S)/fit1$Vmax
  lines(tplot, S)
  
  points(fit2$t, fit2$S, pch=2, col="red")
  S <- seq(from=fit2$S0, to=fit2$St, length.out=100)
  tplot <- (fit2$Km*log(fit2$S0/S)+fit2$S0-S)/fit2$Vmax
  lines(tplot, S, col="red")
}

plot_velocity_comparison <- function(fit1, fit2) {
  S <- seq(from=0, to=fit1$S0, length.out=100)
  v <- fit1$Vmax*S/(fit1$Km+S)
  plot(S, v, xlab="substrate [nM]", ylab="velocity [nM/min]")
}

fit1 <- fit_MM(times, S_HF)
fit2 <- fit_MM(times, S_reg)

fit1l <- fit_MM_linear(times, S_HF)
fit2l <- fit_MM_linear(times, S_reg)



