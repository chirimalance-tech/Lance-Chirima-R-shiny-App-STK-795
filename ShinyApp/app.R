# Attach all required packages explicitly
library(shiny)
library(shinythemes)
library(DT)
library(MASS)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(nlme)
library(AER)
library(gganimate)
library(gifski)
library(COUNT)
library(gamlss)
library(gamlss.dist)
library(VGAM)
library(DDPM)
library(COMPoissonReg)
library(png)
library(magick)

if (requireNamespace("shinythemes", quietly = TRUE)) {
  theme <- shinythemes::shinytheme("superhero")
} else {
  theme <- shiny::shinytheme(NULL)
}



dcnbinom<-function(x, mu, alpha, delta, eta, log=F){
  log_term1 = log(1-delta) + dnbinom(x, size = 1/alpha, mu = mu, log = T)
  log_term2 = log(delta) + dnbinom(x, size = 1/(alpha*eta), mu = mu, log = T)
  d = (exp(log_term1) + exp(log_term2))
  if (log == F) {
    return(d)
  } else {
    return(log(d))
  }
}


ml.ex.cnb <- function(formula, alpha.formula=~1, delta.formula=~1, eta.formula=~1, data, start = NULL, method = "BFGS", reltol=1e-10, maxit=1000, hessian=T) {
  mf <- model.frame(formula, data)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  betaX <- model.matrix(formula, data = data) # design matrix for mean
  alphaU <- model.matrix(alpha.formula, data = data) # design matrix for sigma
  deltaV <- model.matrix(delta.formula, data = data) # design matrix for delta
  etaZ <- model.matrix(eta.formula, data = data)
  nb2.reg.ml <- function(b.hat, X, U, V, Z, y) {
    beta.hat <- b.hat[1:ncol(X)] #coefficients for the mean (mu)
    theta.hat <- b.hat[(ncol(X)+1):(ncol(X)+ncol(U))] #coefficients for alpha
    gamma.hat <-b.hat[(ncol(X)+ncol(U)+1):(ncol(X)+ncol(U)+ncol(V))]
    lambda.hat <- b.hat[(ncol(X)+ncol(U)+ncol(V)+1):(ncol(X)+ncol(U)+ncol(V)+ncol(Z))]

    xb.hat <- X %*% beta.hat # mean regression
    mu.hat <- exp(xb.hat)

    ua.hat <- U %*% theta.hat # alpha regression
    alpha.hat <- exp(ua.hat)

    vg.hat <- V %*% gamma.hat # delta regression
    delta.hat <- exp(vg.hat) / (1 + exp(vg.hat))

    zl.hat <- Z %*% lambda.hat # eta regression
    eta.hat <- exp(zl.hat)

    ll <- sum(dcnbinom	(x = y, mu = mu.hat, alpha = alpha.hat, delta = delta.hat, eta=eta.hat, log=T))
    return(ll)
  }

  if (is.null(start)) {#initial param
    nb <- MASS::glm.nb(formula = formula, data = na.omit(data), control = glm.control(maxit = 5000))
    nb$loglike
    beta.hat <- nb$coefficients
    if (anyNA(beta.hat)) {
      warning("NA detected in glm.nb coefficients. Replacing with small fallback values.")
      beta.hat[is.na(beta.hat)] <- runif(sum(is.na(beta.hat)))
    }
    theta.hat <- rep(log(0.1), ncol(alphaU))
    gamma.hat = rep(log(0.0055/(1-0.0055)),ncol(deltaV))
    lambda.hat = rep(log(0.006), ncol(etaZ))


    start <- c(beta.hat, theta.hat, gamma.hat,lambda.hat)

  }

  print(start)
  summary(start)
  any(!is.finite(start))

  fit <- optim(par = start,
               fn = nb2.reg.ml,
               X = betaX,
               U = alphaU,
               V = deltaV,
               Z = etaZ,
               y = y,
               method = method,
               control = list(fnscale = -1, maxit = maxit, reltol = reltol),
               hessian = hessian
  )
  if (fit$convergence != 0) warning("Optimization may not have converged.")

  beta.hat <- fit$par[1:ncol(betaX)] #coefficients for the mean (mu)
  theta.hat <- fit$par[(ncol(betaX)+1):(ncol(betaX)+ncol(alphaU))] #coefficients for alpha
  gamma.hat <-fit$par[(ncol(betaX)+ncol(alphaU)+1):(ncol(betaX)+ncol(alphaU)+ncol(deltaV))]
  lambda.hat <- fit$par[(ncol(betaX)+ncol(alphaU)+ncol(deltaV)+1):(ncol(betaX)+ncol(alphaU)+ncol(deltaV)+ncol(etaZ))]

  xb.hat <- betaX %*% beta.hat # mean regression
  mu.hat <- exp(xb.hat)

  ua.hat <- alphaU %*% theta.hat # alpha regression
  alpha.hat <- exp(ua.hat)

  vg.hat <- deltaV %*% gamma.hat # delta regression
  delta.hat <- exp(vg.hat) / (1 + exp(vg.hat))

  zl.hat <- etaZ %*% lambda.hat # eta regression
  eta.hat <- exp(zl.hat)


  lc=sum(dcnbinom(y, mu=mu.hat, alpha = alpha.hat, delta=delta.hat, eta = eta.hat, log = T))
  # loglik[i]=lc

  if (hessian==T)
  {cov.mat <- tryCatch(solve(-fit$hessian), error = function(e) MASS::ginv(-fit$hessian))
  std.errors <- sqrt(diag(cov.mat))
  tvalue <- fit$par/std.errors
  pval <- 2*pt(abs(tvalue),df=nrow(data)-length(fit$par),lower.tail = F)
  results <- data.frame(
    Estimate = c(fit$par),
    `Std.Error` = std.errors,
    `t.value` = tvalue,
    `P.t` = round(pval,5)
  )}
  if(hessian==F){
    results <- data.frame(
      Estimate = c(fit$par))
  }
  rownames(results) <- c(paste0("beta_", colnames(betaX)), paste0("theta_", colnames(alphaU)), paste0("gamma_", colnames(deltaV)), paste0("lambda_", colnames(etaZ)))


  AIC <- -2*lc+nrow(results)*2
  BIC <- -2*lc+nrow(results)*log(nrow(data))
  HQIC <- -2*lc +nrow(results)*2*log(log(nrow(data)))

  return(list(
    results = results,
    beta = beta.hat,
    theta = theta.hat,
    gamma = gamma.hat,
    lambda = lambda.hat,
    mu = mu.hat,
    alpha = alpha.hat,
    delta = delta.hat,
    eta = eta.hat,
    X = betaX,
    U = alphaU,
    V = deltaV,
    Z = etaZ,
    y = y,
    loglike = lc,
    AIC = AIC,
    BIC = BIC,
    HQIC = HQIC
  ))
}



ui <- navbarPage(
  id = "navbar",
  title = "Contaminated Negative Binomial Regression for Modeling Count Data",
  theme = theme,

  # Tab 1: Background Info
  # Tab 1: Home
  tabPanel(
    "Home",
    fluidPage(
      tags$head(
        tags$style(HTML("
    #home_background {
     background-image: url('background.jpg');
     background-size: cover;
     background-repeat: no-repeat;
     background-position: center center;
     min-height: 100vh;
     display: flex;
     flex-direction: column;
     align-items: center;
     justify-content: flex-start;
     text-align: center;
     color: white;
     padding-top: 50px;
    }
    #home_background::before {
     content: '';
     position: absolute;
     top: 0; left: 0; right: 0; bottom: 0;
     background: rgba(0, 0, 0, 0.4);
     z-index: 1;
    }
    #home_background > * {
     position: relative;
     z-index: 2;
    }
    .home-btn {
     width: 200px;
     height: 100px;
     font-size: 16px;
     margin: 20px;
     background-color: #2c3e50;
     color: white;
     border-radius: 15px;
     border: 2px solid #ffffff;
     box-shadow: 3px 3px 8px rgba(0,0,0,0.3);
    }
    .home-btn:hover {
     background-color: #34495e;
     cursor: pointer;
    }
    .home-row {
     display: flex;
     justify-content: center;
     flex-wrap: wrap;
    }
   "))
      ),

      div(
        id = "home_background",
        h2("Welcome to the cNB App!"),
        p("Click any of the sections below to navigate:"),
        div(
          class = "home-row",
          actionButton("go_var", "cNB Interactive Visuals", class = "home-btn"),
          # actionButton("go_sens", "Sensitivity Study", class = "home-btn"),
          actionButton("go_data", "Data Application", class = "home-btn"),
          actionButton("go_bg", "Background Info", class = "home-btn")
        ),
        p("In a couple of seconds a dynamic plot will appear below:"),
        br(),
        actionButton("show_animation", "Generate Animation", class = "home-btn"),
        br(),
        imageOutput("animated_plot", width = 800, height = 500)
      )
    )
  )
  ,

  # Tab 2: cNB Vsiuals
  tabPanel("cNB Interactive Visuals",
           sidebarLayout(
             sidebarPanel(
               selectInput("plotType", "Choose Plot Type:",
                           choices = c(
                             "Variance vs μ (varying δ)" = "mu_delta",
                             "Variance vs μ (varying η)" = "mu_eta",
                             "Variance vs α (varying δ)" = "alpha_delta",
                             "Variance vs α (varying η)" = "alpha_eta"
                           )
               ),
               uiOutput("paramControls")
             ),
             mainPanel(
               tabPanel(withMathJax(),
                        h3("Overview"),
                        p("This tab examines the behavior of the contaminated Negative Binomial Distribution (cNB-D) when parameters are allowed to change"),
                        plotOutput("variancePlot", height = "600px")
               )
             )
           )
  ),

 

  tabPanel("Data Application",
           sidebarLayout(
             sidebarPanel(
               h4("Dataset"),
               radioButtons("dataset_choice", "Choose dataset:",
                            choices = c("RecreationDemand (AER)" = "recreation",
                                        "NMES1988 (AER)" = "nmes",
                                        "Upload CSV" = "upload")),
               conditionalPanel(
                 condition = "input.dataset_choice == 'upload'",
                 fileInput("file_upload", "Upload CSV file", accept = ".csv")
               ),
               actionButton("load_data", "Load dataset", class = "btn-success"),
               hr(),
               uiOutput("response_ui"),
               uiOutput("predictors_ui"),
               hr(),
               uiOutput("alpha_ui"),
               uiOutput("delta_ui"),
               uiOutput("eta_ui"),
               hr(),
               actionButton("run_model", "Run model", class = "btn-success"),
               br(), br(),
               textOutput("status")
             ),

             mainPanel(
               tabsetPanel(id = "main_tabs",
                           tabPanel("Preview",
                                    hr(),
                                    h3("DATA APPLICATION INSTRUCTIONS"),
                                    p("1) Select/Upload a Dataset."),
                                    p('    2) Press the "Load dataset" Button.'),
                                    p('        3) Select the covariate/(s) of your choice for each parameter.'),
                                    p('            4) "Run model" and veiw the "Results" and "Model Comparison".'),
                                    hr(),
                                    h4("Formulas used"),
                                    verbatimTextOutput("formula_text"),
                                    hr(),
                                    h4("Head"),
                                    verbatimTextOutput("head_out"),
                                    hr(),
                                    h4("Summary"),
                                    verbatimTextOutput("summary_out"),
                                    hr()
                           ),
                           tabPanel(
                             "Results",
                             hr(),

                             uiOutput("results_table_ui"),

                             tableOutput("info_table"),

                             hr(),

                             p("Some link functions to help with the interpretation of parameter estimates"),

                             # --- Two-column LaTeX layout ---
                             fluidRow(
                               column(
                                 6,
                                 withMathJax(
                                   helpText('
          \\[
          \\begin{align*}
            g_1(\\mu(\\boldsymbol{x};\\boldsymbol{\\beta})) &= \\log(\\mu(\\boldsymbol{x};\\boldsymbol{\\beta})) = \\boldsymbol{\\tilde{x}}^\\top \\boldsymbol{\\beta}, \\\\
            g_2(\\alpha(\\boldsymbol{u};\\boldsymbol{\\theta})) &= \\log(\\alpha(\\boldsymbol{u};\\boldsymbol{\\theta})) = \\boldsymbol{\\tilde{u}}^\\top \\boldsymbol{\\theta}, \\\\
            g_3(\\delta(\\boldsymbol{v};\\boldsymbol{\\gamma})) &= \\text{logit}(\\delta(\\boldsymbol{v};\\boldsymbol{\\gamma})) = \\boldsymbol{\\tilde{v}}^\\top \\boldsymbol{\\gamma}, \\\\
            g_4(\\eta(\\boldsymbol{z};\\boldsymbol{\\lambda})) &= \\log(\\eta(\\boldsymbol{z};\\boldsymbol{\\lambda})) = \\boldsymbol{\\tilde{z}}^\\top \\boldsymbol{\\lambda}.
          \\end{align*}
          \\]
        ')
                                 )
                               ),

                               column(
                                 6,
                                 withMathJax(
                                   helpText('
          \\[
          \\begin{align*}
            \\mu(\\boldsymbol{x};\\boldsymbol{\\beta}) &= g_1^{-1}(\\boldsymbol{\\tilde{x}}^\\top \\boldsymbol{\\beta}) = e^{\\boldsymbol{\\tilde{x}}^\\top \\boldsymbol{\\beta}}, \\\\
            \\alpha(\\boldsymbol{u};\\boldsymbol{\\theta}) &= g_2^{-1}(\\boldsymbol{\\tilde{u}}^\\top \\boldsymbol{\\theta}) = e^{\\boldsymbol{\\tilde{u}}^\\top \\boldsymbol{\\theta}} + 1, \\\\
            \\delta(\\boldsymbol{v};\\boldsymbol{\\gamma}) &= g_3^{-1}(\\boldsymbol{\\tilde{v}}^\\top \\boldsymbol{\\gamma}) = \\frac{e^{\\boldsymbol{\\tilde{v}}^\\top \\boldsymbol{\\gamma}}}{1 + e^{\\boldsymbol{\\tilde{v}}^\\top \\boldsymbol{\\gamma}}}, \\\\
            \\eta(\\boldsymbol{z};\\boldsymbol{\\lambda}) &= g_4^{-1}(\\boldsymbol{\\tilde{z}}^\\top \\boldsymbol{\\lambda}) = e^{\\boldsymbol{\\tilde{z}}^\\top \\boldsymbol{\\lambda}} + 1.
          \\end{align*}
          \\]
        ')
                                 )
                               )
                             ),
                             # --- End LaTeX block ---

                             hr(),

                             downloadButton("download_results", "Download RDS")
                           )
                           ,
                           tabPanel("Model Comparison",
                                    h3("Model Fit Comparison"),
                                    p("This table compares Poisson, Negative Binomial, CNB, ZIP, ZINB, and other models by AIC, BIC, and Log-Likelihood."),
                                    tableOutput("model_comparison_table"),
                                    br(),
                                    h4("Best Model"),
                                    verbatimTextOutput("best_model_text")
                           )
               )
             )
           )
  ),

  # Tab 1: New independent page (now first)
  tabPanel("Background Information",
           withMathJax(),
           fluidRow(
             column(10, offset = 1,
                    h2("Contaminated Negative Binomial Distribution"),

                    p("A random variable \\( Y \\) is said to follow a mean parameterized NB-D if its PMF is given as:"),
                    p("$$\\begin{align}
f(y; \\mu, \\alpha ) =
\\begin{cases}
\\frac{\\Gamma\\left(y + \\frac{1}{\\alpha}\\right)}{\\Gamma(y + 1) \\Gamma\\left(\\frac{1}{\\alpha}\\right)} \\left(\\frac{\\alpha \\mu}{1 + \\alpha \\mu}\\right)^y \\left(\\frac{1}{1 + \\alpha \\mu}\\right)^{\\frac{1}{\\alpha}},\\quad y \\in \\mathbb{N}_{0} \\\\
0,\\quad \\text{otherwise}
\\end{cases}
\\end{align}$$"),

                    p("where the expected value \\( E_{NB}(Y; \\mu) = \\mu > 0 \\) is the mean and \\( \\alpha > 0 \\) is the dispersion parameter. \\( \\Gamma(\\cdot) \\) denotes the gamma function. When a random variable \\( Y \\) follows the distribution specified above, it can be said that \\( Y \\sim NB(\\mu, \\alpha) \\)."),

                    p("When \\( \\alpha \\rightarrow 0^+ \\), the NB-D tends to the Poisson-D with mean \\( \\mu \\) and variance \\( \\text{Var}_{NB}(Y; \\mu, \\alpha) = \\mu + \\alpha\\mu^2 \\). It is worth mentioning that the NB-D and cNB-D are both mean \\( \\mu \\) parameterized. When stated in this manner, the expressions for the variance of both the NB-D and cNB-D indicate that the variation of the data is some linear function of the mean that is greater than or equal to the mean."),

                    p("The Contaminated NB-D (cNB-D) is proposed for handling NB overdispersion. A random variable \\( Y \\) is said to follow a mean parameterized cNB-D if its PMF is given as:"),
                    p("$$\\begin{align}
f(y; \\mu, \\alpha, \\delta, \\eta ) =
\\begin{cases}
(1 - \\delta) \\cdot \\frac{\\Gamma\\left(y + \\frac{1}{\\alpha}\\right)}{\\Gamma(y + 1) \\Gamma\\left(\\frac{1}{\\alpha}\\right)} \\left(\\frac{\\alpha \\mu}{1 + \\alpha \\mu}\\right)^y \\left(\\frac{1}{1 + \\alpha \\mu}\\right)^{\\frac{1}{\\alpha}} \\\\
+ \\delta \\cdot \\frac{\\Gamma\\left(y + \\frac{1}{\\eta\\alpha}\\right)}{\\Gamma(y + 1) \\Gamma\\left(\\frac{1}{\\eta\\alpha}\\right)} \\left(\\frac{\\eta\\alpha \\mu}{1 + \\eta\\alpha \\mu}\\right)^y \\left(\\frac{1}{1 + \\eta\\alpha \\mu}\\right)^{\\frac{1}{\\eta\\alpha}},\\quad y \\in \\mathbb{N}_{0} \\\\
0,\\quad \\text{otherwise}
\\end{cases}
\\end{align}$$"),

                    p("whereby \\( \\delta \\in(0, 1) \\) and \\( \\eta > 1 \\)."),

                    p("As a shorthand, the cNB-D could alternatively be denoted as:"),
                    p("$$f_{cNB}(y; \\mu, \\alpha, \\delta, \\eta) = (1 - \\delta) f_{NB}(y; \\mu, \\alpha) + \\delta f_{NB}(y; \\mu, \\eta \\alpha)$$"),

                    p("When a random variable \\( Y \\) follows the distribution specified above, it can be denoted as \\( Y \\sim cNB(\\mu, \\alpha, \\delta, \\eta) \\)."),

                    p("The parameters \\( \\delta \\) and \\( \\eta \\) have a practical interpretation, whereby \\( \\delta \\) represents the proportion of points from the contaminant distribution, whilst \\( \\eta \\) represents the extent of contamination of the data. The parameter \\( \\eta \\) can be viewed as an inflation parameter: since \\( \\eta > 1 \\), the further \\( \\eta \\) is from 1, the greater the variation of the contaminant distribution will be in comparison to the reference distribution."),

                    p("The expected value \\( E_{cNB}(Y; \\mu) = \\mu > 0 \\) is the mean, and the variance is given by:"),
                    p("$$\\text{Var}_{cNB}(Y; \\mu, \\alpha, \\delta, \\eta) = \\mu + [(1 - \\delta) + \\delta \\eta] \\alpha \\mu^2$$")
             )
           )
  ),

  tags$head(
    tags$style(HTML("
  /* Make DT table text white */
  table.dataTable tbody td {
   color: white !important;
  }
  table.dataTable thead th {
   color: white !important;
  }
  /* Optional: Make table background a bit darker for contrast */
  table.dataTable tbody {
   background-color: #2c3e50 !important;
  }
 "))
  )


)



server <- function(input, output, session) {

  # Add this inside your server function

  observeEvent(input$go_var, {
    updateNavbarPage(session, inputId = "navbar", selected = "cNB Interactive Visuals")
  })

  # observeEvent(input$go_sens, {
  #   updateNavbarPage(session, inputId = "navbar", selected = "Sensitivity Study")
  # })

  observeEvent(input$go_data, {
    updateNavbarPage(session, inputId = "navbar", selected = "Data Application")
  })

  observeEvent(input$go_bg, {
    updateNavbarPage(session, inputId = "navbar", selected = "Background Information")
  })



  # Dynamically render parameter controls based on plot type
  output$paramControls <- renderUI({
    switch(input$plotType,
           "mu_delta" = tagList(
             numericInput("eta", "η (eta):", value = 2, min = 0),
             sliderInput("alpha", "α (alpha):", min = 0.1, max = 5, value = 2, step = 0.1),
             sliderInput("delta1", "δ₁:", min = 0, max = 1, value = 0.04, step = 0.01),
             sliderInput("delta2", "δ₂:", min = 0, max = 1, value = 0.35, step = 0.01),
             sliderInput("delta3", "δ₃:", min = 0, max = 1, value = 0.55, step = 0.01)
           ),
           "mu_eta" = tagList(
             sliderInput("delta", "δ (delta):", min = 0, max = 1, value = 0.04, step = 0.01),
             sliderInput("alpha", "α (alpha):", min = 0.1, max = 5, value = 2, step = 0.1),
             numericInput("eta1", "η₁:", value = 4, min = 0),
             numericInput("eta2", "η₂:", value = 8, min = 0),
             numericInput("eta3", "η₃:", value = 16, min = 0)
           ),
           "alpha_delta" = tagList(
             numericInput("eta", "η (eta):", value = 2, min = 0),
             numericInput("mu", "μ (mu):", value = 2, min = 0),
             sliderInput("delta1", "δ₁:", min = 0, max = 1, value = 0.04, step = 0.01),
             sliderInput("delta2", "δ₂:", min = 0, max = 1, value = 0.35, step = 0.01),
             sliderInput("delta3", "δ₃:", min = 0, max = 1, value = 0.55, step = 0.01)
           ),
           "alpha_eta" = tagList(
             sliderInput("delta", "δ (delta):", min = 0, max = 1, value = 0.04, step = 0.01),
             numericInput("mu", "μ (mu):", value = 2, min = 0),
             numericInput("eta1", "η₁:", value = 4, min = 0),
             numericInput("eta2", "η₂:", value = 8, min = 0),
             numericInput("eta3", "η₃:", value = 16, min = 0)
           )
    )
  })

  # Reactive plot logic
  output$variancePlot <- renderPlot({
    try({
      plotType <- input$plotType

      if (plotType == "mu_delta") {
        mu <- seq(0, 5, by = 0.01)
        eta <- input$eta
        alpha <- input$alpha
        d1 <- input$delta1; d2 <- input$delta2; d3 <- input$delta3

        nb <- mu + alpha * mu^2
        v1 <- mu + ((1 - d1) + d1 * eta) * alpha * mu^2
        v2 <- mu + ((1 - d2) + d2 * eta) * alpha * mu^2
        v3 <- mu + ((1 - d3) + d3 * eta) * alpha * mu^2

        plot(mu, nb, type = "l", lwd = 2, xlab = expression(mu), ylab = "Variance",
             main = "Variance vs μ (varying δ)", col = "#34495E", ylim = range(c(nb, v1, v2, v3)))
        lines(mu, v1, lty = 2, col = "#E74C3C")
        lines(mu, v2, lty = 3, col = "#27AE60")
        lines(mu, v3, lty = 4, col = "#2980B9")
        legend("topleft", legend = c("NB", bquote(delta == .(d1)), bquote(delta == .(d2)), bquote(delta == .(d3))),
               lty = c(1, 2, 3, 4), col = c("#34495E", "#E74C3C", "#27AE60", "#2980B9"), bty = "n")

      } else if (plotType == "mu_eta") {
        mu <- seq(0, 5, by = 0.01)
        delta <- input$delta
        alpha <- input$alpha
        e1 <- input$eta1; e2 <- input$eta2; e3 <- input$eta3

        nb <- mu + alpha * mu^2
        v1 <- mu + ((1 - delta) + delta * e1) * alpha * mu^2
        v2 <- mu + ((1 - delta) + delta * e2) * alpha * mu^2
        v3 <- mu + ((1 - delta) + delta * e3) * alpha * mu^2

        plot(mu, nb, type = "l", lwd = 2, xlab = expression(mu), ylab = "Variance",
             main = "Variance vs μ (varying η)", col = "#34495E", ylim = range(c(nb, v1, v2, v3)))
        lines(mu, v1, lty = 2, col = "#E74C3C")
        lines(mu, v2, lty = 3, col = "#27AE60")
        lines(mu, v3, lty = 4, col = "#2980B9")
        legend("topleft", legend = c("NB", bquote(eta == .(e1)), bquote(eta == .(e2)), bquote(eta == .(e3))),
               lty = c(1, 2, 3, 4), col = c("#34495E", "#E74C3C", "#27AE60", "#2980B9"), bty = "n")

      } else if (plotType == "alpha_delta") {
        alpha <- seq(0, 5, by = 0.01)
        eta <- input$eta
        mu <- input$mu
        d1 <- input$delta1; d2 <- input$delta2; d3 <- input$delta3

        mu_vec <- rep(mu, length(alpha))
        nb <- mu_vec + alpha * mu_vec^2
        v1 <- mu_vec + ((1 - d1) + d1 * eta) * alpha * mu_vec^2
        v2 <- mu_vec + ((1 - d2) + d2 * eta) * alpha * mu_vec^2
        v3 <- mu_vec + ((1 - d3) + d3 * eta) * alpha * mu_vec^2

        plot(alpha, nb, type = "l", lwd = 2, xlab = expression(alpha), ylab = "Variance",
             main = "Variance vs α (varying δ)", col = "#34495E", ylim = range(c(nb, v1, v2, v3)))
        lines(alpha, v1, lty = 2, col = "#E74C3C")
        lines(alpha, v2, lty = 3, col = "#27AE60")
        lines(alpha, v3, lty = 4, col = "#2980B9")
        legend("topleft", legend = c("NB", bquote(delta == .(d1)), bquote(delta == .(d2)), bquote(delta == .(d3))),
               lty = c(1, 2, 3, 4), col = c("#34495E", "#E74C3C", "#27AE60", "#2980B9"), bty = "n")

      } else if (plotType == "alpha_eta") {
        alpha <- seq(0, 5, by = 0.01)
        delta <- input$delta
        mu <- input$mu
        e1 <- input$eta1; e2 <- input$eta2; e3 <- input$eta3

        mu_vec <- rep(mu, length(alpha))
        nb <- mu_vec + alpha * mu_vec^2
        v1 <- mu_vec + ((1 - delta) + delta * e1) * alpha * mu_vec^2
        v2 <- mu_vec + ((1 - delta) + delta * e2) * alpha * mu_vec^2
        v3 <- mu_vec + ((1 - delta) + delta * e3) * alpha * mu_vec^2

        plot(alpha, nb, type = "l", lwd = 2, xlab = expression(alpha), ylab = "Variance",
             main = "Variance vs α (varying η)", col = "#34495E", ylim = range(c(nb, v1, v2, v3)))
        lines(alpha, v1, lty = 2, col = "#E74C3C")
        lines(alpha, v2, lty = 3, col = "#27AE60")
        lines(alpha, v3, lty = 4, col = "#2980B9")
        legend("topleft", legend = c("NB", bquote(eta == .(e1)), bquote(eta == .(e2)), bquote(eta == .(e3))),
               lty = c(1, 2, 3, 4), col = c("#34495E", "#E74C3C", "#27AE60", "#2980B9"), bty = "n")
      }
    }, silent = TRUE) # suppress any errors without showing in UI
  })

  rv <- reactiveValues(data = NULL, name = NULL)

  observeEvent(input$load_data, {
    if (input$dataset_choice == "recreation") {
      data("RecreationDemand", package = "AER")
      rv$data <- RecreationDemand; rv$name <- "RecreationDemand"
      output$status <- renderText("Loaded RecreationDemand (AER)")
    } else if (input$dataset_choice == "nmes") {
      # data("NMES1988", package = "AER")
      # rv$data <- NMES1988; rv$name <- "NMES1988"
      data("NMES1988", package = "AER")

      options(contrasts = c("contr.treatment", "contr.poly"))

      NMES1988_mm <- model.matrix(
        ~ visits + nvisits + ovisits + novisits + emergency + hospital +
          health + chronic + adl + region + afam + gender + married +
          school + income + employed + insurance + medicaid -1,
        data = NMES1988
      )

      rv$data <- as.data.frame(NMES1988_mm)
      rv$name <- "NMES1988"

      output$status <- renderText("Loaded NMES1988 (AER)")
    } else if (input$dataset_choice == "upload" && !is.null(input$file_upload)) {
      rv$data <- tryCatch(read.csv(input$file_upload$datapath, stringsAsFactors = FALSE),
                          error = function(e) NULL)
      rv$name <- input$file_upload$name
      output$status <- renderText(paste("Loaded file:", rv$name))
    }
  })

  varnames <- reactive({ if (is.null(rv$data)) NULL else names(rv$data) })

  output$response_ui <- renderUI({
    if (is.null(varnames())) return(NULL)
    selectInput("response", "Response variable:", choices = varnames())
  })
  output$predictors_ui <- renderUI({
    if (is.null(varnames())) return(NULL)
    checkboxGroupInput("predictors", "Predictors (mean):", choices = varnames())
  })

  output$alpha_ui <- renderUI({
    if (is.null(varnames())) return(NULL)
    checkboxGroupInput("alpha_vars", "Alpha variables:", choices = varnames())
  })
  output$delta_ui <- renderUI({
    if (is.null(varnames())) return(NULL)
    checkboxGroupInput("delta_vars", "Delta variables:", choices = varnames())
  })
  output$eta_ui <- renderUI({
    if (is.null(varnames())) return(NULL)
    checkboxGroupInput("eta_vars", "Eta variables:", choices = varnames())
  })

  output$head_out <- renderPrint({ if (!is.null(rv$data)) head(rv$data) })
  output$summary_out <- renderPrint({ if (!is.null(rv$data)) summary(rv$data) })

  output$formula_text <- renderPrint({
    if (is.null(rv$data) || is.null(input$response)) return("No data")
    resp <- input$response
    preds <- input$predictors
    mean_form <- if (length(preds) == 0) as.formula(paste(resp, "~ 1"))
    else as.formula(paste(resp, "~", paste(preds, collapse = "+")))
    alpha_form <- if (length(input$alpha_vars) == 0) ~1 else as.formula(paste("~", paste(input$alpha_vars, collapse = "+")))
    delta_form <- if (length(input$delta_vars) == 0) ~1 else as.formula(paste("~", paste(input$delta_vars, collapse = "+")))
    eta_form  <- if (length(input$eta_vars) == 0) ~1 else as.formula(paste("~", paste(input$eta_vars, collapse = "+")))
    list(mean = mean_form, alpha = alpha_form, delta = delta_form, eta = eta_form)
  })


  model_fit <- eventReactive(input$run_model, {
    if (is.null(rv$data) || is.null(input$response)) return(NULL)

    withProgress(message = "Running model, please wait...", value = 0, {
      incProgress(0.2)
      resp <- input$response
      preds <- input$predictors
      mean_form <- if (length(preds) == 0) as.formula(paste(resp, "~ 1"))
      else as.formula(paste(resp, "~", paste(preds, collapse = "+")))
      alpha_form <- if (length(input$alpha_vars) == 0) ~1 else as.formula(paste("~", paste(input$alpha_vars, collapse = "+")))
      delta_form <- if (length(input$delta_vars) == 0) ~1 else as.formula(paste("~", paste(input$delta_vars, collapse = "+")))
      eta_form  <- if (length(input$eta_vars) == 0) ~1 else as.formula(paste("~", paste(input$eta_vars, collapse = "+")))

      incProgress(0.5)
      fit <- tryCatch(
        suppressWarnings(
          ml.ex.cnb(formula = mean_form,
                    alpha.formula = alpha_form,
                    delta.formula = delta_form,
                    eta.formula  = eta_form,
                    data = rv$data)
        ),
        error = function(e) NULL
      )
      incProgress(1)
      fit
    })
  })

  observeEvent(model_fit(), {
    if (!is.null(model_fit())) {
      updateTabsetPanel(session, "main_tabs", selected = "Results")
    }
  })

  output$results_table_ui <- renderUI({
    fit <- model_fit()
    if (is.null(fit)) return("No results — check data/model")
    tableOutput("results_table")
  })

  output$results_table <- renderTable({
    fit <- model_fit()
    if (is.null(fit)) return(NULL)
    signif(fit$results, 5)
  }, rownames = TRUE)

  output$info_table <- renderTable({
    fit <- model_fit()
    if (is.null(fit)) return(NULL)
    data.frame(
      Statistic = c("LogLik", "AIC", "BIC", "HQIC"),
      Value = c(fit$loglike, fit$AIC, fit$BIC, fit$HQIC)
    )
  })

  output$optim_info <- renderPrint({
    fit <- model_fit()
    if (is.null(fit)) { cat("No optimization results.") } else { str(fit$optim) }
  })



  output$download_results <- downloadHandler(
    filename = function() paste0("ml_ex_cnb_results_", Sys.Date(), ".rds"),
    content = function(file) saveRDS(model_fit(), file)
  )


  model_comparison <- eventReactive(input$run_model, {
    if (is.null(rv$data) || is.null(input$response)) return(NULL)

    resp <- input$response
    preds <- input$predictors
    mean_form <- if (length(preds) == 0) as.formula(paste(resp, "~ 1"))
    else as.formula(paste(resp, "~", paste(preds, collapse = "+")))

    data <- rv$data

    results_list <- list()

    # Poisson
    poi_fit <- tryCatch(glm(mean_form, data = data, family = poisson), error = function(e) NULL)
    if (!is.null(poi_fit)) {
      results_list$Poisson <- c(
        LogLik = as.numeric(logLik(poi_fit)),
        AIC = AIC(poi_fit),
        BIC = BIC(poi_fit)
      )
    }

    # Negative Binomial (NB2)
    nb2_fit <- tryCatch(MASS::glm.nb(mean_form, data = data), error = function(e) NULL)
    if (!is.null(nb2_fit)) {
      results_list$NegBin2 <- c(
        LogLik = as.numeric(logLik(nb2_fit)),
        AIC = AIC(nb2_fit),
        BIC = BIC(nb2_fit)
      )
    }

    # Zero-inflated Poisson (ZIP)
    zip_fit <- tryCatch(
      gamlss::gamlss(mean_form, family = gamlss.dist::ZIP2, data = data,
                     control = gamlss.control(n.cyc = 200)),
      error = function(e) NULL
    )
    if (!is.null(zip_fit)) {
      results_list$ZIP <- c(
        LogLik = as.numeric(logLik(zip_fit)),
        AIC = AIC(zip_fit),
        BIC = BIC(zip_fit)
      )
    }

    # Zero-inflated NegBin (ZINB)
    zinb_fit <- tryCatch(
      gamlss::gamlss(mean_form, family = gamlss.dist::ZINBI, data = data,
                     control = gamlss.control(n.cyc = 200)),
      error = function(e) NULL
    )
    if (!is.null(zinb_fit)) {
      results_list$ZINB <- c(
        LogLik = as.numeric(logLik(zinb_fit)),
        AIC = AIC(zinb_fit),
        BIC = BIC(zinb_fit)
      )
    }

    # CNB (using your existing model_fit)
    cnb_fit <- model_fit()
    if (!is.null(cnb_fit)) {
      results_list$CNB <- c(
        LogLik = cnb_fit$loglike,
        AIC = cnb_fit$AIC,
        BIC = cnb_fit$BIC
      )
    }

    # Convert to data frame
    comp_df <- do.call(rbind, results_list)
    comp_df <- as.data.frame(comp_df)
    comp_df$Model <- rownames(comp_df)
    rownames(comp_df) <- NULL

    comp_df
  })

  output$model_comparison_table <- renderTable({
    df <- model_comparison()
    if (is.null(df)) return(NULL)
    df <- df[order(df$AIC), ] # sort by AIC ascending
    df
  }, digits = 3)

  output$best_model_text <- renderPrint({
    df <- model_comparison()
    if (is.null(df)) {
      cat("No models fitted yet.")
      return()
    }

    best_aic <- df$Model[which.min(df$AIC)]
    best_bic <- df$Model[which.min(df$BIC)]
    best_ll <- df$Model[which.max(df$LogLik)]

    cat("Best model by AIC:", best_aic, "\n")
    cat("Best model by BIC:", best_bic, "\n")
    cat("Best model by Log-Likelihood:", best_ll, "\n")
  })

  observeEvent(model_comparison(), {
    updateTabsetPanel(session, "main_tabs", selected = "Results")
  })

  output$animated_plot <- renderImage({
    req(input$show_animation)

    outfile <- tempfile(fileext = ".gif")

    x_vals <- 0:20
    mu <- 10
    eta_seq  <- seq(1.2, 3.0, length.out = 20)
    alpha_seq <- seq(0.2, 0.9, length.out = 20)
    delta_seq <- seq(0.05, 0.4, length.out = 20)


    data_eta <- data.table::rbindlist(lapply(seq_along(eta_seq), function(i) {
      data.frame(
        x = x_vals,
        prob = suppressWarnings(pmax(0, dcnbinom(x_vals, mu = mu, alpha = 0.6, delta = 0.2, eta = eta_seq[i]))),
        frame = i,
        phase = "Varying η"
      )
    }))

    data_alpha <- data.table::rbindlist(lapply(seq_along(alpha_seq), function(i) {
      data.frame(
        x = x_vals,
        prob = suppressWarnings(pmax(0, dcnbinom(x_vals, mu = mu, alpha = alpha_seq[i], delta = 0.2, eta = 2.0))),
        frame = i + length(eta_seq),
        phase = "Varying α"
      )
    }))

    data_delta <- data.table::rbindlist(lapply(seq_along(delta_seq), function(i) {
      data.frame(
        x = x_vals,
        prob = suppressWarnings(pmax(0, dcnbinom(x_vals, mu = mu, alpha = 0.6, delta = delta_seq[i], eta = 2.0))),
        frame = i + length(eta_seq) + length(alpha_seq),
        phase = "Varying δ"
      )
    }))

    plot_data <- dplyr::bind_rows(data_eta, data_alpha, data_delta)
    plot_data$phase <- factor(plot_data$phase, levels = c("Varying η", "Varying α", "Varying δ"))

    p <- ggplot(plot_data, aes(x = x, y = prob, color = phase)) +
      geom_line(size = 1.2, show.legend = TRUE) +
      scale_color_manual(
        values = c("Varying η" = "#D55E00",
                   "Varying α" = "#0072B2",
                   "Varying δ" = "#009E73")
      ) +
      labs(
        title = "Contaminated Negative Binomial Distribution",
        subtitle = "Each parameter varies in sequence",
        x = "x (Count)", y = "Probability", color = "Phase"
      ) +
      theme_light(base_size = 14) +
      transition_states(frame, transition_length = 1, state_length = 0) +
      ease_aes("linear") + enter_fade() + exit_fade()

    withProgress(message = "Generating animation...", {
      anim_save(outfile, animate(p, fps = 4, width = 600, height = 400, renderer = gifski_renderer()))
    })

    list(src = outfile, contentType = "image/gif", width = 800, height = 500)
  }, deleteFile = TRUE)



}

shinyApp(ui, server)
