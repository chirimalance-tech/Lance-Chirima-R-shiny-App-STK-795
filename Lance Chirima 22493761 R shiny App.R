
required_packages <- c(
  "shiny", "shinythemes", "DT", "MASS", "dplyr", 
  "ggplot2", "tidyr", "tibble", "nlme", "AER", 
  "gganimate", "gifski", "COUNT", "gamlss",
  "gamlss.dist", "VGAM", "DDPM", "COMPoissonReg",
  "png", "magick"
)

# Function to install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Special handling for ggplot2 version compatibility
if (!requireNamespace("ggplot2", quietly = TRUE) ||
    packageVersion("ggplot2") < "3.5.2") {
  message("Installing ggplot2 ≥ 3.5.2 for gganimate compatibility...")
  install.packages("ggplot2", dependencies = TRUE)
} else {
  message("ggplot2 version is sufficient. Skipping installation.")
}

# Install remaining packages (excluding ggplot2 which is already handled)
for (pkg in setdiff(required_packages, "ggplot2")) {
  install_if_missing(pkg)
}

# Load all required packages
lapply(required_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
})

# Verify all packages loaded successfully
for (pkg in required_packages) {
  if (!suppressWarnings(require(pkg, character.only = TRUE))) {
    stop(paste("Package", pkg, "could not be loaded. Please install it."))
  }
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

rcnbinom_ex_regression <- function(n = 50000, off = 0,
                                   theta = c((1.8 - 1), (1.9 - 1)),
                                   beta = c(2, 0.75, -1.5),
                                   gamma = c(0.05 / (1 - 0.05), 0.04 / (1 - 0.04)),
                                   lambda = c((1.2 - 1), (1.3 - 1)),
                                   contam_perc = 0) {
  dX <- length(beta) - 1 # number of predictors
  dU <- length(theta) - 1 
  dV <- length(gamma) - 1 
  dZ <- length(lambda) - 1 
  
  # X, U, V, Z matrices
  X <- matrix(runif(n * dX, -2, 2), nrow = n, ncol = dX)
  Xstar <- cbind(1, X) 
  
  U <- matrix(runif(n * dU, -2, 2), nrow = n, ncol = dU)
  Ustar <- cbind(1, U) 
  
  V <- matrix(runif(n * dV, -2, 2), nrow = n, ncol = dV)
  Vstar <- cbind(1, V) 
  
  Z <- matrix(runif(n * dZ, -2, 2), nrow = n, ncol = dZ)
  Zstar <- cbind(1, Z) 
  
  # Linear predictors
  lin.pred.mu  <- Xstar %*% beta
  cond.mu    <- exp(lin.pred.mu + off) + 1
  lin.pred.alpha <- Ustar %*% theta
  cond.alpha   <- exp(lin.pred.alpha)
  lin.pred.delta <- Vstar %*% gamma
  cond.delta   <- exp(lin.pred.delta) / (1 + exp(lin.pred.delta))
  lin.pred.eta  <- Zstar %*% lambda
  cond.eta    <- exp(lin.pred.eta)
  
  # Bernoulli mixture
  bernoulli <- rbinom(n, 1, prob = cond.delta)
  
  # Generate counts
  mu_good  <- cond.mu[bernoulli == 1]
  theta_good <- cond.alpha[bernoulli == 1]
  cnby_good <- MASS::rnegbin(length(mu_good), mu = mu_good, theta = theta_good)
  
  if (sum(bernoulli == 0) > 0) {
    mu_bad  <- cond.mu[bernoulli == 0]
    theta_bad <- cond.alpha[bernoulli == 0] * cond.eta[bernoulli == 0]
    cnby_bad <- MASS::rnegbin(length(mu_bad), mu = mu_bad, theta = theta_bad)
    
    cnby <- numeric(n)
    cnby[bernoulli == 1] <- cnby_good
    cnby[bernoulli == 0] <- cnby_bad
  } else {
    cnby <- cnby_good
  }
  
  # Contaminate a percentage of the data if contam_perc > 0
  if (contam_perc > 0) {
    n_replace <- ceiling(n * contam_perc / 100)
    contam_idx <- sample(seq_len(n), n_replace)
    cnby[contam_idx] <- sample(cnby, n_replace, replace = TRUE) + rpois(n_replace, lambda = 20)
  }
  
  out <- data.frame(cbind(cnby, Xstar, Ustar, Vstar, Zstar))
  colnames(out) <- c("cnby",
                     paste0("x", 0:dX),
                     paste0("u", 0:dU),
                     paste0("v", 0:dV),
                     paste0("z", 0:dZ))
  return(out)
}


ml.ex.cnb <- function(formula, alpha.formula=~1, delta.formula=~1, eta.formula=~1, data, start = NULL, method = "Nelder-Mead", reltol=1e-10, maxit=1000, hessian=T) {
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
    mu.hat <- exp(xb.hat) + 1
    
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
    nb <- MASS::glm.nb(formula = formula, data = na.omit(data), control = glm.control(maxit = 1000))
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
  mu.hat <- exp(xb.hat) + 1
  
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
  HIC <- -2*lc +nrow(results)*2*log(log(nrow(data)))
  
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
    HIC = HIC
  ))
}


# Dataset creation
make_dataset <- function(n, beta, p_percent = 0) {
  # Generate dataset with given beta
  reg <- rcnbinom_ex_regression(n = n, beta = beta)
  df <- data.frame(cnby = reg$cnby, x0 = reg$x0, x1 = reg$x1, x2 = reg$x2)
  
  # Apply contamination
  if (p_percent > 0) {
    p <- max(1, floor(n * p_percent / 100))
    minv <- min(df$cnby, na.rm = TRUE)
    maxv <- max(df$cnby, na.rm = TRUE)
    df$cnby[sample(1:n, p)] <- sample(minv:maxv, size = p, replace = TRUE)
  }
  df
}


# Model fitting wrappers
fit_glm_nb <- function(df) {
  tryCatch(
    suppressWarnings(MASS::glm.nb(cnby ~ x1 + x2, data = df, control = glm.control(maxit = 10000))),
    error = function(e) NULL
  )
}

fit_cnb_extended <- function(df) {
  res <- tryCatch({
    suppressWarnings(suppressMessages(
      ml.ex.cnb(formula = cnby ~ x1 + x2, data = df, method = "Nelder-Mead")
    ))
  }, error = function(e) {
    warning("cNB fit failed: ", e$message)
    return(NULL)
  })
  
  # Check for non-finite estimates
  if (!is.null(res) && any(!is.finite(res$beta[1:3]))) return(NULL)
  
  res
}

# Result extractors
extract_glm_nb_results <- function(fit) {
  if (is.null(fit)) return(tibble(term = c("beta0","beta1","beta2"), estimate = NA_real_, std_error = NA_real_))
  s <- summary(fit)$coefficients
  tibble(term = c("beta0","beta1","beta2"),
         estimate = as.numeric(s[c("(Intercept)", "x1", "x2"), "Estimate"]),
         std_error = as.numeric(s[c("(Intercept)", "x1", "x2"), "Std. Error"]))
}

extract_cnb_results <- function(res) {
  if (is.null(res) || !all(c("beta","results") %in% names(res)) ||
      any(!is.finite(res$beta[1:3])) ||
      any(is.na(res$results$Std.Error[1:3]))) {
    return(tibble(term = c("beta0","beta1","beta2"),
                  estimate = NA_real_, std_error = NA_real_))
  }
  
  tibble(term = c("beta0","beta1","beta2"),
         estimate = as.numeric(res$beta[1:3]),
         std_error = as.numeric(res$results$Std.Error[1:3]))
}


# Summarize bias and MSE
process_results_table <- function(results_tbl, true_beta) {
  truth_df <- tibble(term = names(true_beta), truth = as.numeric(true_beta))
  
  results_tbl %>%
    left_join(truth_df, by = "term") %>%
    group_by(term) %>%
    summarise(
      bias  = mean(estimate - truth, na.rm = TRUE),
      mse  = mean((estimate - truth)^2, na.rm = TRUE),
      n_used = sum(!is.na(estimate)),
      .groups = "drop"
    )
}




# MathJax-friendly labels
term_expr <- c(beta0 = expression(beta[0]), beta1 = expression(beta[1]), beta2 = expression(beta[2]))



ui <- navbarPage(
  id = "navbar",
  title = "Contaminated Negative Binomial Regression for Modeling Count Data",
  theme = shinytheme("superhero"),
  
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
          actionButton("go_sens", "Sensitivity Study", class = "home-btn"),
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
  
  # Tab 3: Sensitivity Analysis Simulator (App 2)
  tabPanel("Sensitivity Study",
           sidebarLayout(
             sidebarPanel(
               numericInput("n_datasets", "Number of datasets to generate:", value = 20, min = 10, max = 3000, step = 10),
               numericInput("n_obs", "Observations per dataset:", value = 50, min = 50, max = 10000, step = 50),
               numericInput("seed", "Random seed:", value = 345),
               
               h4("β Parameters"),
               numericInput("beta0", "β₀:", value = 2),
               numericInput("beta1", "β₁:", value = 0.75),
               numericInput("beta2", "β₂:", value = -1.5),
               
               sliderInput("contam_perc", "Contamination (%)", min = 0, max = 100, value = 5, step = 0.5),
               
               actionButton("run_sim", "Run simulation", class = "btn-primary"),
               br(), br(),
               downloadButton("download_csv", "Download combined results (CSV)")
             )
             ,
             mainPanel(
               tabsetPanel(
                 tabPanel("Summary",
                          withMathJax(),
                          br(),
                          hr(),
                          p("Simulation results (Bias & MSE) for NB vs cNB under user-specified contamination."),
                          DTOutput("summary_table"),
                          br(),
                          fluidRow(
                            column(6, plotOutput("plot_bias_nb", height = "300px")),
                            column(6, plotOutput("plot_bias_cnb", height = "300px"))
                          ),
                          br(),
                          fluidRow(
                            column(12, plotOutput("plot_bias_compare", height = "320px"))
                          ),
                          br(),
                          verbatimTextOutput("info_text")
                 ),
                 tabPanel("Raw per-dataset (sample)",
                          h4("NB"), DTOutput("raw_nb"),
                          h4("cNB"), DTOutput("raw_cnb")
                          
                 )
               )
             )
           )
  )
  ,
  
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
                    
                    p("A random variable \\( Y \\) is said to follow a mean parameterized NB-D if its PFM is given as:"),
                    p("$$\\begin{align}
f(y; \\mu, \\alpha ) = 
\\begin{cases}
\\frac{\\Gamma\\left(y + \\frac{1}{\\alpha}\\right)}{\\Gamma(y + 1) \\Gamma\\left(\\frac{1}{\\alpha}\\right)} \\left(\\frac{\\alpha \\mu}{1 + \\alpha \\mu}\\right)^y \\left(\\frac{1}{1 + \\alpha \\mu}\\right)^{\\frac{1}{\\alpha}},\\quad y \\in \\mathbb{N}_{0} \\\\ 
0,\\quad \\text{otherwise}
\\end{cases}
\\end{align}$$"),
                    
                    p("where the expected value \\( E_{NB}(Y; \\mu) = \\mu > 0 \\) is the mean and \\( \\alpha > 0 \\) is the dispersion parameter. \\( \\Gamma(\\cdot) \\) denotes the gamma function. When a random variable \\( Y \\) follows the distribution specified above, it can be said that \\( Y \\sim NB(\\mu, \\alpha) \\)."),
                    
                    p("When \\( \\alpha \\rightarrow 0^+ \\), the NB-D tends to the Poisson-D with mean \\( \\mu \\) and variance \\( \\text{Var}_{NB}(Y; \\mu, \\alpha) = \\mu + \\alpha\\mu^2 \\). It is worth mentioning that the NB-D and cNB-D are both mean \\( \\mu \\) parameterized. When stated in this manner, the expressions for the variance of both the NB-D and cNB-D indicate that the variation of the data is some linear function of the mean that is greater than or equal to the mean."),
                    
                    p("The Contaminated NB-D (cNB-D) is proposed for handling NB overdispersion. A random variable \\( Y \\) is said to follow a mean parameterized cNB-D if its PFM is given as:"),
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
  
  observeEvent(input$go_sens, {
    updateNavbarPage(session, inputId = "navbar", selected = "Sensitivity Study")
  })
  
  observeEvent(input$go_data, {
    updateNavbarPage(session, inputId = "navbar", selected = "Data Application")
  })
  
  observeEvent(input$go_bg, {
    updateNavbarPage(session, inputId = "navbar", selected = "Background Information")
  })
  
  # ---- SERVER: Sensitivity Study ----
  sim_results <- eventReactive(input$run_sim, {
    showModal(modalDialog(
      title = NULL,
      div(style = "text-align:center;", h3("Please wait — running simulations")),
      footer = NULL,
      easyClose = FALSE
    ))
    on.exit(removeModal())
    
    set.seed(input$seed)
    n_datasets <- input$n_datasets
    n_obs   <- input$n_obs
    contam_perc <- input$contam_perc / 100 # convert to proportion
    true_beta <- c(beta0 = input$beta0, beta1 = input$beta1, beta2 = input$beta2)
    
    
    nb_list <- vector("list", n_datasets)
    cnb_list <- vector("list", n_datasets)
    
    withProgress(message = paste0("Running simulations at ", input$contam_perc, "% contamination"), value = 0, {
      for (i in seq_len(n_datasets)) {
        incProgress(1/n_datasets, detail = paste("Dataset", i, "of", n_datasets))
        
        # Generate dataset using your modified function
        d <- make_dataset(n = n_obs, beta = true_beta, p_percent = input$contam_perc)
        
        # Fit models
        nb_fit <- fit_glm_nb(d)
        cnb_fit <- fit_cnb_extended(d)
        
        nb_list[[i]] <- extract_glm_nb_results(nb_fit) %>% mutate(dataset = i)
        cnb_list[[i]] <- extract_cnb_results(cnb_fit)  %>% mutate(dataset = i)
      }
    })
    
    nb_df <- bind_rows(nb_list)
    cnb_df <- bind_rows(cnb_list)
    
    # Summarize results
    true_beta <- c(beta0 = input$beta0,
                   beta1 = input$beta1,
                   beta2 = input$beta2)
    
    nb_summary <- process_results_table(nb_df, true_beta = true_beta) %>%
      rename(Bias_NB = bias, MSE_NB = mse, Nused_NB = n_used)
    
    cnb_summary <- process_results_table(cnb_df, true_beta = true_beta) %>%
      rename(Bias_cNB = bias, MSE_cNB = mse, Nused_cNB = n_used)
    
    
    combined <- nb_summary %>%
      left_join(cnb_summary, by = "term") %>%
      mutate(term_label = case_when(
        term == "b0" ~ "\\(\\beta_0\\)",
        term == "b1" ~ "\\(\\beta_1\\)",
        term == "b2" ~ "\\(\\beta_2\\)",
        TRUE ~ term
      ))
    
    list(combined = combined, nb_df = nb_df, cnb_df = cnb_df)
  })
  
  # ---- SUMMARY TABLE ----
  output$summary_table <- renderDT({
    req(sim_results())
    dat <- sim_results()$combined
    
    # Round numeric columns
    numcols <- sapply(dat, is.numeric)
    dat[numcols] <- lapply(dat[numcols], function(x) round(x, 4))
    
    # Select only relevant columns
    dat_display <- dat %>%
      dplyr::select(term_label, Bias_NB, MSE_NB, Bias_cNB, MSE_cNB)
    
    colnames(dat_display) <- c("Parameter", "Bias NB", "MSE NB", "Bias cNB", "MSE cNB")
    
    datatable(dat_display,
              escape = FALSE,
              options = list(dom='t', columnDefs = list(list(className='dt-center', targets=1:4))),
              class = 'cell-border stripe')
  })
  
  # ---- RAW TABLES ----
  output$raw_nb <- renderDT({ req(sim_results()); sim_results()$nb_df %>% filter(dataset <= 10) %>% datatable() })
  output$raw_cnb <- renderDT({ req(sim_results()); sim_results()$cnb_df %>% filter(dataset <= 10) %>% datatable() })
  
  # ---- INFO TEXT ----
  output$info_text <- renderText({
    df <- sim_results()
    failed_cnb <- df$cnb_df %>%
      group_by(dataset) %>%
      summarise(all_na = all(is.na(estimate)), .groups="drop") %>%
      summarise(fails = sum(all_na)) %>%
      pull(fails)
    paste0("cNB fits failed for ", input$contam_perc, "% contamination: ", failed_cnb)
  })
  
  # ---- DOWNLOAD CSV ----
  output$download_csv <- downloadHandler(
    filename = function() paste0("sensitivity_simulation_results_", Sys.Date(), ".csv"),
    content = function(file) {
      req(sim_results())
      write.csv(sim_results()$combined, file, row.names = FALSE)
    }
  )
  
  # ---- PLOTS ----
  
  # Bias — NB
  output$plot_bias_nb <- renderPlot({
    req(sim_results())
    df <- sim_results()$combined
    
    ggplot(df, aes(x=factor(term, levels=c("b0", "b1", "b2")), y=Bias_NB)) +
      geom_boxplot(fill="#34495E", alpha=0.7, outlier.color="grey40") +
      scale_x_discrete(labels=c(expression(beta[0]), expression(beta[1]), expression(beta[2]))) +
      labs(
        title=paste0("Bias Distribution — NB (", input$contam_perc, "% contamination)"),
        x=NULL, y="Bias"
      ) +
      theme_minimal(base_size=14)
  })
  
  # Bias — cNB
  output$plot_bias_cnb <- renderPlot({
    req(sim_results())
    df <- sim_results()$combined
    
    ggplot(df, aes(x=factor(term, levels=c("b0", "b1", "b2")), y=Bias_cNB)) +
      geom_boxplot(fill="#E74C3C", alpha=0.7, outlier.color="grey40") +
      scale_x_discrete(labels=c(expression(beta[0]), expression(beta[1]), expression(beta[2]))) +
      labs(
        title=paste0("Bias Distribution — cNB (", input$contam_perc, "% contamination)"),
        x=NULL, y="Bias"
      ) +
      theme_minimal(base_size=14)
  })
  
  # Bias comparison — NB vs cNB
  output$plot_bias_compare <- renderPlot({
    req(sim_results())
    df <- sim_results()$combined %>%
      pivot_longer(cols=c(Bias_NB, Bias_cNB), names_to="Model", values_to="Bias")
    
    ggplot(df, aes(
      x=factor(term, levels=c("b0","b1","b2")),
      y=Bias,
      fill=Model
    )) +
      geom_boxplot(alpha=0.7, outlier.color="grey40") +
      scale_x_discrete(labels=c(expression(beta[0]), expression(beta[1]), expression(beta[2]))) +
      scale_fill_manual(
        values=c("Bias_NB"="#34495E", "Bias_cNB"="#E74C3C"),
        labels=c("NB","cNB")
      ) +
      labs(
        title=paste0("Bias Comparison (", input$contam_perc, "% contamination)"),
        x=NULL, y="Bias", fill="Model"
      ) +
      theme_minimal(base_size=14)
  })
  
  # ---- MSE Comparison ----
  output$plot_mse_compare <- renderPlot({
    req(sim_results())
    df <- sim_results()$combined %>%
      pivot_longer(cols=c(MSE_NB, MSE_cNB), names_to="Model", values_to="MSE")
    
    ggplot(df, aes(
      x=factor(term, levels=c("b0","b1","b2")),
      y=MSE,
      fill=Model
    )) +
      geom_boxplot(alpha=0.7, outlier.color="grey40") +
      scale_x_discrete(labels=c(expression(beta[0]), expression(beta[1]), expression(beta[2]))) +
      scale_fill_manual(
        values=c("MSE_NB"="#34495E", "MSE_cNB"="#E74C3C"),
        labels=c("NB","cNB")
      ) +
      labs(
        title=paste0("MSE Comparison (", input$contam_perc, "% contamination)"),
        x=NULL, y="Mean Squared Error", fill="Model"
      ) +
      theme_minimal(base_size=14)
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
      data("NMES1988", package = "AER")
      rv$data <- NMES1988; rv$name <- "NMES1988"
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
      Statistic = c("LogLik", "AIC", "BIC", "HIC"),
      Value = c(fit$loglike, fit$AIC, fit$BIC, fit$HIC)
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
    eta_seq  <- seq(1.2, 3.0, length.out = 60)
    alpha_seq <- seq(0.2, 0.9, length.out = 60)
    delta_seq <- seq(0.05, 0.4, length.out = 60)
    
    data_eta <- lapply(seq_along(eta_seq), function(i) {
      data.frame(
        x = x_vals,
        prob = suppressWarnings(pmax(0, dcnbinom(x_vals, mu = mu, alpha = 0.6, delta = 0.2, eta = eta_seq[i]))),
        frame = i,
        phase = "Varying η"
      )
    }) |> dplyr::bind_rows()
    
    data_alpha <- lapply(seq_along(alpha_seq), function(i) {
      data.frame(
        x = x_vals,
        prob = suppressWarnings(pmax(0, dcnbinom(x_vals, mu = mu, alpha = alpha_seq[i], delta = 0.2, eta = 2.0))),
        frame = i + length(eta_seq),
        phase = "Varying α"
      )
    }) |> dplyr::bind_rows()
    
    data_delta <- lapply(seq_along(delta_seq), function(i) {
      data.frame(
        x = x_vals,
        prob = suppressWarnings(pmax(0, dcnbinom(x_vals, mu = mu, alpha = 0.6, delta = delta_seq[i], eta = 2.0))),
        frame = i + length(eta_seq) + length(alpha_seq),
        phase = "Varying δ"
      )
    }) |> dplyr::bind_rows()
    
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
      theme_minimal(base_size = 16) +
      transition_states(frame, transition_length = 2, state_length = 1) +
      ease_aes("linear") + enter_fade() + exit_fade()
    
    withProgress(message = "Generating animation...", {
      anim_save(outfile, animate(p, fps = 2, width = 800, height = 500, renderer = magick_renderer()))
    })
    
    list(src = outfile, contentType = "image/gif", width = 800, height = 500)
  }, deleteFile = TRUE)
  
  
  
}

shinyApp(ui, server)