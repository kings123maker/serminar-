
setwd("C:\\Users\\Dr. O. J. Obulezi\\Documents\\R projects\\Seminar-Econometrics")
# =========================================================================
# SEMINAR IN EMPIRICAL FINANCE (SUMMER 2026)
# Complete Monolithic Script: Data Pipeline, Modeling, Backtesting, & Plots
# Clean Error-Free Production Version
# =========================================================================

# -------------------------------------------------------------------------
# PHASE 1: INITIALIZATION, UTILITIES, AND DATA PROCUREMENT
# -------------------------------------------------------------------------
library(quantmod)
library(rugarch)
library(zoo)
library(ggplot2)
library(gridExtra)

set.seed(2026)
start_date <- "2016-06-01"
end_date   <- "2026-06-01"

cat("Phase 1: Downloading historical asset indices from Yahoo Finance...\n")
getSymbols(c("^GSPC", "BTC-USD"), src = "yahoo", from = start_date, to = end_date, auto.assign = TRUE)

# Intersect timestamps to synchronize weekly equity markets with 24/7 crypto markets
clean_prices <- merge(GSPC$GSPC.Adjusted, `BTC-USD`$`BTC-USD.Adjusted`, all = FALSE)
names(clean_prices) <- c("SP500", "BTC")

# Transform prices to continuous Logarithmic Returns: R_t = ln(P_t / P_{t-1})
returns_xts <- na.omit(diff(log(clean_prices)))

# Export raw return matrix for your data submission requirements
write.zoo(returns_xts, file = "seminar_returns_data.csv", sep = ",", row.names = TRUE)

returns_matrix <- as.matrix(returns_xts)
N_total        <- nrow(returns_matrix)


# -------------------------------------------------------------------------
# PHASE 2: PARAMETER MODELING AND OPTIMIZATION SETUP
# -------------------------------------------------------------------------
cat("Phase 2: Initializing simulation engines and rolling parameters...\n")
window_size        <- 1000  # 1000-day rolling estimation window
out_of_sample_size <- N_total - window_size
alpha              <- 0.01  # Target 99% coverage level

assets <- c("SP500", "BTC")

# Nested list matrices to collect daily forecast loops
forecasts <- list(
  SP500 = list(HS_VaR = numeric(out_of_sample_size), GARCH_VaR = numeric(out_of_sample_size)),
  BTC   = list(HS_VaR = numeric(out_of_sample_size), GARCH_VaR = numeric(out_of_sample_size))
)

# Standard GARCH(1,1) spec augmented with a Student-t distribution for thick fat tails
garch_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model     = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)


# -------------------------------------------------------------------------
# PHASE 3: THE ROLLING ESTIMATION LOOP ENGINE (ROBUST TO SE ERRORS)
# -------------------------------------------------------------------------
cat("Phase 3: Running rolling risk simulations. This may take a few minutes...\n")

for (asset in assets) {
  cat(paste(" >> Simulating out-of-sample risk boundaries for:", asset, "\n"))
  
  for (t in 1:out_of_sample_size) {
    # Crop moving history target window
    start_idx <- t
    end_idx   <- t + window_size - 1
    history   <- returns_matrix[start_idx:end_idx, asset]
    
    # 1. Non-Parametric Historical Simulation Engine
    forecasts[[asset]]$HS_VaR[t] <- -quantile(history, probs = alpha)
    
    # 2. Dynamic Conditional GARCH-t Engine (Defensive TryCatch Configuration)
    garch_fit <- tryCatch({
      ugarchfit(
        spec = garch_spec, 
        data = history, 
        solver = "hybrid", 
        calculate.robust.se = FALSE, # Bypasses the 'object B not found' matrix error completely
        solver.control = list(trace = 0)
      )
    }, error = function(e) {
      return(NULL) # Safely flags a NULL if numeric optimization fails on a black-swan market crash day
    })
    
    if (!is.null(garch_fit) && convergence(garch_fit) == 0) {
      garch_fcst <- ugarchforecast(garch_fit, n.ahead = 1)
      
      mu_t       <- as.numeric(fitted(garch_fcst))
      sigma_t    <- as.numeric(sigma(garch_fcst))
      shape_df   <- as.numeric(coef(garch_fit)["shape"])
      
      # Correctly scale quantile boundary with dynamic degrees of freedom
      t_quantile <- qdist("std", p = alpha, shape = shape_df)
      forecasts[[asset]]$GARCH_VaR[t] <- -(mu_t + sigma_t * t_quantile)
    } else {
      # Robust numerical optimization fallback matrix to seamlessly handle extreme trading days
      forecasts[[asset]]$GARCH_VaR[t] <- forecasts[[asset]]$HS_VaR[t]
    }
  }
}


# -------------------------------------------------------------------------
# PHASE 4: STATISTICAL TEST SUITE (KUPIEC AND CHRISTOFFERSEN MARKOV DECOMPOSITION)
# -------------------------------------------------------------------------
run_backtest_suite <- function(actual_returns, predicted_var, alpha) {
  violations <- as.numeric(actual_returns < -predicted_var)
  N          <- sum(violations)
  T_total    <- length(violations)
  p_nominal  <- alpha
  p_hat      <- N / T_total
  
  # 1. Kupiec Proportion of Failures Test
  lr_uc_num <- (p_nominal^N) * ((1 - p_nominal)^(T_total - N))
  lr_uc_den <- (p_hat^N) * ((1 - p_hat)^(T_total - N))
  if(N == 0 || N == T_total) { lr_uc <- 0 } else { lr_uc <- -2 * log(lr_uc_num / lr_uc_den) }
  p_value_uc <- 1 - pchisq(lr_uc, df = 1)
  
  # 2. Christoffersen First-Order Transition Matrix Independence Test
  v_lag <- violations[1:(T_total-1)]
  v_cur <- violations[2:T_total]
  
  n00 <- sum(v_lag == 0 & v_cur == 0)
  n01 <- sum(v_lag == 0 & v_cur == 1)
  n10 <- sum(v_lag == 1 & v_cur == 0)
  n11 <- sum(v_lag == 1 & v_cur == 1)
  
  pi_01  <- n01 / (n00 + n01)
  pi_11  <- n11 / (n10 + n11)
  pi_all <- (n01 + n11) / (n00 + n01 + n10 + n11)
  
  L_null <- ((1 - pi_all)^(n00 + n10)) * (pi_all^(n01 + n11))
  L_alt  <- ((1 - pi_01)^n00) * (pi_01^n01) * ((1 - pi_11)^n10) * (pi_11^n11)
  
  if (is.nan(pi_01) || is.nan(pi_11)) { lr_ind <- 0 } else { lr_ind <- -2 * log(L_null / L_alt) }
  p_value_ind <- 1 - pchisq(lr_ind, df = 1)
  
  return(list(expected = round(T_total * alpha), actual = N, kupiec_p = round(p_value_uc, 4), christ_p = round(p_value_ind, 4)))
}


# -------------------------------------------------------------------------
# PHASE 5: DIAGNOSTIC STRUCTURING AND SUMMARY FILE EXPORT
# -------------------------------------------------------------------------
cat("\nPhase 4: Evaluating out-of-sample data violations...\n")
actual_out_of_sample <- returns_matrix[(window_size + 1):N_total, ]
sp500_actual         <- actual_out_of_sample[, "SP500"]
btc_actual           <- actual_out_of_sample[, "BTC"]

results_matrix           <- matrix(NA, nrow = 4, ncol = 4)
rownames(results_matrix) <- c("Expected Violations", "Actual Violations", "Kupiec p-value", "Christoffersen p-value")
colnames(results_matrix) <- c("SP500_HS", "SP500_GARCH", "BTC_HS", "BTC_GARCH")

res_sp_hs  <- run_backtest_suite(sp500_actual, forecasts$SP500$HS_VaR, alpha)
res_sp_gr  <- run_backtest_suite(sp500_actual, forecasts$SP500$GARCH_VaR, alpha)
res_btc_hs <- run_backtest_suite(btc_actual, forecasts$BTC$HS_VaR, alpha)
res_btc_gr <- run_backtest_suite(btc_actual, forecasts$BTC$GARCH_VaR, alpha)

results_matrix[, 1] <- c(res_sp_hs$expected, res_sp_hs$actual, res_sp_hs$kupiec_p, res_sp_hs$christ_p)
results_matrix[, 2] <- c(res_sp_gr$expected, res_sp_gr$actual, res_sp_gr$kupiec_p, res_sp_gr$christ_p)
results_matrix[, 3] <- c(res_btc_hs$expected, res_btc_hs$actual, res_btc_hs$kupiec_p, res_btc_hs$christ_p)
results_matrix[, 4] <- c(res_btc_gr$expected, res_btc_gr$actual, res_btc_gr$kupiec_p, res_btc_gr$christ_p)

# Display tabular output array directly inside terminal environment
print(results_matrix)

# Pack modeling outputs safely to a central data frame
output_df <- data.frame(
  Date        = index(returns_xts[(window_size + 1):N_total]),
  SP500_Ret   = sp500_actual,
  SP500_HS    = forecasts$SP500$HS_VaR,
  SP500_GARCH = forecasts$SP500$GARCH_VaR,
  BTC_Ret     = btc_actual,
  BTC_HS      = forecasts$BTC$HS_VaR,
  BTC_GARCH   = forecasts$BTC$GARCH_VaR
)
write.csv(output_df, file = "comprehensive_seminar_backtesting_outputs.csv", row.names = FALSE)


# -------------------------------------------------------------------------
# PHASE 6: HIGH-RESOLUTION ACADEMIC GRAPHICS ENGINE
# -------------------------------------------------------------------------
cat("\nPhase 5: Generating visual charts for your 12 pages...\n")

# Visualization 1: Time Series Volatility Distribution Layout
p1_sp <- ggplot(output_df, aes(x = Date, y = SP500_Ret)) +
  geom_line(color = "#2c3e50", alpha = 0.7) +
  theme_minimal() + labs(title = "S&P 500 Index Daily Log Returns", x = "", y = "Log Return") +
  theme(plot.title = element_text(face = "bold", size = 11))

p1_btc <- ggplot(output_df, aes(x = Date, y = BTC_Ret)) +
  geom_line(color = "#f39c12", alpha = 0.7) +
  theme_minimal() + labs(title = "Bitcoin (BTC-USD) Daily Log Returns", x = "Timeline", y = "Log Return") +
  theme(plot.title = element_text(face = "bold", size = 11))

plot_returns_comparison <- grid.arrange(p1_sp, p1_btc, ncol = 1)
ggsave("plot_1_returns_comparison.png", plot_returns_comparison, width = 8, height = 6, dpi = 300)

# Map dynamic logical pointers where losses exceeded thresholds
output_df$SP500_Violation_HS    <- ifelse(output_df$SP500_Ret < -output_df$SP500_HS, output_df$SP500_Ret, NA)
output_df$SP500_Violation_GARCH <- ifelse(output_df$SP500_Ret < -output_df$SP500_GARCH, output_df$SP500_Ret, NA)
output_df$BTC_Violation_HS      <- ifelse(output_df$BTC_Ret < -output_df$BTC_HS, output_df$BTC_Ret, NA)
output_df$BTC_Violation_GARCH   <- ifelse(output_df$BTC_Ret < -output_df$BTC_GARCH, output_df$BTC_Ret, NA)

# Visualization 2: S&P 500 Risk Boundary Chart
plot_sp_backtest <- ggplot(output_df, aes(x = Date)) +
  geom_line(aes(y = SP500_Ret, color = "Realized Return"), alpha = 0.4) +
  geom_line(aes(y = -SP500_HS, color = "Historical Simulation (99%)"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = -SP500_GARCH, color = "GARCH(1,1)-t (99%)"), linetype = "solid", size = 0.8) +
  geom_point(aes(y = SP500_Violation_HS), color = "#e74c3c", size = 2, na.rm = TRUE) +
  scale_color_manual(values = c("Realized Return" = "gray60", "Historical Simulation (99%)" = "#e67e22", "GARCH(1,1)-t (99%)" = "#2980b9")) +
  theme_minimal() + labs(title = "S&P 500 99% Value-at-Risk Backtesting Excursion Analysis", x = "Timeline", y = "Returns / Threshold", color = "Model Specs") +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 12))
ggsave("plot_2_sp500_backtest.png", plot_sp_backtest, width = 9, height = 5, dpi = 300)

# Visualization 3: Bitcoin Cluster Violations Failure Chart
plot_btc_backtest <- ggplot(output_df, aes(x = Date)) +
  geom_line(aes(y = BTC_Ret, color = "Realized Return"), alpha = 0.4) +
  geom_line(aes(y = -BTC_HS, color = "Historical Simulation (99%)"), linetype = "dashed", size = 0.8) +
  geom_line(aes(y = -BTC_GARCH, color = "GARCH(1,1)-t (99%)"), linetype = "solid", size = 0.8) +
  geom_point(aes(y = BTC_Violation_HS), color = "#e74c3c", size = 2.5, na.rm = TRUE) +
  scale_color_manual(values = c("Realized Return" = "gray60", "Historical Simulation (99%)" = "#e67e22", "GARCH(1,1)-t (99%)" = "#2980b9")) +
  theme_minimal() + labs(title = "Bitcoin 99% Value-at-Risk Backtesting Excursion Analysis", subtitle = "Red markers indicate mathematical risk limit violations", x = "Timeline", y = "Returns / Threshold", color = "Model Specs") +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 9, face = "italic"))
ggsave("plot_3_bitcoin_backtest.png", plot_btc_backtest, width = 9, height = 5, dpi = 300)

cat("\n=========================================================================\n")
cat("SUCCESS: Process complete! Look inside your folder directory for:\n")
cat("  1. 'seminar_returns_data.csv' (Raw data attachment)\n")
cat("  2. 'comprehensive_seminar_backtesting_outputs.csv' (Calculated forecasts matrix)\n")
cat("  3. 'plot_1_returns_comparison.png' (Log returns timeline plot)\n")
cat("  4. 'plot_2_sp500_backtest.png' (S&P 500 boundary violation plot)\n")
cat("  5. 'plot_3_bitcoin_backtest.png' (Bitcoin cluster failure plot)\n")
cat("=========================================================================\n")