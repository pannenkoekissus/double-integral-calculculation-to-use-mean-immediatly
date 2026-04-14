# --- 1. Base Parameters ---
T_max <- 1000

mu_P_initial <- c(1:15)
S_P <- length(mu_P_initial)
P_t_initial <- rep(100, S_P)
sigma_P <- rep(1, S_P)
r_P <- rep(0.1, S_P)
m_P <- rep(0.01, S_P)
b_P <- rep(1.5, S_P)
M_P0 <- rep(10, S_P)
c_cost <- 0.1
h2_P <- 0.5 # Heritability for Plant traits

mu_A_initial <- c(1:20)
S_A <- length(mu_A_initial)
A_t_initial <- rep(50, S_A)
sigma_A <- rep(1.2, S_A)
b_A <- rep(2.0, S_A)
M_A0 <- rep(15, S_A)
C_A0 <- 0.05
mu_A0 <- 10
h2_A <- 0.5 # Heritability for Animal traits

sigma_c <- 1

# --- 2. Fitness Evaluation Function ---
calc_fitness <- function(mu_P, mu_A, P_t, A_t) {
  # Interaction Matrices
  d_AP <- outer(mu_A, mu_P, "-")
  var_sum_AP <- outer(sigma_A^2, sigma_P^2, "+")
  H_int <- pnorm(d_AP / sqrt(var_sum_AP))

  d_PA <- outer(mu_P, mu_A, "-")
  var_sum_PA <- outer(sigma_P^2, sigma_A^2, "+")
  var_total <- sigma_c^2 + var_sum_PA
  sd_sum_PA <- sqrt(var_sum_PA)
  sd_total <- sqrt(var_total)

  f_int <- pnorm(d_PA / sd_sum_PA) +
    (sigma_c / sd_total) * exp(-(d_PA^2) / (2 * var_total)) *
      pnorm((d_PA * sigma_c) / (sd_sum_PA * sd_total))

  # Plant Math
  denom_H <- as.vector(H_int %*% P_t) + 1e-9
  M_C_A <- sweep(H_int, 1, denom_H, "/")
  M_C_A <- sweep(M_C_A, 1, A_t, "*")
  C_P_val <- c_cost * colSums(M_C_A)

  denom_f <- as.vector(t(f_int) %*% P_t) + 1e-9
  M_P_matrix <- sweep(f_int, 2, denom_f, "/")
  M_P_matrix <- sweep(M_P_matrix, 2, A_t, "*")
  B_P_val <- b_P * (rowSums(M_P_matrix) / (rowSums(M_P_matrix) + M_P0))

  W_P <- exp(r_P + B_P_val - C_P_val - m_P * P_t)

  # Animal Math
  denom_H_anim <- as.vector(t(H_int) %*% A_t) + 1e-9
  M_A_matrix <- sweep(H_int, 2, denom_H_anim, "/")
  M_A_matrix <- sweep(M_A_matrix, 2, P_t, "*")
  B_A_val <- b_A * (rowSums(M_A_matrix) / (rowSums(M_A_matrix) + M_A0))
  C_A_val <- C_A0 * exp((mu_A / mu_A0) + 0.5 * (sigma_A / mu_A0)^2)

  W_A <- exp(B_A_val - C_A_val)

  return(list(W_P = W_P, W_A = W_A))
}

# --- 3. Set Up History Matrices ---
P_history <- matrix(0, nrow = T_max, ncol = S_P)
A_history <- matrix(0, nrow = T_max, ncol = S_A)
mu_P_history <- matrix(0, nrow = T_max, ncol = S_P)
mu_A_history <- matrix(0, nrow = T_max, ncol = S_A)

P_history[1, ] <- P_t_initial
A_history[1, ] <- A_t_initial
mu_P_history[1, ] <- mu_P_initial
mu_A_history[1, ] <- mu_A_initial

# --- 4. Eco-Evolutionary Simulation Loop ---
delta <- 1e-5 # Finite difference step for numeric gradients

for (t in 1:(T_max - 1)) {
  curr_P <- P_history[t, ]
  curr_A <- A_history[t, ]
  curr_mu_P <- mu_P_history[t, ]
  curr_mu_A <- mu_A_history[t, ]

  # Base fitness evaluation
  res <- calc_fitness(curr_mu_P, curr_mu_A, curr_P, curr_A)
  W_P <- res$W_P
  W_A <- res$W_A

  # Calculate selection gradients for Plant traits
  grad_W_P <- numeric(S_P)
  for (i in 1:S_P) {
    pert_mu_P <- curr_mu_P
    pert_mu_P[i] <- pert_mu_P[i] + delta
    W_P_new <- calc_fitness(pert_mu_P, curr_mu_A, curr_P, curr_A)$W_P
    grad_W_P[i] <- (W_P_new[i] - W_P[i]) / delta
  }

  # Calculate selection gradients for Animal traits
  grad_W_A <- numeric(S_A)
  for (j in 1:S_A) {
    pert_mu_A <- curr_mu_A
    pert_mu_A[j] <- pert_mu_A[j] + delta
    W_A_new <- calc_fitness(curr_mu_P, pert_mu_A, curr_P, curr_A)$W_A
    grad_W_A[j] <- (W_A_new[j] - W_A[j]) / delta
  }

  # Update Trait Means using the Quantitative Genetics equation:
  curr_mu_P <- curr_mu_P + h2_P * (sigma_P^2) * (1 / W_P) * grad_W_P
  curr_mu_A <- curr_mu_A + h2_A * (sigma_A^2) * (1 / W_A) * grad_W_A

  # Apply Eco-growth to next gen
  P_history[t + 1, ] <- curr_P * W_P
  A_history[t + 1, ] <- curr_A * W_A
  mu_P_history[t + 1, ] <- curr_mu_P
  mu_A_history[t + 1, ] <- curr_mu_A
}

# --- 5. Export Results ---
png("eco_evo_plot.png", width = 1200, height = 1000, res = 100)
par(mfrow = c(2, 2))

matplot(1:T_max, P_history,
  type = "l", lty = 1, lwd = 2, col = rainbow(S_P),
  xlab = "Generations", ylab = "Abundance", main = "Plant Populations"
)

matplot(1:T_max, mu_P_history,
  type = "l", lty = 1, lwd = 2, col = rainbow(S_P),
  xlab = "Generations", ylab = "Mean Corolla Size", main = "Plant Trait Evolution"
)

matplot(1:T_max, A_history,
  type = "l", lty = 1, lwd = 2, col = rainbow(S_A),
  xlab = "Generations", ylab = "Abundance", main = "Animal Populations"
)

matplot(1:T_max, mu_A_history,
  type = "l", lty = 1, lwd = 2, col = rainbow(S_A),
  xlab = "Generations", ylab = "Mean Proboscis Size", main = "Animal Trait Evolution"
)

dev.off()
print("Simulation complete! Check the new 'eco_evo_plot.png' image file.")
print("Final Plant Traits at T_max:")
print(round(mu_P_history[T_max, ], 3))
print("Final Animal Traits at T_max:")
print(round(mu_A_history[T_max, ], 3))
