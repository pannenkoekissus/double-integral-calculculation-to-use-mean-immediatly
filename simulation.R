# =========================================================
# Ecological Simulation of the Mutualistic Network
# Simulating population dynamics over T_max generations
# =========================================================

# --- 1. Base Parameters ---
T_max <- 100 # Adjust this to simulate more or fewer generations

mu_P <- c(1:15) # Corollas of varied lengths
S_P <- length(mu_P) 
P_t <- rep(100, S_P) # Initial Plant Abundances
sigma_P <- rep(1, S_P)
r_P <- rep(0.1, S_P)
m_P <- rep(0.01, S_P) 
b_P <- rep(1.5, S_P)
M_P0 <- rep(10, S_P)
c_cost <- 0.1

mu_A <- c(1:20) # Proboscises of varied lengths
S_A <- length(mu_A) 
A_t <- rep(50, S_A) # Initial Animal Abundances
sigma_A <- rep(1.2, S_A)
b_A <- rep(2.0, S_A)
M_A0 <- rep(15, S_A)
C_A0 <- 0.05
mu_A0 <- 10

sigma_c <- 1

# --- 2. Calculate Continuous Interaction Matrices ---
# Note: Because traits don't evolve in this specific script, we only 
# need to calculate their overlapping capacities once!
d_AP <- outer(mu_A, mu_P, "-") 
var_sum_AP <- outer(sigma_A^2, sigma_P^2, "+")
sd_sum_AP <- sqrt(var_sum_AP)
H_int <- pnorm(d_AP / sd_sum_AP)

d_PA <- outer(mu_P, mu_A, "-") 
var_sum_PA <- outer(sigma_P^2, sigma_A^2, "+")
var_total <- sigma_c^2 + var_sum_PA
sd_sum_PA <- sqrt(var_sum_PA)
sd_total <- sqrt(var_total)

part1 <- pnorm(d_PA / sd_sum_PA)
part2 <- (sigma_c / sd_total) * exp(-(d_PA^2) / (2 * var_total)) * pnorm((d_PA * sigma_c) / (sd_sum_PA * sd_total))
f_int <- part1 + part2

# --- 3. Set Up History Matrices ---
P_history <- matrix(0, nrow = T_max, ncol = S_P)
A_history <- matrix(0, nrow = T_max, ncol = S_A)

# Inject starting situations
P_history[1, ] <- P_t
A_history[1, ] <- A_t

# --- 4. Simulation Loop ---
for (t in 1:(T_max - 1)) {
  # Grab the populations from the current generation
  curr_P <- P_history[t, ]
  curr_A <- A_history[t, ]
  
  # === Plant Math ===
  denom_H <- as.vector(H_int %*% curr_P) + 1e-9
  M_C_A <- sweep(H_int, 1, denom_H, "/")
  M_C_A <- sweep(M_C_A, 1, curr_A, "*")
  C_P_val <- c_cost * colSums(M_C_A)
  
  denom_f <- as.vector(t(f_int) %*% curr_P) + 1e-9
  M_P_matrix <- sweep(f_int, 2, denom_f, "/")
  M_P_matrix <- sweep(M_P_matrix, 2, curr_A, "*")
  M_P_total <- rowSums(M_P_matrix)
  B_P_val <- b_P * (M_P_total / (M_P_total + M_P0))
  
  W_P <- exp(r_P + B_P_val - C_P_val - m_P * curr_P)
  
  # === Animal Math ===
  denom_H_anim <- as.vector(t(H_int) %*% curr_A) + 1e-9
  M_A_matrix <- sweep(H_int, 2, denom_H_anim, "/")
  M_A_matrix <- sweep(M_A_matrix, 2, curr_P, "*")
  M_A_total <- rowSums(M_A_matrix)
  B_A_val <- b_A * (M_A_total / (M_A_total + M_A0))
  C_A_val <- C_A0 * exp((mu_A / mu_A0) + 0.5 * (sigma_A / mu_A0)^2)
  W_A <- exp(B_A_val - C_A_val)
  
  # === Apply Fitness Growth to Next Generation ===
  P_history[t+1, ] <- curr_P * W_P
  A_history[t+1, ] <- curr_A * W_A
}

# --- 5. Export Results to Plotted Image ---
png("simulation_plot.png", width = 1200, height = 600, res=100)
par(mfrow=c(1,2))

matplot(1:T_max, P_history, type = "l", lty = 1, lwd = 2, col = rainbow(S_P),
        xlab = "Time (Generations)", ylab = "Population Size",
        main = sprintf("Population Dynamics of %d Plant Species", S_P))

matplot(1:T_max, A_history, type = "l", lty = 1, lwd = 2, col = rainbow(S_A),
        xlab = "Time (Generations)", ylab = "Population Size",
        main = sprintf("Population Dynamics of %d Animal Species", S_A))

dev.off()
print("Simulation complete! Check the new 'simulation_plot.png' image file in your folder.")
