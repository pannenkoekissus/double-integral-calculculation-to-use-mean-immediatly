# New exact parameters from the paper text
c_cost <- 0.1 # proportionality constant for cost 'c'
A_t <- 50     # Animal abundance
M_P0 <- 10    # Half-saturation constant for Plant benefits
b_P <- 1.5    # Max benefit for Plant
sigma_c <- 1  # Attenuation of benefit

r_P <- 0.1
m_P <- 0.01
P_t <- 100

mu_P <- 5
sigma_P <- 1
mu_A <- 6
sigma_A <- 1.5

# The Analytical Integration Method

# Variance sums shared across functions
var_sum <- sigma_P^2 + sigma_A^2
sd_sum <- sqrt(var_sum)

var_total <- sigma_c^2 + var_sum
sd_total <- sqrt(var_total)

# 1. Integrated Heaviside H_int replacing the raw step function
d_AP <- mu_A - mu_P
H_int <- pnorm(d_AP / sd_sum)

# 2. Integrated f_int
d_PA <- mu_P - mu_A
part1 <- pnorm(d_PA / sd_sum)
part2 <- (sigma_c / sd_total) * exp(-(d_PA^2) / (2 * var_total)) * pnorm((d_PA * sigma_c) / (sd_sum * sd_total))

f_int <- part1 + part2

# 3. Calculate Benefits and Costs using the analytical smooth terms
# We replace the discontinuous H and f with H_int and f_int to achieve the intermediate model!

denom_f <- P_t * f_int + 1e-9
M_P_val <- (f_int / denom_f) * A_t
B_P_val <- b_P * (M_P_val / (M_P_val + M_P0))

denom_H <- P_t * H_int + 1e-9
M_C_A <- (H_int / denom_H) * A_t
C_P_val <- c_cost * M_C_A

# 4. Final Growth Factor
analytical_growth_factor <- exp(r_P + B_P_val - C_P_val - m_P * P_t)

print(sprintf("Analytical Intermediate Approach Result: %.7f", analytical_growth_factor))
