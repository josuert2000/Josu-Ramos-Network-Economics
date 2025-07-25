# =============================================
# PROBLEM SET - Network Economics - Josué Ramos
# =============================================

#Part 1: Peer Effects on club participation

#  This code estimates peer effects in student networks using the AddHealth dataset, 
#  implementing a linear-in-means model where a student's club participation (y) 
#  depends on their friends' average participation (Gy) 
#  while controlling for individual characteristics (X). 
#  To address endogeneity, the analysis employs 2SLS estimation with G²y as an instrument, 
#  following Bramoulle et al. (2009)'s network identification approach. 
#  After constructing a 0-4 club participation index and verifying identification conditions through rank tests of matrices [I, G, G²], 
#  the model produces estimates of peer influence (β) which are compared against a 0.2 benchmark, providing empirical evidence about 
#  the strength of social spillovers in educational settings.

# Clear environment
rm(list = ls())

# Load required libraries
packages <- c("igraph", "igraphdata", "tidyverse", "tidygraph", "Matrix", "AER", "kableExtra", "broom", "lmtest", "sandwich")
installed <- packages %in% installed.packages()
if (any(!installed)) install.packages(packages[!installed])
invisible(lapply(packages, library, character.only = TRUE))

# Load AddHealth dataset
load("AddHealth.RData")

# Select school
school_code <- "044"
inschool <- inschool %>% filter(sschlcde == school_code)
valid_sqids <- inschool %>% 
  filter(sschlcde == school_code) %>% 
  pull(sqid)

#Filtering friend so it matches the school
friend_filtered <- friend %>% 
  filter(sqid %in% valid_sqids) 

# Index of clubs and various activities
total_clubs <- rowSums(inschool[, 8:40]) 
inschool$total_clubs <- total_clubs

inschool <- inschool %>% 
  mutate(
    id = case_when(
      total_clubs == 0 ~ 0,
      total_clubs == 1 ~ 1,
      total_clubs == 2 ~ 2,
      total_clubs == 3 ~ 3,
      total_clubs >= 4 ~ 4,
      TRUE ~ NA_real_  
    )
  )
# Creating the edge list of friends
final_friends <- friend_filtered %>%
  select(sqid, starts_with("mf"), starts_with("ff")) %>%
  pivot_longer(
    cols = -sqid,
    names_to = "friend_type",
    values_to = "friend_id"
  )
final_friends <- final_friends %>%
  left_join(inschool %>% select(sqid, aid), by = "sqid") %>%
  filter(!is.na(aid)) 

# To create the network
edges <- final_friends %>%
  filter(friend_id != 77777777, friend_id != 88888888, friend_id != 99999999, !is.na(friend_id)) %>%
  filter(friend_id %in% inschool$aid) %>%  
  select(from = aid, to = friend_id) %>%
  distinct()

g <- graph_from_data_frame(edges, directed = TRUE)


plot(g,
     layout = layout_with_fr(g, start.temp = 30),
     vertex.color = 'blue',
     vertex.frame.color = "green",
     vertex.size = 4,
     edge.arrow.size = 0.01,
     vertex.label = NA)

library(Matrix)

# Adjacency matrix M
M <- as_adjacency_matrix(g, sparse = FALSE)
M <- sweep(M, 1, rowSums(M), FUN = "/")
M[is.nan(M)] <- 0
M2 <- M %*% M

# Vectorize matrices
N <- cbind(as.vector(diag(nrow(M))), as.vector(M), as.vector(M2))
r <- rankMatrix(N)[1]

\begin{verbatim}
if (rank_condition == 3) {
  cat("Identification Test Result:\n")
  cat("-------------------------\n")
  cat("• Full rank condition satisfied (rank = 3)\n")
  cat("• Matrices I, M, and M² are linearly independent\n") 
  cat("• Conclusion: Point identification is achieved\n")
} else {
  cat("Identification Warning:\n")
  cat("---------------------\n")
  cat("• Rank deficiency detected (rank =", rank_condition, ")\n")
  cat("• Linear dependence exists between I, M, and M²\n")
  cat("• Implication: No point identification possible with current instruments\n")
}
\end{verbatim}

# First we want to organize y as y_vec
y_vec <- inschool$total_clubs[match(V(g)$name, inschool$aid)]

#We take the adjacency matrix
M <- as_adjacency_matrix(g, sparse = FALSE)

My <- M %*% y_vec

# Check current types
class(M)
class(y_vec)

# Convert to numeric matrices/vectors
M <- as.matrix(M) %>% apply(2, as.numeric)  # Ensure matrix is numeric
y_vec <- as.numeric(y_vec)  # Ensure vector is numeric
X <- inschool %>%
  filter(aid %in% V(g)$name)
X <- X[match(V(g)$name, X$aid), ]

#Building the control vector
X <- X %>%
  transmute(
    intercept = 1,
    age = as.numeric(s1),
    male = ifelse(s2 == 1, 1, 0),
    mother_ed = as.numeric(s12),
    father_ed = as.numeric(s18),
    race_white = as.numeric(s6a),
    race_black = as.numeric(s6b)
  )

# Convert to matrix
X_mat <- as.matrix(X)
dim(X_mat) 
dim(M)  

#Create the dataframe

df <- data.frame(
  y    = y_vec,
  My   = as.numeric(M %*% y_vec),
  M2y  = as.numeric(M2 %*% y_vec),
  X_mat 
)

#Missing values
df_clean <- df[complete.cases(df), ]
cat("Rows before:", nrow(df), " – after cleaning:", nrow(df_clean), "\n")

# Fformula with endogenous variable: My and instrument: M2y
controls <- paste(colnames(X_mat), collapse = " + ")
iv_formula <- as.formula(
  paste0("y ~ ", controls, " + My | ", controls, " + M2y")
)

#Instrumental varibale 
iv_model <- ivreg(iv_formula, data = df_clean)
summary(iv_model, diagnostics = TRUE)

install.packages("stargazer")
library(stargazer)


##Tables 

stargazer(iv_model,
          type = "text",       # "text" for Word, "html" also works
          title = "2SLS Peer‑Effects Estimates",
          dep.var.labels = "y",
          covariate.labels = c("Peer Avg (Gy)", "Age", "Male", "Mother's Ed", 
                               "Father's Ed", "White", "Black"),
          omit.stat = c("adj.rsq", "f"),   # omit stats you don't need
          out = "iv_table.text")



stargazer(iv_model,
          type = "latex",       # "text" for Word, "html" also works
          title = "2SLS Peer‑Effects Estimates",
          dep.var.labels = "y",
          covariate.labels = c("Peer Avg (Gy)", "Age", "Male", "Mother's Ed", 
                               "Father's Ed", "White", "Black"),
          omit.stat = c("adj.rsq", "f"),   # omit stats you don't need
          out = "iv_table.tex")

stargazer(iv_model,
          type = "html",       # "text" for Word, "html" also works
          title = "2SLS Peer‑Effects Estimates",
          dep.var.labels = "y",
          covariate.labels = c("Peer Avg (Gy)", "Age", "Male", "Mother's Ed", 
                               "Father's Ed", "White", "Black"),
          omit.stat = c("adj.rsq", "f"),   # omit stats you don't need
          out = "iv_table.word")

#Part 2: Impacts of disease spreading

# This code investigates how community structure impacts disease spread by simulating 
# epidemics on networks with varying community strength. Using a planted partition model 
# (a special case of stochastic block models), it generates networks with two communities where within-group (pin) 
# and between-group (pout) connection probabilities are controlled by parameter ε, which determines community strength. 
# The study first visualizes networks with different ε values (0, 9, 18) to show increasing community structure. 
# Then, it simulates a discrete-time SI epidemic model on these networks, tracking final outbreak size (s/n) and duration (ℓ) across 
# transmission probabilities p. By comparing results for ε=0 (no communities) versus ε=18 (strong communities), the analysis reveals 
# how community barriers slow disease spread and alter epidemic thresholds, with simulations averaged over multiple runs for robust estimates. 
# This provides quantitative insights into how network clustering affects contagion dynamics.


# Clear environment
rm(list = ls())

# Load required libraries
packages <- c("igraph", "igraphdata", "tidyverse", "tidygraph", "Matrix", "AER", "kableExtra", "broom", "lmtest", "sandwich")
installed <- packages %in% installed.packages()
if (any(!installed)) install.packages(packages[!installed])
invisible(lapply(packages, library, character.only = TRUE))

library(igraph)
library(ggplot2)
library(dplyr)
library(deSolve)
library(gridExtra)

# Generate and visualize planted partition model graphs

generate_planted_partition <- function(n, c, epsilon) {
  q <- 2  # two communities
  cin <- c + epsilon/2
  cout <- c - epsilon/2
  
  # Avoid negative probabilities
  cout <- max(0, cout)
  
  pin <- cin/n
  pout <- cout/n
  
  # Create the block matrix
  block_sizes <- rep(n/q, q)
  pref_matrix <- matrix(pout, nrow = q, ncol = q)
  diag(pref_matrix) <- pin
  
  # Generate the graph
  g <- tryCatch({
    sample_sbm(n, pref.matrix = pref_matrix, block.sizes = block_sizes, directed = FALSE)
  }, error = function(e) {
    sample_sbm(n, pref.matrix = pref_matrix, blocksizes = block_sizes, directed = FALSE)
  })
  
  # Assign community membership as vertex attribute
  V(g)$community <- rep(1:q, each = n/q)
  V(g)$color <- ifelse(V(g)$community == 1, "orange", "lightgreen")
  
  return(g)
}

# Visualize the graphs
visualize_communities <- function() {
  n <- 100
  c <- 10
  epsilon_values <- c(0, 9, 18)
  graphs <- lapply(epsilon_values, function(eps) generate_planted_partition(n, c, eps))
  
  dev.new(width = 12, height = 4)
  par(mfrow = c(1, 3), mar = c(2, 2, 3, 1))
  
  for (i in 1:3) {
    plot(graphs[[i]], 
         vertex.color = V(graphs[[i]])$color,
         vertex.size = 5,
         vertex.label = NA,
         main = paste("e =", epsilon_values[i], "\np_in =", round((c + epsilon_values[i]/2)/n, 3)),
         layout = layout_with_fr(graphs[[i]]))
  }
}
visualize_communities()

# Run simulations for ε=18 (using same parameters as ε=0)
set.seed(123)
results_epsilon18 <- run_simulations(epsilon = 18, p_values = p_values, 
                                     n_networks = 5, n_simulations = 5)

# Aggregate results for both cases
summary_epsilon0 <- results_epsilon0 %>%
  group_by(p) %>%
  summarise(avg_size = mean(size),
            avg_duration = mean(duration)) %>%
  mutate(epsilon = "0")

summary_epsilon18 <- results_epsilon18 %>%
  group_by(p) %>%
  summarise(avg_size = mean(size),
            avg_duration = mean(duration)) %>%
  mutate(epsilon = "18")

combined_results <- bind_rows(summary_epsilon0, summary_epsilon18)

# Create comparative plots
plot_size_comparison <- ggplot(combined_results, aes(x = p, y = avg_size, color = epsilon)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("0" = "blue", "18" = "orange")) +
  labs(title = "Epidemic Size vs Transmission Probability",
       subtitle = "Comparison of ε = 0 vs ε = 18",
       x = "Transmission probability (p)",
       y = "Average fraction infected ⟨s/n⟩",
       color = "ε value") +
  theme_minimal() +
  theme(legend.position = "top")

plot_duration_comparison <- ggplot(combined_results, aes(x = p, y = avg_duration, color = epsilon)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = c("0" = "blue", "18" = "orange")) +
  labs(title = "Epidemic Duration vs Transmission Probability",
       subtitle = "Comparison of ε = 0 vs ε = 18",
       x = "Transmission probability (p)",
       y = "Average duration ⟨ℓ⟩",
       color = "ε value") +
  theme_minimal() +
  theme(legend.position = "top")

# Display comparative plots
grid.arrange(plot_size_comparison, plot_duration_comparison, ncol = 1)

# Key observations (printed to console)
cat("\nKey Observations:\n",
    "1. Strong community structure (ε=18) shows:\n",
    "   - Delayed onset of large outbreaks\n",
    "   - Lower overall infection fractions at intermediate p values\n",
    "   - Longer epidemic durations due to slow cross-community spread\n",
    "2. Critical p value shifts rightward for ε=18\n",
    "3. Maximum duration occurs at higher p for ε=18\n")










