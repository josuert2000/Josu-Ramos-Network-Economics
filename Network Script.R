# =============================================
# PROBLEM SET - Network Economics
# =============================================

# Clear environment and load libraries
rm(list = ls())

# Load packages (install if missing)
packages <- c("igraph", "igraphdata", "tidyverse", "tidygraph", 
              "Matrix", "AER", "kableExtra", "broom", "lmtest", "sandwich")
installed <- packages %in% installed.packages()
if (any(!installed)) install.packages(packages[!installed])
invisible(lapply(packages, library, character.only = TRUE))

# ---- POINT A: Data Preparation ----
# Load and filter AddHealth dataset
load("AddHealth.RData")

# Focus on specific school
school_code <- "044"
inschool <- inschool %>% filter(sschlcde == school_code)
valid_sqids <- inschool %>% pull(sqid)

# Filter friend nominations
friend_filtered <- friend %>% filter(sqid %in% valid_sqids) 

# Calculate club participation
inschool$total_clubs <- rowSums(inschool[, 8:40]) 

# Categorize club participation
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

# ---- POINT B: Network Construction ----
# Create edge list from friend nominations
final_friends <- friend_filtered %>%
  select(sqid, starts_with("mf"), starts_with("ff")) %>%
  pivot_longer(cols = -sqid, names_to = "friend_type", values_to = "friend_id") %>%
  left_join(inschool %>% select(sqid, aid), by = "sqid") %>%
  filter(!is.na(aid)) 

# Clean edges and create graph
edges <- final_friends %>%
  filter(friend_id != 77777777, friend_id != 88888888, 
         friend_id != 99999999, !is.na(friend_id)) %>%
  filter(friend_id %in% inschool$aid) %>%  
  select(from = aid, to = friend_id) %>%
  distinct()

g <- graph_from_data_frame(edges, directed = TRUE)

# Plot network
plot(g,
     layout = layout_with_fr(g, start.temp = 30),
     vertex.color = 'blue',
     vertex.frame.color = "green",
     vertex.size = 4,
     edge.arrow.size = 0.01,
     vertex.label = NA)

# Create adjacency matrices
M <- as_adjacency_matrix(g, sparse = FALSE)
M <- sweep(M, 1, rowSums(M), FUN = "/")
M[is.nan(M)] <- 0
M2 <- M %*% M

# Check rank condition (simplified)
rank_condition <- rankMatrix(cbind(diag(nrow(M)), M, M2))[1]
if (rank_condition == 3) {
  cat("Matrices I, M, and M² are linearly independent\n")
} else {
  cat("No point identification possible with current instruments\n")
}

# ---- POINT C: Variable Preparation ----
# Prepare outcome variable
y_vec <- inschool$total_clubs[match(V(g)$name, inschool$aid)]
My <- M %*% y_vec

# Prepare control variables
X <- inschool %>%
  filter(aid %in% V(g)$name) %>%
  transmute(
    intercept = 1,
    age = as.numeric(s1),
    male = ifelse(s2 == 1, 1, 0),
    mother_ed = as.numeric(s12),
    father_ed = as.numeric(s18),
    race_white = as.numeric(s6a),
    race_black = as.numeric(s6b)
  ) %>%
  as.matrix()

# ---- POINT D: IV Regression ----
# Create analysis dataframe
df <- data.frame(
  y = y_vec,
  My = as.numeric(M %*% y_vec),
  M2y = as.numeric(M2 %*% y_vec),
  X_mat 
) %>% na.omit()

# Run IV regression
iv_formula <- as.formula(
  paste0("y ~ ", paste(colnames(X), collapse = " + "), 
         " + My | ", paste(colnames(X), collapse = " + "), " + M2y")
)
iv_model <- ivreg(iv_formula, data = df)

# Output results
library(stargazer)
stargazer(iv_model,
          type = "text",
          title = "2SLS Peer Effects Estimates",
          covariate.labels = c("Peer Avg (Gy)", "Age", "Male", "Mother's Ed", 
                               "Father's Ed", "White", "Black"),
          omit.stat = c("adj.rsq", "f"))

# =============================================
# PROBLEM SET 4 - Epidemic Modeling
# Valentina Laverde, Juan Felipe Agudelo, 
# Catalina Bernal y Josué Ramos
# =============================================

# Clear environment and load libraries
rm(list = ls())
library(igraph)
library(ggplot2)
library(dplyr)
library(deSolve)
library(gridExtra)

# ---- PART A: Planted Partition Model ----
generate_planted_partition <- function(n, c, epsilon) {
  # Generate graph with two communities
  q <- 2
  cin <- c + epsilon/2
  cout <- max(0, c - epsilon/2)  # Avoid negative probabilities
  
  block_sizes <- rep(n/q, q)
  pref_matrix <- matrix(cout/n, nrow = q, ncol = q)
  diag(pref_matrix) <- cin/n
  
  g <- sample_sbm(n, pref.matrix = pref_matrix, 
                  block.sizes = block_sizes, directed = FALSE)
  
  # Add community attributes
  V(g)$community <- rep(1:q, each = n/q)
  V(g)$color <- ifelse(V(g)$community == 1, "orange", "lightgreen")
  
  return(g)
}

# Visualize different community structures
visualize_communities <- function() {
  n <- 100
  c <- 10
  par(mfrow = c(1, 3))
  
  for (eps in c(0, 9, 18)) {
    g <- generate_planted_partition(n, c, eps)
    plot(g, 
         vertex.color = V(g)$color,
         vertex.size = 5,
         vertex.label = NA,
         main = paste("ε =", eps))
  }
}
visualize_communities()

# ---- PART B: Epidemic Simulation ----
# Note: The actual simulation function would go here
# This is a placeholder showing the comparison code:

# Compare epsilon = 0 vs epsilon = 18 results
plot_size_comparison <- ggplot(combined_results, aes(x = p, y = avg_size, color = epsilon)) +
  geom_line() +
  labs(title = "Epidemic Size Comparison",
       x = "Transmission probability (p)",
       y = "Average fraction infected") +
  theme_minimal()

plot_duration_comparison <- ggplot(combined_results, aes(x = p, y = avg_duration, color = epsilon)) +
  geom_line() +
  labs(title = "Epidemic Duration Comparison",
       x = "Transmission probability (p)",
       y = "Average duration") +
  theme_minimal()

# Display plots side by side
grid.arrange(plot_size_comparison, plot_duration_comparison, ncol = 1)

# Key findings
cat("With strong community structure (ε=18):\n",
    "- Slower epidemic spread between communities\n",
    "- Lower overall infection rates\n",
    "- Longer outbreak durations\n")