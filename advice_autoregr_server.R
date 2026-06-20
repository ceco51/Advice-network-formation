.libPaths(c("~/R/libs", .libPaths()))

# Install anything missing on the fly
pkgs <- c("codetools", "iterators", "foreach", "doParallel", 
          "igraph", "data.table", "writexl", "readr")
missing <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(missing) > 0) {
  install.packages(missing, lib="~/R/libs", 
                   repos="https://cloud.r-project.org", 
                   dependencies=TRUE)
}

library(igraph)
#library(netUtils)
library(writexl)
library(data.table)
library(parallel)
library(doParallel)
library(foreach)

# Function to run a single simulation
run_simulation <- function(
    sim_id, 
    n, 
    delta, 
    sigma_zero, 
    sigma_T1,
    sigma_T2,
    sigma_T3,
    xi_in, 
    xi_out, 
    p, 
    q, 
    r, 
    theta, 
    alpha_zero_adj, 
    beta_zero_adj
) {
  
  # Create empty network
  adv_net_sim <- make_empty_graph(n)
  
  # Initialize status and skills variables
  alpha <- alpha_zero_adj
  beta <- beta_zero_adj
  
  # Initialize metrics dataframe
  # Plus 1 for the initial state (iteration 0)
  # Saving data only every 10 iterations, out of a total of 1000
  metric_rows <- ceiling(1000/10) + 1
  metrics <- data.frame(
    iteration = numeric(metric_rows),
    n_gurus = numeric(metric_rows),
    n_partners = numeric(metric_rows),
    transitivity_net = numeric(metric_rows),
    centr_in = numeric(metric_rows),
    dens_net = numeric(metric_rows),
    centr_out = numeric(metric_rows),
    edges_net = numeric(metric_rows),
    modularity_net = numeric(metric_rows)
  )
  
  # Initialize first row with iteration 0 data
  metrics[1, ] <- c(0, 0, 0)
  
  # Track the current row in the metrics dataframe
  metric_row_index <- 2
  
  # Pre-compute distance matrix for homophily using outer vector product
  # Convert logical to numeric with * 1
  age_diff <- outer(theta$age_group, theta$age_group, "!=") * 1  
  gender_diff <- outer(theta$gender, theta$gender, "!=") * 1
  office_diff <- outer(theta$office, theta$office, "!=") * 1
  
  # Combine the differences and apply sqrt
  homophily_distance <- sqrt(age_diff + gender_diff + office_diff)
  
  # Set diagonal to 0 (no self-comparisons)
  diag(homophily_distance) <- 0
  
  # Main loop
  for(i in 1:1000) {
    
    # Select a random node to update its connections
    j <- sample(1:n, 1)
    
    # Compute differences between j and all other nodes (vectorized)
    alpha_diffs <- alpha - alpha[j]  # alpha[k] - alpha[j] for all k
    beta_diffs <- beta - beta[j]     # beta[k] - beta[j] for all k
    
    # Extract the homophily distances for node j
    homophily_j <- homophily_distance[j, ]
    
    # Calculate all probabilities at once for j
    # In other words, here we compute p_{jk} for each k != j, k = 1, ..., n
    P <- 1/(1 + exp(-(sigma_zero + 
                        sigma_T1 * alpha_diffs + 
                        sigma_T2 * beta_diffs + 
                        sigma_T3 * homophily_j)))
    
    # Set j's self-probability to 0
    P[j] <- 0
    
    # Update all connections for node j efficiently
    for(k in 1:n) {
      if(k != j) {
        if(!are_adjacent(adv_net_sim, j, k) && P[k] >= delta) {
          adv_net_sim <- add_edges(adv_net_sim, c(j, k))
        } else if(are_adjacent(adv_net_sim, j, k) && P[k] < delta) {
          adv_net_sim <- delete_edges(
                                     adv_net_sim, 
                                     get.edge.ids(adv_net_sim, c(j, k))
                                     )
        }
      }
    }
    
    # Calculate normalized degrees once (vectorized)
    out_degree <- degree(adv_net_sim, V(adv_net_sim), "out", normalized = TRUE)
    in_degree <- degree(adv_net_sim, V(adv_net_sim), "in", normalized = TRUE)
    
    # Update status and expertise levels 
    tmp_gamma <- alpha^p #status loss
    tmp_phi <- alpha^q #status gain
    tmp_lambda <- beta^r #up-skilling
    
    # Update alpha and beta efficiently
    alpha <- pmin(pmax((alpha - out_degree * tmp_gamma + in_degree * tmp_phi), 0), 1)
    beta <- pmin(pmax((beta + out_degree * tmp_lambda), 0), 1)
    
    # Compute network metrics
    # Only compute network metrics every 10 steps
    if(i %% 10 == 0 || i == 1) {
      
      edges_net <- ecount(adv_net_sim)
      dens_net <- edge_density(adv_net_sim)
      transitivity_net <- transitivity(adv_net_sim)
      
      modularity_net <- tryCatch({
        if(edges_net > 0) modularity(cluster_walktrap(adv_net_sim)) else NA
      }, error = function(e) NA)
      
      centr_in <- tryCatch({
        if(edges_net > 0) centr_degree(adv_net_sim, mode = "in")$centralization else NA
      }, error = function(e) NA)
      
      centr_out <- tryCatch({
        if(edges_net > 0) centr_degree(adv_net_sim, mode = "out")$centralization else NA
      }, error = function(e) NA)
      
      # Count gurus
      indeg_condition <- in_degree > xi_in
      outdeg_condition <- out_degree < xi_out
      n_gurus <- sum(indeg_condition & outdeg_condition)
      guru_indices <- which(indeg_condition & outdeg_condition)
      if(length(guru_indices) == 0) {
        n_partners <- 0
      } else {
        n_partners <- sum(theta$status[guru_indices] == 1, na.rm = TRUE)
      }
      #n_partners <- table(theta$status[guru_indices])[1]

      # Update metrics dataframe only when we calculate metrics
      metrics[metric_row_index, ] <- c(i, n_gurus, n_partners,
                                       transitivity_net, centr_in, 
                                       dens_net, centr_out,
                                       edges_net, modularity_net)
      
      # Increment the row index for next time
      metric_row_index <- metric_row_index + 1
    }
  }
  
  # Remove unused rows in the metrics dataframe
  metrics <- metrics[1:(metric_row_index-1), ]
  
  # Add simulation ID and parameters
  metrics$sim <- sim_id
  metrics$p <- p
  metrics$q <- q
  metrics$r <- r
  metrics$sigma_zero <- sigma_zero
  metrics$sigma_T1 <- sigma_T1
  metrics$sigma_T2 <- sigma_T2
  metrics$sigma_T3 <- sigma_T3
  
  return(metrics)
}

# Main execution script
main <- function() {
  # Uploading the file from Lazega (from hardcoded path)
  adv_data <- read.table("/home/users/NETSIM/AutoAdvice/ELadv.dat", 
                         header = FALSE)
  colnames(adv_data) <- as.character(c(1:ncol(adv_data)))
  adv_data <- as.matrix(adv_data)
  
  adv_att <- read.table("/home/users/NETSIM/AutoAdvice/ELattr.dat",
                        header = FALSE)
  colnames(adv_att) <- c("senority", "status", "gender", "office", 
                         "years with the firm", "age", "practice", "law school")
  
  # Create the network from Lazega data
  adv_net_lazega <- graph_from_adjacency_matrix(adv_data, "directed")
  
  # Set parameters
  n <- vcount(adv_net_lazega)
  d_max <- edge_density(adv_net_lazega)
  theta <- adv_att[, c(2,3,4,5,6,8)] # status, gender, office, years with the firm, age, law school
  
  # Create age groups 
  theta$age_group <- cut(theta$age, 
                         breaks = c(0, 35, 45, 55, Inf), 
                         labels = 1:4,
                         right = FALSE)
  theta$age_group <- as.numeric(as.character(theta$age_group))
  
  # Generate initial reputation and expertise levels
  alpha_zero <- runif(n, min = 0, max = 1)
  beta_zero <- runif(n, min = 0, max = 1)
  
  # Order in descending order
  alpha_zero <- alpha_zero[order(alpha_zero, decreasing = TRUE)]
  beta_zero <- beta_zero[order(beta_zero, decreasing = TRUE)]
  
  # Adjust by characteristics
  adjustment_factor1 <- 0.8 
  #adjustment_factor2 <- 0.1
  adjustment_factor3 <- 1.8
  
  alpha_zero_adj <- as.vector(
    alpha_zero + 
      #adjustment_factor1*scale(1/theta$status) + 
      adjustment_factor1*scale(1/theta$`law school`)
  )
  beta_zero_adj <- as.vector(
    beta_zero + 
      #adjustment_factor2*scale(1/theta$status) + 
      adjustment_factor3*scale(as.numeric(theta$`years with the firm`))
  )
  
  # Normalize to [0,1]
  alpha_zero_adj <- (alpha_zero_adj - min(alpha_zero_adj))/(max(alpha_zero_adj) - min(alpha_zero_adj))
  beta_zero_adj <- (beta_zero_adj - min(beta_zero_adj))/(max(beta_zero_adj) - min(beta_zero_adj))
  
  # Check correlation
  cat("Correlation between alpha and beta:", cor(alpha_zero_adj, beta_zero_adj), "\n")
  
  # Should be 0.60 out from (0.8, 1.8) as adjustment factors 

  # Set simulation parameters
   #dislike for heterophily --> Penalty 
  xi_in <- 0.1
  xi_out <- 0.1
  delta <- 0.5
  #p <- 1.2 # Status loss
  #r <- 0.2 #upskilling fixed
  #q <- 0.5 #Status gain
  # Define parameter grid
  param_grid <- expand.grid(
    sigma_T1 = c(0.5, 1, 1.5, 2, 2.5),
    sigma_T2 = c(0.5, 1, 1.5, 2, 2.5),
    sigma_zero = c(-1.5, -1, -0.5, -0.3, -0.1),
    sigma_T3 = c(-1.5, -1, -0.5, -0.3, -0.1),
    q = c(0.05, 0.25, 0.50, 0.75, 0.95),
    r = c(0.2, 0.8)  
  )
  
  
  # Setup parallel processing
  num_cores <- 40
  cl <- makeCluster(num_cores)
  on.exit(stopCluster(cl))
  registerDoParallel(cl)
  
  cat("Running simulations on", num_cores, "cores\n")
  
  # Number of parallel simulations per parameter set
  num_sims <- 20
  
  # Run parallel simulations
  results <- foreach(i = 1:nrow(param_grid), 
                     .packages = c("igraph", "data.table"),
                     .combine = rbind,
                     .export = c("run_simulation")) %:%
    foreach(sim = 1:num_sims, 
            .combine = rbind,
            .packages = c("igraph", "data.table")) %dopar% {
              
              # Extract parameter combination from grid
              # All the others are fixed to the values defined above
              sigma_T1_val <- param_grid$sigma_T1[i]
              sigma_T2_val <- param_grid$sigma_T2[i]
              sigma_zero_val <- param_grid$sigma_zero[i]
              sigma_T3_val <- param_grid$sigma_T3[i]
              q_val <- param_grid$q[i]
              r_val <- param_grid$r[i]

              run_simulation(
                sim_id = paste0("sim_", i, "_", sim),
                n = n,
                delta = delta,
                sigma_zero = sigma_zero_val,
                sigma_T1 = sigma_T1_val, 
                sigma_T2 = sigma_T2_val, 
                sigma_T3 = sigma_T3_val,
                xi_in = xi_in,
                xi_out = xi_out,
                q = q_val,
                p = (1 / q_val),
                r = r_val,
                theta = theta,
                alpha_zero_adj = alpha_zero_adj,
                beta_zero_adj = beta_zero_adj
              )
            }
  
  # Stop the cluster
  stopCluster(cl)
  
  # Convert to data.table for efficient analysis
  metrics_all <- as.data.table(results)
  
  # Calculate means
  metrics_mean <- metrics_all[, lapply(.SD, function(x) mean(x, na.rm = TRUE)), 
                              by = .(iteration, p, q, r, sigma_zero, sigma_T1, sigma_T2, sigma_T3), 
                              .SDcols = c("n_gurus", "n_partners", 
                                          "transitivity_net", "centr_in", 
                                          "dens_net", "centr_out",
                                          "edges_net", "modularity_net")]
  
  # Rename mean columns
  setnames(metrics_mean, 
           old = c("n_gurus", "n_partners", 
                    "transitivity_net", "centr_in", 
                    "dens_net", "centr_out",
                    "edges_net", "modularity_net"),
           new = paste0("mean_", c("n_gurus", "n_partners", 
                                   "transitivity_net", "centr_in", 
                                   "dens_net", "centr_out",
                                   "edges_net", "modularity_net")))
  
  # Calculate standard deviations
  metrics_sd <- metrics_all[, lapply(.SD, function(x) sd(x, na.rm = TRUE)), 
                            by = .(iteration, p, q, r, sigma_zero, sigma_T1, sigma_T2, sigma_T3), 
                            .SDcols = c("n_gurus", "n_partners", 
                                        "transitivity_net", "centr_in", 
                                        "dens_net", "centr_out",
                                        "edges_net", "modularity_net")]
  
  # Rename sd columns
  setnames(metrics_sd, 
           old = c("n_gurus", "n_partners", 
                   "transitivity_net", "centr_in", 
                   "dens_net", "centr_out",
                   "edges_net", "modularity_net"),
           new = paste0("sd_", c("n_gurus", "n_partners", 
                                 "transitivity_net", "centr_in", 
                                 "dens_net", "centr_out",
                                 "edges_net", "modularity_net")))
  
  # Merge the two data tables
  metrics_summary <- merge(metrics_mean, metrics_sd, 
                           by = c("iteration", "p", "q", "r", "sigma_zero", "sigma_T1", "sigma_T2", "sigma_T3"))
  # Save results
  write_xlsx(metrics_summary, 
             "/home/users/NETSIM/AutoAdvice/metrics_summary_parallel_new_zeta.xlsx")
  
  cat("Simulations complete! Results saved to Excel files.\n")
  return(metrics_summary)
}

# Execute the main function
results <- main()
