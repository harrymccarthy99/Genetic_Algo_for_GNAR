# Loading Packages

library(dplyr)
library(ggplot2)
library(reshape2)
library(lubridate)
library(tidyverse)
library(geodist)
library(countrycode)
library(geodist)
library(dplyr)
library(maps)
library(ggmap)
library(igraph)
library(geosphere)
library(GNAR)
library(MCMCpack)
library(glue)
library(GA)
library(here)

# Read function definitions
source(here("R", "ga_functions.R"))

# Load data and networks
unemployment_ts <- read.csv(here("data", "unemployment_ts.csv"))
unemployment_vts <- ts(subset(unemployment_ts, select = -Month_Year))
load(here('data', 'MST_Distance_Networks.RData'))

# Plotting and saving networks
plot_from_matrix(MST_37, 'Distance Based Minimum Spanning Tree Network', 'MST_37_network.png')
plot_from_matrix(distance_matrix_600, 'Distance Inferred Network (600km)', 'distance_600_network.png')

# Finding the unique additional  edges in distance_matrix_600 that are not in MST_37
MST_37_edges <- adj_matrix_to_edge_list(MST_37)
distance_matrix_600_edges <- adj_matrix_to_edge_list(distance_matrix_600)
unique_additional_edges <- length(find_unique_edges(MST_37_edges, distance_matrix_600_edges))

# Merge and plot full network
merged_adj_matrix <- merge_adjacency_matrices_igraph(distance_matrix_600, MST_37)
plot_from_matrix(merged_adj_matrix, 'MST & Proximity Inferred Network', 'MST_proximity.png')

 
# GA Parameters
population_size <- 100
number_of_bits <- 630
crossover_prob <- 0.8
mutation_prob <- 0.01
max_generations <- 100


# Binary
