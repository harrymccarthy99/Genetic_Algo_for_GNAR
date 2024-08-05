# Genetic_Algo_for_GNAR
Code to implement a Genetic Algorithm to search for optimal network structures to maximise Generalised Network Autoregressive model predictive accuracy

GNARs are used to model nodal time series, providing a parsimonious model that improves forecasting accuracy in many scenarios versus competing models eg. VAR

Crucially, these models rely on an inputted network structure indicating the relevant links between nodes. This Genetic Algorithm treats these networks as chromosomes and searches for the network that maximises within-sample predictive accuracy through evolution of the population. 
