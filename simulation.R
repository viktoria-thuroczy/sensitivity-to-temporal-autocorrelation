rm(list = ls())
# load packages
library(popdemo)

library(stats)
library(dplyr)

setwd("C:/Users/vikto/OneDrive/Documents/stochastic simulation paper/code/data")

load("comadre.RData")
load("compadre.RData")
LHT <- read.csv("LifeHistoryTraits.csv")

### Prep data for simulation ----

#comadre
comadre_df <- as.data.frame(comadre)
comadre_df$lambda <- NA

for(i in 1:nrow(comadre_df)) {
  matA  <- comadre_df$mat[[i]]@matA
  
  l <- eigs(matA, what = "lambda")
  
  comadre_df$lambda[i] <- l
}

comadre_df$matA <- lapply(comadre_df$mat, function(x) x@matA)

comadre_df <- comadre_df %>%
  select(SpeciesAccepted, MatrixPopulation, MatrixDimension, mat, matA, lambda, MatrixStartYear, MatrixEndYear) %>%
  group_by(SpeciesAccepted) %>%
  mutate(Pops = n_distinct(MatrixPopulation)) %>%
  mutate(MPMs = n()) %>%
  group_by(MatrixPopulation, SpeciesAccepted) %>%
  mutate(MPMs_Pop = n()) %>%
  ungroup(MatrixPopulation)


#compadre 
compadre_df <- as.data.frame(compadre)
compadre_df$lambda <- NA

for(i in 1:nrow(compadre_df)){
  matA <- compadre_df$mat[[i]]@matA
  
  l <- eigs(matA, what = "lambda")
  
  compadre_df$lambda[i] <- l
}

compadre_df$matA <- lapply(compadre_df$mat, function(x) x@matA)

compadre_df <- compadre_df %>%
  select(SpeciesAccepted, MatrixPopulation, MatrixDimension, mat, matA, lambda, MatrixStartYear, MatrixEndYear) %>%
  group_by(SpeciesAccepted) %>%
  mutate(Pops = n_distinct(MatrixPopulation)) %>% 
  mutate(MPMs = n()) %>%
  group_by(MatrixPopulation, SpeciesAccepted) %>% 
  mutate(MPMs_Pop = n()) %>%
  ungroup(MatrixPopulation)


#combine two dataframes
sim_data <- rbind(comadre_df, compadre_df)

#filter out species which have MPMs with multiple dimensions (if one dimension is found 5 or more times
#keep those MPMs)
sim_data <- sim_data %>%
  group_by(SpeciesAccepted, MatrixDimension) %>%
  summarize(num_rows = n()) %>%
  ungroup() %>%
  filter(num_rows > 4) %>%
  group_by(SpeciesAccepted) %>%
  slice(which.max(num_rows)) %>%
  ungroup() %>%
  left_join(sim_data, by = c("SpeciesAccepted", "MatrixDimension"))

#save data
save(sim_data, file = "sim_data.RData")


### Load functions ----
build_stochastic_matrix_continuous <- function(dimensions, autocorrelation, lambda){
  
  mat <- matrix(0, dimensions, dimensions)
  
  lambda_range <- max(lambda) - min(lambda)
  
  
  for(column in 1:dimensions){
    
    for(row in 1:dimensions){
      
      if(column == row){
        
        mat[row, column] <- ((1/dimensions) + autocorrelation * (1/dimensions)) / (2/dimensions)
        
      } else if(column != row){
        
        m <- (-2*autocorrelation) / dimensions
        x <- (abs(lambda[column] - lambda[row])) / lambda_range
        b <- ((1/dimensions) + autocorrelation * (1/dimensions))
        
        mat[row, column] <- (m * x + b) / (2/dimensions)
        
        
      }
    }
  }
  
  
  # Now make the columns sum to 1
  
  for(column in 1:dimensions){
    
    mat[, column] <- mat[, column] / sum(mat[, column])
    
  }
  
  
  
  if(sum(mat < 0) > 0){
    
    cat("Negative values were created :(")
    
    return()
    
  }
  
  return(mat)
  
}

project_states <- function(stochastic_matrix, initial_state = 1, timesteps = 100){
  
  states <- c(1:ncol(stochastic_matrix))
  
  state_vector <- c()
  
  state_vector[1] <- initial_state
  
  for(timestep in 1:(timesteps-1)){
    
    probabilities <- stochastic_matrix[,state_vector[timestep]]
    
    state_vector[timestep+1] <- sample(states, size = 1, prob = probabilities)
    
    
  }
  
  
  
  return(state_vector)
  
}




### Simulation ----
iterations <- 1000
discard <- 100
autos <- c(-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

# Extract unique species names from sim_data
species_names <- unique(sim_data$SpeciesAccepted)

# Create an empty dataframe to store results
species_lambda_df <- data.frame(SpeciesAccepted = character(), Autocorrelation = double(), Stoch_GR = double())

# Iterate over each species
for (species in species_names) {
  # Load species data
  species_data <- sim_data %>%
    filter(SpeciesAccepted == species) %>%
    arrange(lambda)
  
  lambda_values <- species_data$lambda
  matA_values <- species_data$matA
  dim_value <- length(lambda_values)
  
  # Iterate over each auto value
  for (auto in autos) {
    # Make dataframe where stochastic pop growth rate will be stored
    sgr_df <- data.frame(sgr = double())
    
    # Loop to run simulation multiple times
    for (i in 1:50) {
      continuous_stochastic_matrix <- build_stochastic_matrix_continuous(dim_value, auto, lambda_values)
      continuous_projection <- project_states(continuous_stochastic_matrix, initial_state = sample(dim_value, 1), timesteps = iterations)
      
      vector <- runif(dim(matA_values[[1]])[1])
      vector <- vector / sum(vector)
      
      pr <- project(A = matA_values, vector = vector, time = iterations, Aseq = continuous_projection,
                    standard.A = FALSE, standard.vec = FALSE, return.vec = TRUE, PREcheck = TRUE)
      
      gr <- pr[(discard:iterations) + 1] / pr[discard:iterations]
      gr <- gr[is.finite(gr) & gr != 0]
      sgr <- mean(log(gr))
      
      sgr_df <- rbind(sgr_df, data.frame(sgr = sgr))
    }
    
    # Calculate mean sgr value
    mean_lambda <- mean(sgr_df$sgr)
    
    # Append results to species_lambda_df
    species_lambda_df <- rbind(species_lambda_df, data.frame(SpeciesAccepted = species, Autocorrelation = auto, Stoch_GR = mean_lambda))
  }
}

#make adjustments so it fits phylogenetic analysis: 
species_lambda_df$SpeciesAccepted[species_lambda_df$SpeciesAccepted=="Chamaecrista lineata var. keyensis"]="Chamaecrista lineata"
species_lambda_df$SpeciesAccepted[species_lambda_df$SpeciesAccepted=="Echinospartum ibericum subsp. algibicum"]="Echinospartum ibericum"
species_lambda_df$SpeciesAccepted[species_lambda_df$SpeciesAccepted=="Eriogonum longifolium var. gnaphalifolium"]="Eriogonum longifolium"
species_lambda_df$SpeciesAccepted[species_lambda_df$SpeciesAccepted=="Oenothera coloradensis subsp. coloradensis"]="Oenothera coloradensis"
species_lambda_df$SpeciesAccepted[species_lambda_df$SpeciesAccepted=="Petrocoptis pyrenaica subsp. pseudoviscosa"]="Petrocoptis pyrenaica"

merged_data <- merge(species_lambda_df, LHT, by = "SpeciesAccepted", all.x =  TRUE)
write.csv(merged_data, "LHT_SPGR.csv", row.names = FALSE)


### Calculate Sensitivity ----
species_names <- unique(merged_data$SpeciesAccepted)

# Create an empty dataframe to store results
species_sensitivity <- data.frame(SpeciesAccepted = character(), Sensitivity = double())

# Iterate over each species
for (species in species_names) {
  # Load species data
  species_data <- merged_data %>%
    filter(SpeciesAccepted == species) 
  
  stoch_gr <- species_data$Stoch_GR
  auto <- species_data$Autocorrelation
  
  mod <- lm(stoch_gr ~ auto)
  sen <- abs(coef(mod)[[2]])
  
  species_sensitivity <- rbind(species_sensitivity, data.frame(SpeciesAccepted = species, Sensitivity = sen))
}

sens_data <- merge(species_sensitivity, LHT, by = "SpeciesAccepted", all.x =  TRUE)
write.table(sens_data, "LHT_sensitivity.csv", row.names = FALSE, sep = ";")



