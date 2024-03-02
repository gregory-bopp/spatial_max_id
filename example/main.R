library(stablemix)

# Function for setting up simulation
source('./simulation_fuction.R')
# Simulation settings
simset <- list(
            alpha = 0.1,               # alpha-stable index
            theta = 0.1,               # exponential tilting parameter
            gvar = 25,                 # GP Var
            gscl = 0.5,                # GP Scale
            n = 30,                    # number of replicates (years)
            nloc = 100,                # number of spatial locations
            L = 15                     # number of random effects per year
            )

# Run example simulation (~10 minutes)
mh <- example_sim(1000, simset)

# Save posterior samples
dir.create('./data')
save(mh, file = './data/mcmc.Rdata')