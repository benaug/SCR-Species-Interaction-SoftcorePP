# SCR-Species-Interaction-SoftcorePP
A 2 species SCR model where species 2 avoids species 1 using a soft-core point process approach

This repository contains a basic implementation of a soft-core point process for 2-species interaction models. 
These models use an "interaction function", "h", for species 1's effect on the placement of species 2's activity centers.
The function used is:

dists <- e2dist(s1,dSS) #distance between species 1 activity centers and all state space locations 

kern <- exp(-dists^2/(2*omega^2)) 

log.h <- colSums(log(1-kern)) 

h <- exp(log.h)

This is what is used in the base model "Softcore". "Softcore baseline" includes a second parameter in the interaction function.

kern <- h0*exp(-dists^2/(2*omega^2))

This baseline interaction parameter allows less extreme repulsion. Note, that without h0, whenever a species 1 activity center is located
exactly at a state space location, h for that cell is 0 and species 2 cannot have an activity center in this cell. This happens very frequently
when using a discrete state space as used in this repo. With a continuous state space, as species 1 activity center approaches the state space
centroid, we get the same effect. But now, the arbitrary movement of species 1 activity centers inside a state space cell, perhaps assumed to be
uniformly distributed, changes the magnitude of repulsion of species 2. Also, changing the state space cell sizes also arbitrarily controls how often
complete exclusion from a cell happens (this seems less than ideal behavior, so I prefer a discrete state space.
for this model). 

Including h0 allows for less extreme repulsion, allowing species 2 to potentially have activity centers in cells housing
a species 1 activity center. But, two interaction function parameters are less identifiable than one.

"Softcore baseline shareSig" attempts to remedy the reduced identifiability from introducing the h0 parameter by assuming the spatial scale of repulsion
is the same as the detection spatial scale of species 1. This might be interpreted as species 1 excluding species 2 at the home range scale, but this is not
the exact interpretation because a species 2 activity center placed at the edge of a species 1 home range will still have implied space use inside of the species 1
home range, governed by species 2's detection sigma.

kern <- h0*exp(-dists^2/(2*sigma1^2))

I do not expect these models to work well with most data sets that can be collected in practice. The interaction function parameters and density parameters
are often very weakly identifiable. The geometry of the landscape and trap configurations can have a large impact on identifiability. Also, identifiability is better
when densities are higher. The data simulators produce useful plots for determining how well a scenario might work, including a plot of the interaction function.

The testscripts are set up for a scenario where both species respond strongly to the same habitat covariate, so they prefer the same parts of the landscape. But species 2 avoids species 1.
I was able to get nearly unbiased estimates with roughly nominal 95% coverage from the "Softcore" model with 225 traps, and average of 30 species 1 individuals and 150 species 2 individuals. A similar scenario for
"Softcore baseline shareSig" yielded nearly unbiased estimates, but coverage was low for some parameters (down to about 80%). Both the testscripts and simulation investigation assume detection probabilities are higher than is usually realistic in practice.

Generally, I think the "hard exclusion" assumption of h0=1 goes a long way to providing identifiability of process model parameters. But maybe there are smarter scenarios.


Notes on MCMC: I use a "pragmatic" reversible jump MCMC approach for both species where the process model latent variables remain in the model, but the individuals in the observation model are turned on and off. This is less efficient that regular RJMCMC, but can be set up relatively easily in nimble and allows an N ~ Poisson(lambda) assumption.
To speed up computation, I use the approach of Herliansyah et al. (2024) in the custom N/z and activity center updates.

Also, these models are still very slow to run. The species 1 "data augmentation" (still have to give nimble an upper limit) has a large effect on run time due to species 2 process model parameters depending on species 1 latent variables.
Trimming of the process and observation model computations could futher improve efficiency.


https://link.springer.com/article/10.1007/s13253-023-00598-3


Two references for softcore point processes in spatial capture-recapture are:

Diana, Alex, et al. "A vector of point processes for modeling interactions between and within species using capture‐recapture data." Environmetrics 33.8 (2022): e2781.
https://onlinelibrary.wiley.com/doi/full/10.1002/env.2781

Gaya, Heather E., and Richard B. Chandler. "Individual‐level biotic interactions and species distribution models." Journal of Biogeography 51.11 (2024): 2071-2083.
https://onlinelibrary.wiley.com/doi/full/10.1111/jbi.14972

