# Models that we want (to be continued)
## Brownian motion
- single trait, single sample (HQ)
- multiple traits, single sample (HQ)
- multiple traits, multiple samples (HQ)
- multiple correlated traits (HQ)

## Ornstein-Uhlenbeck
- single trait, single sample (PL)
- multiple traits, single sample (PL)
- multiple traits, multiple samples (PL)
- multiple correlated traits (PL)

## State dependence in models
- Don't include it for now

# Calculating the likelihood
## Likelhood type
- REML
- Standard <-- we use this

## Likelihood algorithm
- Generalized least squares
- Pruning (FitzJohn 2012, see supplement page 8 https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fj.2041-210X.2012.00234.x&file=mee3234_sm_diversitree-suppl.pdf )

We implement both for validation, use pruning for the real analysis.

# Data structure
- Treat the tree as known, or already provided
- Use the ape format for representing the tree
- Consider using the `tidytree` format for the data (https://yulab-smu.top/treedata-book/chapter2.html )
- Other options: data.frame, NEXUS file

# Milestones
1. Implement the likelihood function for all of the models
2. Validate that the likelihood functions are correct
3. Hill climber algorithm to find the maximum-likelihood parameter estimates
4. Simple Bayesian inference (or use another package for this, for example TESS)
5. Write a vignette with examples on how to use the package

