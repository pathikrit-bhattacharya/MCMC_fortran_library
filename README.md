# MCMC_fortran_library
Suite of adaptive MCMC routines in Fortran:
The MCMC algorithm used here was based on the Metropolis-Hastings (MH)
algorithm [Metropolis et al., 1953; Hastings, 1970] and makes use of the
strategy of sampling the posterior using Bayes' theorem. The code is
designed for fast convergence on correlated parameter spaces by
adaptively learning the parameter covariance of the sampled region. The
proposal jumps are then planned based on this covariance structure with
step size prefixed according to some linear scaling of the covariance
matrix. The scaling used is exact for large dimensional spaces as long
as the moves of the proposals are diffusive. The document 'Description
of using the MCMC code for a rate-state inverse problem.pdf' (which is really
the Supplementary Material from Bhattacharya et al. [2015]) explains
many of the working details of the code. The code also has options for
setting a-priori constraints on any subset of the parameters. There are
also options for running non-adaptive chains with no covariance
structures and those with pre-defined covariance structures.
