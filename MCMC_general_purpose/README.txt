The general structure of the code is as fllows:

!-------------------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------------------!
A) The main file is bayesian_inv_xxxxxxx...f90
This file contains the main subroutine fixedmcmc
with call structure:

fixedmcmc(p,func,cov,uni_range,iterate,se,burnin,skip,psw,param,proposal,prior,dir_scaling,acc_req,covar,adapt_covar, &
threshold,thresh_ind,thresh_vec)

1) p: input parameter list, N-tuple vector
2) func: Handle for the misfit calculation function, whatever the name of the cost function calculator is.

3) cov: Initial value of the covariance matrix between parameters, to be used for the centr-shifting N-d Gaussian proposal.

4) iterate: Total number of proposals to try, the user can terminate the run at any time; the output till that time is recorded.

5) se: scalar norm of data error. Assumes that all data points are uncorrelated and equally error prone.

6) burnin: no. of initial members to be disregarded for statistics calculation, around 2 to 10 times N is usually fine. Important choice for adpative runs.

7) skip: no. of members to be skipped after burnin for statistics calculation, again around 2 to 10 N usually fine.
psw: non-zero value for small world proposal of Guan et al. [2006], zero for non-small world.

8) param: Integer arguments for parameter scaling definition: 1-'Linear' (the parameter space is searched in linear increments), 2-'Logarithmic' 
(the parameter space is searched in logarithmic increments).

9) proposal: Integer argument for form of proposal: 1-'Uniform' (proposal is correlated uniform), 2-'Gauss' (proposal is correlated Gaussian).

10) prior: Form of the prior: 1-'Uniform', 2-'Exponential'
dir_scaling: Yan-Bai, 2009 suggestion to tune covariance of proposal by scaling the size of jumps along principal axes directions of the covariance
learnt from the sampled region and based on a required acceptance rate.

11) acc_req: Required accepance rate for convergence, has a theoretical limit of around 25% for large N (Roberts and Rosenthal, 2006).

12) covar: Learn covariance for proposal with a fixed number of samples after discarding burnin (assigned as covmat_len) from an initial chain. Logical variable, true or false

13) adapt_covar: Optional full Roberts_Rosenthal (2006) adaptive proposal, keeps updating the proposal-covariance even after burnin+covmat_len samples, true or false if present.

14) threshold: Constrain on prior, if ||p_prop||/||p_init|| > threshold => Prior = 0, optional argument; where p_prop is the current proposal, p_init is initial. Ensures that
the chain does not sample too far from initial point.

15) thresh_ind: Vector of indices or positions of the parameters on which the prior constraint is to be applied e.g. if constraints on 1st and 3rd params thresh_ind = [1,3];
used to impose absolute constraints on parameter values.

16) thresh_vec: matrix of thresholds, same number of cols as thresh_ind. 2 rows: 1st row-Upper Bound, 2nd row-Lower Bound.

!-------------------------------------------------------------------------------------------------------------------------------------!
!-------------------------------------------------------------------------------------------------------------------------------------!

B) matutils.f90 is a module containing various matrix processing utilities using LAPACK. If the user has MKL and threaded LAPACK, most of these routines can be parallelized. A list with brief descriptions follows:

1) function covmat(x, w): calculates covariance matrix of x,
 with optional (frequency) weights in w,
 returns covariance matrix in covmat.
2) function emp_covmat(x): calculates the empirical covariance matrix given n-samples of a m-dimensional vector as an nxm matrix
. Formula from Yan-Bai, 2009 (An Adaptive Directional Metropolis-within-Gibbs algorithm), corrected Yan_Bai's
 formula to accommodate the full covariance matrix.
3) function update_covmat(x,x_mean,n,std): Updates the existing covariance matrix iteratively by adding the most recently sampled observation
 to implement Roberts and Rosenthals adaptive approach.
4) subroutine covtor(cmat,R,info): Calculates upper diagonal Cholesky transformation
 of symmetric matrix 'cmat'
.
5) function covtor_svd(cmat,info): Singular Value Decomposition of dense matrix, need to supply matrix in column major order. Of course, covariance matrix is symmetric. Need to compute matrix square root.
6) function invertmat_sym(x,info0): Invert symmetric matrix using Cholesky decomposition.
7) function det_sym(A): Determinant of a symmetric matrix A.
8) function mahalanobis(x,mu,cmat): Mahalanobis distance of rows of x from mu given a covariance matrix cmat.
9) subroutine copylower(x): Copy the lower diagonal from upper diagonal for symmetric x.
10) subroutine sqrtm(cmat,R,n,INFO): Square Root of a positive semidefinite square matrix by diagonalization.
11) subroutine linsolve(Amat,solvec, info): Solve a dense, square, linear system using dgesv from LAPACK.