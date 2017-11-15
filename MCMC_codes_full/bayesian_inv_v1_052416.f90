module bayesian_inv

use global_fric
use randn_gen
use random_gen
use matutils
use utilities
use util_mcmc

implicit none


contains
!--------------------------------------------------------------------------------------------------------------------------------------------------------
subroutine fixedmcmc(p,cov,uni_range,iterate,se,burnin,skip,psw,param,proposal,prior,dir_scaling,acc_req,covar,adapt_covar,threshold,thresh_ind)

! p: input parameter list, N-tuple vector
! cov: covariance matrix of the original proposal
! iterate: no. of members to be chosen for MCMC
! se: Square Error, scalar
! burnin: no. of initial members to be disregarded for statistics calculation
! skip: no. of members to be skipped after burnin for each chosen for statistics calculation
! psw: non-zero value for small world proposal, zero for non-small world
! param: Integer arguments for parameter scale definition: 1-'Linear', 2-'Logarithmic'
! proposal: Integer argument for form of proposal: 1-'Uniform', 2-'Gauss'
! prior: Form of the prior: 1-'Uniform', 2-'Exponential'
! dir_scaling: Yan-Bai, 2009 suggestion to tune covariance of proposal by scaling the size of jumps along PA directions based on a required acc. rate
! acc_req: Required accepance rate 
! covar: Learn covariance for proposal with a fixed number of samples after discarding burnin (assigned as covmat_len) from an initial chain, logical variable, true/false
! adapt_covar : Optional full Roberts_Rosenthal (2006) adaptive proposal, keeps updating the proposal-covariance even after burnin+covmat_len samples, true/false if present
! threshold: Constrain on prior, if ||p_prop||/||p_init|| > threshold => Prior = 0, optional argument
! thresh_ind : Vector of indices or positions of the parameters on which the prior constraint is to be applied e.g. if constrain on 1st and 3rd params thresh_ind = [1,3]

	real(8), parameter                          :: pi = 4.d0*atan(1.d0)
	real(8), parameter                          :: mix_factor = 0.9d0 ! If covar = .true., how often to use the updated covariance to propose
	integer, parameter                          :: covmat_len = 3*NP ! How many accepted proposals to be used for covmat
	real(8), dimension(:), intent(inout)        :: p
	real(8), dimension(:,:), intent(inout)      :: cov
	real(8), dimension(:), intent(in)           :: uni_range
	real(8), intent(in)                         :: psw, se, acc_req	
	logical, intent(in)                         :: covar
	logical, intent(in), optional               :: adapt_covar
	integer, intent(in)                         :: iterate, burnin, skip, proposal, param, prior, dir_scaling
	real(8), intent(in), optional               :: threshold
	integer, dimension(:), intent(in), optional :: thresh_ind
	real(8), dimension(:,:),allocatable         :: covar_inp, cov_nor, cov_temp, cov_safety, cov_init, std, scal_mat, I_NP
	real(8), dimension(:), allocatable          :: p_low, p_upp, p_prop, rand_p, p_loc, p_temp, p_check, p_max_ap, p_init, mu_p, lambda, mu_nor
	integer, dimension(:), allocatable          :: counter, count_vec
	integer                                     :: i, hist_num, bincnt, j, keep, k, lim, swp, nc
	integer                                     :: cov_flag, write_cov, cnt, reshaper, info, keep_sum, keep_choose
	real(8), dimension(2)                       :: rand_dum
	real(8), dimension(1)                       :: temp
	real(8)                                     :: likelihoodq, likelihoodp, rat_prior, swu, muloc(NP), dyad_vec(NP,1)
	real(8)                                     :: rmsep,chi2p,rmseq,chi2q,unirand, p_max, p_min, rebin, lebin, bin
	real(8)                                     :: alpha, keep_ind, fac, c_init, Vc_init, c_loc, Vc_loc
	real(8), dimension(:,:), allocatable        :: p_mat, p_keep, hist_mat, std_init, std_nor
	character(len=2)	                    :: strnp
	character(len=50)	                    :: fmat
	logical                                     :: constr, std_update
	
	unirand  = 0.d0
	swu      = 0.d0
	temp     = 1.d0
	swp      = 0
	cov_flag = 0 ! Counter for covariance update
	std_update = .false. ! To ensure that till the time std is updated (if at all), proposals use the initial std
	write_cov = 0 ! When covariance matrix is written for the first time
	allocate(covar_inp(iterate,NP))
	allocate(p_low(NP),p_upp(NP),p_prop(NP),rand_p(NP), p_temp(NP), p_max_ap(NP), p_check(NP), p_init(NP),mu_nor(NP))
	allocate(p_mat(NP+2,iterate),cov_nor(NP,NP),mu_p(NP),lambda(NP),cov_temp(NP,NP),cov_safety(NP,NP),cov_init(NP,NP),std(NP,NP))
	allocate(std_nor(NP,NP),std_init(NP,NP),scal_mat(NP,NP),I_NP(NP,NP))
	!-----------------------------------------------------------------------
	! Construct the NP rank identity matrix
	!-----------------------------------------------------------------------
	I_NP = 0.d0
	do nc = 1,NP
		I_NP(nc,nc) = 1.d0
	enddo
	!-----------------------------------------------------------------------
	!-----------------------------------------------------------------------
	cov_temp    = 0.d0
	cov_safety  = 0.d0	
	scal_mat    = I_NP
	cov_init    = cov
	call sqrtm(cov_init,std_init,NP,info)
	write(*,*) std_init
	std         = std_init
	hist_num    = 150
	rmsep 	    = func(p)
	if (present(threshold)) p_check = p ! Constraint on parameter values
	if ((present(threshold)) .and. (law=="K") .and. (choose_vc.eq.2)) p_check(NP) = 10**p(NP) ! Constraint on Vc? Store initial value of Vc
	if ((present(threshold)) .and. ((strain_rate_dep /= 0).or.((fric_law_form == 3).and.(BC_refvel >= 3)))) & 
	p_check(NP) = 10**p(NP) ! Constraint on Vstar? Store initial value of Vstar
	chi2p 	    = (rmsep**2.d0)/(se**2.d0)
	likelihoodp = -0.5d0*(NMcMc-1)*chi2p !! Log likelihood
	keep         = 1
	keep_ind     = 1.d0
	call init_random_seed()
	!-------------------------------------------------------
	! Covariances output file opened
	!-------------------------------------------------------
	if (covar) then
		open(unit=77,file=covmatout,STATUS='REPLACE')
		write(strnp,'(I2)') NP*NP
		fmat = trim('(')//trim(adjustl(strnp))//trim('e19.6)')
	endif
	!-------------------------------------------------------
		
	if (param==1) then
		if ((law.eq."K").and.(choose_vc.eq.2)) then
			p(NP) = 10.d0**(p(NP))
		elseif ((strain_rate_dep /= 0).or.((fric_law_form == 3).and.(BC_refvel >= 3))) then
			p(NP) = 10.d0**(p(NP))
		endif
		p_low       = 10.d0**(floor(log10(p))-uni_range)
		p_upp       = 10.d0**(ceiling(log10(p))+uni_range)
		p_init      = p
	elseif (param==2) then
		if ((law.eq."K").and.(choose_vc.eq.2)) then
			p(NP) = 10.d0**(p(NP))
		elseif ((strain_rate_dep /= 0).or.((fric_law_form == 3).and.(BC_refvel >= 3))) then
			p(NP) = 10.d0**(p(NP))
		endif
		p_low       = floor(log10(p)) - uni_range
		p_upp       = ceiling(log10(p)) + uni_range
		p           = log10(p)
		p_init      = p
	else
		write(*,*) 'Wrong value of param, values can be either 1 or 2'
		call ABORT
	endif
	lambda = 1.d0/p_init
	p_mat(:,1)  = (/p_init,keep_ind,rmsep/)
!------------------------------------------------------------------------------------------------!
	do i = 2,iterate ! Constructing the Markov chain
!------------------------------------------------------------------------------------------------!
		keep_ind    = 0.d0
2222		mu_p = p
		! Prior Selection and proposal of next state
		select case(proposal)
			case(1)			
				call random_number(rand_p)
				p_prop = (p_upp - p_low)*rand_p + p_low
			case(2)
				if ((proposal_std_constr==0).or.((proposal_std_constr==1).and.(.not.std_update))) then
					if (param==1) then
						p_prop = -1.d0
						do while (any(p_prop<=0.d0))
							call random_number(p_prop)
							call sqrtm(cov,std,NP,info)
							call nrmrnd(p_prop,mu_p,matmul(scal_mat,std),'Basic') ! The 'Polar' generation is less robust
						enddo
					else
						call random_number(p_prop)
						call sqrtm(cov,std,NP,info)
						call nrmrnd(p_prop,mu_p,matmul(scal_mat,std),'Basic') ! The 'Polar' generation is less robust
					endif
				 elseif (((proposal_std_constr==1).and.std_update).or.(proposal_std_constr==2)) then
		        	 	if (param==1) then
						p_prop = -1.d0
						do while (any(p_prop<=0.d0))
							call random_number(p_prop)
							call sqrtm(cov,std,NP,info)
							p_prop = mix_gauss2(p_prop,mu_p,matmul(scal_mat,std),mu_p,std_init,mix_factor) ! Proposal from mixture
						enddo
					else
							call random_number(p_prop)
							call sqrtm(cov,std,NP,info)
							p_prop = mix_gauss2(p_prop,mu_p,matmul(scal_mat,std),mu_p,std_init,mix_factor) ! Proposal from mixture
					endif
				endif
				
			case(3)
				mu_nor = 0.d0
 				std_nor = I_NP
 				call random_number(p_prop)
 				call sqrtm(cov_nor,std_nor,NP,info)
				call nrmrnd(p_prop,mu_nor,std_nor,'Basic')
				p_prop = corr_uni_rand(p_prop,matmul(scal_mat,std))
				p_prop = (p_upp - p_low)*p_prop + p_low	
			case default
				call ABORT
		end select
		
!----------------------------------------------------------------------------------------------------------------!
! Impose prior constraints: Option 2 ==> Better option if the forward model slows down too much for wild samples
!----------------------------------------------------------------------------------------------------------------!
	if (present(threshold)) then
		p_temp = p_check  !! Just put all values equal to p_check to start with
		if (present(thresh_ind)) then
			if (param==1) p_temp(thresh_ind) = p_prop(thresh_ind)
			if (param==2) p_temp(thresh_ind) = 10**p_prop(thresh_ind)
		elseif (.not.present(thresh_ind)) then
			if (param==1) p_temp = p_prop
			if (param==2) p_temp = 10**p_prop
		endif		
		if (any(p_temp/p_check>=threshold).or.any(p_temp/p_check<1.d-2)) then
			write(*,*) 'Prior constraint violated!! Resampling'
			goto 2222
		endif
	endif	        	 		
!---------------------------------------------------------------------------------------------------------------------!
!---------------------------------------------------------------------------------------------------------------------!		        	
		
		! Small World Proposal to facilitate mixing
		if (psw>0.d0.and.proposal==2) then
			call random_number(swu)
			if (swu<=psw) then
				call random_number(rand_p)
				p_prop = (p_upp - p_low)*rand_p + p_low
				swp = swp + 1
			endif			
		endif
			
		! Calculating rmse for proposed parameters, constrained if the value is too different from initial pt. by threshold
		if (param==1) then
			rmseq 	    = func(p_prop)
		else
			if ((law=='K').and.(choose_vc==2)) then
				rmseq       = func((/10.d0**p_prop(1:NP-1),p_prop(NP)/))
			elseif ((strain_rate_dep /= 0).or.((fric_law_form == 3).and.(BC_refvel >= 3))) then
				rmseq      = func((/10.d0**p_prop(1:NP-1),p_prop(NP)/))
			else
				rmseq       = func(10.d0**(p_prop))
			endif			
		endif
		if(isnan(rmseq)) rmseq = huge(rmseq)
		chi2q 	    = (rmseq**2.d0)/(se**2.d0)
		likelihoodq = -0.5d0*(NMcMc-1)*chi2q !! Log likelihood
		if (prior==1) then
			alpha       = min(1.d0,exp(likelihoodq-likelihoodp))
		elseif (prior==2) then
			rat_prior = exp(sum(-lambda*(p_prop-p)))
			alpha     = min(1.d0,exp(likelihoodq-likelihoodp)*rat_prior)
		endif
		
		fac = 1.d0 ! Setting default value for fac = 1.d0, alpha is not changed as a result

!----------------------------------------------------------------------------------------------------------------!
! Impose prior constraints: Option 1 ==> Not very good if the forward model slows down too much for wild samples
!----------------------------------------------------------------------------------------------------------------!		
		! Impose constraint on c, if beyond limits, set alpha -ve
		!if ((present(threshold)) .and. (law=="N")) then
		!	if (param==1) c_loc = p(NP)
		!	if (param==2) c_loc = 10**p(NP)
		!	constr     = (c_loc/c_init.ge.threshold).or.(c_loc/c_init.le.1.d-2)! Evaluating whether Nagata constraint is satisfied
		!elseif ((present(threshold)) .and. (law=="K") .and. (choose_vc==2)) then
		!	Vc_loc = 10**p(NP)
		!	constr = (Vc_loc/Vc_init.ge.threshold).or.(Vc_loc/Vc_init.le.1.d-3)! Evaluating whether Kato constraint is satisfied
		!endif
		!if (constr) fac = -1.d0 ! Constraint imposed for linear search
		
		!alpha = alpha*fac ! If c is beyond constraints, alpha is -ve hence never chosen
!---------------------------------------------------------------------------------------------------------------------!
!---------------------------------------------------------------------------------------------------------------------!
		
		call random_number(unirand)
		if (unirand .le. alpha) then
			keep 	    = keep + 1
			keep_ind    = 1.d0
			p   	    = p_prop
			likelihoodp = likelihoodq     
		endif
		p_mat(:,i) = (/p_prop,keep_ind,rmseq/)
		
!----------------------------------------------------------------------------------------------------------------!
! Another policy of choosing a different sample once edge is reached. Essentially restarts the chain.
!----------------------------------------------------------------------------------------------------------------!
		
!		keep_sum   = 0		
!		if ((fac==-1.d0).and.(i>1)) then ! If you encounter the edge of the prior support restart the chain from one of the accepted samples
!			keep_choose = int(keep*unirand) ! Choose the starting point randomly from the accepted samples
!			nc = 1
!			do while(keep_sum<keep_choose)
!				if (p_mat(NP+1,nc)==1.d0) keep_sum = keep_sum + 1
!				nc = nc + 1
!			enddo
!			p  = p_mat(1:NP,nc-1) ! Actually choosing the starting point
!			write(*,*) 'Edge of prior support, sampling from accepted samples!!'
!		endif
			
			
		
		write(*,'(A,i10,A,f7.3,A,i1,A,f16.10,A,i6)') 'Iteration: ', i, ', Acc. Rate: ', real(keep)/real(i)*100.d0, & 
		', Acc. Ind.: ', int(keep_ind), ', alpha: ', alpha ,', Small World: ', swp
!-------------------------------------------------------------------------------------------!
	! Adaptive proposal part, update std if required for next iteration of the chain                       
!-------------------------------------------------------------------------------------------!		
		! Covariance matrix for proposal formulation
		if (covar.and.(keep_ind==1.d0)) then
			cov_flag   = cov_flag + 1
			covar_inp(cov_flag,:) = p
			if ((cov_flag-burnin==covmat_len).or.(i==iterate)) then				
				cov_temp = emp_covmat(covar_inp(burnin+1:cov_flag,:)) ! Covariance if adaptive proposal reqd.
				muloc = sum(covar_inp(burnin+1:cov_flag,:))/(cov_flag-burnin)
				cov = (2.38**2.d0)**2*cov_temp/dble(NP) ! Roberts and Rosenthal scaling
				cov_safety = cov ! For use with scale_mat when running adaptation of variance is not used
				write(77,fmat) reshape(cov,[NP*NP, 1])
				write(*,*) 'Covar updated'
				std_update = .true.
				if (.not.present(adapt_covar).or.(.not.adapt_covar)) EXIT
		      !-------------------------------------------------------------------------------------------!
		      ! Update std at every iteration beyond covmat_len: Not very robust                          !
		      !-------------------------------------------------------------------------------------------!				
			elseif ((cov_flag-burnin>covmat_len).and.(present(adapt_covar).and.(adapt_covar))) then ! update std on the fly
				rewind(77)
				muloc = (dble(cov_flag-burnin-1)*muloc + p)/(cov_flag-burnin)
				cov_temp = emp_covmat(covar_inp(burnin+1:cov_flag,:)) !update_covmat(reshape(p,[NP,1]),reshape(muloc,[NP,1]),cov_flag-
				! burnin,cov_temp)
				cov = (2.38**2.d0)**2*cov_temp/dble(NP)
				write(77,fmat) reshape(cov,[NP*NP, 1])
		     !----------------------------------------------------------------------------------------------------!
		     !                                                                                                    !
		     !----------------------------------------------------------------------------------------------------!
			endif
		endif
!------------------------------------------------------------------------------------------------------------------!
! Adaptive directional scaling: Based on Yan-Bai's ADMG suggestion, Yan-Bai, 2009 and PhD thesis
!------------------------------------------------------------------------------------------------------------------!	
		if ((dir_scaling==1).and.((keep>burnin).or.(i>1000))) then  ! Do not make std scaling continually smaller at the edges 
		! of prior support
			do nc = 1,NP
				scal_mat(nc,nc) = exp(2.d0*NP*(real(keep)/real(i) - acc_req)) + 0.01 ! Define a variable scaling dependent on acc. rate history
			enddo
		else
			scal_mat = I_NP
		endif

		
	enddo ! End of MCMC chain
!------------------------------------------------------------------------------------------------------------------!
	if (covar) close(77)
	! Writing out to file results
	if (param==2) p_mat(1:NP,:) = 10.d0**(p_mat(1:NP,:))
	open(unit=75,file=statfile,STATUS='REPLACE') 
	write(strnp,'(I1)') NP+2
	fmat = trim('(')//trim(strnp)//trim('e19.6)')
	write(75,fmat) p_mat
	close(75)
	! Statistics Calculation
	allocate(p_keep(NP+1,keep))
	allocate(hist_mat(2*NP,hist_num))
	allocate(p_loc(keep),counter(keep),count_vec(hist_num))
	p_loc = 0.d0
	j = 1
	do i = 1,iterate
		if (p_mat(NP+1,i)==1.d0) then
			p_keep(:,j) = p_mat(1:NP,i)
			j = j + 1
		endif
	enddo
	lim          = ((keep-(burnin+1))/skip)+1
	do i = 1,NP
		p_loc(1:lim) = p_keep(i,burnin+1:keep:skip)
		p_max        = maxval(p_loc(1:lim))
		p_min        = minval(p_loc(1:lim))
		bin   = (p_max - p_min)/hist_num
		do j = 1,hist_num
			counter = 0
			bincnt = 0
			lebin  = p_min + (j-1)*bin
			rebin  = p_min + j*bin
			where ((p_loc > lebin) .and. (p_loc <= rebin)) counter = 1
			bincnt = sum(counter)
			hist_mat(2*i-1,j) = 0.5d0*(lebin + rebin)
			hist_mat(2*i,j) = bincnt
		enddo
		count_vec   = hist_mat(2*i,:)
		p_max_ap(i) = hist_mat(2*i-1,maxloc(count_vec,1))
	enddo	
	p = p_max_ap
	!call combined_model(p_keep(1:NP,burnin+1:keep:skip))
	! Writing out to file histogram results
	open(unit=76,file=histfile,STATUS='REPLACE')	
	write(strnp,'(I2)') 2*NP
	fmat = trim('(')//trim(strnp)//trim('e19.6)')  
	write(76,fmat) hist_mat
	close(76)
	deallocate(p_low, p_upp, p_prop, rand_p, p_loc, p_max_ap,counter,count_vec)
	deallocate(p_mat)
	deallocate(p_keep)
	deallocate(hist_mat)
	
	end subroutine fixedmcmc
!-----------------------------------------------------------------------------------------------

 end module bayesian_inv
