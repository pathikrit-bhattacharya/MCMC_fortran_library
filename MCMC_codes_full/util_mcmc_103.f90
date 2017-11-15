module util_mcmc

use global_fric
use nrtype
use nrutil
use nr_bsstep
use utilities
use matutils
use parameters
use randn_gen
use random_gen

implicit none


contains 

subroutine combined_model(p_acc)

	  real(8), dimension(:,:), intent(in) :: p_acc
	  real(8), dimension(:), allocatable  :: x_f, x_lp, t, v_lp, tau, slipv, v_f
	  real(8), dimension(:), allocatable  :: state, sliprate, rand, state_av, fmu, sliprate_av, slipv_av, tau_av, p_waste
	  real(8)                             :: t1
	  real(8)                             :: aa,bb,dc,Vc,stiffk,vo,c
	  character                           :: ln
	  character, dimension(29)            :: fm
	  integer                             :: i, Nav, sz, ind
	  integer, dimension(:), allocatable  :: collect
	  character(len=100)                  :: mcmcavout,fitdata
	  
	  allocate(p_waste(NP))
	  do ind = 1,2
	  	select case(ind)
			case(1)
				mcmcavout = mcmcavout1
				fitdata   = fitdata1
				Nrange    = Nr1
			case(2)
				mcmcavout = mcmcavout2
				fitdata   = fitdata2
				Nrange    = Nr2
		end select
			  
	  	allocate(x_f(Nrange),tau(Nrange),x_lp(Nrange),t(Nrange),v_lp(Nrange))
	  	allocate(slipv(Nrange),sliprate(Nrange),state(Nrange),fmu(Nrange),v_f(Nrange))
	  	allocate(slipv_av(Nrange),sliprate_av(Nrange),state_av(Nrange),tau_av(Nrange))
	  	Nav = 100
	  	allocate(rand(Nav))
	  	allocate(collect(Nav))
	  	call init_random_seed()
		call random_number(rand)
	  	collect=1+int(dble(size(p_acc,2))*rand)
		open(unit=112,file=mcmcavout,STATUS='REPLACE')
		open(unit=113,file=fitdata,STATUS='OLD') 
		do i = 1,Nrange
			read(113,'(7e20.9)') x_lp(i), t(i), fmu(i), v_lp(i), x_f(i), v_f(i), t1
		enddo	
		close(113)
		
		slipv_av     = 0.d0
		sliprate_av = 0.d0
		state_av    = 0.d0
		tau_av      = 0.d0
		
		do i = 1, Nav
			call getparams(p_acc(:,collect(i)),p_waste,ind,aa,bb,dc,Vc,stiffk,c)
			call fmodel(p_acc(:,collect(i)), tau, slipv, sliprate, state, x_f, x_lp, t, v_lp, ind)
			slipv_av     = slipv_av + slipv
			sliprate_av = sliprate_av + sliprate
			state_av    = state_av + state
			tau_av      = tau_av + tau
		enddo
		slipv_av     = slipv_av/Nav
		sliprate_av = sliprate_av/Nav
		state_av    = state_av/Nav
		tau_av      = tau_av/Nav
		do i = 1,Nrange
			write(112,'(2f14.5,e14.5,3f14.5,e14.5,i6,3e14.5)')t(i),x_f(i), &
		  	sliprate_av(i),slipv_av(i),tau_av(i),fmu(i),sliprate(i)*state_av(i)/dc,i,v_f(i), t(i), x_lp(i)
		enddo
	deallocate(x_f,tau,x_lp,t,v_lp)
	deallocate(slipv,sliprate,state,fmu,v_f)
	deallocate(slipv_av,sliprate_av,state_av,tau_av)
	deallocate(rand)
	deallocate(collect)
	enddo
end subroutine combined_model

!------------------------------------------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------------------------------------------!

subroutine fmodel(rnstry, tau, slipv, sliprate, state, x_f, x_lp, t, v_lp, runind)

! This version has no sn! because output is in terms of f.
! stiffk/sn has been changed to stiffk
! The velocity used to update is the change in l.p. disp. from the previous line 
!   divided by the time since the previous line.
!
	implicit none
      	  real(DP), dimension(:), intent(in) :: rnstry
	  real(8), dimension(:), intent(in)  :: x_f, x_lp, t, v_lp
	  integer, intent(in)                :: runind
	  real(8)                            :: aa,bb,dc,Vc,stiffk,vo,c
	  real(8), dimension(:), intent(inout) :: tau, slipv, sliprate, state
      	  real(8), dimension(:), allocatable :: rnstry_true
	  integer, parameter                 :: nmax=100000, npmax=10
	  integer, parameter                 :: neqs=2
	  real(8), parameter                 :: eps = 1.e-30
	  real(8)                            :: theta
	  real(8)                            :: fric,slip
	  real(8), dimension(neqs)           :: yt, dydt ,yscal, ytsv, dydtsv
	  real(8)                            :: stiff_temp, fitmis
	  real(8)                            :: fric_init, stiff_old
	  real(8)                            :: dx,dt_try, temp1, temp2
	  real(8)                            :: dt_next, acc, dt_did
	  real(8)                            :: time, timesv, slipsv
	  real(8), dimension(:), allocatable :: stiff_vec, par
	  integer                            :: i, ioerr, ne, Ncpts, stiff_err,fnloc
	  integer, dimension(:), allocatable :: IndCpts, n_stiff
	  character                          :: ln
	  character, dimension(29)           :: fm
	      
	     
	    allocate(rnstry_true(size(rnstry)))
	    allocate(par(7))
	    call getparams(rnstry,rnstry_true,runind,aa,bb,dc,Vc,stiffk,c) 

!-------------------------------------------------------------------------------!
! Initialize
!-------------------------------------------------------------------------------!

      fitmis      = 0.d0
      time        = t(1)
      slipv(1)     = x_f(1)
      theta       = dc/v_lp(1)
      fric_init   = aa*log(v_lp(1))+bb*log(theta)
      yt(2)       = v_lp(1)
      yt(1)       = theta      
      acc         = 1.e-12
      dt_try      = t(2) - t(1)
      stiff_err   = 0
      stiff_old   = stiffk
      tau(1)      = fric_init
      sliprate(1) = yt(2)
      state(1)    = yt(1)
      
!-------------------------------------------------------------------------------!
! Time loop
!-------------------------------------------------------------------------------!
            do i = 2, Nrange
		  vo = v_lp(i)
		  par = [aa,bb,dc,c,Vc,stiffk,vo]
		  do while (abs(time-t(i)) >= acc)
		  	  call derivs(time,yt,dydt,par)
		  	  yscal(:)=abs(yt(:))+abs(dt_try*dydt(:))
		  	  ytsv = yt
		  	  dydtsv = dydt
		  	  timesv = time
		  	  slipsv = slip
		  	  call rkqs(yt,dydt,time,dt_try,acc,yscal,dt_did,dt_next,derivs,par,fnloc)
		  	  if (time.gt.t(i)) then
		  	  	yt   = ytsv
		  	  	dydt = dydtsv
		  	  	time = timesv
		  	  	dt_try = t(i) - time  
		  	  else
		  	  	slipv(i)= slipv(i-1) + yt(2)*dt_did
		  	  	dt_try = dt_next
		  	  endif
          	  end do
		  state(i)     = yt(1)
		  sliprate(i)  = yt(2)
		  tau(i)       = aa*log(sliprate(i))+bb*log(state(i))
      	   enddo
      tau = tau - fric_init
   end subroutine fmodel
   !------------------------------------------------------------------------------------------------!
   !------------------------------------------------------------------------------------------------!
   
   function gauss(x,mu,std)
   
   real(8)                              :: gauss
   real(8), dimension(:), intent(in)    :: x, mu
   real(8), dimension(:,:), intent(in)  :: std
   real(8), dimension(:,:), allocatable :: cinv
   real(8)                              :: arg, N, denom, det_s
   real(8), parameter                   :: pi = 4.d0*atan(1.d0)
   
   allocate(cinv(size(std,1),size(std,2)))
   cinv=invertmat_sym(std)
   ! Mahalanobis distance for x
   arg   = -0.5*dot_product(matmul(cinv,x-mu),(x-mu))
   N     = dble(size(x))
   det_s = det_sym(std)
   denom = ((sqrt(2.d0*pi))**N)*det_s
   gauss = exp(arg)/denom
   
   end function gauss
   
   !---------------------------------------------------------------------------------------------------!
   !---------------------------------------------------------------------------------------------------!
   function corr_uni_rand(x,std) result(outvec)
   
   real(8), parameter                   :: pi = 4.d0*atan(1.d0)
   real(8), dimension(:), intent(in)    :: x
   real(8), dimension(:,:), intent(in)  :: std
   real(8), dimension(:,:), allocatable :: corr_BP,corr_SPE,corr_SPEUD
   real(8), dimension(:), allocatable   :: outvec
   integer                              :: i,j,info,Ns
 Ns = size(std,1)
 ! x must be normally distributed random numbers
 allocate(corr_SPE(Ns,Ns),corr_BP(Ns,Ns),corr_SPEUD(Ns,Ns),outvec(Ns))

 ! Convert STD matrix to corrcoef
 do i = 1,Ns
 	do j = i,Ns
		corr_BP(i,j) = std(i,j)/sqrt(std(i,i)*std(j,j))
		corr_BP(j,i) = corr_BP(i,j)
	enddo
 enddo
 corr_SPE = corr_BP
 ! adjust correlations for uniforms
 do i = 1,Ns
     do j = max(Ns-1,i),Ns
         if (i /= j) then
             corr_SPE(i, j) = 2.d0 * sin(pi * corr_BP(i, j) / 6.d0);
             corr_SPE(j, i) = 2.d0 * sin(pi * corr_BP(j, i) / 6.d0);
         endif
     enddo
 enddo

 ! induce correlation
 call covtor(corr_SPE,corr_SPEUD,info)
 outvec = matmul(corr_SPEUD,x);
 
 ! convert to uniform rand
 outvec = 0.5d0*(1.d0+erf(outvec/sqrt(2.d0)))
 
 end function corr_uni_rand
 !---------------------------------------------------------------------------------------------------!
 !---------------------------------------------------------------------------------------------------!
 function mix_gauss2(p,mu1,std1,mu2,std2,mix_factor) result(pout)
 !!! Function for generating a parameter proposal from the mixture of two Gaussians with different
 !!! means and std-s according to a predefined probability, mix_factor, to choose either.
 implicit none
 
 real(8), dimension(:), intent(in)       :: p
 real(8), dimension(:), intent(in)       :: mu1, mu2 
 real(8), dimension(:,:), intent(in)     :: std1, std2
 real(8), intent(in)                     :: mix_factor 
 real(8), dimension(size(p))		 :: pout
 real(8), dimension(size(mu1))           :: mu
 real(8), dimension(size(mu1),size(mu1)) :: std
 real(8)                                 :: unirand
 
   pout = p
   call random_number(unirand)
   if (unirand.le.mix_factor) then
	std = std1
	mu  = mu1
   else
	std = std2
	mu  = mu2
   endif
   call nrmrnd(pout,mu,std,'Basic') ! The 'Polar' generation is less robust        
   end function mix_gauss2
 !---------------------------------------------------------------------------------------------------!
 !---------------------------------------------------------------------------------------------------!
   end module util_mcmc
   
