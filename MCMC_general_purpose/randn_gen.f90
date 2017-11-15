module randn_gen

use random_gen

implicit none


contains

!--------------------------------------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------------------------------------
	 subroutine nrmrnd_scl(random,mu,sigma,method)
	 ! This subroutine uses the Box-Muller transform to compute univariate normal deviates
	 ! random : A scalar uniform random number, transformed into a normal random number on output
	 ! mu: Scalar mean
	 ! sigma: Scalar Standard Deviation
	 ! method: A string signifying the method used: 1) 'Polar' - Uses points inside the 1 circle, use with a fast uniform random number generator
	 !						2) 'Basic' - Uses trigonometric functions, for use in platforms with fast trig functions
		 real(8), parameter                    :: pi = 4.d0*atan(1.d0)
		 real(8), dimension(:), intent(inout)  :: random
		 real(8), intent(in)                   :: mu
		 real(8), intent(in)                   :: sigma
		 real(8), dimension(2)                 :: r_prox
		 character(len=5), intent(in)          :: method
		 integer                               :: nr, i
	 
		 nr = size(random)
		 i  = 1
		 if (modulo(nr,2)==0) then
		 	do while (i <= nr-1) 
				 call stndnrmrnd(random(i:i+1),method)
				 i = i + 2
			enddo
		else
			do while (i <= nr-2) 
				 call stndnrmrnd(random(i:i+1),method)
				 i = i + 2
			enddo
			call random_number(r_prox)
			call stndnrmrnd(r_prox,method)
			random(nr) = r_prox(1)
		endif
		random = sigma*random + mu
		 
	 end subroutine nrmrnd_scl
!--------------------------------------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------------------------------------------------------------
	 subroutine nrmrnd(random,mu,sigma,method)
	 ! This subroutine uses the Box-Muller transform to compute univariate normal deviates
	 ! random : A n-tuple vector of uniform random numbers, transformed into n normal random numbers on output
	 ! mu: Mean (can be a vector or scalar)
	 ! sigma: Standard Deviation (can be a vector or scalar)
	 ! method: A string signifying the method used: 1) 'Polar' - Uses points inside the 1 circle, use with a fast uniform random number generator
	 !						2) 'Basic' - Uses trigonometric functions, for use in platforms with fast trig functions
		 real(8), parameter                    :: pi = 4.d0*atan(1.d0)
		 real(8), dimension(:), intent(inout)  :: random
		 real(8), dimension(:), intent(in)     :: mu
		 real(8), dimension(:,:), intent(in)   :: sigma
		 real(8), dimension(2)                 :: r_prox
		 character(len=5), intent(in)          :: method
		 integer                               :: nr, i
	 
		 nr = size(random)
		 i  = 1
		 if (modulo(nr,2)==0) then
		 	do while (i <= nr-1) 
				 call stndnrmrnd(random(i:i+1),method)
				 i = i + 2
			enddo
		else
			do while (i <= nr-2) 
				 call stndnrmrnd(random(i:i+1),method)
				 i = i + 2
			enddo
			call random_number(r_prox)
			call stndnrmrnd(r_prox,method)
			random(nr) = r_prox(1)
		endif
		random = matmul(sigma,random) + mu
		 
	 end subroutine nrmrnd
!--------------------------------------------------------------------------------------------------------------------------------------------------------
	subroutine stndnrmrnd(rand,method)
	
	! This subroutine uses the Box-Muller transform to compute univariate standard normal deviates
	 ! rand  : A 2-tuple vector of uniform random numbers, transformed into 2 standard normal random numbers on output
	 ! method: A string signifying the method used: 1) 'Polar' - Uses points inside the 1 circle, use with a fast uniform random number generator
	 !						2) 'Basic' - Uses trigonometric functions, for use in platforms with fast trig functions
		 real(8), parameter                    :: pi = 4.d0*atan(1.d0)
		 character(len=5), intent(in)          :: method
		 real(8), dimension(2), intent(inout)  :: rand
		 real(8)                               :: R, Theta
	
				select case(method)
				 	case("Polar")
				 		rand = -1.d0 + 2.d0*rand
				 		R = sum(rand**2.d0)
				 		do while (R>=1.d0)
							call random_number(rand)
							rand = -1.d0 + 2.d0*rand
							R = sum(rand**2.d0)
						enddo
						rand = sqrt(-2.d0*log(R)/R)*rand
					case("Basic")
						R = sqrt(-2.d0*log(rand(1)))
						Theta = 2.d0*pi*rand(2)
						rand(1) = R*cos(Theta)
						rand(2) = R*sin(Theta)
					case default
						write(*,*) 'Unknown method, set either Polar or Basic'
						call ABORT
				end select
				
	end subroutine stndnrmrnd

end module randn_gen
