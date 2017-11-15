module matutils

  implicit none

contains

!!! we need these matrix multiplication utilities to use the Cholesky
!!! factor provided by LAPACK routine, where C=R'R, R upper diagonal
!!! and the lower diagonal entries are arbitrary.

!!!
!!! calculates covariance matrix of x,
!!! with optional (frequency) weights in w,
!!! returns covmat
!!! 
  function covmat(x, w)
    implicit none
    
    real(8), intent(in) :: x(:,:)
    real(8), intent(in), optional :: w(:)

    ! local variables
    real(8) :: covmat(size(x,2),size(x,2))
    real(8) :: xmean(size(x,2))
    integer :: n,p,i,j
    real(8) :: wsum2, w2, w3

    n=size(x,1)
    p=size(x,2)

    !! optional weight in w
    if (present(w) .and. size(w) == n) then
       w2 = -1.d0
       wsum2 = sum(w)
    else if (present(w) .and. size(w) == 1) then
       w2 = w(1)
       wsum2 = dble(n)*w2
    else if (.not. present(w)) then
       w2 = 1.d0
       wsum2 = dble(n)
    else
       stop 'covmat: invalid weight in w'
    end if

   do i=1,p
          if (w2 == -1.d0) then
             xmean(i) = sum(x(:,i)*w)/wsum2
          else
             xmean(i) = sum(x(:,i)*w2)/wsum2
          end if
   end do
       do i=1,p
          do j=1,i
             if ( w2 == -1.d0 ) then
                covmat(i,j) = sum((x(:,i)-xmean(i))*(x(:,j)-xmean(j))*w)/(wsum2 - 1.d0)
             else
                covmat(i,j) = sum((x(:,i)-xmean(i))*(x(:,j)-xmean(j))*w2)/(wsum2 - 1.d0)
             end if
             if (i/=j) covmat(j,i) = covmat(i,j)
          end do
       end do

  end function covmat
!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!
!!! Calculates the empirical covariance matrix given n-samples of a m-dimensional vector as an nxm matrix
!!! Formula from Yan-Bai, 2009 (An Adaptive Directional Metropolis-within-Gibbs algorithm), corrected Yan_Bai's
!!! formula to accommodate the full covariance matrix.

function emp_covmat(x)   result(covmat)
    implicit none
    
    real(8), intent(in) :: x(:,:)

    ! local variables
    real(8) :: covmat(size(x,2),size(x,2))
    real(8) :: xmean(size(x,2),1),x_temp(size(x,2),size(x,1)),x_loc(size(x,2),1)
    integer :: n,npar,i
    
    n = size(x,1) ! Number of samples
    npar = size(x,2) ! Number of parameters
    
    x_temp      = transpose(x)
    xmean(:,1)  = sum(x_temp,2)/n
    covmat      = 0.d0
    do i = 1, n
    	x_loc(:,1)  = x_temp(:,i) - xmean(:,1)
    	covmat = covmat + matmul(x_loc,transpose(x_loc))   ! Eq. 2.5 from Bai (2009)
    enddo
    covmat = covmat/(n-1)   

end function emp_covmat
!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!
! Update the existing covariance matrix by adding the most recently sampled observation
! Iterative update of covariance matrix to implement Roberts and Rosenthals adaptive approach
function update_covmat(x,x_mean,n,std)  result(stdout)

 real(8), intent(in)    :: x(:,:), x_mean(:,:) ! x - newest sampled vector of size [NP,1], x_mean - updated mean over all samples, size [NP,1]
 real(8), intent(in)    :: std(:,:) ! Old covariance matrix, size [NP,NP]
 real(8)                :: stdout(size(std,1),size(std,2))
 integer, intent(in)    :: n ! Number of samples drawn
 
 stdout = (n-2)/(n-1)*std + matmul(x-x_mean,transpose(x-x_mean))/(n-1)

end function update_covmat

!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!
!!! Calculates upper diagonal Cholesky transformation
!!! of symmetric matrix 'cmat'
  subroutine covtor(cmat,R,info)
    implicit none
    real(8) :: cmat(:,:)
    real(8), optional :: R(:,:)
    integer, intent(out), optional :: info

    integer :: info2, n

    n = size(cmat,1)
    if (size(cmat,2) /= n) then
       stop 'size error in covtor'
    end if

    if (present(R)) then
       if (size(R,1) /= n .or. size(R,1) /= size(R,2)) then
          stop 'size error in R in covtor'
       end if
       R = cmat
       call dpotrf('l',n,R,n,info2)
    else
       call dpotrf('l',n,cmat,n,info2)
    end if

    if (present(info)) then
       info=info2
    else  if (info2 .ne. 0) then
       stop 'error in Chol'
    end if

  end subroutine covtor
  
!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!  
!!! svd of covariance matrix
!!!
  function covtor_svd(cmat,info) result(y)
    implicit none
    real(8),intent(in) :: cmat(:,:)
    integer, intent(out), optional :: info
    real(8), allocatable :: y(:)
    integer :: info2, n, i
    real(8), allocatable :: work(:)
    real(8), allocatable :: u(:,:), s(:), R(:,:)
    integer :: lwork

    real(8) :: tol, cond

    n = size(cmat,1)
    if (size(cmat,2) /= n) then
       stop 'size error in covtor'
    end if

    lwork = 5*n
    allocate(work(lwork),u(n,n),R(n,n),s(n),y(n),stat=info2)
    if (info2 /= 0) then
       stop 'memory allocation in svd'
    end if

    if (size(R,1) /= n .or. size(R,1) /= size(R,2)) then
       stop 'size error in R in covtor'
    end if
    R = cmat

!!! call dgesvd(jobu,jobvt,m,n,A,lda,s,u,ldu,vt,ldvt,work,lwork,info)
    call dgesvd('A','N',n,n,R,n,s,u,n,u,n,work,lwork,info2)
!   write(*,*) 'S       : ', real(s(max(1,n-3):n))
    if (s(1) == 0.0d0) then
       write(*,*) 'SVD: s(1)=0'
       if (present(info)) then
          info=n
       end if
       return
    end if
    if (info2 > 0) then
       write(*,*) 'svd info = ', info2
    end if
    y = s

    deallocate(work,u,s)

  end function covtor_svd
!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!
!!! invert symmetric matrix using Cholesky decomposition
!!!
  function invertmat_sym(x,info0) result(y)
    implicit none
    real(8), intent(in) :: x(:,:)
    integer, optional, intent(inout) :: info0
    real(8) :: y(size(x,1),size(x,2))

    integer :: npar, info

    npar = size(x,1)
    if (npar .ne. size(x,2)) then
       write(*,*) 'ERROR: x must be square'
       stop
    end if
    y = x
    call dpotrf('u',npar,y,npar,info)
    if (info .ne. 0) then
       if (present(info0)) then
          info0 = info
          return
       else
          write(*,*) 'ERROR: cannot factor cmat, info = ', info
          stop
       end if
    end if
    call dpotri('u',npar,y,npar,info)
    if (info .ne. 0) then
       if (present(info0)) then
          info0 = info
          return
       else
          write(*,*) 'ERROR: cannot invert cmat, info = ', info
          stop
       end if
    end if

    !! fill the lower diagonal with the upper
    call copylower(y)
    if (present(info0)) info0=0

  end function invertmat_sym

!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!
!!! Determinant of a symmetric matrix 
  function det_sym(A) result(y)
  implicit none
  real(8),dimension(:,:)              :: A
  real(8)                             :: y
  real(8), dimension(:), allocatable  :: eig
  
  allocate(eig(size(A,1)))
  eig = covtor_svd(A)
  y = product(eig)
  
  end function det_sym
  
!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------! 
  
  	
!!! Mahalanobis distance of rows of x from mu given cmat
  function mahalanobis(x,mu,cmat) result(d)
    implicit none
    real(8) :: x(:,:), mu(:), cmat(:,:), d(size(x,1))
    real(8) :: cinv(size(cmat,1),size(cmat,2)), x0(size(x,1),size(x,2))
    integer :: i, n

    n = size(x,1)
    ! remove mean
    do i=1,n
       x0(i,:) = x(i,:)-mu
    end do
    ! invert cmat
    cinv=invertmat_sym(cmat)
    ! Mahalanobis distance for each row in x
    d = sum(matmul(x0,cinv)*x0,2)

  end function mahalanobis
  
!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!
!!! copy the lower diagonal from upper diagonal for symmetric x
!!!
  subroutine copylower(x)
    implicit none
    real(8), intent(inout) :: x(:,:)
    integer i,m,n
    m = size(x,1)
    n = size(x,2)
    do i=1,min(n,m)-1
       x((i+1):n,i) = x(i,(i+1):n)
    end do
  end subroutine copylower
!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!
!!! Square Root of a positive semidefinite square matrix by diagonalization
!!!
subroutine sqrtm(cmat,R,n,INFO)

	implicit none
	real(8), intent(in)       :: cmat(:,:)
	real(8), intent(inout)    :: R(:,:)	
	integer, intent(in)	  :: n
	real(8)                   :: cin(n,n)
	integer, intent(inout)    :: INFO
	real(8)                   :: abstol
	integer                   :: i
	! Variables specific to the LAPACK routine
	integer	                  :: LWORK, LIWORK ! Default minimum values as outlined in LAPACK docs
	integer                   :: M
	real(8)                   :: W(26*n), Z(n,n), WORK(26*n), lambda(n,n)
	integer                   :: ISUPPZ(2*n), IWORK(10*n)
	
	
! First use LAPACK libraries to define abstol for use in diagonalization routine
	!abstol = DLAMCH('S') ! ABSTOL set to 'safe minimum' as noted in LAPACK library documentation
	
	if (size(cmat,2) /= size(cmat,1)) then
       		stop 'size error in sqrtm' ! Check for square matrix
    	end if
    	
    	cin = cmat
    	
    	! Call the eigenvalue computation LAPACK subroutine	
	CALL DSYEVR('V','A','U',n,cin,n,0.d0,0.d0,0,0,1.d-12,M,W,Z,n,ISUPPZ,WORK,50*n,IWORK,20*n,INFO)
	
	! Initiate the sqrt of eigenvalue matrix
	lambda = 0.d0
	do i = 1,n
		lambda(i,i) = sqrt(W(i))
	enddo
	! Perform similarity transform to compute square-root
	R = matmul(matmul(Z,lambda),transpose(Z))
		
end subroutine sqrtm	

!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!
!!! Solve a square linear system using dgesv from LAPACK
!!!
 subroutine linsolve(Amat,solvec, info)

 real(8), dimension(:,:), intent(inout)   :: Amat ! Input matrix, replaced on output by LU factors
 real(8), dimension(:), intent(inout)     :: solvec ! Input vector b in Ax=b, replaced on output by solution x
 integer, intent(inout)                   :: info
 integer, dimension(:), allocatable       :: ipiv
 integer :: i, n

	if (size(Amat,2) /= size(Amat,1)) then
       		stop 'size error in sqrtm' ! Check for square matrix
    	end if    
	n = size(Amat,1) ! Amat is square
	allocate(ipiv(n))
	call dgesv(n, 1, Amat, n, ipiv, solvec, n, info)
	!     Check for the exact singularity.
        if ( info.gt.0 ) then
         	write(*,*)'Warning! Matrix is singular to working precision, not computed'
		solvec = 0.d0
         	! stop
        end if
 end subroutine linsolve

!-----------------------------------------------------------------------------------------------------!
!-----------------------------------------------------------------------------------------------------!

end module matutils
