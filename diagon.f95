module numeric_kind
  integer, parameter :: sp = kind(1.0),&
                        dp = selected_real_kind(2 * precision(1.0_sp)), &
                        qp = selected_real_kind(2 * precision(1.0_dp))
  real(kind=dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
end module numeric_kind
! ============================================================
module zeigens
  use numeric_kind
  !- eigenvalue calculator module
  contains
  subroutine eigens_z(mtx,eigns,dim)
    !==============================================================
    !- Subroutine to compute eigenvalues of a complex square matrix
    !
    !    Date        Programmer            Comment
    !  =======       ==========            =======
    !  6/12/2012     Amin Ahmadi           version f95 using Lapack
    !
    !- input :  mtx   : matrix to compute the eigenvelues
    !           dim   : dimension of the matrix
    !
    !- output:  eigns : eigenvalues assoociated with mtx matrix
    !==============================================
    implicit none 
    
    integer, intent(in) :: dim
    complex(kind=dp), intent(in), dimension(dim,dim) :: mtx
    complex(kind=dp), dimension(dim) :: eigns
    complex(kind=dp),dimension(dim,dim) :: mtx_aux           ! copy of mtx 
    complex(kind=dp) :: vl(1,1), vr(1,1)
    complex(kind=dp), dimension(2*dim) :: work
    real(kind=dp) , dimension(2*dim) :: rwork
    integer ::  lwork
    integer :: lda, ldvl, ldvr
    Integer :: info = 0
    ! the eigenvectors are not computed, 'N' and 'N'
    character(len=1), parameter :: option_l = 'N', option_r = 'N' 
    ! change to 'V' to compute
    
    ! ======================================================================
    ! - general usage of zgeev subroutine
    ! subroutine zgeev (JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, 
    !                                        LDVR, WORK, LWORK, RWORK, INFO) 
    ! ZGEEV computes the eigenvalues and, optionally, 
    ! the left and/or right eigenvectors for GE matrices
    ! ======================================================================
    
    
    mtx_aux = mtx                  !- copy of the Hamiltonian
    lda = dim                      !- leading dimension max(1,mlat)
    ldvl = 1; ldvr = 1             !- since the eigenvectors aren't computed 
    lwork = 4 * dim
    call ZGEEV(option_l, option_r, dim, mtx_aux, lda, eigns, vl, ldvl, vr, ldvr,&
         work, lwork, rwork, info)
    
    
    
    if (info .ne. 0) then
       if (info .lt. 0) write(*,*) info, "argument has illegal value!"
       if (info .gt. 0) write(*,*) info, "failed to compute"
    end if
    
    
  end subroutine eigens_z
  
end module zeigens
!===========================================================================
!- test program to test the zeigens module and eigens_z subroutine
!- to test if the diagon subroutine works fine
program zeigens_test
  use numeric_kind
  use zeigens
  implicit none
  integer, parameter :: dim = 2
  complex(kind=dp), dimension(dim,dim) :: a_mtx
  complex(kind=dp), dimension(dim*dim) :: a_arr
  complex(kind=dp), dimension(dim) :: eig_v
  character(len=30) :: fmt
  a_arr = (/(1.0,0.0), (2.0,0.0), (3.0,0.0), (4.0,0.0)/)
  a_mtx = reshape(a_arr, (/2,2/))
  call eigens_z(a_mtx,eig_v,dim)
  print *, "============================================================="
  print *, "         This is a test for eigenvalues subroutine           "
  print *, "============================================================="
  print *, "        matrix A = [[1,2],                                   "
  print *, "                    [3,4]]"
  print *, "        ===================                                  "
  print *, "        Eigenvalues are =                                    "
  fmt = '(10x,f7.3,2x,f7.3)'
  print fmt, eig_v(1), eig_v(2)
  print *, "============================================================="
  print *,""
  print *,""
   
  ! notice that the arguments of subroutines zgeev and dgeev are different
  ! add a test file for a complex matrix as well
end program zeigens_test
