module numeric_kind
  integer, parameter :: sp = kind(1.0),&
                        dp = selected_real_kind(2* precision(1.0_sp)), &
                        qp = selected_real_kind(2* precision(1.0_dp))
  real(kind=dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
end module numeric_kind
module eigens
  use numeric_kind
!==============================================
!- Subroutine to compute eigenvalues
!
!    Date        Programmer            Comment
!  =======       ==========            =======
!  6/12/2012       Amin                version f95 using Lapack
!
!- Pass the dimension of the Hamiltonian matrix to this subroutine
!==============================================
  contains
  subroutine diagon_d(mtx,eigens_r,dim)
    ! -inputs : mtx : the matrix to find the eigenvalues
    !           dim : the size of matrix
    ! -output : eigens : eigenvalues
    !===================================================
    
    Implicit None
    
    integer, intent(in) :: dim                            
    real(kind=dp), intent(in), dimension(dim,dim) :: mtx
    real(kind=dp),dimension(dim,dim) :: mtx_aux                 ! copy of Hamiltonian 
    real(kind=dp), dimension(dim,dim) :: vl, vr
    real(kind=dp), dimension(dim) :: eigens_r, eigens_i
    real(kind=dp), dimension(2*dim) :: work
    ! real(kind=dp) , dimension(2*dim) :: rwork 
    integer :: lda, ldvl, ldvr, lwork
    Integer :: info = 0
    character(len=1), parameter :: option_l = 'N', option_r = 'N' ! the eigenvectors are not computed
    ! change to 'V' to compute
    
    ! subroutine zgeev (JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO)
    !   ZGEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices
    
    
    mtx_aux = mtx                  !- copy of the Hamiltonian
    lda = dim                      !- leading dimension max(1,mlat)
    ldvl = 1; ldvr = 1             !- since the eigenvectors aren't computed 
    lwork = 4 * dim
    call dgeev(option_l, option_r, dim, mtx_aux, lda, eigens_r, eigens_i, vl, ldvl, vr, ldvr,&
         work, lwork, info)
    
    
    
    if (info .ne. 0) then
       if (info .lt. 0) write(*,*) info, "argument has illegal value!"
       if (info .gt. 0) write(*,*) info, "failed to compute"
    end if
    
    
  end subroutine diagon_d
end module eigens
!!========================================================================================================
!- to test if the diagon subroutine works fine
program eigen_test
  use numeric_kind
  use eigens
  implicit none
  integer, parameter :: dim = 2
  real(kind=dp), dimension(dim,dim) :: a_mtx
  real(kind=dp), dimension(dim*dim) :: a_arr
  real(kind=dp), dimension(dim) :: eig_v
  character(len=30) :: fmt
  a_arr = (/1,2,3,4/)
  a_mtx = reshape(a_arr, (/2,2/))
  call diagon_d(a_mtx,eig_v,dim)
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

end program eigen_test
