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
!=================================================================================
module hamiltonian
  use numeric_kind
contains
  !===============================================
  ! construct the Hamiltonian and U-matrix for a
  ! square lattice with width mlat. The structure 
  ! has a form of
  !=============================
  !        8               nlat*mlat
  !  4 0-----0-----0-----0-----0-----4
  !    |     |     |     |     |     
  !    |   7 |     |     |     |  PBC 
  !  3 0-----0-----0-----0-----0-----3
  !    |     |     |     |     |
  !    |   6 |     |     |     |
  !  2 0-----0-----0-----0-----0-----2
  !    |     |     |     |     |   
  !    |   5 |     |     |     |  PBC 
  !  1 0-----0-----0-----0-----0-----1
  !=============================
  ! two columns must be considered as a unit cell to 
  ! take into account the Bloch translation prefactor.
  ! ================================================
  subroutine sqr_hu(hm, mlat, phi)
    !=======================================================
    ! - input - mlat : size of one column of square lattice
    !           phi  : the wave number "k"
    !
    ! - output - hm :  the Hamiltonian matrix
    !=======================================================
    use numeric_kind
    implicit none
    integer, intent(in) :: mlat                  ! size of one column
    complex(kind=dp), intent(out), dimension(2*mlat,2*mlat) :: hm
    real(kind=dp), parameter :: t = 1.0_dp       ! hopping paramater
    integer :: i,j
    real(kind=dp), intent(in) :: phi             ! wave number 'k'

    !- initializing
    hm = (0.d0, 0.d0)
    !==========================
    !- Constructing Hamiltonian
    !==========================
    do i = 1, mlat-1                       ! vertical connection
       hm(i,i+1) = -t*(1.d0,0.d0)
       hm(i+1,i) = conjg(hm(i,i+1))
       j = i + mlat
       hm(j,j+1) = -t*(1.d0,0.d0)
       hm(j+1,j) = conjg(hm(j,j+1))
    end do 

    do i = 1, mlat                         !- horizontal connection
       hm(i,i+mlat) = -t*(1.d0,0.d0)
       hm(i+mlat,i) = conjg(hm(i,i+mlat))
    end do
    !  no PBC
    !===============================================
    !- Constructing translational part of Hamiltonian
    !===============================================

    do i = 1, mlat
       hm(i,i+mlat) = hm(i,i+mlat)  - t*exp(phi*(0.d0,1.d0))
       hm(i+mlat,i) = conjg(hm(i,i+mlat))
    end do
  end subroutine sqr_hu
end module hamiltonian
!===============================================================
program sqr_band
  use numeric_kind
  use hamiltonian
  use zeigens
  implicit  none

  integer, parameter :: mlat = 4
  real(kind=dp) :: phi = 1.d0
  complex(kind=dp), dimension(2*mlat,2*mlat) :: hm
  integer :: i, ik, num = 100!,j
  complex(kind=dp), dimension(2*mlat) :: eigns
  character(len=30) :: fmt

  fmt = '(30(2x,f5.2))'
  do ik = 0, num
     phi = -pi + (2_dp*pi/num) * ik
     call sqr_hu(hm, mlat, phi)
     call eigens_z(hm,eigns,2*mlat)
     ! do i = 1 , mlat
     !    print fmt, (real(hm(i,j)), j=1,mlat)
     ! end do
     print fmt, phi,  (real(eigns(i)), i = 1, 2*mlat)
  end do
end program sqr_band
