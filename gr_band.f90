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
  subroutine sqr_hu(hm, n_sites, phi)
    !=======================================================
    ! - input - mlat : size of one column of square lattice
    !           phi  : the wave number "k"
    !
    ! - output - hm :  the Hamiltonian matrix
    !=======================================================
    implicit none
    integer, intent(in) :: n_sites                 ! unitcell size
    complex(kind=dp), intent(out), dimension(n_sites,n_sites) :: hm
    real(kind=dp), parameter :: t = 1.0_dp       ! hopping paramater
    integer :: i,j, mlat
    real(kind=dp), intent(in) :: phi             ! wave number 'k'

    !- initializing
    hm = (0.d0, 0.d0)
    mlat = int(n_sites/2)
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
  !==========================================================
  ! - To construct the Hamiltonian and U-matrix for an armachair 
  ! - graphene. The input parameter is Mlat. The structure is defined 
  ! - as follows:  
  !   Mlat = 2  ===> n_sites = 4*Mlat + 2 = 10 # of site in the unit cell
  
  !      9*---^10 *---^
  !       |   |   |   |
  !      7^   *8--^   *
  !       |   |   |   |
  !      5*---^6  *---^
  !       |   |   |   |
  !      3^   *4--^3 4*
  !       |   |   |   |
  !      1*---^2  *---^
  !               1   2
  !=====================================================================
  
  subroutine ac_graph_h(hm,n_sites,phi)
    implicit none
    
    integer, intent(in) :: n_sites
    real(kind=dp), intent(in) :: phi
    complex(kind=dp), intent(out), dimension(n_sites,n_sites) :: hm
    real(kind=dp), parameter :: t = 1.d0
    integer :: i,j, mlat

    !- Initializing 
    hm = (0.d0,0.d0)
    mlat = int((n_sites - 2)/4)
    !==========================
    !- Constructing Hamiltonian
    !==========================
    do i = 0, mlat-1 
       j = i*4
       hm(j+1,j+2) = -t*(1.d0,0.d0)
       hm(j+2,j+1) = conjg(hm(j+1,j+2))
       
       hm(j+1,j+3) = -t*(1.d0,0.d0)
       hm(j+3,j+1) = conjg(hm(j+1,j+3))
       
       hm(j+2,j+4) = -t*(1.d0,0.d0)
       hm(j+4,j+2) = conjg(hm(j+2,j+4))
       !===================================
       if (i .ne. 0) then
          hm(j+3,j+5) = -t*(1.d0,0.d0)
          hm(j+5,j+3) = conjg(hm(j+3,j+5))
          
          hm(j+4,j+6) = -t*(1.d0,0.d0)
          hm(j+6,j+4) = conjg(hm(j+4,j+6))
       end if
       if (i .eq. mlat-1) then
          hm(j+5,j+6) = -t*(1.d0,0.d0)
          hm(j+6,j+5) = conjg(hm(j+5,j+6))
       end if
    end do
    !=================
    !- Constructing UR
    !=================
    
    do i = 4, n_sites, 4
       hm(i,i-1) = hm(i,i-1) - t*exp(phi*(0.d0,1.d0))
       hm(i-1,i) = conjg(hm(i,i-1))
    end do
    
  end subroutine ac_graph_h
  
end module hamiltonian
!===============================================================
program sqr_band
  use numeric_kind
  use hamiltonian
  use zeigens
  implicit  none

  integer, parameter :: mlat = 2
  real(kind=dp) :: phi = 1.d0
  complex(kind=dp), allocatable, dimension(:,:) :: hm
  integer :: i, ik, num = 100, opt = 0, n_sites!,j
  complex(kind=dp), allocatable, dimension(:) :: eigns
  character(len=30) :: fmt

  fmt = '(30(2x,f5.2))'
  print *, "Enter 0 for square lattice and 1 for graphene"
  read *,opt
  if (opt == 0) then
     n_sites = 2*mlat
     allocate(hm(n_sites,n_sites), eigns(n_sites))
     do ik = 0, num
        phi = -pi + (2_dp*pi/num) * ik
        call sqr_hu(hm, n_sites, phi)
        call eigens_z(hm,eigns,n_sites)
        print fmt, phi,  (real(eigns(i)), i = 1, n_sites)
     end do
  elseif (opt == 1) then
     n_sites = 4*mlat + 2
     allocate(hm(n_sites,n_sites), eigns(n_sites))
     do ik = 0, num
        phi = -pi + (2_dp*pi/num) * ik
        call ac_graph_h(hm, n_sites, phi)
        call eigens_z(hm,eigns,n_sites)
        print fmt, phi,  (real(eigns(i)), i = 1, n_sites)
     end do
  end if
  deallocate(hm,eigns)
end program sqr_band
