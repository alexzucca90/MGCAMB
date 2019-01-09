    ! Module to analyze mock WL shear datasets
    ! Based on the standard wl.f90 module and on 1210.2194
    ! It can also be used to generate mock datasets (just set do_mock=T in your .ini file)
    ! Last update: 04/02/2016

    module wl_mock
    use settings
    use CosmologyTypes
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    use WL_MockMaker
    implicit none
    private

    type, extends(TCosmoCalcLikelihood)          :: WLmockLikelihood

        integer                                  :: num_z_bins,num_cross,ell_max,ell_min,num_z_p
        real(mcp), allocatable, dimension(:)     :: z_bins
        real(mcp), allocatable, dimension(:)     :: dl_obs
        real(mcp), allocatable, dimension(:,:)   :: p, data_auto, data_cross, data_noise
        real(mcp), allocatable, dimension(:,:,:) :: data_covmat
        real(mcp), allocatable, dimension(:)     :: z_p,ell                    
        real(mcp)                                :: fsky
        logical                                  :: use_non_linear                        ! Whether to use non-linear corrections
        logical                                  :: use_weyl                              ! Whether to use Weyl potential or matter P(k)

    contains
    procedure          :: LogLike => WLmock_LnLike
    procedure          :: ReadIni => WLmock_ReadIni
    procedure, private :: get_convergence !computes the shear C_l

    end type WLmockLikelihood


    logical :: do_mock
    integer :: mock_type,bin_type
    character(100) :: mock_root, sys_root, sys_dir

    logical :: WL_debugging=.false.

    ! Integration accuracy parameters
    real(mcp), parameter :: const_pi = 3.1415926535897932384626433832795

    public WLmockLikelihood, WLmockLikelihood_Add
    contains


    subroutine WLmockLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(WLmockLikelihood), pointer :: this
    Type(TSettingIni) :: DataSets
    integer i
    logical :: nonlinear, useweyl
    real(mcp) :: kmax_temp


    if (WL_debugging) then
       write(0,*)'-------------------------------------------'
       write(0,*)'                  WARNING                  '
       write(0,*)'you are using the module in debugging mode '
       write(0,*)'the MCMC will do one iteration, compute    '
       write(0,*)'some test files and then stop              '
       write(0,*)'To run normally change the WL_debugging    '
       write(0,*)'flag in wl_mock.f90                        '
       write(0,*)'-------------------------------------------'
    end if


    if (Ini%Read_Logical('use_WL_mock',.false.)) then
        nonlinear = Ini%Read_Logical('wl_use_non_linear')!,.true.)
        useweyl = Ini%Read_Logical('wl_use_weyl')!,.true.)

        do_mock = Ini%Read_Logical('do_mock',.false.)
        if (do_mock) then
            mock_type = Ini%Read_Int('mock_type')
            bin_type = Ini%Read_Int('bin_type')
            mock_root = Ini%ReadFileName('mock_root')
        end if

        allocate(this)
        this%needs_nonlinear_pk = nonlinear
        this%use_non_linear = nonlinear
        this%use_weyl = useweyl
        call this%ReadDatasetFile(Ini%ReadFileName('wl_mock_dataset'))
        this%LikelihoodType = 'WL_mock'
        this%needs_powerspectra = .true.
        this%needs_weylpower = useweyl
        call LikeList%Add(this)

        if (Feedback>1) write(*,*) 'read WL data sets'
    end if

    end subroutine WLmockLikelihood_Add


    subroutine WLmock_ReadIni(this, Ini)
    use MatrixUtils
    class(WLmockLikelihood) this
    class(TSettingIni) :: Ini
    character(LEN=:), allocatable :: measurements_autofile, window_file, measurements_crossfile, noise_file, covmat_file
    character(2) :: indbin
    Type(TTextFile) :: F
    Type(WLspec) :: SP
    real(mcp) :: dummy1,dummy2,pnorm
    integer :: il,i,iz,it,ib,jb,j,k,nt,izl,izh

    if (.not.do_mock) then

        if (Feedback > 0) write (*,*) 'reading WL data set: '//trim(this%name)

        this%num_z = 100 !Look better at this!!! WITHOUT THIS EVERYTHING IS FUCKED UP... WHY???????
        this%max_z = Ini%Read_Double('max_z')
        this%ell_max = Ini%Read_Int('ell_max')
        this%ell_min = Ini%Read_Int('ell_min')
        this%kmax = Ini%Read_Double('kmax')
        this%fsky = Ini%Read_Double('fsky')
        this%num_z_bins = Ini%Read_Int('num_z_bins')
        this%num_cross = (this%num_z_bins*(this%num_z_bins+1)/2)-this%num_z_bins
        this%num_z_p = Ini%Read_Int('num_z_window')

        allocate(this%ell(this%ell_max))
        allocate(this%z_bins(this%num_z_bins+1))
        allocate(this%data_covmat(this%ell_max,this%num_z_bins,this%num_z_bins))
        allocate(this%dl_obs(this%ell_max))
        allocate(this%data_noise(this%ell_max,this%num_z_bins))
        allocate(this%z_p(this%num_z_p))
        allocate(this%p(this%num_z_p,this%num_z_bins))


        do iz=1,this%num_z_bins+1
           write(indbin,'(I2)') iz
           this%z_bins(iz)=Ini%Read_Double('zbin('//trim(adjustl(indbin))//')')
        end do


        covmat_file            = Ini%ReadFileName('WL_covmat')        
        noise_file             = Ini%ReadFileName('WL_noise')
        window_file            = Ini%ReadFileName('window_file')

        open(42,file=trim(window_file))
        do iz=1,this%num_z_p
           read (42,*) this%z_p(iz),this%p(iz,:)
        end do
        close(42)


        call F%Open(covmat_file)
        do il=this%ell_min,this%ell_max
           do ib = 1, this%num_z_bins
              read (F%unit,*) (this%data_covmat(il,ib,jb),jb=1,this%num_z_bins)
           end do
        end do
        call F%Close()


        call F%Open(noise_file)
        do il=this%ell_min,this%ell_max
           read (F%unit,*) this%ell(il), (this%data_noise(il,ib),ib=1,this%num_z_bins)
        end do
        call F%Close()

        !adding noise to covmat
        do il=this%ell_min,this%ell_max
           do ib=1,this%num_z_bins
              this%data_covmat(il,ib,ib) = this%data_covmat(il,ib,ib) + this%data_noise(il,ib)
           end do
        end do


        do il=this%ell_min,this%ell_max
           this%dl_obs(il) = exp(MatrixSym_LogDet(this%data_covmat(il,:,:)))
        end do


    else

        write(0,*) 'Doing mock WL dataset'
        call GetSpecs(mock_type,mock_root,bin_type,SP)


        this%num_z = Ini%Read_Int('nz_wl',100)!SP%num_z_temp
        this%max_z = SP%zmax
        this%ell_max = SP%lmax
        this%ell_min = SP%lmin
        this%kmax = SP%kmax
        this%num_z_bins = SP%numbins
        this%num_cross = SP%numcross
        this%num_z_p = SP%num_z_temp

        allocate(this%ell(this%ell_max))
        allocate(this%z_bins(this%num_z_bins+1))
        allocate(this%z_p(this%num_z_p))
        allocate(this%p(this%num_z_p,this%num_z_bins))


        do i=1,this%num_z_bins+1
            this%z_bins(i)=SP%zbins(i)
        end do

        open(42,file='mock_output/WL_mock/'//trim(mock_root)//'_window.dat')
        do iz=1,this%num_z_p
           read (42,*) this%z_p(iz),this%p(iz,:)
        end do
        close(42)

        deallocate(SP%zbins,SP%frac)

    end if

    !Normalize window functions p so \int p(z) dz = 1
    do ib=1,this%num_z_bins
        pnorm = 0
        do iz=2,this%num_z_p
            pnorm = pnorm + 0.5d0*(this%p(iz-1,ib)+this%p(iz,ib))*(this%z_p(iz)-this%z_p(iz-1))
        end do
        this%p(:,ib) = this%p(:,ib)/pnorm
    end do


    if (WL_debugging) then
       open(77,file='wldebug_files/normalized_windows.dat')
       do iz=1,this%num_z_p
          write(77,55) this%z_p(iz),(this%p(iz,i),i=1,this%num_z_bins)
       end do
       close(77)
    end if

    55 FORMAT(100(1X,e15.9))

    end subroutine WLmock_ReadIni







    function WLmock_LnLike(this, CMB, Theory, DataParams)
    use MatrixUtils
    Class(WLmockLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) WLmock_LnLike, sum, term
    real(mcp), allocatable, dimension(:,:) :: TH_auto,TH_cross,TH_auto2,TH_cross2,dmix_temp
    real(mcp), allocatable, dimension(:,:,:) :: TH_matrix
    real(mcp), allocatable, dimension(:,:) :: mixmats,invmat,outmat
    real(mcp), allocatable, dimension(:) :: dl_th, dl_mix
    real(mcp), allocatable, dimension(:) :: ell
    real(mcp) :: DMGT,temp,detp,detm,ppart,mpart,trace
    integer :: l,bin,ib,jb,k,t,i

    allocate(TH_auto(this%ell_max,this%num_z_bins))
    allocate(TH_cross(this%ell_max,this%num_cross))
    allocate(TH_matrix(this%ell_max,this%num_z_bins,this%num_z_bins))
    allocate(dl_th(this%ell_max),dl_mix(this%ell_max))
    allocate(mixmats(this%num_z_bins,this%num_z_bins),invmat(this%num_z_bins,this%num_z_bins),outmat(this%num_z_bins,this%num_z_bins))
    allocate(ell(this%ell_max))
    allocate(dmix_temp(this%num_z_bins,this%ell_max))

    call this%get_convergence(CMB,Theory, TH_auto, TH_cross, TH_matrix)

    if (do_mock) then
       call wl_mock_maker(mock_root, mock_type, sys_root, sys_dir, TH_auto, TH_cross, TH_matrix, &
                          this%ell_min,this%ell_max,this%num_z_bins, this%num_cross)

       open(42,file='mock_output/'//trim(mock_root)//'_params.dat',status='unknown')
       write(42,'(A,1F10.2)')'ombh2=',CMB%ombh2
       write(42,'(A,1F10.2)')'omch2=',CMB%omch2
       write(42,'(A,1F10.2)')'tau=',CMB%tau
       write(42,'(A,1F10.2)')'ns=',CMB%InitPower(ns_index)
       write(42,'(A,1F10.2)')'log[10^10 A_s]=',CMB%InitPower(As_index)
       write(42,'(A,1F10.2)')'Omega_m=',CMB%omdm+CMB%omb
       write(42,'(A,1F10.2)')'Omega_l=',1-(CMB%omdm+CMB%omb)
       write(42,'(A,1F10.2)')'H_0=',CMB%H0
       close(42)


       write(0,*) 'WL SHEAR MOCK DONE!'
       write(0,*) 'you can find it in your mock_output folder'
       write(0,*) 'So long and thanks for all the fish!'
       stop
    end if


    do ib=1,this%num_z_bins
       TH_matrix(:,ib,ib) = TH_matrix(:,ib,ib) + this%data_noise(:,ib)
    end do


!old method for dl_mix... will be removed. Keeping it just to be safe
!    !Computing determinants
!    dl_mix(:) = 0.d0
!    do l=this%ell_min,this%ell_max
!       ell(l) = l
!       call diatest(this%num_z_bins,TH_matrix(l,:,:),dl_th(l))
!       dl_th(l) = exp(MatrixSym_LogDet(TH_matrix(l,:,:)))
!       do ib=1,this%num_z_bins
!          mixmats(:,:) = TH_matrix(l,:,:)
!          mixmats(:,ib) = this%data_covmat(l,:,ib)
!          dl_mix(l) = dl_mix(l) + exp(MatrixSym_LogDet(mixmats))!FindDet(mixmats(ib,l,:,:),this%num_z_bins)
!       end do
!    end do
!stop

    !computing chi square
    sum = 0.d0
    WLmock_LnLike=0.d0
    do l=this%ell_min,this%ell_max
       trace = 0.d0
       ell(l) = l
       dl_th(l) = exp(MatrixSym_LogDet(TH_matrix(l,:,:)))
       invmat(:,:) = TH_matrix(l,:,:)
       call Matrix_Inverse(invmat)
       call Matrix_Mult(invmat,this%data_covmat(l,:,:),outmat)
       do ib=1,this%num_z_bins
          trace = trace + outmat(ib,ib)
       end do
!       if (ell(l).ge.this%ell_min) then
!          if (abs(log(dl_th(l)/this%dl_obs(l))).lt.1e-5) then !MMcheck: this condition should be made with some sense
!             term = (2*ell(l)+1)*this%fsky*((dl_mix(l)/dl_th(l))-this%num_z_bins)
!          else
!             term = (2*ell(l)+1)*this%fsky*((dl_mix(l)/dl_th(l))+log(dl_th(l)/this%dl_obs(l))-this%num_z_bins)
             term = (2*ell(l)+1)*this%fsky*(trace+log(dl_th(l)/this%dl_obs(l))-this%num_z_bins)
!          end if
!       else
!          term = 0.d0
!       end if
       sum = sum + term
    end do

    WLmock_LnLike = sum/2.d0

    if(feedback>1) write(*,*) trim(this%name)//' WL likelihood = ', WLmock_LnLike

    if (WL_debugging) then
      open(45,file='wldebug_files/mock_auto.dat',status='unknown')
      do l=this%ell_min, this%ell_max
         write(45,35)ell(l),(TH_auto(l,bin)*(ell(l)*ell(l)+1)/6.28,bin=1,this%num_z_bins)
      end do
      close(45)

      open(46,file='wldebug_files/mock_cross.dat',status='unknown')
      do l=this%ell_min, this%ell_max
         write(46,35)ell(l),(TH_cross(l,bin)*(ell(l)*ell(l)+1)/6.28,bin=1,this%num_cross)
      end do
      close(46)

      open(47,file='wldebug_files/determinants.dat',status='unknown')
      do l=this%ell_min, this%ell_max
         write(47,35)ell(l),dl_th(l),this%dl_obs(l),dl_mix(l)
      end do
      close(47)

      open(47,file='wldebug_files/thcovmat.dat',status='unknown')
      do l=this%ell_min, this%ell_max
         do ib=1,this%num_z_bins
            write(47,35) (TH_matrix(l,ib,jb),jb=1,this%num_z_bins)
         end do
      end do
      close(47)

      open(47,file='wldebug_files/obscovmat.dat',status='unknown')
      do l=this%ell_min, this%ell_max
         do ib=1,this%num_z_bins
            write(47,35) (this%data_covmat(l,ib,jb),jb=1,this%num_z_bins)
         end do
      end do
      close(47)



!fai stampare anche la covmat


      open(42,file='wldebug_files/used_params.dat',status='unknown')
      write(42,'(A,1F10.2)')'ombh2=',CMB%ombh2
      write(42,'(A,1F10.2)')'omch2=',CMB%omch2
      write(42,'(A,1F10.2)')'tau=',CMB%tau
      write(42,'(A,1F10.2)')'theta='
      write(42,'(A,1F10.2)')'ns=',CMB%InitPower(ns_index)
      write(42,'(A,1F10.2)')'log[10^10 A_s]=',CMB%InitPower(As_index)
      write(42,'(A,1F10.2)')'Omega_m=',CMB%omdm+CMB%omb
      write(42,'(A,1F10.2)')'Omega_l=',1-(CMB%omdm+CMB%omb)
      write(42,'(A,1F10.2)')'H_0=',CMB%H0
!      write(42,'(A,1F10.2)')'sigma8=',Theory%sigma8
      close(42)

       write(0,*)'-------------------------------------------'
       write(0,*)'           EVERYTHING WENT FINE!           '
       write(0,*)'the single iteration run had no problem... '
       write(0,*)'some files can be found in wldebug_files/  '
       write(0,*)'Now it is time to run this for real!       '
       write(0,*)'-------------------------------------------'

      stop
    end if

    deallocate(TH_auto,TH_cross,ell,TH_matrix)
    deallocate(dl_th,dl_mix,dmix_temp)
    deallocate(mixmats)




    35 FORMAT(100(1X,e15.9))


   end function WLmock_LnLike


   subroutine diatest(n,m0,deter)

   implicit none
   integer, intent(in) :: n
   real(mcp), dimension(n,n), intent(in) :: m0
   integer :: i
   real(mcp), dimension(n,n) :: m1,m2
   real(mcp), dimension(n) :: eig
   real(mcp) :: deter

   m1(:,:)=m0(:,:)

   call diasym(m1,eig,n)
   deter=1.d0
   do i=1,n
      deter = deter*eig(i)
   end do
   
   m2=matmul(transpose(m1),m0)
!   m0=matmul(m2,m1)

   end subroutine diatest


   !Calls the LAPACK diagonalization subroutine DSYEV        
   !input:  a(n,n) = real symmetric matrix to be diagonalized
   !            n  = size of a                               
   !output: a(n,n) = orthonormal eigenvectors of a           
   !        eig(n) = eigenvalues of a in ascending order     

   subroutine diasym(a,eig,n)
   implicit none

   integer n,l,inf
   real*8  a(n,n),eig(n),work(n*(3+n/2))

   l=n*(3+n/2)
   call dsyev('V','U',n,a,n,eig,work,l,inf)

   end subroutine diasym



   real(kind=8) function det(N, mat)

   implicit none
   integer, intent(in) :: N 
   real(kind=8), intent(inout), dimension(:,:) :: mat
   integer(kind=8) :: i, info
   integer, allocatable :: ipiv(:)

   real(kind=8) :: sgn

   allocate(ipiv(N))

   ipiv = 0
   call zgetrf(N, N, mat, N, ipiv, info)

   det = 1.d0!ONE

   do i = 1, N
      det = det*mat(i, i)
   end do

   sgn = 1.d0!ONE

   do i = 1, N
      if(ipiv(i) /= i) then
         sgn = -sgn
      end if
   end do

   det = sgn*det   

   end function det


!Author : Louisda16th a.k.a Ashwith J. Rego
!Description: The subroutine is based on two key points:
!1] A determinant is unaltered when row operations are performed: Hence, using this principle,
!row operations (column operations would work as well) are used
!to convert the matrix into upper traingular form
!2]The determinant of a triangular matrix is obtained by finding the product of the diagonal elements
!
    FUNCTION FindDet(matrix, n)
    IMPLICIT NONE
    REAL*8, DIMENSION(n,n) :: matrix
    INTEGER, INTENT(IN) :: n
    REAL*8 :: m, temp,FindDet
    INTEGER :: i, j, k, l
    LOGICAL :: DetExists = .TRUE.
    l = 1
    !Convert to upper triangular form
    DO k = 1, n-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0
                return
            END IF
        ENDIF
        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
   
    !Calculate determinant by finding product of diagonal elements
    FindDet = l
    DO i = 1, n
        FindDet = FindDet * matrix(i,i)
    END DO
   
END FUNCTION FindDet








    subroutine get_convergence(this,CMB,Theory,TH_auto,TH_cross,Cl)
    use Interpolation
    Class(WLmockLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    type(TCosmoTheoryPK), pointer :: PK
    type(TCubicSpline) :: r_z, dzodr_z, P_z
    type(TCubicSpline),  allocatable :: C_l(:,:), xi1_theta(:,:),xi2_theta(:,:)
    real(mcp) :: h,z,kh,k
    real(mcp), allocatable :: r(:),dzodr(:)
    real(mcp), allocatable :: rbin(:),gbin(:,:)
    real(mcp), allocatable :: ll(:),PP(:)
    real(mcp), allocatable :: integrand(:)
    real(mcp), allocatable :: Cl(:,:,:)
    !Systematics variables----------------------------------------------------
    real(mcp), allocatable                 :: m(:), mulsys(:,:),addsys(:,:,:)
    real(mcp)                              :: term01, term11, term02, term12, long1, long2, lambda, lameta, term01bis, term02bis
    real(mcp)                              :: A1, A2, A3, A4 
    real(mcp), allocatable, dimension(:)   :: cf1,cf2,sysa1,sysa2
    real(mcp), allocatable, dimension(:,:) :: xi01, xi11, xi02, xi12
    real(mcp), allocatable :: mean_bin(:)
    !-------------------------------------------------------------------------
    real(mcp), dimension(this%ell_max,this%num_z_bins),intent(out) :: TH_auto
!    real(mcp), dimension(this%ell_max,this%num_z_bins,this%num_z_bins),intent(out) :: TH_matrix
    real(mcp), dimension(this%ell_max,this%num_cross),intent(out) :: TH_cross
    real(mcp) :: khmin, khmax, lmin, lmax, xmin, xmax, x, lll
    real(mcp) :: i1, i2, lp
    real(mcp) :: a2r
    integer :: i,ib,jb,il,it,iz,nr,nrs,izl,izh,j
    integer :: num_z, ntheta, min_iz, nlmax, aa
    real(mcp), parameter :: pai = 3.1415926535897932384626433832795

    nlmax=this%ell_max


    if (this%use_non_linear) then
        if (this%use_weyl) then
            PK => Theory%NL_MPK_WEYL
        else
            PK => Theory%NL_MPK
        end if
    else
        if (this%use_weyl) then
            PK => Theory%MPK_WEYL
        else
            PK => Theory%MPK
        end if
    end if

    h = CMB%H0/100
    num_z = PK%ny
    khmin = exp(PK%x(1))
    khmax = exp(PK%x(PK%nx))
    a2r = pi/(180._mcp*60._mcp)

    !-----------------------------------------------------------------------
    ! Compute comoving distance r and dz/dr
    !-----------------------------------------------------------------------

    allocate(r(num_z),dzodr(num_z))
    min_iz=1
    do iz=1,num_z
        z = PK%y(iz)
        if (z==0) min_iz=iz+1
        r(iz) = this%Calculator%ComovingRadialDistance(z)
        dzodr(iz) = this%Calculator%Hofz(z)
    end do
    call r_z%Init(PK%y,r,n=num_z)
    call dzodr_z%Init(PK%y,dzodr,n=num_z)

    !-----------------------------------------------------------------------
    ! Compute lensing efficiency
    !-----------------------------------------------------------------------

    allocate(rbin(this%num_z_p),gbin(this%num_z_p,this%num_z_bins))
    rbin=0
    gbin=0
    do iz=1,this%num_z_p
        rbin(iz) = r_z%Value(this%z_p(iz))
    end do
    do ib=1,this%num_z_bins
        do nr=2,this%num_z_p-1
            do nrs=nr+1,this%num_z_p
                gbin(nr,ib)=gbin(nr,ib)+0.5*(dzodr_z%Value(this%z_p(nrs))*this%p(nrs,ib)*(rbin(nrs)-rbin(nr))/rbin(nrs) &
                    + dzodr_z%Value(this%z_p(nrs-1))*this%p(nrs-1,ib)*(rbin(nrs-1)-rbin(nr))/rbin(nrs-1))*(rbin(nrs)-rbin(nrs-1))
            end do
        end do
    end do

!open(666,file='window.dat')
!    do nr=1,this%num_z_p
!       gbin(nr,:) = (3./2.)*(100*h/299790.d0)**2.*(CMB%omdm+CMB%omb)*(1+this%z_p(nr))*rbin(nr)*gbin(nr,:)
!       write(666,*)this%z_p(nr),(3./2.)*(100*h/299790.d0)**2.*(CMB%omdm+CMB%omb)*(1+this%z_p(nr))*rbin(nr)*gbin(nr,1),(3./2.)*(100*h/299790.d0)**2.*(CMB%omdm+CMB%omb)*(1+this%z_p(nr))*rbin(nr)*gbin(nr,3)
!    end do
!close(666)

!    gbin(:,:)=2*gbin(:,:) !devpoint

    !-----------------------------------------------------------------------
    ! Find convergence power spectrum using Limber approximation
    !-----------------------------------------------------------------------

    allocate(ll(nlmax),PP(num_z))
    allocate(integrand(this%num_z_p))
!    allocate(Cl(nlmax,this%num_z_bins,this%num_z_bins))
    Cl = 0
    do il=1,nlmax
        ll(il)=il!1.*exp(dlnl*(il-1._mcp)) !devpoint
        PP=0
        do iz=min_iz,num_z
            k = (ll(il)+0.5)/r(iz)
            kh = k/h ! CAMB wants k/h values
            z = PK%y(iz)
            if ((kh .le. khmin) .or. (kh .ge. khmax)) then
                PP(iz)=0.0d0
            else
                PP(iz)= PK%PowerAt(kh,z)
            end if
        end do

        ! Compute integrand over comoving distance
        call P_z%Init(r,PP,n=num_z)
        do ib=1,this%num_z_bins
            do jb=1,this%num_z_bins
                integrand = 0
                do nr=1,this%num_z_p
                    integrand(nr) = gbin(nr,ib)*gbin(nr,jb)*P_z%Value(rbin(nr))
                    if (.not. this%use_weyl) integrand(nr) = integrand(nr) *(1.0+this%z_p(nr))**2.0
                end do
                do nr=2,this%num_z_p
                    Cl(il,ib,jb)=Cl(il,ib,jb)+0.5d0*(integrand(nr)+integrand(nr-1))*(rbin(nr)-rbin(nr-1))
                end do
            end do
        end do

        if (.not. this%use_weyl) then
            Cl(il,:,:) = (Cl(il,:,:)/h**3.0*9._mcp/4._mcp*(h*1e5_mcp/const_c)**4.0*(CMB%omdm+CMB%omb)**2)!*(h*100)**2./6.
        end if
    end do

!    do il=1,nlmax
!       Cl(il,:,:) = ll(il)*(ll(il)+1)*Cl(il,:,:)/6.28
!    end do


    !Packing C_l
    aa=0
    do ib=1,this%num_z_bins
       TH_auto(:,ib)=Cl(:,ib,ib)
       do jb=2,this%num_z_bins
          if (jb.gt.ib) then
             aa=aa+1
             TH_cross(:,aa)=Cl(:,ib,jb)
          end if
       end do
    end do

    if (aa.ne.this%num_cross) then
       write(0,*) 'ERROR!!!',aa,this%num_cross
    end if



    end subroutine get_convergence



    end module wl_mock
