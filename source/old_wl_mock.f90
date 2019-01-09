    ! Module to analyze mock WL shear datasets
    ! Based on the standard wl.f90 module
    ! It can also be used to generate mock datasets (just set do_mock=T in your .ini file)
    ! Last update: 20/08/2015

    module wl_mock
    use settings
    use CosmologyTypes
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    use WL_MockMaker
    implicit none
    private

    type, extends(TCosmoCalcLikelihood) :: WLmockLikelihood
        real(mcp), allocatable, dimension(:,:) :: wl_invcov,wl_cov
        integer :: num_z_bins,num_cross,ell_max,ell_min
        real(mcp), allocatable, dimension(:) :: z_bins
        integer :: num_z_p ! Source galaxy distribution p(z,bin)
        real(mcp), allocatable, dimension(:,:) :: p,data_auto,data_cross,data_autoerror,data_crosserror
        real(mcp), allocatable, dimension(:) :: z_p,ell
        logical :: use_non_linear ! Whether to use non-linear corrections
        logical :: use_weyl ! Whether to use Weyl potential or matter P(k)
    contains
    procedure :: LogLike => WLmock_LnLike
    procedure :: ReadIni => WLmock_ReadIni
    procedure, private :: get_convergence
    end type WLmockLikelihood


    logical :: do_mock
    integer :: mock_type,bin_type
    logical :: systematics
    character(100) :: mock_root, sys_root, sys_dir

    logical :: WL_debugging=.false.

    ! Integration accuracy parameters
    integer, parameter :: nlmax = 65
    real(mcp), parameter :: dlnl = 0.2d0
    real(mcp) :: dlntheta = 0.25d0
    real(mcp), parameter :: dx = 0.02d0
    real(mcp), parameter :: xstop = 200.0d0
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

    if (Ini%Read_Logical('use_WL_mock',.false.)) then
        nonlinear = Ini%Read_Logical('wl_use_non_linear')!,.true.)
        useweyl = Ini%Read_Logical('wl_use_weyl')!,.true.)

        do_mock = Ini%Read_Logical('do_mock',.false.)
        if (do_mock) then
            mock_type = Ini%Read_Int('mock_type')
            bin_type = Ini%Read_Int('bin_type')
            mock_root = Ini%ReadFileName('mock_root')
        end if
        systematics = Ini%Read_Logical('systematics',.false.)
        if (systematics) then
            sys_dir = Ini%ReadFileName('sys_dir')
            sys_root = Ini%ReadFileName('sys_file')
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
    character(LEN=:), allocatable :: measurements_autofile, window_file, measurements_crossfile
    character(2) :: indbin
    Type(TTextFile) :: F
    Type(WLspec) :: SP
    real(mcp) :: dummy1,dummy2,pnorm
    integer i,iz,it,ib,j,k,nt,izl,izh

    if (.not.do_mock) then

        if (Feedback > 0) write (*,*) 'reading WL data set: '//trim(this%name)

        this%num_z = 100 !Look better at this!!!
        this%max_z = Ini%Read_Double('max_z')
        this%ell_max = Ini%Read_Int('ell_max')
        this%ell_min = Ini%Read_Int('ell_min')
        this%kmax = Ini%Read_Double('kmax')
        this%num_z_bins = Ini%Read_Int('num_z_bins')
        this%num_cross = (this%num_z_bins*(this%num_z_bins+1)/2)-this%num_z_bins
        this%num_z_p = Ini%Read_Int('num_z_p')

        allocate(this%ell(this%ell_max))
        allocate(this%z_bins(this%num_z_bins+1))
        allocate(this%data_auto(this%ell_max,this%num_z_bins))
        allocate(this%data_cross(this%ell_max,this%num_cross))
        allocate(this%data_autoerror(this%ell_max,this%num_z_bins))
        allocate(this%data_crosserror(this%ell_max,this%num_cross))
        allocate(this%z_p(this%num_z_p))
        allocate(this%p(this%num_z_p,this%num_z_bins))


        do i=1,this%num_z_bins+1
           write(indbin,'(I2)') i
           this%z_bins(i)=Ini%Read_Double('zbin('//trim(adjustl(indbin))//')')
        end do


        measurements_autofile  = Ini%ReadFileName('WL_auto_file')
        measurements_crossfile = Ini%ReadFileName('WL_cross_file')

        window_file  = Ini%ReadFileName('window_file')

!        do ib=1,this%num_z_bins
!           call F%Open(FormatString(window_file,ib,allow_unused=.true.))
!           do iz=1,this%num_z_p
!              read (F%unit,*) this%z_p(iz),this%p(iz,ib)
!           end do
!           call F%Close()
!        end do


        open(42,file=trim(window_file))
        do iz=1,this%num_z_p
           read (42,*) this%z_p(iz),this%p(iz,:)
        end do
        close(42)


        call F%Open(measurements_autofile)
        do iz=this%ell_min,this%ell_max
           read (F%unit,*) this%ell(iz), (this%data_auto(iz,ib),ib=1,this%num_z_bins), &
                           (this%data_autoerror(iz,ib),ib=1,this%num_z_bins)
        end do
        call F%Close()

        call F%Open(measurements_crossfile)
        do iz=this%ell_min,this%ell_max
           read (F%unit,*) this%ell(iz), (this%data_cross(iz,ib),ib=1,this%num_cross), &
                           (this%data_crosserror(iz,ib),ib=1,this%num_cross)
        end do
        call F%Close()
        do iz=this%ell_min,this%ell_max
           this%data_auto(iz,:)=this%data_auto(iz,:)*2*const_pi/(this%ell(iz)*(this%ell(iz)+1.))
           this%data_cross(iz,:)=this%data_cross(iz,:)*2*const_pi/(this%ell(iz)*(this%ell(iz)+1.))
           this%data_autoerror(iz,:)=this%data_autoerror(iz,:)*2*const_pi/(this%ell(iz)*(this%ell(iz)+1.))
           this%data_crosserror(iz,:)=this%data_crosserror(iz,:)*2*const_pi/(this%ell(iz)*(this%ell(iz)+1.))
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

        open(42,file='mock_output/'//trim(mock_root)//'_window.dat')
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

    end subroutine WLmock_ReadIni

    function WLmock_LnLike(this, CMB, Theory, DataParams)
    use MatrixUtils
    Class(WLmockLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) WLmock_LnLike
    real(mcp), allocatable, dimension(:,:) :: TH_auto,TH_cross,TH_auto2,TH_cross2
    real(mcp), allocatable, dimension(:) :: ell
    integer :: l,bin

    allocate(TH_auto(this%ell_max,this%num_z_bins))
    allocate(TH_cross(this%ell_max,this%num_cross))
    allocate(ell(this%ell_max))

    call this%get_convergence(CMB,Theory, TH_auto, TH_cross)

   if (do_mock) then
      call wl_mock_maker(mock_root, mock_type, sys_root, sys_dir, TH_auto, TH_cross, this%ell_max,this%num_z_bins, this%num_cross, systematics)

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


   WLmock_LnLike=0.d0

   do l=1, this%ell_max
      ell(l)=l
      if (ell(l).ge.this%ell_min) then  
          do bin=1,this%num_z_bins
             WLmock_LnLike = WLmock_LnLike + ((TH_auto(l,bin)-this%data_auto(l,bin))**2)/(this%data_autoerror(l,bin)**2)
          end do
          if (this%num_z_bins.gt.1) then
             do bin=1,this%num_cross
                WLmock_LnLike = WLmock_LnLike + ((TH_cross(l,bin)-this%data_cross(l,bin))**2)/(this%data_crosserror(l,bin)**2)
             end do
          end if
      end if
   end do
   WLmock_LnLike = WLmock_LnLike/2.d0

   if (WL_debugging) then
      open(45,file='out_spectra/mock_auto.dat',status='unknown')
      do l=1, this%ell_max
         write(45,35)ell(l),(TH_auto(l,bin)*(ell(l)*ell(l)+1)/6.28,bin=1,this%num_z_bins)
      end do
      close(45)

      open(46,file='out_spectra/mock_cross.dat',status='unknown')
      do l=1, this%ell_max
         write(46,35)ell(l),(TH_cross(l,bin)*(ell(l)*ell(l)+1)/6.28,bin=1,this%num_cross)
      end do
      close(46)

      open(42,file='out_spectra/mock_params.dat',status='unknown')
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

      write(0,*) 'EVERYTHING WENT FINE!'
   end if

   deallocate(TH_auto,TH_cross,ell)

   35 FORMAT(100(1X,e15.9))
   if(feedback>1) write(*,*) trim(this%name)//' WL likelihood = ', WLmock_LnLike

   if (WL_debugging) then
      stop
   end if

    end function WLmock_LnLike


    subroutine get_convergence(this,CMB,Theory,TH_auto,TH_cross)
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
    real(mcp), allocatable :: m(:), mulsys(:,:),addsys(:,:,:)
    real(mcp), allocatable :: sf01(:), sf02(:), sf11(:), sf12(:), bF01(:,:), bF11(:,:), bF02(:,:), bF12(:,:)
    real(mcp), allocatable :: lambda(:),lambdaP(:),cf1(:),cf2(:),bP01(:),bP02(:),bP11(:),bP12(:)
    real(mcp), allocatable :: mean_bin(:)
    real(mcp), dimension(this%ell_max,this%num_z_bins),intent(out) :: TH_auto
    real(mcp), dimension(this%ell_max,this%num_cross),intent(out) :: TH_cross
    real(mcp) :: khmin, khmax, lmin, lmax, xmin, xmax, x, lll
    real(mcp) :: i1, i2, lp
    real(mcp) :: a2r
    integer :: i,ib,jb,il,it,iz,nr,nrs,izl,izh,j
    integer :: num_z, ntheta, min_iz, nlmax, aa

    nlmax=this%ell_max
write(0,*)'entro'

    if (this%use_non_linear) then
        if (this%use_weyl) then
            PK => Theory%NL_MPK_WEYL
        else
write(0,*)'test?'
            PK => Theory%NL_MPK
        end if
    else
        if (this%use_weyl) then
            PK => Theory%MPK_WEYL
        else
write(0,*)'test?'
            PK => Theory%MPK
        end if
    end if

    h = CMB%H0/100
    num_z = PK%ny
write(0,*)'passed'
stop
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

!    gbin(:,:)=2*gbin(:,:) !devpoint

    !-----------------------------------------------------------------------
    ! Find convergence power spectrum using Limber approximation
    !-----------------------------------------------------------------------

    allocate(ll(nlmax),PP(num_z))
    allocate(integrand(this%num_z_p))
    allocate(Cl(nlmax,this%num_z_bins,this%num_z_bins))
    Cl = 0

    do il=1,nlmax
        ll(il)=il!1.*exp(dlnl*(il-1._mcp)) !devpoint
        PP=0
        do iz=min_iz,num_z
            k = ll(il)/r(iz)
            kh = k/h ! CAMB wants k/h values
            z = PK%y(iz)
            if ((kh .le. khmin) .or. (kh .ge. khmax)) then
                PP(iz)=0.0d0
            else
                PP(iz)= PK%PowerAt(kh,z)
!if (iz.eq.min_iz) write(42,*) kh,PK%PowerAt(kh,z),z
                ! Testing
!if ((z.gt.0.99).and.(z.lt.1.03)) then
!                write(42,'(10E15.5)') k,PK%PowerAt(kh,z)!,9.0/(8.0*pi**2.0)&
!end if
                !*PK%PowerAt(kh,z)/(h**3.0)*(h*1e5_mcp/const_c)**4.0*(CMB%omdm+CMB%omb)**2*(1+z)**2.0
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
!            Cl(il,:,:) = Cl(il,:,:)/h**3.0*9._mcp/16._mcp*(h*1e5_mcp/const_c)**4.0*(CMB%omdm+CMB%omb)**2
!            Cl(il,:,:) = Cl(il,:,:)*(CMB%h/2997.9)**4*(9./4.)*(CMB%omb+CMB%omdm)**2
        end if
    end do

!    do il=1,nlmax
!       Cl(il,:,:) = ll(il)*(ll(il)+1)*Cl(il,:,:)/6.28
!    end do


!ADDING SYSTEMATICS--------------------------------------------------------------
    allocate(m(this%num_z_bins),mulsys(this%num_z_bins,this%num_z_bins))
    allocate(addsys(nlmax,this%num_z_bins,this%num_z_bins))
    allocate(sf01(this%num_z_bins),sf11(this%num_z_bins),sf02(this%num_z_bins),sf12(this%num_z_bins))
    allocate(bF01(this%num_z_bins,this%num_z_bins),bF11(this%num_z_bins,this%num_z_bins),bF02(this%num_z_bins,this%num_z_bins),bF12(this%num_z_bins,this%num_z_bins))
    allocate(cf1(this%num_z_bins),cf2(this%num_z_bins))
    allocate(lambda(nlmax),lambdaP(nlmax),bP01(nlmax),bP02(nlmax),bP11(nlmax),bP12(nlmax))
    allocate(mean_bin(this%num_z_bins))

    do ib=2,this%num_z_bins+1
       mean_bin(ib-1)=(this%z_bins(ib)+this%z_bins(ib-1))/2
    end do


    !multiplicative sys
    do ib=1,this%num_z_bins
       m(ib)=CMB%m0+CMB%m1*mean_bin(ib)/(1+mean_bin(ib))
    end do

    do ib=1,this%num_z_bins
       do jb=1,this%num_z_bins
          mulsys(ib,jb)=m(ib)+m(jb)+m(ib)*m(jb)
       end do
    end do


    !additive sys
    do ib=1,this%num_z_bins
       cf1(ib)=CMB%c01+CMB%c11*mean_bin(ib)/(1+mean_bin(ib))
       cf2(ib)=CMB%c01+CMB%c12*mean_bin(ib)/(1+mean_bin(ib))
       sf01(ib)=CMB%f0_01+CMB%f1_01*mean_bin(ib)/(1+mean_bin(ib))
       sf11(ib)=CMB%f0_11+CMB%f1_11*mean_bin(ib)/(1+mean_bin(ib))
       sf02(ib)=CMB%f0_02+CMB%f1_02*mean_bin(ib)/(1+mean_bin(ib))
       sf12(ib)=CMB%f0_12+CMB%f1_12*mean_bin(ib)/(1+mean_bin(ib))
    end do

    do ib=1,this%num_z_bins
       do jb=1,this%num_z_bins
          bF01(ib,jb)=(2*sf01(ib)*(sf01(jb)+sf11(jb)*CMB%etas)+CMB%etas*sf11(ib)*(2*sf01(jb)+3*sf11(jb))*CMB%etas)*(CMB%etas**2./8)
          bF02(ib,jb)=(2*sf02(ib)*(sf02(jb)+sf12(jb)*CMB%etas)+CMB%etas*sf12(ib)*(2*sf02(jb)+3*sf12(jb))*CMB%etas)*(CMB%etas**2./8)
          bF11(ib,jb)=sf01(ib)*sf11(jb)+sf11(ib)*sf01(jb)*CMB%etas
          bF12(ib,jb)=sf02(ib)*sf12(jb)+sf12(ib)*sf02(jb)*CMB%etas
       end do
    end do

    do il=1,nlmax
       lambda(il) = ll(il)*CMB%thetas*CMB%etas
       lambdaP(il) = sqrt(1+lambda(il)**2.)
       bP01(il)=-2*(1+lambda(il)**2.)**2.*((lambdaP(il)-4)*lambda(il)**2.+6*(lambdaP(il)-1))
       bP11(il)=-4*(1+lambda(il)**2.)**2.*CMB%etas*(2*lambda(il)**4.+(-6*lambdaP(il)+9)*lambda(il)**2.-6*(lambdaP(il)-1))
       bP02(il)=(1+lambda(il)**2.)*((2*lambdaP(il)-7)*lambda(il)**4.+2*lambda(il)**2.*(7*lambdaP(il)-10)+12*(lambdaP(il)-1))
       bP12(il)=CMB%etas*(7*lambda(il)**6.+lambda(il)**4.*(-24*lambdaP(il)+46)+lambda(il)**2.*(-48*lambdaP(il)+60)-24*(lambdaP(il)-1))
    end do

    do il=1,nlmax
       do ib=1,this%num_z_bins
          do jb=1,this%num_z_bins
             addsys(il,ib,jb)=(4*3.14*CMB%thetas**4.*CMB%etas**2./(lambda(il)**4.*(lambda(il)**2.+1)**(5./2.)))* &
                              (cf1(ib)*cf1(jb)*(bF01(ib,jb)*bP01(il)+bF11(ib,jb)*bP11(il))+ &
                              cf2(ib)*cf2(jb)*(bF02(ib,jb)*bP02(il)+bF12(ib,jb)*bP12(il)))
          end do
       end do
    end do

    !adding systematics to the spectra

    do il=1,nlmax
       Cl(il,:,:) = (1+mulsys(:,:))*Cl(il,:,:)+addsys(il,:,:)
    end do
!--------------------------------------------------------------------------------------

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


    deallocate(m,mulsys)
    deallocate(addsys)
    deallocate(sf01,sf11,sf02,sf12)
    deallocate(bF01,bF11,bF02,bF12)
    deallocate(cf1,cf2)
    deallocate(lambda,lambdaP,bP01,bP02,bP11,bP12)
    deallocate(mean_bin)

    end subroutine get_convergence

    end module wl_mock


