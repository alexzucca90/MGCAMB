module redshiftdrift
    use MatrixUtils
    use settings
    use CosmologyTypes
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    use RD_mockmaker
    use IniObjects
    implicit none

    private

    type, extends(TCosmoCalcLikelihood) :: RDLikelihood

    integer :: num_red_rd
    real(mcp), allocatable, dimension(:) :: rd_z, rd_obs, rd_err,rd_deltat
!    character(LEN=100) :: rd_file

    contains
    procedure :: LogLike => RD_LnLike
    procedure :: ReadIni => RD_ReadIni
    end type RDLikelihood


    type, extends(RDLikelihood) :: ELTLikelihood
        real(mcp) :: temp
    contains
    procedure :: LogLike => RD_ELT_Lnlike
    end type ELTLikelihood

    type, extends(RDLikelihood) :: SKALikelihood
        real(mcp) :: temp
    contains
    procedure :: LogLike => RD_SKA_Lnlike
    end type SKALikelihood

    type, extends(RDLikelihood) :: CHIMELikelihood
        real(mcp) :: temp
    contains
    procedure :: LogLike => RD_CHIME_Lnlike
    end type CHIMELikelihood


    logical :: do_mock
    integer :: mock_type
    character(100) :: mock_root

    public RDLikelihood, RDLikelihood_Add
    contains

    subroutine RDLikelihood_Add(LikeList, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    class(RDLikelihood), pointer :: this
    Type(TSettingIni) :: DataSets, OverrideSettings
    integer i

    if (.not. Ini%Read_Logical('use_RD',.false.)) return

    call Ini%TagValuesForName('RD_dataset', DataSets, filename=.true.)
    if (DataSets%Count==0) call MpiStop('Use_RD but no RD_dataset[NAMETAG] defined')

    do_mock = Ini%Read_Logical('do_mock',.false.)
    if (do_mock) then
        call Ini%SettingValuesForTagName('RD_dataset',DataSets%Name(1),OverrideSettings)
        mock_type = Ini%Read_Int('mock_type')
        mock_root = Ini%ReadFileName('mock_root')
        write(0,*) 'Doing redshift drift mock'
        allocate(RDLikelihood::this)
        call this%ReadDatasetFile(Datasets%Value(1),OverrideSettings)
        this%needs_background_functions = .true.
        call LikeList%Add(this)
    else
        write(0,*)'test',DataSets%Count

        do i= 1, DataSets%Count
           call Ini%SettingValuesForTagName('RD_dataset',DataSets%Name(i),OverrideSettings)
           if (Datasets%Name(i)=='ELT') then
               allocate(ELTLikelihood::this)
           else if (Datasets%Name(i)=='SKA') then
               allocate(SKALikelihood::this)
           else if (Datasets%Name(i)=='CHIME') then
               allocate(CHIMELikelihood::this)
           else
               allocate(RDLikelihood::this)
           end if
           call this%ReadDatasetFile(Datasets%Value(i),OverrideSettings)
           this%tag = Datasets%Name(i)
write(0,*)'test',Datasets%Name(i)
           this%LikelihoodType = 'RD'
           this%needs_background_functions = .true.
           call LikeList%Add(this)
        end do

!    allocate(this)
!    call this%ReadDatasetFile(Ini%ReadFileName('RD_dataset'))
!    this%LikelihoodType = 'RD'
!    this%needs_background_functions = .true.
!        this%num_red_rd=Ini%Read_Int('num_red_rd',0)
!        this%RD_debugging=Ini%Read_Logical('RD_debugging',.false.)
!        this%rd_file = Ini%ReadFileName('RD_file')
!write(0,*)this%num_red_rd,this%RD_debugging
!    call LikeList%Add(this)
        if (Feedback>1) write(*,*) 'read RD data sets'
    end if

    end subroutine RDLikelihood_Add

    subroutine RD_ReadIni(this, Ini)
    class(RDLikelihood) this
    class(TSettingIni) :: Ini
    character(LEN=100) :: rd_file
    integer i,iopb
    Type(TTextFile) :: F
    Type(RDspec) :: SP


    if (do_mock) then
       call GetSpecs(mock_type,mock_root,SP)
       this%num_red_rd=SP%num_red
       allocate(this%rd_z(this%num_red_rd))
       allocate(this%rd_deltat(this%num_red_rd))
       do i=1,this%num_red_rd
          this%rd_z(i)=SP%zbins(i)
          this%rd_deltat(i)=SP%deltat_years
       end do
       deallocate(SP%zbins)
    else
       if (Feedback > 0) write (*,*) 'reading RD data set: '//trim(this%name)
       this%num_red_rd=Ini%Read_Int('num_red_rd',0)

 
       allocate(this%rd_z(this%num_red_rd))
       allocate(this%rd_obs(this%num_red_rd))
       allocate(this%rd_err(this%num_red_rd))
       allocate(this%rd_deltat(this%num_red_rd))

       rd_file = Ini%ReadFileName('RD_file')
       call F%Open(trim(adjustl(rd_file)))
       do i=1,this%num_red_rd
            read (F%unit,*, iostat=iopb) this%rd_z(i),this%rd_obs(i),this%rd_err(i),this%rd_deltat(i)
       end do
       call F%Close()

    end if

    end subroutine RD_ReadIni


    function RD_LnLike(this, CMB, Theory, DataParams)
    Class(RDLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory !useless
    real(mcp), allocatable :: hofz
    real(mcp) :: DataParams(:) !useless
    integer j,k,i,num_z
    real(mcp) RD_LnLike,z
    real(mcp), allocatable :: RD_theory(:)
    real(mcp) :: deltat,acca,omega_m,E
    real(mcp), parameter :: conv_time=3.15e7, conv_mpc=1/(3.09e19), MPC_in_sec = 1.029272d14, clight=2.998e10 !cm/s
    character(100) derphi


    RD_LnLike=0
    allocate(RD_theory(this%num_red_rd))
    deltat=this%rd_deltat(1)*conv_time
    acca=CMB%H0*conv_mpc
    omega_m=(CMB%ombh2+CMB%omch2)/((CMB%H0/100)**2)

    do j=1, this%num_red_rd
       z=this%rd_z(j)
       E= this%Calculator%Hofz_Hunit(z)/CMB%H0
       RD_theory(j)=clight*acca*deltat*(1.-(E/(1+this%rd_z(j))))
    end do

    if (do_mock) then
       write(0,*) 'redshift drift computed for target'
       call rd_mock_maker(mock_root, mock_type,RD_theory)

       open(42,file='mock_output/'//trim(mock_root)//'_params.dat',status='unknown')
       write(42,'(A,1F10.2)')'ombh2=',CMB%ombh2
       write(42,'(A,1F10.2)')'omch2=',CMB%omch2
       write(42,'(A,1F10.2)')'tau=',CMB%tau
       write(42,'(A,1F10.2)')'ns=',CMB%InitPower(ns_index)
       write(42,'(A,1F10.2)')'A_s=',CMB%InitPower(As_index)
       write(42,'(A,1F10.2)')'Omega_m=',CMB%omdm+CMB%omb
       write(42,'(A,1F10.2)')'Omega_l=',1-(CMB%omdm+CMB%omb)
       write(42,'(A,1F10.2)')'H_0=',CMB%H0
       close(42)


       write(0,*) 'REDSHIFT DRIFT MOCK DONE!'
       write(0,*) 'you can find it in your mock_output folder'
       write(0,*) 'So long and thanks for all the fish!'

       stop
    else

       do j=1, this%num_red_rd
          RD_LnLike = RD_LnLike + ((RD_theory(j)-this%rd_obs(j))**2)/(this%rd_err(j)**2)
       end do
       RD_LnLike = RD_LnLike/2.d0


       if(feedback>1) write(*,*) trim(this%name)//' RD likelihood = ', RD_LnLike

    end if

    deallocate(RD_theory)


35 FORMAT(100(1X,e15.9))

    end function RD_LnLike

    function RD_ELT_LnLike(this, CMB, Theory, DataParams)
    Class(ELTLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory !useless
    real(mcp), allocatable :: hofz
    real(mcp) :: DataParams(:) !useless
    integer j,k,i,num_z
    real(mcp) RD_ELT_LnLike,z
    real(mcp), allocatable :: RD_theory(:)
    real(mcp) :: deltat,acca,omega_m,E
    real(mcp), parameter :: conv_time=3.15e7, conv_mpc=1/(3.09e19), MPC_in_sec = 1.029272d14, clight=2.99792458e10 !cm/s
    character(100) derphi

    real(mcp) :: firstder,seconder
    real(mcp),allocatable :: decpar
!    real(mcp), parameter :: MPC_in_sec = 1.029272d14 in SI units
    real(mcp), parameter :: c = 2.99792458e8


    RD_ELT_LnLike=0
    allocate(RD_theory(this%num_red_rd))

    deltat=this%rd_deltat(1)*conv_time
    acca=CMB%H0*conv_mpc
    omega_m=0.3!(CMB%ombh2+CMB%omch2)/((CMB%H0/100)**2)


    do j=1, this%num_red_rd
       z=this%rd_z(j)
       E= this%Calculator%Hofz_Hunit(z)/CMB%H0
       RD_theory(j)=clight*acca*deltat*(1.-(E/(1+this%rd_z(j))))
    end do

    do j=1, this%num_red_rd
       RD_ELT_LnLike = RD_ELT_LnLike + ((RD_theory(j)-this%rd_obs(j))**2)/(this%rd_err(j)**2)
    end do
    RD_ELT_LnLike = RD_ELT_LnLike/2.d0

    deallocate(RD_theory)

35 FORMAT(100(1X,e15.9))
    if(feedback>1) write(*,*) trim(this%name)//' ELT likelihood = ', RD_ELT_LnLike
    end function RD_ELT_LnLike


    function RD_SKA_LnLike(this, CMB, Theory, DataParams)
    Class(SKALikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory !useless
    real(mcp), allocatable :: hofz
    real(mcp) :: DataParams(:) !useless
    integer j,k,i,num_z
    real(mcp) RD_SKA_LnLike,z
    real(mcp), allocatable :: RD_theory(:)
    real(mcp) :: deltat,acca,omega_m,E
    real(mcp), parameter :: conv_time=3.15e7, conv_mpc=1/(3.09e19), MPC_in_sec = 1.029272d14, clight=2.99792458e10 !cm/s
    character(100) derphi

    RD_SKA_LnLike=0
    allocate(RD_theory(this%num_red_rd))

    deltat=this%rd_deltat(1)*conv_time
    acca=CMB%H0*conv_mpc
    omega_m=(CMB%ombh2+CMB%omch2)/((CMB%H0/100)**2)


    do j=1, this%num_red_rd
       z=this%rd_z(j)
       E= this%Calculator%Hofz_Hunit(z)/CMB%H0 
       RD_theory(j)=clight*acca*deltat*(1.-(E/(1+this%rd_z(j))))
    end do

    do j=1, this%num_red_rd
       RD_SKA_LnLike = RD_SKA_LnLike + ((RD_theory(j)-this%rd_obs(j))**2)/(this%rd_err(j)**2)
    end do
    RD_SKA_LnLike = RD_SKA_LnLike/2.d0

    deallocate(RD_theory)

35 FORMAT(100(1X,e15.9))
    if(feedback>1) write(*,*) trim(this%name)//' SKA likelihood = ', RD_SKA_LnLike
    end function RD_SKA_LnLike



    function RD_CHIME_LnLike(this, CMB, Theory, DataParams)
    Class(CHIMELikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory !useless
    real(mcp), allocatable :: hofz
    real(mcp) :: DataParams(:) !useless
    integer j,k,i,num_z
    real(mcp) RD_CHIME_LnLike,z
    real(mcp), allocatable :: RD_theory(:)
    real(mcp) :: deltat,acca,omega_m,E
    real(mcp), parameter :: conv_time=3.15e7, conv_mpc=1/(3.09e19), MPC_in_sec = 1.029272d14, clight=2.99792458e10 !cm/s
    character(100) derphi

    RD_CHIME_LnLike=0
    allocate(RD_theory(this%num_red_rd))

    deltat=this%rd_deltat(1)*conv_time
    acca=CMB%H0*conv_mpc
    omega_m=(CMB%ombh2+CMB%omch2)/((CMB%H0/100)**2)


    do j=1, this%num_red_rd
       z=this%rd_z(j)
       E= this%Calculator%Hofz_Hunit(z)/CMB%H0 
       RD_theory(j)=clight*acca*deltat*(1.-(E/(1+this%rd_z(j))))
    end do

    do j=1, this%num_red_rd
       RD_CHIME_LnLike = RD_CHIME_LnLike + ((RD_theory(j)-this%rd_obs(j))**2)/(this%rd_err(j)**2)
    end do
    RD_CHIME_LnLike = RD_CHIME_LnLike/2.d0

    deallocate(RD_theory)

35 FORMAT(100(1X,e15.9))
    if(feedback>1) write(*,*) trim(this%name)//' CHIME likelihood = ', RD_CHIME_LnLike
    end function RD_CHIME_LnLike




      SUBROUTINE splint(xa,ya,y2a,n,x,y)
      INTEGER n
      real xa(n),y2a(n),ya(n)
      real*8 :: x,y
      INTEGER k,khi,klo
      REAL a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      END SUBROUTINE splint

      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
      INTEGER n,NMAX
      real x(n),y(n),y2(n)
      real*8 :: yp1,ypn
      PARAMETER (NMAX=200000)
      INTEGER i,k
      REAL p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      END SUBROUTINE spline







end module redshiftdrift
