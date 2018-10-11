!>


module MGCAMB_MuGamma

    use Precision
    use IniFile
    use AMLutils
    use MGCAMB_cache
    use MGCAMB_def

    implicit none

    private

    public MGCAMB_model_MuGamma

    type, extends ( MGCAMB_abstract_model ) :: MGCAMB_model_MuGamma

        !> some flags
        integer :: MGCAMBMuGammaModel

        class( parametrized_function_1D ), allocatable :: MGCAMBMuGamma_Mu     !< the mu function
        class( parametrized_function_1D ), allocatable :: MGCAMBMuGamma_Gamma  !< the gamma function

    contains

        ! initialization of the model:
        procedure :: read_model_selection            => MGCAMBMuGammaReadModelSelectionFromFile  !< subroutine that reads the parameters of the model from file
        procedure :: allocate_model_selection        => MGCAMBMuGammaAllocateModelSelection      !< subroutine that allocates the model selection.
        procedure :: init_model_parameters           => MGCAMBMuGammaInitModelParameters         !< subroutine that initializes the model parameters based on the values found in an input array.
        procedure :: init_model_parameters_from_file => MGCAMBMuGammaInitModelParametersFromFile !< subroutine that reads the parameters of the model from file.
        ! utility functions:
        procedure :: compute_param_number  => MGCAMBMuGammaComputeParametersNumber    !< subroutine that computes the number of parameters of the model.
        procedure :: feedback              => MGCAMBMuGammaFeedback                   !< subroutine that prints on the screen feedback information about the model.
        procedure :: parameter_names       => MGCAMBMuGammaParameterNames             !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_names_latex => MGCAMBMuGammaParameterNamesLatex        !< subroutine that returns the i-th parameter name of the model.
        procedure :: parameter_values      => MGCAMBMuGammaParameterValues            !< subroutine that returns the i-th parameter value.
        ! CAMB related procedures:
        procedure :: compute_background_MG_functions  => MGCAMBMuGammaBackgroundEFTFunctions   !< subroutine that computes the value of the background EFT functions at a given time.
        procedure ::

    end type MGCAMB_model_MuGamma

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the parameters of the model from file.
    subroutine MGCAMBMuGammaReadModelSelectionFromFile( self, Ini )

        implicit none

        class( MGCAMB_model_MuGamma ) :: self   !< the base class
        type(TIniFile)                :: Ini    !< Input ini file

        ! read model selection flags:
        self%MGCAMBMuGammaModel = Ini_Read_Int_File( Ini, 'model'  , 0 )

    end subroutine MGCAMBMuGammaReadModelSelectionFromFile


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the model selection.
    subroutine MGCAMBMuGammaAllocateModelSelection( self )

        implicit none

        class( MGCAMB_model_MuGamma )                       :: self              !< the base class

    ! allocate Mu and Gamma according to the model
    if ( allocated(self%PureEFTOmega) ) deallocate(self%PureEFTOmega)
    select case ( self%PureEFTmodelOmega )
    case(0)
    allocate( zero_parametrization_1D::self%PureEFTOmega )
    case(1)
    allocate( constant_parametrization_1D::self%PureEFTOmega )
    case(2)
    allocate( linear_parametrization_1D::self%PureEFTOmega )
    case(3)
    allocate( power_law_parametrization_1D::self%PureEFTOmega )
    call self%PureEFTOmega%set_param_names( ['EFTOmega0  ', 'EFTOmegaExp'], ['\Omega_0^{\rm EFT}', 'n^{\rm EFT}       '] )
    case(4)
    allocate( exponential_parametrization_1D::self%PureEFTOmega )
    call self%PureEFTOmega%set_param_names( ['EFTOmega0  ', 'EFTOmegaExp'], ['\Omega_0^{\rm EFT}', 'n^{\rm EFT}       '] )
    !--------
    ! GBD MOD: start
    !> adding a new case here
    case(5)
    allocate(interpolated_function_1D::self%PureEFTOmega)
    call self%PureEFTOmega%set_param_names(['Omega_filename'])
    ! GBD MOD: end
    !--------
    case default
    write(*,'(a,I3)') 'No model corresponding to PureEFTmodelOmega =', self%PureEFTmodelOmega
    write(*,'(a)')    'Please select an appropriate model.'
    end select
    ! allocate wDE:
    if ( allocated(self%PureEFTwDE) ) deallocate(self%PureEFTwDE)
    select case ( self%EFTwDE )
    case(0)
    allocate( wDE_LCDM_parametrization_1D::self%PureEFTwDE )
    case(1)
    allocate( constant_parametrization_1D::self%PureEFTwDE )
    case(2)
    allocate( CPL_parametrization_1D::self%PureEFTwDE )
    call self%PureEFTwDE%set_param_names( ['EFTw0', 'EFTwa'], ['w_0', 'w_a'] )
    case(3)
    allocate( JBP_parametrization_1D::self%PureEFTwDE )
    call self%PureEFTwDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTwn'], [ 'w_0', 'w_a', 'n  ' ] )
    case(4)
    allocate( turning_point_parametrization_1D::self%PureEFTwDE )
    call self%PureEFTwDE%set_param_names( ['EFTw0 ', 'EFTwa ', 'EFTwat'], ['w_0', 'w_a', 'a_t'] )
    case(5)
    allocate( taylor_parametrization_1D::self%PureEFTwDE )
    call self%PureEFTwDE%set_param_names( ['EFTw0', 'EFTwa', 'EFTw2', 'EFTw3'], ['w_0', 'w_a', 'w_2', 'w_3'] )
    !------
    ! GBD MOD: start
    !> adding a new case here: interpolated function
    case(6)
    allocate( interpolated_function_1D::self%PureEFTwDE )
    call self%PureEFTwDE%set_param_names(['wDE_filename  '])
    ! Nothing else to add
    ! GBD MOD: end
    !------
    case default
    write(*,'(a,I3)') 'No model corresponding to EFTwDE =', self%EFTwDE
    write(*,'(a)')    'Please select an appropriate model.'
    end select
    ! allocate Gamma1:
    if ( allocated(self%PureEFTGamma1) ) deallocate(self%PureEFTGamma1)
    select case ( self%PureEFTmodelGamma1 )
    case(0)
    allocate( zero_parametrization_1D::self%PureEFTGamma1 )
    case(1)
    allocate( constant_parametrization_1D::self%PureEFTGamma1 )
    case(2)
    allocate( linear_parametrization_1D::self%PureEFTGamma1 )
    case(3)
    allocate( power_law_parametrization_1D::self%PureEFTGamma1 )
    call self%PureEFTGamma1%set_param_names( ['EFTGamma10  ', 'EFTGamma1Exp'], ['\gamma_0^{(1) {\rm EFT}}        ', '\gamma_{\rm exp}^{(1) {\rm EFT}}'] )
    case(4)
    allocate( exponential_parametrization_1D::self%PureEFTGamma1 )
    call self%PureEFTGamma1%set_param_names( ['EFTGamma10  ', 'EFTGamma1Exp'], ['\gamma_0^{(1) {\rm EFT}}        ', '\gamma_{\rm exp}^{(1) {\rm EFT}}'] )
    case default
    write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma1 =', self%PureEFTmodelGamma1
    write(*,'(a)')    'Please select an appropriate model.'
    end select
    ! allocate Gamma2:
    if ( allocated(self%PureEFTGamma2) ) deallocate(self%PureEFTGamma2)
    select case ( self%PureEFTmodelGamma2 )
    case(0)
    allocate( zero_parametrization_1D::self%PureEFTGamma2 )
    case(1)
    allocate( constant_parametrization_1D::self%PureEFTGamma2 )
    case(2)
    allocate( linear_parametrization_1D::self%PureEFTGamma2 )
    case(3)
    allocate( power_law_parametrization_1D::self%PureEFTGamma2 )
    call self%PureEFTGamma2%set_param_names( ['EFTGamma20  ', 'EFTGamma2Exp'], ['\gamma_0^{(2) {\rm EFT}}        ', '\gamma_{\rm exp}^{(2) {\rm EFT}}'] )
    case(4)
    allocate( exponential_parametrization_1D::self%PureEFTGamma2 )
    call self%PureEFTGamma2%set_param_names( ['EFTGamma20  ', 'EFTGamma2Exp'], ['\gamma_0^{(2) {\rm EFT}}        ', '\gamma_{\rm exp}^{(2) {\rm EFT}}'] )
    case default
    write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma2 =', self%PureEFTmodelGamma2
    write(*,'(a)')    'Please select an appropriate model.'
    end select
    ! allocate Gamma3:
    if ( allocated(self%PureEFTGamma3) ) deallocate(self%PureEFTGamma3)
    select case ( self%PureEFTmodelGamma3 )
    case(0)
    allocate( zero_parametrization_1D::self%PureEFTGamma3 )
    case(1)
    allocate( constant_parametrization_1D::self%PureEFTGamma3 )
    case(2)
    allocate( linear_parametrization_1D::self%PureEFTGamma3 )
    case(3)
    allocate( power_law_parametrization_1D::self%PureEFTGamma3 )
    call self%PureEFTGamma3%set_param_names( ['EFTGamma30  ', 'EFTGamma3Exp'], ['\gamma_0^{(3) {\rm EFT}}        ', '\gamma_{\rm exp}^{(3) {\rm EFT}}'] )
    case(4)
    allocate( exponential_parametrization_1D::self%PureEFTGamma3 )
    call self%PureEFTGamma3%set_param_names( ['EFTGamma30  ', 'EFTGamma3Exp'], ['\gamma_0^{(3) {\rm EFT}}        ', '\gamma_{\rm exp}^{(3) {\rm EFT}}'] )
    case default
    write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma3 =', self%PureEFTmodelGamma3
    write(*,'(a)')    'Please select an appropriate model.'
    end select
    ! allocate the other functions only if not Horndeski:
    if ( .not. self%PureEFTHorndeski ) then
    ! allocate Gamma4:
    if ( allocated(self%PureEFTGamma4) ) deallocate(self%PureEFTGamma4)
    select case ( self%PureEFTmodelGamma4 )
    case(0)
    allocate( zero_parametrization_1D::self%PureEFTGamma4 )
    case(1)
    allocate( constant_parametrization_1D::self%PureEFTGamma4 )
    case(2)
    allocate( linear_parametrization_1D::self%PureEFTGamma4 )
    case(3)
    allocate( power_law_parametrization_1D::self%PureEFTGamma4 )
    call self%PureEFTGamma4%set_param_names( ['EFTGamma40  ', 'EFTGamma4Exp'], ['\gamma_0^{(4) {\rm EFT}}        ', '\gamma_{\rm exp}^{(4) {\rm EFT}}'] )
    case(4)
    allocate( exponential_parametrization_1D::self%PureEFTGamma4 )
    call self%PureEFTGamma4%set_param_names( ['EFTGamma40  ', 'EFTGamma4Exp'], ['\gamma_0^{(4) {\rm EFT}}        ', '\gamma_{\rm exp}^{(4) {\rm EFT}}'] )
    case default
    write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma4 =', self%PureEFTmodelGamma4
    write(*,'(a)')    'Please select an appropriate model.'
    end select
    ! allocate Gamma5:
    if ( allocated(self%PureEFTGamma5) ) deallocate(self%PureEFTGamma5)
    select case ( self%PureEFTmodelGamma5 )
    case(0)
    allocate( zero_parametrization_1D::self%PureEFTGamma5 )
    case(1)
    allocate( constant_parametrization_1D::self%PureEFTGamma5 )
    case(2)
    allocate( linear_parametrization_1D::self%PureEFTGamma5 )
    case(3)
    allocate( power_law_parametrization_1D::self%PureEFTGamma5 )
    call self%PureEFTGamma5%set_param_names( ['EFTGamma50  ', 'EFTGamma5Exp'], ['\gamma_0^{(5) {\rm EFT}}        ', '\gamma_{\rm exp}^{(5) {\rm EFT}}'] )
    case(4)
    allocate( exponential_parametrization_1D::self%PureEFTGamma5 )
    call self%PureEFTGamma5%set_param_names( ['EFTGamma50  ', 'EFTGamma5Exp'], ['\gamma_0^{(5) {\rm EFT}}        ', '\gamma_{\rm exp}^{(5) {\rm EFT}}'] )
    case default
    write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma5 =', self%PureEFTmodelGamma5
    write(*,'(a)')    'Please select an appropriate model.'
    end select
    ! allocate Gamma6:
    if ( allocated(self%PureEFTGamma6) ) deallocate(self%PureEFTGamma6)
    select case ( self%PureEFTmodelGamma6 )
    case(0)
    allocate( zero_parametrization_1D::self%PureEFTGamma6 )
    case(1)
    allocate( constant_parametrization_1D::self%PureEFTGamma6 )
    case(2)
    allocate( linear_parametrization_1D::self%PureEFTGamma6 )
    case(3)
    allocate( power_law_parametrization_1D::self%PureEFTGamma6 )
    call self%PureEFTGamma6%set_param_names( ['EFTGamma60  ', 'EFTGamma6Exp'], ['\gamma_0^{(6) {\rm EFT}}        ', '\gamma_{\rm exp}^{(6) {\rm EFT}}'] )
    case(4)
    allocate( exponential_parametrization_1D::self%PureEFTGamma6 )
    call self%PureEFTGamma6%set_param_names( ['EFTGamma60  ', 'EFTGamma6Exp'], ['\gamma_0^{(6) {\rm EFT}}        ', '\gamma_{\rm exp}^{(6) {\rm EFT}}'] )
    case default
    write(*,'(a,I3)') 'No model corresponding to PureEFTmodelGamma6 =', self%PureEFTmodelGamma6
    write(*,'(a)')    'Please select an appropriate model.'
    end select
    end if

    ! initialize the names:
    call self%PureEFTOmega%set_name ( 'EFTOmega' , '\Omega'       )
    call self%PureEFTwDE%set_name   ( 'EFTw'     , 'w'            )
    call self%PureEFTGamma1%set_name( 'EFTGamma1', '\gamma^{(1)}' )
    call self%PureEFTGamma2%set_name( 'EFTGamma2', '\gamma^{(2)}' )
    call self%PureEFTGamma3%set_name( 'EFTGamma3', '\gamma^{(3)}' )

    if ( .not. self%PureEFTHorndeski ) then
    call self%PureEFTGamma4%set_name( 'EFTGamma4', '\gamma^{(4)}' )
    call self%PureEFTGamma5%set_name( 'EFTGamma5', '\gamma^{(5)}' )
    call self%PureEFTGamma6%set_name( 'EFTGamma6', '\gamma^{(6)}' )
    end if

    end subroutine MGCAMBMuGammaAllocateModelSelection




end module MGCAMB_MuGamma
