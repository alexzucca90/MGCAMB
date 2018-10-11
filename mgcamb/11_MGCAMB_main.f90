!----------------------------------------------------------------------------------------
!
! This file is part of MGCAMB.
!
! Copyright (C) 2013-2017 by the MGCAMB authors
!
! The MGCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file mgcamb/LICENSE at
! the top level of the MGCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 11_MGCAMB_main.f90
!! This file contains the general MGCAMB driver. We have a type that encapsulate all
!! MGCAMB parameters and stuff and this is owned by CAMB to do all MG calculations.


!----------------------------------------------------------------------------------------
!> This module contains the general MGCAMB driver. We have a type that encapsulate all
!! MGCAMB parameters and stuff and this is owned by CAMB to do all MG calculations.

!> @author Bin Hu, Marco Raveri, Simone Peirone

module MGCAMB_main

    use precision
    use IniFile
    use AMLutils
    use MG_def
    use MGCAMB_abstract_model
    use MGCAMB_abstract_model_full
    use MGCAMB_abstract_model_designer
    use MGCAMB_pure_MG_std
    use MGCAMB_Reparametrized_Horndeski
    use MGCAMB_designer_fR
    use MGCAMB_designer_mc_quintessence
    use MGCAMB_LE_Horava
    implicit none

    private

    public MGCAMB

    !----------------------------------------------------------------------------------------
    !> This is the main object for MGCAMB and contains all the necessary ingredients to
    !! perform the MG calculations. For the physical details of the calculations implemented
    !! please refer to the Numerical Notes.
    type MGCAMB

        ! MGCAMB output root:
        character(LEN=:), allocatable :: outroot !< The root for auxiliary MGCAMB output.

        ! MGCAMB model selection flags:
        integer   :: MGflag              !< Main MGCAMB model selection flag. Decides one of the four modes to run MGCAMB.
        integer   :: PureMGmodel         !< Model selection flag for pure MG models.
        integer   :: AltParMGmodel       !< Model selection flag for alternative MG parametrizations.
        integer   :: DesignerMGmodel     !< Model selection flag for designer mapping MG models.
        integer   :: FullMappingMGmodel  !< Model selection flag for full mapping MG models.

        ! MGCAMB stability flags:
        logical   :: MG_mathematical_stability  !< Flag that extablishes wether to use mathematical stability.
        logical   :: MG_physical_stability      !< Flag that extablishes wether to use physical stability.
        logical   :: MG_additional_priors       !< Flag that extablishes wether to use some additional priors that are related to the specific model.

        ! MGCAMB model:
        class(MGCAMB_model), allocatable :: model !< This is the MGCAMB model in the main class.

        ! MGCAMB working flags:
        integer   :: MGCAMB_feedback_level      !< Amount of feedback that is printed to screen.
        real(dl)  :: MGCAMB_turn_on_time        !< Scale factor at which MGCAMB becomes active. Default set to MGturnonpiInitial in 01_MG_def.f90.
        logical   :: MGCAMB_model_is_designer   !< Logical flag that establishes whether the model is designer or not.

    contains

        ! utility functions:
        procedure :: MGCAMB_init_from_file        => read_MGCAMB_flags_from_file  !< subroutine that initializes MGCAMB from an INI file.
        procedure :: MGCAMB_init_model_from_file  => init_MGCAMB_model_from_file  !< subroutine that initializes the selected model from file.
        procedure :: MGCAMB_print_header          => print_MGCAMB_header          !< subroutine that prints to screen the MGCAMB header.
        procedure :: MGCAMB_print_CosmoMC_header  => print_MGCosmoMC_header       !< subroutine that prints to screen the MGCosmoMC header.
        procedure :: MGCAMB_print_model_feedback  => print_MGCAMB_flags           !< subroutine that prints to screen the model flags and parameters.
        ! model allocation:
        procedure :: MGCAMB_allocate_model            => allocate_MGCAMB_model            !< subroutine that, based on the model selection flags allocates the MGCAMB model.
        procedure :: MGCAMB_read_model_selection      => read_MGCAMB_model_selection      !< subroutine that reads the model selection parameters. Just a wrapper to the model specific subroutine.
        procedure :: MGCAMB_allocate_model_functions  => allocate_MGCAMB_model_functions  !< subroutine that, based on the model specific selection flags allocates the MGCAMB model functions.
        procedure :: MGCAMB_read_model_parameters     => read_MGCAMB_model_parameters     !< subroutine that reads the model parameters. Just a wrapper to the model specific subroutine.

    end type MGCAMB

    ! ---------------------------------------------------------------------------------------------

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads from file the values of the MGCAMB flags.
    subroutine read_MGCAMB_flags_from_file( self, Ini )

        implicit none

        class(MGCAMB)      :: self       !< the base class
        type(TIniFile)      :: Ini        !< Input ini file

        ! read from the INI file the main MG flag:
        self%MGflag              = Ini_Read_Int_File( Ini, 'MGflag'            , 0 )
        ! read the model selection flags:
        self%PureMGmodel         = Ini_Read_Int_File( Ini, 'PureMGmodel'       , 0 )
        self%AltParMGmodel       = Ini_Read_Int_File( Ini, 'AltParMGmodel'     , 0 )
        self%DesignerMGmodel     = Ini_Read_Int_File( Ini, 'DesignerMGmodel'   , 0 )
        self%FullMappingMGmodel  = Ini_Read_Int_File( Ini, 'FullMappingMGmodel', 0 )

        ! MGCAMB working stuff:
        self%MGCAMB_feedback_level     = Ini_Read_Int_File( Ini, 'feedback_level', 1 )
        self%MGCAMB_turn_on_time       = Ini_Read_Double_File( Ini, 'MGCAMB_turn_on_time', MGturnonpiInitial )

    end subroutine read_MGCAMB_flags_from_file

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the selected model from file.
    subroutine init_MGCAMB_model_from_file( self, Ini )

        implicit none

        class(MGCAMB)      :: self       !< the base class
        type(TIniFile)      :: Ini        !< Input ini file

        ! allocate model:
        call self%MGCAMB_allocate_model()
        ! read the parameters defining the model from file:
        call self%MGCAMB_read_model_selection( Ini )
        ! allocate model functions and parameters:
        call self%MGCAMB_allocate_model_functions( )
        ! read model parameters from file:
        call self%MGCAMB_read_model_parameters( Ini )
        ! compute model number of parameters:
        call self%model%compute_param_number()

    end subroutine init_MGCAMB_model_from_file

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the MGCAMB header.
    subroutine print_MGCAMB_header( self )

        implicit none

        class(MGCAMB) :: self       !< the base class

        ! check feedback level:
        if ( .not. self%MGCAMB_feedback_level > 0 ) return
        ! if GR return:
        if ( self%MGflag == 0 ) return
        ! print the header:
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') "     __   _______  ________   __  ______  "
        write(*,'(a)') "    /  \/  / ___/ / ___/ _ | /  |/  / _ ) "
        write(*,'(a)') "   / /\_/ / /_// / /__/ __ |/ /|_/ / _  | "
        write(*,'(a)') "  /_/  /_/____/  \___/_/ |_/_/  /_/____/  "//" "//MGCAMB_version
        write(*,'(a)') "  "
        write(*,'(a)') "***************************************************************"

    end subroutine print_MGCAMB_header

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the MGCosmoMC header.
    subroutine print_MGCosmoMC_header( self )

        implicit none

        class(MGCAMB) :: self       !< the base class

        ! check feedback level:
        if ( .not. self%MGCAMB_feedback_level > 0 ) return
        ! if GR return:
        if ( self%MGflag == 0 ) return
        ! print the header:
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') "     ___________________                    __  ________"
        write(*,'(a)') "    / __/ __/_  __/ ___/__  ___ __ _  ___  /  |/  / ___/"
        write(*,'(a)') "   / _// _/  / / / /__/ _ \(_-</  ' \/ _ \/ /|_/ / /__  "
        write(*,'(a)') "  /___/_/   /_/  \___/\___/___/_/_/_/\___/_/  /_/\___/  "
        write(*,'(a)') "  "
        write(*,'(a)') "  "//MGCAMB_version
        write(*,'(a)') "***************************************************************"

    end subroutine print_MGCosmoMC_header

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen informations about the model and the model parameters.
    subroutine print_MGCAMB_flags( self, print_params )

        implicit none

        class(MGCAMB)      :: self         !< the base class
        logical, optional   :: print_params !< optional flag that decised whether to print numerical values
                                            !! of the parameters.

        character(len=500)  :: temp_name
        real(dl)            :: temp_value
        integer             :: i

        ! check the allocation of the model:
        if ( .not. allocated(self%model) ) then
            write(*,*) 'MGCAMB WARNING: trying to call MGCAMB model feedback without allocating the model'
            call MpiStop('MGCAMB error')
        end if

        ! check feedback level:
        if ( .not. self%MGCAMB_feedback_level > 0 ) return
        ! if GR return:
        if ( self%MGflag == 0 ) return
        ! print feedback flag:
        write(*,*)
        write(*,'(a, I3)') ' MGCAMB feedback level  =', self%MGCAMB_feedback_level

        ! print stability flags:
        write(*,*)
        write(*,*) 'MGCAMB stability flags:'

        write(*,*) ' Mathematical stability = ', self%MG_mathematical_stability
        write(*,*) ' Physical stability     = ', self%MG_physical_stability
        write(*,*) ' Additional priors      = ', self%MG_additional_priors
        write(*,*)
        ! print model selection flags:
        write(*,*)              'MGCAMB model flags:'
        write(*,"(A24,I3)")     '   MGflag             =', self%MGflag
        if ( self%MGflag == 1 ) &
            write(*,"(A24,I3)") '   PureMGmodel        =', self%PureMGmodel
        if ( self%MGflag == 2 ) &
            write(*,"(A24,I3)") '   AltParMGmodel      =', self%AltParMGmodel
        if ( self%MGflag == 3 ) &
            write(*,"(A24,I3)") '   DesignerMGmodel    =', self%DesignerMGmodel
        if ( self%MGflag == 4 ) &
            write(*,"(A24,I3)") '   FullMappingMGmodel =', self%FullMappingMGmodel
        ! print model informations:
        call self%model%feedback( print_params )
        ! leave one white line:
        write(*,*)

    end subroutine print_MGCAMB_flags

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that allocates the MGCAMB model based on the model selection flags.
    !! If implementing a new model this is the place to allocate it.
    subroutine allocate_MGCAMB_model( self )

        implicit none

        class(MGCAMB)      :: self       !< the base class

        ! check the allocation of the model:
        if ( allocated(self%model) ) deallocate(self%model)

        ! do the allocation:
        select case ( self%MGflag )

            case (0)     ! GR: no need to allocate

            case (1)     ! Pure MG:

                select case ( self%PureMGmodel )
                    case(1)
                        allocate( MGCAMB_std_pure_MG::self%model )
                        call self%model%init( 'Standard Pure MG', 'Standard Pure MG' )
                    case(2)
                        allocate( MGCAMB_mod_pure_MG::self%model)
                        call self%model%init( 'Modified Pure MG', 'Modified Pure MG' )
                    case default
                        write(*,'(a,I3)') 'No model corresponding to MGFlag =', self%MGflag
                        write(*,'(a,I3)') 'and PureMGmodel =', self%PureMGmodel
                        write(*,'(a)')    'Please select an appropriate model:'
                        write(*,'(a)')    'PureMGmodel=1  standard Pure MG'
                        call MpiStop('MGCAMB error')
                end select

            case (2)     ! Alternative MG:

                select case ( self%AltParMGmodel )
                    case(1)
                        allocate( MGCAMB_RPH::self%model )
                        call self%model%init( 'RPH', 'RPH' )
                    case default
                        write(*,'(a,I3)') 'No model corresponding to MGFlag =', self%MGflag
                        write(*,'(a,I3)') 'and AltParMGmodel =', self%AltParMGmodel
                        write(*,'(a)')    'Please select an appropriate model:'
                        write(*,'(a)')    'AltParMGmodel=1  reparametrized Horndeski'
                        call MpiStop('MGCAMB error')
                end select

            case (3)     ! Designer mapping MG:

                select case ( self%DesignerMGmodel )
                    case(1)
                        allocate( MGCAMB_fR_designer::self%model )
                        call self%model%init( 'Designer f(R)', 'Designer f(R)' )
                    case(2)
                        allocate( MGCAMB_des_mc_quint::self%model )
                        call self%model%init( 'Designer minimally coupled quintessence', 'Designer minimally coupled quintessence' )
                    !> GBD mod: adding the GBD designer case
                    case(3)
                        allocate( MGCAMB_GBD_designer::self%model)
                        call self%model%init( 'Designer GBD', 'Designer GBD' )
                    case(4)
                        allocate( MGCAMB_GBD_designer_2::self%model )
                        call self%model%init( 'Designer GBD', 'Designer GBD' )
                    case(5)
                        allocate( MGCAMB_fR_designer_mod::self%model )
                        call self%model%init( 'Designer f(R)', 'Designer f(R)' )
                    case default
                        write(*,'(a,I3)') 'No model corresponding to MGFlag =', self%MGflag
                        write(*,'(a,I3)') 'and DesignerMGmodel =', self%DesignerMGmodel
                        write(*,'(a)')    'Please select an appropriate model:'
                        write(*,'(a)')    'DesignerMGmodel=1  designer f(R)'
                        write(*,'(a)')    'DesignerMGmodel=2  designer minimally coupled quintessence'
                        write(*,'(a)')    'DesignerMGmodel=3  designer GBD'
                        write(*,'(a)')    'DesignerMGmodel=3  designer GBD 2 (using DE density instead of wDE)'
                    !> GBD mod end.
                        call MpiStop('MGCAMB error')
                end select

            case (4)     ! Full mapping MG:

                select case ( self%FullMappingMGmodel )
                    case(1)
                        allocate( MGCAMB_Horava::self%model )
                        call self%model%init( 'Horava', 'Horava' )
                    case default
                        write(*,'(a,I3)') 'No model corresponding to MGFlag =', self%MGflag
                        write(*,'(a,I3)') 'and FullMappingMGmodel =', self%FullMappingMGmodel
                        write(*,'(a)')    'Please select an appropriate model:'
                        write(*,'(a)')    'FullMappingMGmodel=1  Horava gravity'
                        call MpiStop('MGCAMB error')
                end select

            case default ! not found:

                write(*,'(a,I3)') 'No model corresponding to MGFlag =', self%MGflag
                write(*,'(a)') 'Please select an appropriate model:'
                write(*,'(a)') 'MGFlag=0  GR code'
                write(*,'(a)') 'MGFlag=1  Pure MG'
                write(*,'(a)') 'MGFlag=2  MG alternative parametrizations'
                write(*,'(a)') 'MGFlag=3  designer mapping MG'
                write(*,'(a)') 'MGFlag=4  full mapping MG'
                call MpiStop('MGCAMB error')

        end select

        ! now store the designer flag:
        select type ( model => self%model )
            class is ( MGCAMB_full_model )
            self%MGCAMB_model_is_designer = .False.
            class is ( MGCAMB_designer_model )
            self%MGCAMB_model_is_designer = .True.
        end select

    end subroutine allocate_MGCAMB_model

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the model selection parameters. Just a wrapper to the model specific subroutine.
    subroutine read_MGCAMB_model_selection( self, Ini )

        implicit none

        class(MGCAMB)      :: self       !< the base class
        type(TIniFile)      :: Ini        !< Input ini file

        ! check the allocation of the model:
        if ( .not. allocated(self%model) ) then
            write(*,*) 'MGCAMB WARNING: trying to call MGCAMB model read_model_selection'
            write(*,*) ' without allocating the model'
            call MpiStop('MGCAMB error')
        end if

        ! call the model specific read parameters:
        call self%model%read_model_selection( Ini )

    end subroutine read_MGCAMB_model_selection

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that, based on the model specific selection flags allocates the MGCAMB model functions.
    !! Just a wrapper to the model specific subroutine.
    subroutine allocate_MGCAMB_model_functions( self )

        implicit none

        class(MGCAMB)      :: self       !< the base class

        ! check the allocation of the model:
        if ( .not. allocated(self%model) ) then
            write(*,*) 'MGCAMB WARNING: trying to call MGCAMB model allocate_model_selection'
            write(*,*) ' without allocating the model'
            call MpiStop('MGCAMB error')
        end if

        ! call the model specific read parameters:
        call self%model%allocate_model_selection( )

    end subroutine allocate_MGCAMB_model_functions

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the model parameters. Just a wrapper to the model specific subroutine.
    subroutine read_MGCAMB_model_parameters( self, Ini )

        implicit none

        class(MGCAMB)      :: self       !< the base class
        type(TIniFile)      :: Ini        !< Input ini file

        ! check the allocation of the model:
        if ( .not. allocated(self%model) ) then
            write(*,*) 'MGCAMB WARNING: trying to call MGCAMB model read_model_parameters_from_file'
            write(*,*) ' without allocating the model'
            call MpiStop('MGCAMB error')
        end if

        ! call the model specific read parameters:
        call self%model%init_model_parameters_from_file( Ini )

    end subroutine read_MGCAMB_model_parameters

    ! ---------------------------------------------------------------------------------------------

end module MGCAMB_main

!----------------------------------------------------------------------------------------
