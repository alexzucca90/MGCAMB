!----------------------------------------------------------------------------------------
!
! This file is part of MGCAMB.
!
! Copyright (C) 2007-201 by the MGCAMB authors
!
! The MGCAMB code is free software;
! You can use it, redistribute it, and/or modify it under the terms
! of the GNU General Public License as published by the Free Software Foundation;
! either version 3 of the License, or (at your option) any later version.
! The full text of the license can be found in the file MGcamb/LICENSE at
! the top level of the MGCAMB distribution.
!
!----------------------------------------------------------------------------------------

!> @file 06_abstract_MGCAMB_model.f90
!! This file contains the abstract definition of all the places where MGCAMB interacts
!! with CAMB. All MGCAMB models should inherit from this class or the two derived
!! classes contained in 06p2_abstract_MGCAMB_full.f90 or 06p3_abstract_MGCAMB_designer.f90


!----------------------------------------------------------------------------------------
!> This module contains the abstract definition of all the places where MGCAMB interacts
!! with CAMB. All MGCAMB models should inherit from this class or the two derived
!! classes contained in 06p2_abstract_MGCAMB_full.f90 or 06p3_abstract_MGCAMB_designer.f90

!> @author Alex Zucca

module MGCAMB_abstract_model

    use precision
    use IniFile
    use MGCAMB_cache

    implicit none

    private

    public MGCAMB_model

    !----------------------------------------------------------------------------------------
    !> This is the abstract type for MGCAMB models. As a rule, when there is a
    !! new model it should be declared as a class inheriting from MGCAMB_model.
    !! This guarantees maximum performances as well as maximum flexibility.
    type, abstract :: MGCAMB_model

        integer                       :: parameter_number !< number of parameters of the model.
        character(len=:), allocatable :: name             !< name of the model.
        character(len=:), allocatable :: name_latex       !< latex name of the model.

    contains

        ! initialization of the model:
        procedure :: init                        => MGCAMBModelInitialize                                !< subroutine that initializes the name and latex name of the model.
        procedure(MGCAMBModelReadModelSelectionFromFile  ), deferred :: read_model_selection             !< subroutine that reads the parameters of the model from file.
        procedure(MGCAMBModelAllocateModelSelection      ), deferred :: allocate_model_selection         !< subroutine that allocates the model selection.
        procedure(MGCAMBModelInitModelParameters         ), deferred :: init_model_parameters            !< subroutine taht initializes the model parameters based on the values found in an input array.
        procedure(MGCAMBModelInitModelParametersFromFile ), deferred :: init_model_parameters_from_file  !< subroutine that reads the parameters of the model from file.

        ! utility functions:
        procedure(MGCAMBModelComputeParametersNumber     ), deferred :: compute_param_number             !< subroutine that computes the number of parameters of the model.
        procedure(MGCAMBModelFeedback                    ), deferred :: feedback                         !< subroutine that prints on the screen feedback information about the model.
        procedure(MGCAMBModelParameterNames              ), deferred :: parameter_names                  !< subroutine that returns the i-th parameter name of the model.
        procedure(MGCAMBModelParameterNamesLatex         ), deferred :: parameter_names_latex            !< subroutine that returns the i-th parameter name of the model.
        procedure(MGCAMBModelParameterValues             ), deferred :: parameter_values                 !< subroutine that returns the i-th parameter value.

        ! background initialization functions:
        procedure :: initialize_background       => MGCAMBModelInitBackground                            !< subroutine that initializes the background of the model, if needed.

        ! CAMB related procedures:
        procedure(MGCAMBModelBackgroundMGFunctions ), deferred :: compute_background_MG_functions      !< subroutine that computes the value of the background MG functions at a given time.

        procedure(MGCAMBModelComputeSigma), deferred :: compute_sigma   !< subroutine that computes the shear perturbation sigma
        procedure(MGCAMBModelComputeZ),     deferred :: compute_z       !< subroutine that computes the expansion perturbation Z

        !> in case we modify the background
        procedure(MGCAMBModelComputeDtauda          ), deferred :: compute_dtauda                        !< function that computes dtauda = 1/sqrt(a^2H^2).
        procedure(MGCAMBModelComputeAdotoa          ), deferred :: compute_adotoa                        !< subroutine that computes adotoa = H.
        procedure(MGCAMBModelComputeHubbleDer       ), deferred :: compute_H_derivs                      !< subroutine that computes the two derivatives wrt conformal time of H.


    end type MGCAMB_model

    ! ---------------------------------------------------------------------------------------------
    ! MGCAMB abstract interfaces: these are all the model procedures that the user HAS to override
    ! when writing its own model.
    ! ---------------------------------------------------------------------------------------------

    abstract interface

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that computes the number of parameters of the model.
        subroutine MGCAMBModelComputeParametersNumber( self )
            import MGCAMB_model
            implicit none
            class(MGCAMB_model)       :: self   !< the base class
        end subroutine MGCAMBModelComputeParametersNumber

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that reads the parameters of the model from file.
        subroutine MGCAMBModelReadModelSelectionFromFile( self, Ini )
            use IniFile
            import MGCAMB_model
            implicit none
            class(MGCAMB_model)  :: self   !< the base class
            type(TIniFile)        :: Ini    !< Input ini file
        end subroutine MGCAMBModelReadModelSelectionFromFile

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that reads the parameters of the model from file.
        subroutine MGCAMBModelAllocateModelSelection( self )
            import MGCAMB_model
            implicit none
            class(MGCAMB_model)  :: self   !< the base class
        end subroutine MGCAMBModelAllocateModelSelection

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that reads the parameters of the model from file.
        subroutine MGCAMBModelInitModelParameters( self, array )
            use precision
            import MGCAMB_model
            implicit none
            class(MGCAMB_model)                                   :: self   !< the base class
            real(dl), dimension(self%parameter_number), intent(in) :: array  !< input array with the values of the parameters.
        end subroutine MGCAMBModelInitModelParameters

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that reads the parameters of the model from file.
        subroutine MGCAMBModelInitModelParametersFromFile( self, Ini )
            use IniFile
            import MGCAMB_model
            implicit none
            class(MGCAMB_model)  :: self   !< the base class
            type(TIniFile)        :: Ini    !< Input ini file
        end subroutine MGCAMBModelInitModelParametersFromFile

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that prints on the screen feedback information about the model.
        subroutine MGCAMBModelFeedback( self, print_params )
            import MGCAMB_model
            implicit none
            class(MGCAMB_model)  :: self         !< the base class
            logical, optional     :: print_params !< optional flag that decised whether to print numerical values
                                                  !! of the parameters.
        end subroutine MGCAMBModelFeedback

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that returns the i-th parameter name of the model.
        subroutine MGCAMBModelParameterNames( self, i, name )
            import MGCAMB_model
            implicit none
            class(MGCAMB_model)      :: self   !< the base class
            integer     , intent(in)  :: i      !< the index of the parameter
            character(*), intent(out) :: name   !< the output name of the i-th parameter
        end subroutine MGCAMBModelParameterNames

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that returns the i-th parameter name of the model.
        subroutine MGCAMBModelParameterNamesLatex( self, i, latexname )
            import MGCAMB_model
            implicit none
            class(MGCAMB_model)       :: self        !< the base class
            integer     , intent(in)  :: i          !< The index of the parameter
            character(*), intent(out) :: latexname  !< the output latex name of the i-th parameter
        end subroutine MGCAMBModelParameterNamesLatex

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that returns the i-th parameter name of the model.
        subroutine MGCAMBModelParameterValues( self, i, value )
            use precision
            import MGCAMB_model
            implicit none
            class(MGCAMB_model)  :: self   !< the base class
            integer , intent(in)  :: i      !< The index of the parameter
            real(dl), intent(out) :: value  !< the output value of the i-th parameter
        end subroutine MGCAMBModelParameterValues

        ! ---------------------------------------------------------------------------------------------
        !> Subroutine that computes the value of the background MG functions at a given time.
        subroutine MGCAMBModelBackgroundMGFunctions( self, a, MG_par_cache, MG_cache )
            use precision
            use MGCAMB_cache
            import MGCAMB_model
            implicit none
            class(MGCAMB_model)                         :: self          !< the base class.
            real(dl), intent(in)                        :: a             !< the input scale factor.
            type(MGCAMB_parameter_cache), intent(inout) :: MG_par_cache !< the MGCAMB parameter cache that contains all the physical parameters.
            type(MGCAMB_timestep_cache ), intent(inout) :: MG_cache     !< the MGCAMB timestep cache that contains all the physical values.
        end subroutine MGCAMBModelBackgroundMGFunctions


        subroutine MGCAMBModelComputeSigma( self, a, MG_par_cache, MG_cache )
            use precision
            use MGCAMB_cache
            import MGCAMB_model
            implicit none
            class(MGCAMB_model)                         :: self          !< the base class.
            real(dl), intent(in)                        :: a             !< the input scale factor.
            type(MGCAMB_parameter_cache), intent(inout) :: MG_par_cache !< the MGCAMB parameter cache that contains all the physical parameters.
            type(MGCAMB_timestep_cache ), intent(inout) :: MG_cache     !< the MGCAMB timestep cache that contains all the physical values.
        subroutine MGCAMBModelComputeSigma


        subroutine MGCAMBModelComputeZ( self, a, MG_par_cache, MG_cache )
            use precision
            use MGCAMB_cache
            import MGCAMB_model
            implicit none
            class(MGCAMB_model)                         :: self          !< the base class.
            real(dl), intent(in)                        :: a             !< the input scale factor.
            type(MGCAMB_parameter_cache), intent(inout) :: MG_par_cache !< the MGCAMB parameter cache that contains all the physical parameters.
            type(MGCAMB_timestep_cache ), intent(inout) :: MG_cache     !< the MGCAMB timestep cache that contains all the physical values.
        subroutine MGCAMBModelComputeZ


    ! ---------------------------------------------------------------------------------------------

    end interface

contains

    ! ---------------------------------------------------------------------------------------------
    ! MGCAMB abstract model implementation: the following are all the procedures that can be
    ! be safely implemented for the abstract class and are not harmful if not overritten.
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the name and latex name of the model.
    subroutine MGCAMBModelInitialize( self, name, latexname )

        implicit none

        class(MGCAMB_model)      :: self      !< the base class
        character(*), intent(in) :: name      !< the name of the function
        character(*), intent(in) :: latexname !< the latex name of the function

        self%name       = TRIM(name)
        self%name_latex = TRIM(latexname)

    end subroutine MGCAMBModelInitialize


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes the background of the model, if needed.
    subroutine MGCAMBModelInitBackground( self, params_cache, feedback_level, success )

        implicit none

        class(MGCAMB_model)                         :: self           !< the base class
        type(MGCAMB_parameter_cache), intent(in)    :: params_cache   !< a MGCAMB parameter cache containing cosmological parameters
        integer                     , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                     , intent(out)   :: success        !< wether the background initialization succeded or not

        success = .True.

    end subroutine MGCAMBModelInitBackground



end module MGCAMB_abstract_model

!----------------------------------------------------------------------------------------
