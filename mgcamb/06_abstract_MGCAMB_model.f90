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

        !> this might be useless
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
        integer                      , intent(in)    :: feedback_level !< level of feedback from the background code. 0=none; 1=some; 2=chatty.
        logical                      , intent(out)   :: success        !< wether the background initialization succeded or not

        success = .True.

    end subroutine MGCAMBModelInitBackground

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that computes \rho_Q and P_Q. For details refer to the numerical notes.
    subroutine EFTCAMBModelComputeRhoQPQ( self, a, eft_par_cache, eft_cache )

    implicit none

    class(EFTCAMB_model)                         :: self          !< the base class
    real(dl), intent(in)                         :: a             !< the input scale factor.
    type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
    type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

    real(dl) :: a2, adotoa2, aomegaP

    ! precompute some parts:
    a2      = a*a
    adotoa2 = eft_cache%adotoa**2
    aomegaP = a*eft_cache%EFTOmegaP

    ! do the computations:
    eft_cache%grhoq     = 2._dl*eft_cache%EFTc -eft_cache%EFTLambda -3._dl*adotoa2*aomegaP
    eft_cache%gpresq    = eft_cache%EFTLambda + a2*adotoa2*eft_cache%EFTOmegaPP +aomegaP*(eft_cache%Hdot+2._dl*adotoa2)
    eft_cache%grhodotq  = 3._dl*eft_cache%adotoa*(-eft_cache%grhoq-eft_cache%gpresq+adotoa2*aomegaP )
    eft_cache%gpresdotq = eft_cache%EFTLambdadot &
    & +adotoa2*eft_cache%adotoa*(a*a2*eft_cache%EFTOmegaPPP-2._dl*aomegaP+2._dl*a2*eft_cache%EFTOmegaPP) &
    & +aomegaP*eft_cache%Hdotdot &
    & +3._dl*eft_cache%adotoa*eft_cache%Hdot*( aomegaP+a2*eft_cache%EFTOmegaPP )

    end subroutine EFTCAMBModelComputeRhoQPQ

! ---------------------------------------------------------------------------------------------
!> Subroutine that computes the Einstein equations factors. For details refer to the numerical notes.
!! The implementation might look a bit intricated because this is a crucial part performance wise.
subroutine EFTCAMBModelComputeEinsteinFactors( self, a, eft_par_cache, eft_cache )

implicit none

class(EFTCAMB_model)                         :: self          !< the base class
real(dl), intent(in)                         :: a             !< the input scale factor.
type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

real(dl) :: one_plus_omega, adotoa2, aomegaP, k2, a2, hdot_m_adotoa2, pidot_p_H_pi

! precompute some common parts:
one_plus_omega = 1._dl+eft_cache%EFTOmegaV
adotoa2        = eft_cache%adotoa**2
hdot_m_adotoa2 = eft_cache%Hdot-adotoa2
aomegaP        = a*eft_cache%EFTOmegaP
k2             = eft_cache%k**2
a2             = a**2
pidot_p_H_pi   = eft_cache%pidot+eft_cache%adotoa*eft_cache%pi
!write(*,*) "pi, pidot:", eft_cache%pi,eft_cache%pidot
! compute the coefficients:
eft_cache%EFTeomF     = 1.5_dl/(eft_cache%k*one_plus_omega)*( (eft_cache%grhoq+eft_cache%gpresq)*eft_cache%pi &
& + (aomegaP*eft_cache%adotoa+a*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V)*pidot_p_H_pi &
& + eft_cache%pi*( k2*(eft_cache%EFTGamma3V+eft_cache%EFTGamma4V) -(3._dl*eft_cache%EFTGamma3V+eft_cache%EFTGamma4V)*hdot_m_adotoa2  ) )
eft_cache%EFTeomG     = +1._dl +0.5/one_plus_omega*( aomegaP &
& +a*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V/eft_cache%adotoa +3._dl*eft_cache%EFTGamma3V +eft_cache%EFTGamma4V )
eft_cache%EFTeomL     = +0.5_dl/one_plus_omega*( (2._dl*eft_cache%EFTc*pidot_p_H_pi +eft_cache%grhodotq*eft_cache%pi)/eft_cache%adotoa &
& -3._dl*aomegaP*( (3._dl*adotoa2-eft_cache%Hdot+k2/3._dl)*eft_cache%pi +eft_cache%adotoa*eft_cache%pidot) &
& +4._dl*a2*eft_par_cache%h0_mpc**2*eft_cache%EFTGamma1V/eft_cache%adotoa*pidot_p_H_pi &
& +3._dl*a*eft_par_cache%h0_mpc*eft_cache%EFTGamma2V*( (eft_cache%Hdot-2._dl*adotoa2-k2/3._dl)*eft_cache%pi/eft_cache%adotoa -eft_cache%pidot ) &                                                                                  ! Gamma2
& -eft_cache%pi*( k2-3._dl*hdot_m_adotoa2)*(3._dl*eft_cache%EFTGamma3V+eft_cache%EFTGamma4V) +8._dl*eft_cache%EFTGamma6V*k2*pidot_p_H_pi/eft_cache%adotoa )&
&+2._dl*eft_cache%EFTGamma5V/one_plus_omega*eft_cache%k**2*eft_cache%pi
eft_cache%EFTeomM     = eft_cache%gpresdotq*eft_cache%pi +(eft_cache%grhoq+eft_cache%gpresq+a2*adotoa2*eft_cache%EFTOmegaPP)*pidot_p_H_pi &
& +aomegaP*eft_cache%adotoa*( eft_cache%pidotdot +(eft_cache%Hdot+4._dl*adotoa2)*eft_cache%pidot/eft_cache%adotoa +2._dl*(eft_cache%Hdot+3._dl*adotoa2+k2/3._dl)*eft_cache%pi ) &
& +a*eft_par_cache%h0_mpc*( eft_cache%EFTGamma2V*eft_cache%pidotdot &
& +(4._dl*eft_cache%EFTGamma2V +a*eft_cache%EFTGamma2P)*eft_cache%adotoa*eft_cache%pidot +(3._dl*adotoa2*eft_cache%EFTGamma2V &
& +eft_cache%Hdot*eft_cache%EFTGamma2V +a*adotoa2*eft_cache%EFTGamma2P)*eft_cache%pi) &
& -(hdot_m_adotoa2-k2/3._dl)*( (3._dl*eft_cache%EFTGamma3V+eft_cache%EFTGamma4V)*eft_cache%pidot &
& +2._dl*eft_cache%adotoa*(+3._dl*eft_cache%EFTGamma3V+eft_cache%EFTGamma4V+1.5_dl*a*eft_cache%EFTGamma3P+0.5_dl*a*eft_cache%EFTGamma4P)*eft_cache%pi ) &
& -(3._dl*eft_cache%EFTGamma3V+eft_cache%EFTGamma4V)*(eft_cache%Hdotdot-2._dl*eft_cache%adotoa*eft_cache%Hdot)*eft_cache%pi&
& -4._dl*eft_cache%EFTGamma5V*k2*pidot_p_H_pi/3._dl
eft_cache%EFTeomN     = eft_cache%k/one_plus_omega*( eft_cache%adotoa*eft_cache%pi*(-aomegaP+2._dl*eft_cache%EFTGamma4V+a*eft_cache%EFTGamma4P) &
& +eft_cache%EFTGamma4V*eft_cache%pidot +2._dl*eft_cache%EFTGamma5V*pidot_p_H_pi )
eft_cache%EFTeomNdot  =  eft_cache%k/one_plus_omega*( -eft_cache%Hdot*aomegaP*eft_cache%pi &
& -eft_cache%adotoa*aomegaP*eft_cache%pidot &
& -adotoa2*(aomegaP+a2*eft_cache%EFTOmegaPP-aomegaP**2/one_plus_omega)*eft_cache%pi &
& +eft_cache%EFTGamma4V*eft_cache%pidotdot +a*eft_cache%adotoa*eft_cache%pidot*&
&( +eft_cache%EFTGamma4P -eft_cache%EFTGamma4V*eft_cache%EFTOmegaP/one_plus_omega)&
& +2._dl*(eft_cache%EFTGamma4V +0.5_dl*a*eft_cache%EFTGamma4P)*( eft_cache%Hdot*eft_cache%pi +eft_cache%adotoa*eft_cache%pidot)&
& +2._dl*a*adotoa2*eft_cache%pi*(+0.5_dl*a*eft_cache%EFTGamma4PP +1.5_dl*eft_cache%EFTGamma4P&
& -eft_cache%EFTOmegaP/one_plus_omega*(eft_cache%EFTGamma4V +0.5_dl*a*eft_cache%EFTGamma4P))&
& +2._dl*eft_cache%EFTGamma5V*( eft_cache%pidotdot+eft_cache%adotoa*eft_cache%pidot+eft_cache%Hdot*eft_cache%pi)&
& +2._dl*eft_cache%adotoa*pidot_p_H_pi*( +a*eft_cache%EFTGamma5P-eft_cache%EFTGamma5V*aomegaP/one_plus_omega) )
eft_cache%EFTeomU     = 1._dl +(+1.5_dl*eft_cache%EFTGamma3V+0.5_dl*eft_cache%EFTGamma4V)/one_plus_omega
eft_cache%EFTeomV     = +0.5_dl/one_plus_omega*( aomegaP -2._dl*eft_cache%EFTGamma4V -a*eft_cache%EFTGamma4P )
eft_cache%EFTeomVdot  = 0.5_dl*eft_cache%adotoa/one_plus_omega*( aomegaP-3._dl*a*eft_cache%EFTGamma4P &
& +a2*(eft_cache%EFTOmegaPP-eft_cache%EFTGamma4PP) +aomegaP/one_plus_omega*(-aomegaP+2._dl*eft_cache%EFTGamma4V+a*eft_cache%EFTGamma4P))
eft_cache%EFTeomX     = 1._dl -eft_cache%EFTGamma4V/one_plus_omega
eft_cache%EFTeomXdot  = -a*eft_cache%adotoa/one_plus_omega*( +eft_cache%EFTGamma4P &
& -eft_cache%EFTGamma4V*eft_cache%EFTOmegaP/one_plus_omega)
eft_cache%EFTeomY     = +0.5_dl/one_plus_omega*( aomegaP &
& +3._dl*eft_cache%EFTGamma3V+eft_cache%EFTGamma4V   &
& +0.5_dl*a*(3._dl*eft_cache%EFTGamma3P+eft_cache%EFTGamma4P) )
eft_cache%EFTeomQ =1._dl+2._dl*eft_cache%EFTGamma5V/one_plus_omega

end subroutine EFTCAMBModelComputeEinsteinFactors


! ---------------------------------------------------------------------------------------------
!> Subroutine that computes the factors for the tensor propagation equation. For details refer to the numerical notes.
subroutine EFTCAMBModelComputeTensorFactors( self, a, eft_par_cache, eft_cache )

implicit none

class(EFTCAMB_model)                         :: self          !< the base class
real(dl), intent(in)                         :: a             !< the input scale factor.
type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

real(dl) :: one_plus_omega

! precompute some common parts:
one_plus_omega = 1._dl+eft_cache%EFTOmegaV

! compute the coefficients:
eft_cache%EFTAT = one_plus_omega -eft_cache%EFTGamma4V
eft_cache%EFTBT = 2._dl*eft_cache%adotoa*( one_plus_omega -eft_cache%EFTGamma4V +0.5_dl*a*eft_cache%EFTOmegaP -0.5_dl*a*eft_cache%EFTGamma4P )
eft_cache%EFTDT = one_plus_omega

end subroutine EFTCAMBModelComputeTensorFactors

! ---------------------------------------------------------------------------------------------
!> Subroutine that computes the kinetic and gradient terms. For details refer to the numerical notes.
subroutine EFTCAMBModelComputeStabilityFactors( self, a, eft_par_cache, eft_cache )

implicit none

class(EFTCAMB_model)                         :: self          !< the base class
real(dl), intent(in)                         :: a             !< the input scale factor.
type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.


!write(*,*) "Stability factors."

eft_cache%EFT_kinetic  = 9._dl*( 1._dl +eft_cache%EFTOmegaV -eft_cache%EFTGamma4V )*( 4._dl*eft_cache%EFTc*( 1._dl +eft_cache%EFTOmegaV -eft_cache%EFTGamma4V ) &
& +3._dl*eft_cache%adotoa**2*eft_cache%EFTOmegaP**2*a**2 + a**2*eft_par_cache%h0_Mpc*( eft_par_cache%h0_Mpc*( 3._dl*eft_cache%EFTGamma2V**2 +8._dl*eft_cache%EFTGamma1V* &
&( 1._dl +eft_cache%EFTOmegaV -eft_cache%EFTGamma4V ) +6._dl*eft_cache%adotoa*eft_cache%EFTGamma2V*eft_cache%EFTOmegaP ) ) )
!
eft_cache%EFT_gradient = 9._dl*(8._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5P - 16._dl*eft_cache%adotoa**2*eft_cache%EFTGamma5V**2 + 16._dl*eft_cache%EFTc*eft_cache%EFTGamma5V**2 &
&- 2._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma4P*eft_cache%EFTOmegaP + 4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma5P*eft_cache%EFTOmegaP -&
&4._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP - 4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma4P*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP &
&- 8._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5V**2*eft_cache%EFTOmegaP + 3._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTOmegaP**2 +&
&4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP**2 + 16._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5P*eft_cache%EFTOmegaV &
&+ 16._dl*eft_cache%EFTc*eft_cache%EFTGamma5V*eft_cache%EFTOmegaV - 16._dl*eft_cache%adotoa**2*eft_cache%EFTGamma5V**2*eft_cache%EFTOmegaV -&
&2._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma4P*eft_cache%EFTOmegaP*eft_cache%EFTOmegaV + 4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma5P*eft_cache%EFTOmegaP*eft_cache%EFTOmegaV &
&- 4._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP*eft_cache%EFTOmegaV +3._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTOmegaP**2*eft_cache%EFTOmegaV &
&+ 8._dl*a*eft_cache%adotoa**2*eft_cache%EFTGamma5P*eft_cache%EFTOmegaV**2 - a**2*eft_cache%EFTGamma2V**2*eft_par_cache%h0_mpc**2*(1 + eft_cache%EFTOmegaV) +&
&4._dl*a**2*eft_cache%adotoa**2*eft_cache%EFTGamma5V*eft_cache%EFTOmegaPP*(1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV) + 4._dl*eft_cache%EFTc*(4._dl*eft_cache%EFTGamma5V &
&+ (1._dl + eft_cache%EFTOmegaV)**2) -2._dl*a*eft_cache%adotoa*eft_par_cache%h0_Mpc*(a*eft_cache%EFTGamma2P*(1._dl - eft_cache%EFTGamma4V + eft_cache%EFTOmegaV)*(1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV) +&
&eft_cache%EFTGamma2V*(eft_cache%EFTGamma4V*(-1._dl + 2._dl*a*eft_cache%EFTGamma5P + 2._dl*eft_cache%EFTGamma5V + a*eft_cache%EFTOmegaP - eft_cache%EFTOmegaV) +&
&(1._dl + eft_cache%EFTOmegaV)*(1._dl + a*(eft_cache%EFTGamma4P - 2._dl*eft_cache%EFTGamma5P - eft_cache%EFTOmegaP) + eft_cache%EFTOmegaV) - 2._dl*eft_cache%EFTGamma5V*(1._dl - a*eft_cache%EFTGamma4P &
&+ a*eft_cache%EFTOmegaP + eft_cache%EFTOmegaV))) +8._dl*eft_cache%EFTGamma5V*eft_cache%Hdot + 16._dl*eft_cache%EFTGamma5V**2*eft_cache%Hdot + 4._dl*a*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP*eft_cache%Hdot &
&+ 8._dl*a*eft_cache%EFTGamma5V**2*eft_cache%EFTOmegaP*eft_cache%Hdot + 16._dl*eft_cache%EFTGamma5V*eft_cache%EFTOmegaV*eft_cache%Hdot +&
&16._dl*eft_cache%EFTGamma5V**2*eft_cache%EFTOmegaV*eft_cache%Hdot + 4._dl*a*eft_cache%EFTGamma5V*eft_cache%EFTOmegaP*eft_cache%EFTOmegaV*eft_cache%Hdot + 8._dl*eft_cache%EFTGamma5V*eft_cache%EFTOmegaV**2*eft_cache%Hdot +&
&4._dl*eft_cache%EFTGamma4V**2*(eft_cache%adotoa**2*(1._dl + 2._dl*a*eft_cache%EFTGamma5P + 4._dl*eft_cache%EFTGamma5V + a*eft_cache%EFTOmegaP + eft_cache%EFTOmegaV) - (1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV)*eft_cache%Hdot) +&
&2._dl*eft_cache%EFTGamma4V*(eft_cache%adotoa**2*(-(a**2*eft_cache%EFTOmegaP**2) + a**2*eft_cache%EFTOmegaPP*(1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV) -&
&4._dl*(1._dl + eft_cache%EFTOmegaV)*(1._dl + 2._dl*a*eft_cache%EFTGamma5P + 4._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV) - a*eft_cache%EFTOmegaP*(3._dl + 2._dl*a*eft_cache%EFTGamma5P &
&+ 2._dl*eft_cache%EFTGamma5V + 3._dl*eft_cache%EFTOmegaV)) +(1._dl + 2._dl*eft_cache%EFTGamma5V + eft_cache%EFTOmegaV)*(4._dl + a*eft_cache%EFTOmegaP + 4._dl*eft_cache%EFTOmegaV)*eft_cache%Hdot))

!write(*,*) "a,K,G:",a,eft_cache%EFT_kinetic,eft_cache%EFT_gradient

end subroutine EFTCAMBModelComputeStabilityFactors

! ---------------------------------------------------------------------------------------------
!> Function that computes model specific stability requirements.
function EFTCAMBModelAdditionalModelStability( self, a, eft_par_cache, eft_cache )

implicit none

class(EFTCAMB_model)                         :: self          !< the base class
real(dl), intent(in)                         :: a             !< the input scale factor.
type(EFTCAMB_parameter_cache), intent(inout) :: eft_par_cache !< the EFTCAMB parameter cache that contains all the physical parameters.
type(EFTCAMB_timestep_cache ), intent(inout) :: eft_cache     !< the EFTCAMB timestep cache that contains all the physical values.

logical :: EFTCAMBModelAdditionalModelStability               !< the return value of the stability computation. True if the model specific stability criteria are met, false otherwise.

EFTCAMBModelAdditionalModelStability = .True.

end function EFTCAMBModelAdditionalModelStability



! ---------------------------------------------------------------------------------------------
!> Subroutine that reads the parameters of the model for sampling
subroutine EFTCAMBModelInitModelParametersSampling( self, array, sampling_params )
use precision
use EFT_sampler
!import EFTCAMB_model
implicit none
class(EFTCAMB_model)  :: self   !< the base class
real(dl), dimension(self%parameter_number), intent(in) :: array
type(EFTSamplingParameters) :: sampling_params
end subroutine EFTCAMBModelInitModelParametersSampling

! ---------------------------------------------------------------------------------------------





end module MGCAMB_abstract_model

!----------------------------------------------------------------------------------------
