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

!> @file 03_abstract_MGCAMB_cache.f90
!! This file contains the definition of the MGCAMB caches.
!! These are used to store parameters that can be used by MGCAMB, in MGCAMB_parameter_cache
!! or are used to store partial results when solving the time evolution of perturbations,
!! in MGCAMB_timestep_cache.


!----------------------------------------------------------------------------------------
!> This module contains the definition of the MGCAMB caches.
!! These are used to store parameters that can be used by MGCAMB, in MGCAMB_parameter_cache
!! or are used to store partial results when solving the time evolution of perturbations,
!! in MGCAMB_timestep_cache.

!> @author Alex Zucca azucca@sfu.ca

module MGCAMB_cache

    use precision
    use IniFile
    use AMLutils

    implicit none

    private

    public MGCAMB_parameter_cache, MGCAMB_timestep_cache

    ! some settings:
    character(*), parameter :: cache_output_format = 'e18.10'

    !----------------------------------------------------------------------------------------
    !> This is the type that defines the MGCAMB parameter cache. The idea is to copy in here
    !! the cosmological parameters that we need from CAMB and then use this for the interfaces.
    !! This also contains a wrapper for Massive Neutrinos quantities that can be used in
    !! the background code.
    type :: MGCAMB_parameter_cache

        ! 1) relative densities:
        real(dl) :: omegac         !< the value of \f$ \Omega_{\rm CDM}^0 \f$.
        real(dl) :: omegab         !< the value of \f$ \Omega_{\rm b}^0 \f$.
        real(dl) :: omegav         !< the value of \f$ \Omega_{\Lambda}^0 \f$.
        real(dl) :: omegak         !< the value of \f$ \Omega_{\rm K}^0 \f$.
        real(dl) :: omegan         !< the value of \f$ \Omega_{m\nu}^0 \f$ of massive neutrinos.
        real(dl) :: omegag         !< the value of \f$ \Omega_{ \gamma }^0 \f$.
        real(dl) :: omegar         !< the value of \f$ \Omega_{ \nu }^0 \f$ of massless neutrinos.
        ! 2) Hubble constant:
        real(dl) :: h0             !< reduced Hubble constant \f$ H_0/100 \f$.
        real(dl) :: h0_Mpc         !< the Hubble constant in MegaParsec \f$ 10^3 \cdot H_0/c \f$.
        ! 3) densities:
        real(dl) :: grhog          !< the value of \f$ 8 \pi G_{N} \rho_{\gamma}(t_0) \f$.
        real(dl) :: grhornomass    !< the value of \f$ 8 \pi G_{N} \rho_{\nu}(t_0) \f$.
        real(dl) :: grhoc          !< the value of \f$ 8 \pi G_{N} \rho_{\rm CDM}(t_0) \f$.
        real(dl) :: grhob          !< the value of \f$ 8 \pi G_{N} \rho_{\rm b}(t_0) \f$.
        real(dl) :: grhov          !< the value of \f$ 8 \pi G_{N} \rho_{\Lambda}(t_0) \f$.
        real(dl) :: grhok          !< the value of \f$ 8 \pi G_{N} \rho_{\rm K}(t_0) \f$.
        ! 4) massive neutrinos:
        integer  :: Num_Nu_Massive                        !< number of massive neutrinos
        integer  :: Nu_mass_eigenstates                   !< number of mass eigenstates
        real(dl), allocatable, dimension(:) :: grhormass  !< densities of neutrinos in each mass eigenstate \f$ 8 \pi G_{N} \rho_{ m\nu }(t_0) \f$
        real(dl), allocatable, dimension(:) :: nu_masses  !< neutrino masses
        ! 5) massive neutrinos wrapper:
        procedure( Nu_background_Wrapper    ), pointer, nopass :: Nu_background    => null()  !< wrapper to the subroutine that computes the background massive neutrinos density and pressure.
        procedure( Nu_rho_Wrapper           ), pointer, nopass :: Nu_rho           => null()  !< wrapper to the subroutine that computes the background massive neutrinos density.
        procedure( Nu_drho_Wrapper          ), pointer, nopass :: Nu_drho          => null()  !< wrapper to the subroutine that computes the time derivative of the background massive neutrinos density.
        procedure( Nu_pidot_Wrapper         ), pointer, nopass :: Nu_pidot         => null()  !< wrapper to the function that computes the background massive neutrinos time derivative of pressure.
        procedure( Nu_pidotdot_Wrapper      ), pointer, nopass :: Nu_pidotdot      => null()  !< wrapper to the function that computes the background massive neutrinos second time derivative of pressure.
        procedure( Nu_pidotdotdot_Wrapper   ), pointer, nopass :: Nu_pidotdotdot   => null()  !< wrapper to the function that computes the background massive neutrinos third time derivative of pressure.

    contains

        procedure :: initialize => MGCAMBParameterCacheInit  !< subroutine that initializes to zero all the elements of the parameter cache.
        procedure :: is_nan     => MGCAMBParameterCacheIsNan !< Subroutine that check if an element of the MGCAMB_parameter_cache is Nan.
        procedure :: print      => MGCAMBParameterCachePrint !< subroutine that prints the MGCAMB parameters cache to screen.

    end type MGCAMB_parameter_cache

    !----------------------------------------------------------------------------------------
    ! Interface containing the wrapper to massive neutrinos stuff.
    interface
        !----------------------------------------------------------------------------------------
        !> Wrapper to the subroutine that computes the background massive neutrinos density
        !! and pressure.
        subroutine Nu_background_Wrapper( am, rhonu, pnu )
            use precision
            implicit none
            real(dl), intent(in)  :: am     !< input scale factor times the neutrino mass
            real(dl), intent(out) :: rhonu  !< output neutrino density \f$ \frac{\rho_{\nu} a^2}{m_0^2} \f$
            real(dl), intent(out) :: pnu    !< output neutrino pressure \f$ \frac{ P_{\nu} a^2}{m_0^2} \f$
        end subroutine Nu_background_Wrapper
        !----------------------------------------------------------------------------------------
        !> Wrapper to the subroutine that computes the background massive neutrinos density.
        subroutine Nu_rho_Wrapper( am, rhonu )
            use precision
            implicit none
            real(dl), intent(in)  :: am     !< input scale factor times the neutrino mass
            real(dl), intent(out) :: rhonu  !< output neutrino density \f$ \frac{\rho_{\nu} a^2}{m_0^2} \f$
        end subroutine Nu_rho_Wrapper
        !----------------------------------------------------------------------------------------
        !> Wrapper to the subroutine that computes the time derivative of the background
        !! massive neutrinos density.
        function Nu_drho_Wrapper( am, adotoa, rhonu )
            use precision
            implicit none
            real(dl), intent(in)  :: am     !< input scale factor times the neutrino mass
            real(dl)              :: adotoa !< input conformal Hubble
            real(dl)              :: rhonu  !< input neutrino density
            real(dl)              :: Nu_drho_Wrapper  !< output value of the time derivative of neutrino density
        end function Nu_drho_Wrapper
        !----------------------------------------------------------------------------------------
        !> Wrapper to the function that computes the background massive neutrinos
        !! time derivative of pressure.
        function Nu_pidot_Wrapper( am, adotoa, presnu )
            use precision
            implicit none
            real(dl), intent(in)  :: am     !< input scale factor times the neutrino mass
            real(dl)              :: adotoa !< input conformal Hubble
            real(dl)              :: presnu !< input neutrino pressure
            real(dl)              :: Nu_pidot_Wrapper !< output value of the time derivative of neutrino pressure
        end function Nu_pidot_Wrapper
        !----------------------------------------------------------------------------------------
        !> Wwrapper to the function that computes the background massive neutrinos
        !! second time derivative of pressure.
        function Nu_pidotdot_Wrapper( am, adotoa, Hdot, presnu, presnudot )
            use precision
            implicit none
            real(dl), intent(in)  :: am        !< input scale factor times the neutrino mass
            real(dl)              :: adotoa    !< input conformal Hubble
            real(dl)              :: Hdot      !< input time derivative of conformal Hubble
            real(dl)              :: presnu    !< input neutrino pressure
            real(dl)              :: presnudot !< input time derivative of neutrino pressure
            real(dl)              :: Nu_pidotdot_Wrapper !< output value of the second time derivative of neutrino pressure
        end function Nu_pidotdot_Wrapper
        !----------------------------------------------------------------------------------------
        !> Wrapper to the function that computes the background massive neutrinos
        !! third time derivative of pressure.
        function Nu_pidotdotdot_Wrapper( am, adotoa, Hdot,Hdotdot, presnu, presnudot, presnudotdot )
            use precision
            implicit none
            real(dl), intent(in)  :: am           !< input scale factor times the neutrino mass
            real(dl)              :: adotoa       !< input conformal Hubble
            real(dl)              :: Hdot         !< input time derivative of conformal Hubble
            real(dl)              :: Hdotdot      !< input time derivative of conformal Hubble
            real(dl)              :: presnu       !< input neutrino pressure
            real(dl)              :: presnudot    !< input time derivative of neutrino pressure
            real(dl)              :: presnudotdot !< input second time derivative of neutrino pressure
            real(dl)              :: Nu_pidotdotdot_Wrapper !< output value of the second time derivative of neutrino pressure
        end function Nu_pidotdotdot_Wrapper
        !----------------------------------------------------------------------------------------
    end interface


    !----------------------------------------------------------------------------------------
    !> This is the type that defines the MGCAMB time step cache.
    type :: MGCAMB_timestep_cache

        ! 1) time and k:
        real(dl) :: a             !< the value of the scale factor at which the cache is being used.
        real(dl) :: tau           !< the value of conformal time at the given scale factor.
        real(dl) :: k             !< the scale that is being solved for. In \f$ \mbox{Mpc}^{-1} \f$.
        ! 2) total matter densities:
        real(dl) :: grhoa2        !< the input value of \f$ \sum_m\rho_m / a^2 m_0^2 \f$.
        real(dl) :: grhom_t       !< the value of \f$ \sum_m\rho_m a^2 /m_0^2 \f$.
        real(dl) :: gpresm_t      !< the value of \f$ \sum_m P_m a^2 /m_0^2 \f$.
        real(dl) :: gpresdotm_t   !< the value of \f$ \sum_m\dot{P}_m a^2 /m_0^2 \f$.
        ! 3) densities and pressure of the various species:
        real(dl) :: grhob_t       !< the value of \f$ \rho_b a^2 / m_0^2 \f$
        real(dl) :: grhoc_t       !< the value of \f$ \rho_{cdm} a^2 / m_0^2 \f$
        real(dl) :: grhor_t       !< the value of \f$ \rho_{\nu} a^2 / m_0^2 \f$
        real(dl) :: grhog_t       !< the value of \f$ \rho_{\gamma} a^2 / m_0^2 \f$
        real(dl) :: grhov_t       !< the value of \f$ \rho_{\Lambda} a^2 / m_0^2 \f$. Used if needed, especially by designer models.
        real(dl) :: gpiv_t        !< the value of \f$ \P_{\Lambda} a^2 / m_0^2 \f$. Used if needed, especially by designer models.
        real(dl) :: grhonu_tot    !< the value of \f$ \sum_\nu \rho_{m\nu} a^2 / m_0^2 \f$
        real(dl) :: gpinu_tot     !< the value of \f$ \sum_\nu P_{m\nu} a^2 / m_0^2 \f$
        real(dl) :: grhonudot_tot !< the value of \f$ \sum_\nu \dot{\rho}_{m\nu} a^2 / m_0^2 \f$
        real(dl) :: gpinudot_tot  !< the value of \f$ \sum_\nu \dot{P}_{m\nu} a^2 / m_0^2 \f$
        ! 4) expansion history:
        real(dl) :: adotoa        !< the value of \f$ \mathcal{H} \f$ at the given scale factor.
        real(dl) :: Hdot          !< the value of \f$ d\mathcal{H} /d \tau \f$ at the given scale factor.
        real(dl) :: Hdotdot       !< the value of \f$ d^2 \mathcal{H} / d \tau^2 \f$ at the given scale factor.

        ! 5) MG functions:
        real(dl) :: MGMuV       !< the value of mu(a,k)
        real(dl) :: MGMuP       !< the value of the derivative w.r.t the scale factor of mu(a,k)
        real(dl) :: MGMuPP      !< the value of the second derivative w.r.t the scale factor of mu(a,k)
        real(dl) :: MGGammaV    !< the value of gamma(a,k)
        real(dl) :: MGGammaP    !< the value of the derivative w.r.t the scale factor of gamma(a,k)
        real(dl) :: MGGammaPP   !< the value of the second derivative w.r.t the scale factor of gamma(a,k)

        !> Alternative Parametrizations
        real(dl) :: MGQV        !< the value of Q(a,k)
        real(dl) :: MGQP        !< the value of the derivative w.r.t the scale factor of Q(a,k)
        real(dl) :: MGQPP       !< the value of the second derivative w.r.t the scale factor of Q(a,k)
        real(dl) :: MGRV        !< the value of R(a,k)
        real(dl) :: MGRP        !< the value of the derivative w.r.t the scale factor of R(a,k)
        real(dl) :: MGRPP       !< the value of the second derivative w.r.t the scale factor of R(a,k)

        ! 6) other background quantities:
        !real(dl) :: grhoq         !< the value of the effective density of the Q field. Refer to the Numerical Notes for the definition.
        !real(dl) :: gpresq        !< the value of the effective pressure of the Q field. Refer to the Numerical Notes for the definition.
        !real(dl) :: grhodotq      !< the value of the time derivative of the effective density of the Q field. Refer to the Numerical Notes for the definition.
        !real(dl) :: gpresdotq     !< the value of the time derivative of the effective pressure of the Q field. Refer to the Numerical Notes for the definition.
        ! 7) the Einstein equations coefficients:
        !real(dl) :: MGeomF       !< the value of the Einstein equations coefficient F. Refer to the Numerical Notes for the definition.
        !real(dl) :: MGeomN       !< the value of the Einstein equations coefficient N. Refer to the Numerical Notes for the definition.
        !real(dl) :: MGeomNdot    !< the value of the Einstein equations coefficient dN/dtau. Refer to the Numerical Notes for the definition.
        !real(dl) :: MGeomX       !< the value of the Einstein equations coefficient X. Refer to the Numerical Notes for the definition.
        !real(dl) :: MGeomXdot    !< the value of the Einstein equations coefficient dX/dtau. Refer to the Numerical Notes for the definition.
        !real(dl) :: MGeomY       !< the value of the Einstein equations coefficient Y. Refer to the Numerical Notes for the definition.
        !real(dl) :: MGeomG       !< the value of the Einstein equations coefficient G. Refer to the Numerical Notes for the definition.
        !real(dl) :: MGeomU       !< the value of the Einstein equations coefficient U. Refer to the Numerical Notes for the definition.
        !real(dl) :: MGeomL       !< the value of the Einstein equations coefficient L. Refer to the Numerical Notes for the definition.
        !real(dl) :: MGeomM       !< the value of the Einstein equations coefficient M. Refer to the Numerical Notes for the definition.
        !real(dl) :: MGeomV       !< the value of the Einstein equations coefficient V. Refer to the Numerical Notes for the definition.
        !real(dl) :: MGeomVdot    !< the value of the Einstein equations coefficient dV/dtau. Refer to the Numerical Notes for the definition.
        !real(dl) :: MGeomQ       !< the value of the Einstein equations coefficient Q. Refer to the Numerical Notes for the definition.\

        ! 8) scalar perturbations quantities:
        real(dl) :: z             !< Syncronous gauge Z perturbation.
        real(dl) :: dz            !< Syncronous gauge dot Z perturbation. This is used to store the non-RSA result.
        real(dl) :: sigma         !< Syncronous gauge sigma perturbation. This is used to store the non-RSA result.
        real(dl) :: sigmadot      !< Syncronous gauge dot sigma perturbation. This is used to store the non-RSA result.
        real(dl) :: clxc          !< Syncronous gauge cdm density perturbation.
        real(dl) :: clxb          !< Syncronous gauge baryon density perturbation.
        real(dl) :: clxg          !< Syncronous gauge radiation density perturbation.
        real(dl) :: clxr          !< Syncronous gauge massless neutrinos density perturbation.
        real(dl) :: vb            !< Syncronous gauge baryon velocity.
        real(dl) :: dgpnu         !< Syncronous gauge massive neutrinos pressure perturbation.
        real(dl) :: dgrho         !< Syncronous gauge total density perturbation.
        real(dl) :: dgq           !< Syncronous gauge total velocity perturbation.
        ! 9) tensor perturbations quantities:
        real(dl) :: MGAT          !< the value of the tensor equation coefficient A. Refer to the Numerical Notes for the definition.
        real(dl) :: MGBT          !< the value of the tensor equation coefficient B. Refer to the Numerical Notes for the definition.
        real(dl) :: MGDT          !< the value of the tensor equation coefficient D. Refer to the Numerical Notes for the definition.
        ! 10) other quantities useful for debug purposes:
        real(dl) :: MGISW         !< Source for ISW effect.
        real(dl) :: MGLensing     !< Source for lensing effect.
        real(dl) :: CMBTSource    !< Full source of CMB temperature fluctuations.
        real(dl) :: Psi           !< Perturbation in the 00 component of the metric in Conformal Newtonian gauge.
        real(dl) :: Phi           !< Perturbation in the space component of the metric in Conformal Newtonian gauge.
        real(dl) :: mu            !< Effective perturbation gravitational constant. \f$ k^2\Psi/\Delta_{m} $\f
        real(dl) :: gamma         !< Ratio between the two gravitational potentials. \f$ \gamma = \Phi/\Psi $\f

    contains

        procedure :: initialize        => MGCAMBTimestepCacheInit      !< subroutine that initializes to zero all the elements of the cache.
        procedure :: is_nan            => MGCAMBTimestepCacheIsNan     !< Subroutine that check if an element of the MGCAMB_timestep_cache is Nan.
        procedure :: open_cache_files  => MGCAMBTimestepCacheOpenFile  !< subroutine that opens the files to dump the cache to file.
        procedure :: dump_cache_files  => MGCAMBTimestepCacheDumpFile  !< subroutine that dumps the cache to the files.
        procedure :: close_cache_files => MGCAMBTimestepCacheCloseFile !< subroutine that closes the files where the cache has benn dumped.

    end type MGCAMB_timestep_cache

contains

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes to zero all the elements of the cache.
    subroutine MGCAMBTimestepCacheInit( self )

        implicit none

        class(MGCAMB_timestep_cache)  :: self !< the base class.

        ! initialize all class members to zero:
        ! 1) time and k:
        self%a             = 0._dl
        self%tau           = 0._dl
        self%k             = 0._dl
        ! 2) expansion history:
        self%adotoa        = 0._dl
        self%Hdot          = 0._dl
        self%Hdotdot       = 0._dl
        ! 3) total matter densities:
        self%grhoa2        = 0._dl
        self%grhom_t       = 0._dl
        self%gpresm_t      = 0._dl
        self%gpresdotm_t   = 0._dl
        ! 4) densities and pressure of the various species:
        self%grhob_t       = 0._dl
        self%grhoc_t       = 0._dl
        self%grhor_t       = 0._dl
        self%grhog_t       = 0._dl
        self%grhov_t       = 0._dl
        self%gpiv_t        = 0._dl
        self%grhonu_tot    = 0._dl
        self%gpinu_tot     = 0._dl
        self%grhonudot_tot = 0._dl
        self%gpinudot_tot  = 0._dl
        ! 5) MG functions:
        self%MGMuV          = 0._dl
        selff%MGMuP         = 0._dl
        self%MGMuPP         = 0._dl
        self%MGGammaV       = 0._dl
        selff%MGGammaP      = 0._dl
        self%MGGammaPP      = 0._dl
        self%MGQV           = 0._dl
        selff%MGQP          = 0._dl
        self%MGQPP          = 0._dl
        self%MGRV           = 0._dl
        selff%MGRP          = 0._dl
        self%MGRPP          = 0._dl
        ! 6) other background quantities:
        !self%grhoq         = 0._dl
        !self%gpresq        = 0._dl
        !self%grhodotq      = 0._dl
        !self%gpresdotq     = 0._dl
        ! 7) the Einstein equations coefficients:
        !self%MGeomF       = 0._dl
        !self%MGeomN       = 0._dl
        !self%MGeomNdot    = 0._dl
        !self%MGeomX       = 0._dl
        !self%MGeomXdot    = 0._dl
        !self%MGeomY       = 0._dl
        !self%MGeomG       = 0._dl
        !self%MGeomU       = 0._dl
        !self%MGeomL       = 0._dl
        !self%MGeomM       = 0._dl
        !self%MGeomV       = 0._dl
        !self%MGeomVdot    = 0._dl
        !self%MGeomQ       = 0._dl
        ! 10) perturbations quantities:
        self%z             = 0._dl
        self%dz            = 0._dl
        self%sigma         = 0._dl
        self%sigmadot      = 0._dl
        self%clxc          = 0._dl
        self%clxb          = 0._dl
        self%clxg          = 0._dl
        self%clxr          = 0._dl
        self%vb            = 0._dl
        self%dgpnu         = 0._dl
        self%dgrho         = 0._dl
        self%dgq           = 0._dl
        ! 11) tensor perturbations quantities:
        self%MGAT         = 0._dl
        self%MGBT         = 0._dl
        self%MGDT         = 0._dl
        ! 13) other quantities usefull for debug purposes:
        self%MGISW        = 0._dl
        self%MGLensing    = 0._dl
        self%CMBTSource    = 0._dl
        self%Psi           = 0._dl
        self%Phi           = 0._dl
        self%mu            = 0._dl
        self%gamma         = 0._dl

    end subroutine MGCAMBTimestepCacheInit

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that check if an element of the MGCAMB_timestep_cache is Nan.
    subroutine MGCAMBTimestepCacheIsNan( self, HaveNan )

        implicit none

        class(MGCAMB_timestep_cache), intent(in)  :: self    !< The base class.
        logical, intent(out)                       :: HaveNan !< Logical variable which describes the presence of a Nan variable.
                                                              !< If an element of the MGCAMB_timestep_cache is Nan, you get HaveNan=.True.

        HaveNan = .False.
        HaveNan = HaveNan.or.IsNaN(self%a)
        HaveNan = HaveNan.or.IsNaN(self%tau)
        HaveNan = HaveNan.or.IsNaN(self%k)
        HaveNan = HaveNan.or.IsNaN(self%adotoa)
        HaveNan = HaveNan.or.IsNaN(self%Hdot)
        HaveNan = HaveNan.or.IsNaN(self%Hdotdot)
        HaveNan = HaveNan.or.IsNaN(self%grhoa2)
        HaveNan = HaveNan.or.IsNaN(self%grhom_t)
        HaveNan = HaveNan.or.IsNaN(self%gpresm_t)
        HaveNan = HaveNan.or.IsNaN(self%gpresdotm_t)
        HaveNan = HaveNan.or.IsNaN(self%grhob_t)
        HaveNan = HaveNan.or.IsNaN(self%grhoc_t)
        HaveNan = HaveNan.or.IsNaN(self%grhor_t)
        HaveNan = HaveNan.or.IsNaN(self%grhog_t)
        HaveNan = HaveNan.or.IsNaN(self%grhov_t)
        HaveNan = HaveNan.or.IsNaN(self%gpiv_t)
        HaveNan = HaveNan.or.IsNaN(self%grhonu_tot)
        HaveNan = HaveNan.or.IsNaN(self%gpinu_tot)
        HaveNan = HaveNan.or.IsNaN(self%grhonudot_tot)
        HaveNan = HaveNan.or.IsNaN(self%gpinudot_tot)

        !> MG functions
        HaveNan = HaveNan.or.IsNaN(self%MGMuV)
        HaveNan = HaveNan.or.IsNaN(self%MGMuP)
        HaveNan = HaveNan.or.IsNaN(self%MGMuPP)
        HaveNan = HaveNan.or.IsNaN(self%MGGammaV)
        HaveNan = HaveNan.or.IsNaN(self%MGGammaP)
        HaveNan = HaveNan.or.IsNaN(self%MGGammaPP)
        HaveNan = HaveNan.or.IsNaN(self%MGQV)
        HaveNan = HaveNan.or.IsNaN(self%MGQP)
        HaveNan = HaveNan.or.IsNaN(self%MGQPP)
        HaveNan = HaveNan.or.IsNaN(self%MGRV)
        HaveNan = HaveNan.or.IsNaN(self%MGRP)
        HaveNan = HaveNan.or.IsNaN(self%MGRPP)

        !> other quantities
        !HaveNan = HaveNan.or.IsNaN(self%grhoq)
        !HaveNan = HaveNan.or.IsNaN(self%gpresq)
        !HaveNan = HaveNan.or.IsNaN(self%grhodotq)
        !HaveNan = HaveNan.or.IsNaN(self%gpresdotq)

        !HaveNan = HaveNan.or.IsNaN(self%MGeomF)
        !HaveNan = HaveNan.or.IsNaN(self%MGeomN)
        !HaveNan = HaveNan.or.IsNaN(self%MGeomNdot)
        !HaveNan = HaveNan.or.IsNaN(self%MGeomX)
        !HaveNan = HaveNan.or.IsNaN(self%MGeomXdot)
        !HaveNan = HaveNan.or.IsNaN(self%MGeomY)
        !HaveNan = HaveNan.or.IsNaN(self%MGeomG)
        !HaveNan = HaveNan.or.IsNaN(self%MGeomU)
        !HaveNan = HaveNan.or.IsNaN(self%MGeomL)
        !HaveNan = HaveNan.or.IsNaN(self%MGeomM)
        !HaveNan = HaveNan.or.IsNaN(self%MGeomV)
        !HaveNan = HaveNan.or.IsNaN(self%MGeomVdot)
        !HaveNan = HaveNan.or.IsNaN(self%MGeomQ)
        !> perturbation quantities
        HaveNan = HaveNan.or.IsNaN(self%z)
        HaveNan = HaveNan.or.IsNaN(self%dz)
        HaveNan = HaveNan.or.IsNaN(self%sigma)
        HaveNan = HaveNan.or.IsNaN(self%sigmadot)
        HaveNan = HaveNan.or.IsNaN(self%clxc)
        HaveNan = HaveNan.or.IsNaN(self%clxb)
        HaveNan = HaveNan.or.IsNaN(self%clxg)
        HaveNan = HaveNan.or.IsNaN(self%clxr)
        HaveNan = HaveNan.or.IsNaN(self%vb)
        HaveNan = HaveNan.or.IsNaN(self%dgpnu)
        HaveNan = HaveNan.or.IsNaN(self%dgrho)
        HaveNan = HaveNan.or.IsNaN(self%dgq)

        !> Tensor perturbations
        !HaveNan = HaveNan.or.IsNaN(self%MGAT)
        !HaveNan = HaveNan.or.IsNaN(self%MGBT)
        !HaveNan = HaveNan.or.IsNaN(self%MGDT)

        !> Debug quantities
        HaveNan = HaveNan.or.IsNaN(self%MGISW)
        HaveNan = HaveNan.or.IsNaN(self%MGLensing)
        HaveNan = HaveNan.or.IsNaN(self%CMBTSource)
        HaveNan = HaveNan.or.IsNaN(self%Psi)
        HaveNan = HaveNan.or.IsNaN(self%Phi)
        HaveNan = HaveNan.or.IsNaN(self%mu)
        HaveNan = HaveNan.or.IsNaN(self%gamma)

    end subroutine MGCAMBTimestepCacheIsNan

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that opens the files to dump the cache to file.
    subroutine MGCAMBTimestepCacheOpenFile( self, outroot )

        implicit none

        class(MGCAMB_timestep_cache)  :: self    !< the base class.
        character(len=*), intent(in)   :: outroot !< output root of the file.

        logical :: is_open

        ! print some feedback:
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') 'MGCAMB cache opening print files in:'
        write(*,'(a)')  outroot
        write(*,'(a)') "***************************************************************"

        ! test whether the files are already open:
        call test_open( 111  )
        call test_open( 222  )
        call test_open( 333  )
        call test_open( 444  )
        call test_open( 555  )
        call test_open( 666  )
        call test_open( 777  )
        call test_open( 888  )
        call test_open( 999  )
        call test_open( 1111 )
        call test_open( 2222 )
        call test_open( 3333 )
        call test_open( 4444 )
        call test_open( 5555 )

        ! open the files:
        call CreateTxtFile( TRIM(outroot)//'cache_FRW.dat'           ,111 )
        call CreateTxtFile( TRIM(outroot)//'cache_BDens.dat'         ,222 )
        call CreateTxtFile( TRIM(outroot)//'cache_BPres.dat'         ,333 )
        call CreateTxtFile( TRIM(outroot)//'cache_BOmegas.dat'       ,444 )
        call CreateTxtFile( TRIM(outroot)//'cache_BackgroundMG.dat'  ,555 )
        !call CreateTxtFile( TRIM(outroot)//'cache_BackgroundQ.dat'   ,777 )
        !call CreateTxtFile( TRIM(outroot)//'cache_EinsteinCoeff.dat' ,888 )
        call CreateTxtFile( TRIM(outroot)//'cache_EinsteinSol.dat'   ,2222)
        call CreateTxtFile( TRIM(outroot)//'cache_TensorCoeff.dat'   ,3333)
        call CreateTxtFile( TRIM(outroot)//'cache_Sources.dat'       ,4444)
        call CreateTxtFile( TRIM(outroot)//'cache_MetricMG.dat'      ,5555)

        ! write the headers:
        write (111 ,'(12a)')  '# ', 'a ', 'tau ', 'k ', 'adotoa ', 'Hdot ', 'Hdotdot '
        write (222 ,'(12a)')  '# ', 'a ', 'tau ', 'k ', 'grhom_t ', 'grhob_t ', 'grhoc_t ', 'grhor_t ', 'grhog_t ', 'grhov_t ', 'grhonu_tot ', 'grhonudot_tot '
        write (333 ,'(12a)')  '# ', 'a ', 'tau ', 'k ', 'gpresm_t ', 'gpresdotm_t ', 'gpiv_t ', 'gpinu_tot ', 'gpinudot_tot '
        write (444 ,'(12a)')  '# ', 'a ', 'tau ', 'k ', 'omegam ', 'omegab ', 'omegac ', 'omegar ', 'omegag ', 'omegav ', 'omeganu_tot '
        write (555 ,'(20a)')  '# ', 'a ', 'tau ', 'k ', 'MGMuV', 'MGMuP', 'MGMuPP', 'MGGammaV', 'MGGammaP', 'MGGammaPP', 'MGQV', 'MGQP','MGQPP', 'MGRV', 'MGRP', 'MGRPP'
        !write (777 ,'(12a)')  '# ', 'a ', 'tau ', 'k ', 'grhoq ', 'gpresq ', 'grhodotq ', 'gpresdotq '
        !write (888 ,'(18a)')  '# ', 'a ', 'tau ', 'k ', 'MGeomF ', 'MGeomN ', 'MGeomNdot ', 'MGeomX ', 'MGeomXdot ', 'MGeomY ', 'MGeomG ', 'MGeomU ', 'MGeomL ', 'MGeomM ', 'MGeomV ', 'MGeomVdot ', 'MGeomQ'
        write (2222,'(20a)')  '# ', 'a ', 'tau ', 'k ', 'z ', 'sigma', 'clxc ', 'clxb ', 'clxg ', 'clxr ', 'vb ', 'dgpnu ', 'dgrho ', 'dgq '
        write (3333,'(12a)')  '# ', 'a ', 'tau ', 'k ', 'MGAT ', 'MGBT ', 'MGDT '
        write (4444,'(12a)')  '# ', 'a ', 'tau ', 'k ', 'MGISW ', 'MGLensing ', 'CMBTSource '
        write (5555,'(12a)')  '# ', 'a ', 'tau ', 'k ', 'Psi ', 'Phi ', 'mu ', 'gamma '

    contains

        ! Temporary subroutine that tests wether a cache file is open.
        subroutine test_open( number )

            implicit none

            integer, intent(in) :: number
            logical             :: is_open

            inquire( unit=number, opened=is_open )
            if ( is_open ) then
                write(*,*) 'MGCAMB ERROR: Oputput unit', number, 'is already open.'
                write(*,*) 'MGCAMB cannot use it and cannot proceed.'
                call MpiStop('MGCAMB error')
            end if

        end subroutine test_open

    end subroutine MGCAMBTimestepCacheOpenFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that closes the files where the cache has benn dumped.
    subroutine MGCAMBTimestepCacheCloseFile( self )

        implicit none

        class(MGCAMB_timestep_cache)  :: self !< the base class.

        ! test wether the files can be closed:
        call test_close( 111  )
        call test_close( 222  )
        call test_close( 333  )
        call test_close( 444  )
        call test_close( 555  )
        call test_close( 666  )
        call test_close( 777  )
        call test_close( 888  )
        call test_close( 999  )
        call test_close( 1111 )
        call test_close( 2222 )
        call test_close( 3333 )
        call test_close( 4444 )
        call test_close( 5555 )

        ! close the files:
        close( 111  )
        close( 222  )
        close( 333  )
        close( 444  )
        close( 555  )
        close( 666  )
        close( 777  )
        close( 888  )
        close( 999  )
        close( 1111 )
        close( 2222 )
        close( 3333 )
        close( 4444 )
        close( 5555 )

        ! print some feedback:
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') 'MGCAMB cache printing done.'
        write(*,'(a)') "***************************************************************"

    contains

        ! Temporary subroutine that tests wether a cahce file is open.
        subroutine test_close( number )

            implicit none

            integer, intent(in) :: number
            logical             :: is_open

            inquire( unit=number, opened=is_open )
            if ( .not. is_open ) then
                write(*,*) 'MGCAMB ERROR: Oputput unit', number, 'is not open.'
                write(*,*) 'MGCAMB is trying to close it and cannot proceed.'
                call MpiStop('MGCAMB error')
            end if

        end subroutine test_close

    end subroutine MGCAMBTimestepCacheCloseFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that dumps the cache to the files.
    subroutine MGCAMBTimestepCacheDumpFile( self )

        implicit none

        class(MGCAMB_timestep_cache)  :: self !< the base class.

        ! write the background expansion history:
        write (111 ,'(12'//cache_output_format//')')  self%a, self%tau, self%k, self%adotoa, self%Hdot, self%Hdotdot
        ! write the background densities:
        write (222 ,'(14'//cache_output_format//')')  self%a, self%tau, self%k, self%grhom_t, self%grhob_t, self%grhoc_t, self%grhor_t, self%grhog_t, self%grhov_t, self%grhonu_tot, self%grhonudot_tot
        ! write the background pressure:
        write (333 ,'(12'//cache_output_format//')')  self%a, self%tau, self%k, self%gpresm_t, self%gpresdotm_t, self%gpiv_t, self%gpinu_tot, self%gpinudot_tot
        ! write the background omegas:
        write (444 ,'(14'//cache_output_format//')')  self%a, self%tau, self%k, self%grhom_t/(3._dl*self%adotoa**2), self%grhob_t/(3._dl*self%adotoa**2), self%grhoc_t/(3._dl*self%adotoa**2), self%grhor_t/(3._dl*self%adotoa**2), self%grhog_t/(3._dl*self%adotoa**2), self%grhov_t/(3._dl*self%adotoa**2), self%grhonu_tot/(3._dl*self%adotoa**2)
        ! write background MG functions:
        write (555 ,'(14'//cache_output_format//')')  self%a, self%tau, self%k, self%MGMuV, self%MGMuP, self%MGMuPP, self%MGGammaV, self%MGGammaP, self%MGGammaPP, self%MGQV, self%MGQP, self%MGQPP, self%MGRV, self%MGRP, self%MGRPP
        ! write some perturbations:
        write (2222,'(20'//cache_output_format//')')  self%a, self%tau, self%k, self%z, self%sigma, self%clxc, self%clxb, self%clxg, self%clxr, self%dgpnu, self%dgrho, self%dgq
        ! write tensor coefficients:
        write (3333,'(12'//cache_output_format//')')  self%a, self%tau, self%k, self%MGAT, self%MGBT, self%MGDT
        ! write sources:
        write (4444,'(12'//cache_output_format//')')  self%a, self%tau, self%k, self%MGISW, self%MGLensing, self%CMBTSource
        ! write metric potentials:
        write (5555,'(12'//cache_output_format//')')  self%a, self%tau, self%k, self%Psi, self%Phi, self%mu, self%gamma

    end subroutine MGCAMBTimestepCacheDumpFile

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that initializes to zero all the elements of the parameter cache.
    subroutine MGCAMBParameterCacheInit( self )

        implicit none

        class(MGCAMB_parameter_cache)  :: self !< the base class.

        ! initialize all class members to zero:
        self%omegac      = 0._dl
        self%omegab      = 0._dl
        self%omegav      = 0._dl
        self%omegak      = 0._dl
        self%omegan      = 0._dl
        self%omegag      = 0._dl
        self%omegar      = 0._dl
        self%h0          = 0._dl
        self%h0_Mpc      = 0._dl
        self%grhog       = 0._dl
        self%grhornomass = 0._dl
        self%grhoc       = 0._dl
        self%grhob       = 0._dl
        self%grhov       = 0._dl
        self%grhok       = 0._dl
        self%Num_Nu_Massive       = 0
        self%Nu_mass_eigenstates  = 0
        if ( allocated(self%grhormass)         ) deallocate(self%grhormass)
        if ( allocated(self%nu_masses)         ) deallocate(self%nu_masses)
        if ( associated(self%Nu_background)    ) nullify(self%Nu_background)
        if ( associated(self%Nu_rho)           ) nullify(self%Nu_rho)
        if ( associated(self%Nu_pidot)         ) nullify(self%Nu_pidot)
        if ( associated(self%Nu_pidotdot)      ) nullify(self%Nu_pidotdot)
        if ( associated(self%Nu_pidotdotdot)   ) nullify(self%Nu_pidotdotdot)

    end subroutine MGCAMBParameterCacheInit

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints the MGCAMB parameters cache to screen.
    subroutine MGCAMBParameterCachePrint( self )

        implicit none

        class(MGCAMB_parameter_cache)  :: self !< the base class.

        integer  :: i

        ! print to screen the parameter cache:
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') 'MGCAMB parameters cache content:'
        write(*,'(a)') "***************************************************************"
        write(*,'(a14,E13.6)') ' Omega_CDM  : ', self%omegac
        write(*,'(a14,E13.6)') ' Omega_b    : ', self%omegab
        write(*,'(a14,E13.6)') ' Omega_v    : ', self%omegav
        write(*,'(a14,E13.6)') ' Omega_k    : ', self%omegak
        write(*,'(a14,E13.6)') ' Omega_n    : ', self%omegan
        write(*,'(a14,E13.6)') ' Omega_g    : ', self%omegag
        write(*,'(a14,E13.6)') ' Omega_r    : ', self%omegar
        write(*,'(a14,F12.6)') ' h          : ', self%h0
        write(*,'(a14,E13.6)') ' h_Mpc      : ', self%h0_Mpc
        write(*,'(a14,E13.6)') ' grhog      : ', self%grhog
        write(*,'(a14,E13.6)') ' grnonomass : ', self%grhornomass
        write(*,'(a14,E13.6)') ' grhoc      : ', self%grhoc
        write(*,'(a14,E13.6)') ' grhob      : ', self%grhob
        write(*,'(a14,E13.6)') ' grhov      : ', self%grhov
        write(*,'(a14,E13.6)') ' grhok      : ', self%grhok
        write(*,'(a22,I10)') ' Num_Nu_Massive      : ', self%Num_Nu_Massive
        write(*,'(a22,I10)') ' Nu_mass_eigenstates : ', self%Nu_mass_eigenstates
        do i=1, self%Nu_mass_eigenstates
            write(*,'(a11,I3,a9,E13.6)') ' grhormass(',i,')      : ', self%grhormass(i)
            write(*,'(a11,I3,a9,E13.6)') ' nu_masses(',i,')      : ', self%nu_masses(i)
        end do
        write(*,'(a)') "***************************************************************"

    end subroutine MGCAMBParameterCachePrint

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that check if an element of the MGCAMB_parameter_cache is Nan.
    subroutine MGCAMBParameterCacheIsNan( self, HaveNan )

        implicit none

        class(MGCAMB_parameter_cache), intent(in)  :: self    !< The base class.
        logical, intent(out)                        :: HaveNan !< Logical variable which describes the presence of a Nan variable.
                                                               !< If an element of the MGCAMB_parameter_cache is Nan, you get HaveNan=.True.

        integer                                     :: i

        HaveNan = .False.
        HaveNan = HaveNan.or.IsNaN(self%omegac)
        HaveNan = HaveNan.or.IsNaN(self%omegab)
        HaveNan = HaveNan.or.IsNaN(self%omegav)
        HaveNan = HaveNan.or.IsNaN(self%omegak)
        HaveNan = HaveNan.or.IsNaN(self%omegan)
        HaveNan = HaveNan.or.IsNaN(self%omegag)
        HaveNan = HaveNan.or.IsNaN(self%omegar)
        HaveNan = HaveNan.or.IsNaN(self%h0)
        HaveNan = HaveNan.or.IsNaN(self%h0_Mpc)
        HaveNan = HaveNan.or.IsNaN(self%grhog)
        HaveNan = HaveNan.or.IsNaN(self%grhornomass)
        HaveNan = HaveNan.or.IsNaN(self%grhoc)
        HaveNan = HaveNan.or.IsNaN(self%grhob)
        HaveNan = HaveNan.or.IsNaN(self%grhov)
        HaveNan = HaveNan.or.IsNaN(self%grhok)
        HaveNan = HaveNan.or.IsNaN(self%Num_Nu_Massive*1.0_dl)
        HaveNan = HaveNan.or.IsNaN(self%Nu_mass_eigenstates*1.0)

        do i=1, self%Nu_mass_eigenstates
            HaveNan = HaveNan.or.IsNaN(self%grhormass(i)).or.IsNaN(self%nu_masses(i))
        end do

    end subroutine MGCAMBParameterCacheIsNan

    ! ---------------------------------------------------------------------------------------------

end module MGCAMB_cache

!----------------------------------------------------------------------------------------
