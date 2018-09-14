!> @file 01_MG_def.f90
!! This file contains the EFTCAMB compile time options.
!> @author: Alex Zucca azucca@sfu.ca

module MG_def

    use Precision

    implicit none

    character(LEN=*), parameter :: MGCAMB_version = 'V3.0 Mar18'

    real(dl), parameter :: MGbackgroundcutoff = 0._dl !< Smallest scale factor that the code should
    !!    consider when copmputing the background. Set to zero, change background at all times.


    real(dl), parameter :: MGtoGR = 1.d-8  !< Return to GR flag:
    !real(dl), parameter :: MGtoGR = 1.d-2 !< Return to GR flag:
    !!    This is the threshold at which a theory is considered to be exactly GR.

    integer , parameter :: MG_names_max_length       = 20    !< maximum length of names for MG functions and parameters.
    integer , parameter :: MG_names_latex_max_length = 40    !< maximum length of latex names for MG functions and parameters.

    integer , parameter :: MG_RGR_num_points   = 1000        !< number of points to sample logaritmically the time in the return to GR module.

    #ifdef DEBUG
    logical , parameter :: DebugMGCAMB = .true.              !< MGCAMB debug flag.This will turn on printing of many things to aid debugging the code.
    #else
    logical , parameter :: DebugMGCAMB = .false.             !< MGCAMB debug flag.This will turn on printing of many things to aid debugging the code.
    #endif

end module MG_def
