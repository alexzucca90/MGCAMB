module RD_MockMaker

!----------------------------------------------!
!MOCK ESTIMATOR FOR REDSHIFT DRIFT (MERD) v1.0 !
!                  29/07/2015                  !
!----------------------------------------------!


use settings

implicit none


type RDspec
     integer :: num_red
     integer :: Nsource !number of sources in each redshift bin
     real(mcp) :: SN_ratio
     real(mcp) :: deltat_years
     real(mcp), allocatable, dimension(:) :: zbins
end type RDspec


!FIND ANOTHER SOLUTION
integer :: intra_num_red
integer :: intra_Nsource
real(mcp) :: intra_SN_ratio
real(mcp) :: intra_deltat_years
real(mcp), allocatable, dimension(:) :: intra_zbins


contains

subroutine GetSpecs(mock_type,mock_root,SP)

implicit none
integer :: mock_type
character(100) :: mock_root
integer :: nz,ib
Type(RDspec) :: SP
real(mcp) :: zinf,zsup,sig_z
real(mcp), allocatable, dimension(:) :: win_z, win_tot
real(mcp), allocatable, dimension(:,:) :: win


if (mock_type.eq.1) then
   write(0,*) 'Using E-ELT specifications'

   !SPECIFICATIONS
   SP%num_red = 5
   SP%SN_ratio= 3000
   SP%Nsource = 6
   SP%deltat_years= 30

   allocate(SP%zbins(SP%num_red),intra_zbins(SP%num_red))


   SP%zbins=[2.0, 2.8, 3.5, 4.2, 5.0]


else if (mock_type.eq.2) then
   write(0,*) 'Using SKA specifications'
   write(0,*)'Not available yet, sorry! Try again later'
   stop
   !SPECIFICATIONS
   SP%num_red = 5
   SP%SN_ratio= 3000
   SP%Nsource = 6
   SP%deltat_years= 30

   allocate(SP%zbins(SP%num_red),intra_zbins(SP%num_red))


   SP%zbins=[2.0, 2.8, 3.5, 4.2, 5.0]

else if (mock_type.eq.3) then
   write(0,*) 'Using CHIME specifications'
   write(0,*)'Not available yet, sorry! Try again later'
   stop
   !SPECIFICATIONS
   SP%num_red = 5
   SP%SN_ratio= 3000
   SP%Nsource = 6
   SP%deltat_years= 30

   allocate(SP%zbins(SP%num_red),intra_zbins(SP%num_red))


   SP%zbins=[2.0, 2.8, 3.5, 4.2, 5.0]

else if (mock_type.eq.4) then
   write(0,*) 'Using SKA second derivative specifications'
   write(0,*)'Not available yet, sorry! Try again later'

stop
   !SPECIFICATIONS
   SP%num_red = 8
   SP%SN_ratio= 3000
   SP%Nsource = 6
   SP%deltat_years= 30

   allocate(SP%zbins(SP%num_red),intra_zbins(SP%num_red))


   SP%zbins=[0.5,1.0,1.5,2.0, 2.5, 3.0, 3.5,4.0]

else
   write(0,*) 'No option for this, what do you want?'
   stop
end if

   intra_num_red = SP%num_red
   intra_SN_ratio = SP%SN_ratio
   intra_Nsource = SP%Nsource
   intra_deltat_years = SP%deltat_years

   intra_zbins=SP%zbins




10 format(500(1e12.6,1x))

end subroutine GetSpecs

subroutine rd_mock_maker(mock_root, mock_type, RD_theory)

implicit none
character(100) :: mock_root
integer :: mock_type
real(mcp), dimension(intra_num_red) :: RD_theory,sigma
integer :: i
real(mcp), parameter :: pi=3.14159


!compute errors -----------------------------------------------
if (mock_type.eq.1) then
   write(0,*)'Computing E-ELT errors...'
   do i=1,intra_num_red
      if (intra_zbins(i).le.4.0) then
         sigma(i)=1.35*(2370.d0/intra_SN_ratio)*sqrt(30.d0/intra_Nsource)*(5/(1+intra_zbins(i)))**1.7d0!RD errors E-ELT
      else
         sigma(i)=1.35*(2370.d0/intra_SN_ratio)*sqrt(30.d0/intra_Nsource)*(5/(1+intra_zbins(i)))**0.9d0
      end if
   end do
else if (mock_type.eq.2) then
   write(0,*)'Computing SKA errors...'
   write(0,*)'Not available yet, sorry! Try again later'
   stop
else if (mock_type.eq.3) then
   write(0,*)'Computing CHIME errors...'
   write(0,*)'Not available yet, sorry! Try again later'
   stop
else if (mock_type.eq.4) then
   write(0,*)'Computing SKA second derivative errors...'
   write(0,*)'Not available yet, sorry! Try again later'
   stop
else
    write(0,*) 'No option for this... how did you even get here?'
    stop
end if


write(0,*) '...done!'
!--------------------------------------------------------------




!writing mock and dataset files--------------------------------
write(0,*) 'Writing output files'
open(42,file='mock_output/'//trim(mock_root)//'.dataset')
write(42,'(A,A)') 'name = ',trim(mock_root)
write(42,'(A,I5)')'num_red_rd =',intra_num_red
write(42,'(A,A)')'RD_file =','%DATASETDIR%'//trim(mock_root)//'_mock.dat'
close(42)

open(42,file='mock_output/'//trim(mock_root)//'_mock.dat')
do i=1,intra_num_red
   write(42,10) intra_zbins(i),RD_theory(i),sigma(i),intra_deltat_years
end do
close(42)

!---------------------------------------------------------------

   write(0,*) 'FINISHED!'
   write(0,*) 'produced files in output folder'
   write(0,*) trim(mock_root)//'.dataset'
   write(0,*) trim(mock_root)//'_mock.dat'






write(0,*)'BELLA PE TUTTI!'



10 format(500(1e12.6,1x))

deallocate(intra_zbins)
end subroutine rd_mock_maker


end module RD_MockMaker
