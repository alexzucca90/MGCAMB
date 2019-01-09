module WL_MockMaker

!----------------------------------------------!
!MOCK ESTIMATOR FOR SHEAR SPECTRA (MESS) v3.0  !
!                  07/02/2016                  !
!----------------------------------------------!


use settings

implicit none

integer,parameter :: num_z_win=1000

type WLspec
     integer :: num_z_temp
     integer :: lmax
     integer :: lmin
     integer :: numbins
     integer :: numcross
     real(mcp) :: fsky
     real(mcp) :: d2
     real(mcp) :: ng
     real(mcp) :: z0
     real(mcp) :: zmax
     real(mcp) :: zmin
     real(mcp) :: kmax
     real(mcp) :: sig_z0

     real(mcp), allocatable, dimension(:) :: zbins,frac
end type WLspec


!FIND ANOTHER SOLUTION
integer :: intra_num_z_temp
integer :: intra_lmax
integer :: intra_lmin
integer :: intra_numbins
integer :: intra_numcross
real(mcp) :: intra_fsky
real(mcp) :: intra_d2
real(mcp) :: intra_ng
real(mcp) :: intra_z0
real(mcp) :: intra_zmax
real(mcp) :: intra_zmin
real(mcp) :: intra_kmax
real(mcp) :: intra_sig_z0
real(mcp), allocatable, dimension(:) :: intra_zbins,intra_frac
real(mcp), allocatable, dimension(:) :: win_z, win_tot, totwin2
real(mcp), allocatable, dimension(:,:) :: win

contains

subroutine GetSpecs(mock_type,mock_root,bin_type,SP)

implicit none
integer :: mock_type,bin_type
character(100) :: mock_root
integer :: nz,ib,nz2
Type(WLspec) :: SP
real(mcp) :: zinf,zsup,sig_z,integral,totwin
real(mcp), allocatable, dimension(:,:) :: win_temp
!real(mcp), dimension(num_z_win) :: diff
real(mcp) :: diff




if (mock_type.eq.1) then
   write(0,*) 'Using Euclid specifications'

   !SPECIFICATIONS
   SP%num_z_temp=num_z_win
   SP%numbins = 10
   SP%numcross = (SP%numbins*(SP%numbins+1)/2)-SP%numbins
   SP%lmin = 30
   SP%lmax = 2000
   SP%fsky=0.38
   SP%d2=2*0.0484
   SP%ng=30.
   SP%z0=0.9
   SP%zmax=2.d0
   SP%zmin=0.d0
   SP%kmax=20.d0
   SP%sig_z0=0.05


else if (mock_type.eq.2) then

   write(0,*) 'Using Code Comparison specifications'

   !SPECIFICATIONS
   SP%num_z_temp=num_z_win
   SP%numbins = 3
   SP%numcross = (SP%numbins*(SP%numbins+1)/2)-SP%numbins
   SP%lmin = 10
   SP%lmax = 3000
   SP%fsky=0.363
   SP%d2=2*0.0484
   SP%ng=30.
   SP%z0=0.9
   SP%zmax=2.5d0
   SP%zmin=0.d0
   SP%kmax=20.d0
   SP%sig_z0=0.001!0.05


else if (mock_type.eq.3) then
   write(0,*) 'Using LSST specifications'
   write(0,*) 'This will be done eventually'
   write(0,*) 'now it is just the same Euclid thing but with Santiago specs'



   !SPECIFICATIONS
!   SP%num_z_temp=num_z_win
!   SP%numbins = 10
!   SP%numcross = (SP%numbins*(SP%numbins+1)/2)-SP%numbins
!   SP%lmin = 2
!   SP%lmax = 2000
!   SP%fsky= 20000/(4*pi*(180/pi)**2.) 
!   SP%d2=2*0.0484
!   SP%ng=35.
!   SP%z0=0.9
!   SP%zmax=3.0
!   SP%zmin=0.001d0
!   SP%kmax=1.d0
!   SP%sig_z0=0.05
   SP%num_z_temp=num_z_win
   SP%numbins = 10
   SP%numcross = (SP%numbins*(SP%numbins+1)/2)-SP%numbins
   SP%lmin = 10
   SP%lmax = 3000
   SP%fsky=0.38
   SP%d2=2*0.0484
   SP%ng=30.
   SP%z0=0.9
   SP%zmax=2.5d0
   SP%zmin=0.d0
   SP%kmax=1.d0
   SP%sig_z0=0.05


else
   write(0,*) 'Do not be stupid!'
   stop
end if

   intra_num_z_temp=SP%num_z_temp
   intra_lmax=SP%lmax
   intra_lmin=SP%lmin
   intra_numbins=SP%numbins
   intra_numcross=SP%numcross
   intra_fsky=SP%fsky
   intra_d2=SP%d2
   intra_ng=SP%ng
   intra_z0=SP%z0
   intra_zmax=SP%zmax
   intra_zmin=SP%zmin
   intra_kmax=SP%kmax
   intra_sig_z0=SP%sig_z0


   allocate(SP%zbins(SP%numbins+1),SP%frac(SP%numbins),intra_zbins(SP%numbins+1),intra_frac(SP%numbins))




   !COMPUTING WINDOW FUNCTIONS
   allocate(win_z(num_z_win),win(num_z_win,SP%numbins),win_tot(num_z_win),win_temp(num_z_win,SP%numbins),totwin2(SP%numbins))
   do nz=1,num_z_win
      win_z(nz)=SP%zmin+(SP%zmax-SP%zmin)/(num_z_win)*(nz-1)
      if (mock_type.eq.2) then
         win_tot(nz)=(1.50/0.20)*exp(-(win_z(nz)-0.7)**2./0.32**2.)+exp(-(win_z(nz)-1.2)**2./0.46**2.)
      else
         win_tot(nz)=(win_z(nz)**2*exp(-(win_z(nz)/(SP%z0/1.412))**1.5))
      end if
   end do

   totwin=0
   do nz=2,num_z_win
      if ((win_z(nz-1).ge.SP%zmin).and.(win_z(nz).le.SP%zmax)) then
         totwin=totwin+0.5*(win_tot(nz-1)+win_tot(nz))*(win_z(nz)-win_z(nz-1))
      end if
   end do

   open(42,file='mock_file_forplot/'//trim(mock_root)//'_window_total.dat')
   do nz=1,num_z_win
      write(42,10) win_z(nz),win_tot(nz)
   end do
   close(42)

!   win_tot(:) = win_tot(:)
!write(0,*)'total='totwin,totwin/SP%numbins

   !HERE BIN COMPUTATION
   if (bin_type.eq.1) then       !equipopulated bins
      write(0,*) 'Computing equipopulated redshift bins'
      SP%zbins(1)=SP%zmin
      SP%zbins(SP%numbins+1)=SP%zmax
      do ib=1,SP%numbins
         zinf=SP%zbins(ib)
         if (ib.lt.SP%numbins) then
            diff=1000
            do nz2=2,num_z_win
               zsup=win_z(nz2)
!               diff=1000
               do nz=1,num_z_win
                  sig_z=(1+win_z(nz))*SP%sig_z0
                  win_temp(nz,ib)=win_tot(nz)!*0.5*(erf((win_z(nz)-zinf)/(sqrt(2.)*sig_z))-erf((win_z(nz)-zsup)/(sqrt(2.)*sig_z)))
               end do
               integral=0.d0
               do nz=1,num_z_win-1
                  if ((win_z(nz).ge.zinf).and.(win_z(nz).lt.zsup)) then
                     integral=integral+0.5*(win_temp(nz,ib)+win_temp(nz+1,ib))*(win_z(nz+1)-win_z(nz))
                  end if
               end do
               if ((abs(integral/totwin-1.d0/real(SP%numbins))).lt.diff) then
                  win(:,ib) = win_temp(:,ib)*0.5*(erf((win_z(:)-zinf)/(sqrt(2.)*sig_z))-erf((win_z(:)-zsup)/(sqrt(2.)*sig_z)))
                  diff = abs(integral/totwin-1.D0/real(SP%numbins))
                  SP%zbins(ib+1)=zsup
               end if
            end do
         else
            zsup = SP%zmax
            do nz=1,num_z_win
               sig_z=(1+win_z(nz))*SP%sig_z0
               win(nz,ib)=win_tot(nz)*0.5*(erf((win_z(nz)-zinf)/(sqrt(2.)*sig_z))-erf((win_z(nz)-zsup)/(sqrt(2.)*sig_z)))
            end do
            integral=0.d0
            do nz=1,num_z_win-1
               if ((win_z(nz).ge.zinf).and.(win_z(nz).lt.zsup)) then
                  integral=integral+0.5*(win_tot(nz)+win_tot(nz+1))*(win_z(nz+1)-win_z(nz))
               end if
            end do
         end if
         write(0,*) 'redshift limit',ib+1,'=',SP%zbins(ib+1)
!         write(0,*) 'Fraction of galaxies in bin',ib,'is',integral/totwin
      end do
     write(0,*)'COMPUTED EQUIPOPULATED REDSHIFT BINS'
     do ib=1,SP%numbins
        write(0,*)ib,' bin = [',SP%zbins(ib),':',SP%zbins(ib+1),']'
     end do

   else if (bin_type.eq.2) then  !equispaced bins
      write(0,*) 'Computing equispaced redshift bins'
      do ib=1,SP%numbins+1
            SP%zbins(ib)=(ib-1)*(SP%zmax-SP%zmin)/SP%numbins
      end do
      do ib=1,SP%numbins
         zinf=SP%zbins(ib)
         zsup=SP%zbins(ib+1)
         do nz=1,num_z_win
            sig_z=(1+win_z(nz))*SP%sig_z0
            win(nz,ib)=win_tot(nz)*0.5*(erf((win_z(nz)-zinf)/(sqrt(2.)*sig_z))-erf((win_z(nz)-zsup)/(sqrt(2.)*sig_z)))
         end do
      end do
   else if (bin_type.eq.3) then !user defined
      write(0,*) 'Computing user defined redshift bins'
      SP%zbins(1)=0.43
      SP%zbins(2)=0.54!0.418
      SP%zbins(3)=0.89!0.560
      SP%zbins(4)=0.99!0.678
!      SP%zbins(5)=0.15518248E+01!0.789
!      SP%zbins(6)=0.17313063E+01!0.900
!      SP%zbins(7)=0.18985152E+01!1.019
!      SP%zbins(8)=0.20570165E+01!1.155
!      SP%zbins(9)=0.22091091E+01!1.324
!      SP%zbins(10)=0.23563838E+01!1.576
!      SP%zbins(11)=0.25000000E+01!5.d0
      do ib=1,SP%numbins
         zinf=SP%zbins(ib)
         zsup=SP%zbins(ib+1)
         do nz=1,num_z_win
            sig_z=(1+win_z(nz))*SP%sig_z0
!            win(nz,ib)=1.50*exp(-(win_z(nz)-0.7)**2./0.32**2.)+0.20*exp(-(win_z(nz)-1.2)**2./0.46**2.)
            win(nz,ib)=win_tot(nz)*0.5*(erf((win_z(nz)-zinf)/(sqrt(2.)*sig_z))-erf((win_z(nz)-zsup)/(sqrt(2.)*sig_z)))
         end do
      end do
   else
      write(0,*) 'Why are you doing this?'
      stop
   end if

   !Normalization
   do ib=1,SP%numbins
      totwin2(ib)=0
      do nz=2,num_z_win
         if ((win_z(nz-1).ge.SP%zmin).and.(win_z(nz).le.SP%zmax)) then
            totwin2(ib)=totwin2(ib)+0.5*(win(nz-1,ib)+win(nz,ib))*(win_z(nz)-win_z(nz-1))
         end if
      end do
      !win(:,ib) = win(:,ib)/totwin2
   end do
   !writing windows to file
   open(42,file='mock_output/WL_mock/'//trim(mock_root)//'_window.dat')
   do nz=1,num_z_win
      write(42,10) win_z(nz),(win(nz,ib)/totwin2(ib),ib=1,SP%numbins)
   end do
   close(42)

   intra_frac=0
   do ib=1,SP%numbins
      do nz=1,num_z_win-1
         intra_frac(ib)=intra_frac(ib)+0.5*(win(nz,ib)+win(nz+1,ib))*(win_z(nz+1)-win_z(nz))
      end do
      intra_frac(ib)=intra_frac(ib)/totwin
      write(0,*)'fraction of galaxies in bin',ib,'=',intra_frac(ib)
   end do

   intra_zbins=SP%zbins


   deallocate(win_temp)



10 format(500(1e12.6,1x))

end subroutine GetSpecs

subroutine wl_mock_maker(mock_root, mock_type, sys_root, sys_dir, TH_auto, TH_cross, TH_matrix, THlmin, THlmax, THbin, THncr)

implicit none

!input stuff
integer, intent(in) :: THlmin, THlmax, THbin, THncr, mock_type
character(100),intent(in) :: mock_root, sys_root, sys_dir
real(mcp), dimension(THlmax,THbin) :: TH_auto
real(mcp), dimension(THlmax,THncr) :: TH_cross
real(mcp), dimension(THlmax,THbin,THbin) :: TH_matrix


!output stuff
real(mcp), dimension(THlmax,THbin) :: noise


!other stuff
logical, parameter :: chatty=.true.
real(mcp), dimension(THlmax,THbin) :: sig_auto
real(mcp), dimension(THlmax,THncr) :: sig_cross


real(mcp), allocatable, dimension(:) :: ng_r
integer :: i,l,bin,j,aa
character(len=2) :: indbin
real(mcp), parameter :: pi=3.14159
character(len=2) :: bchari, bcharj

real(mcp), dimension(intra_lmax,intra_numbins,intra_numbins) :: add_sys
real(mcp), dimension(intra_numbins) :: mul_sys
real(mcp), dimension(intra_numbins,intra_numbins) :: mul_fac
real(mcp) :: fake

!compute errors -----------------------------------------------
if (chatty) then
   write(0,*)  'computing errors...'
end if

allocate(ng_r(intra_numbins))

ng_r(:)=(intra_ng*3600.*(180./pi)**2.)*intra_frac(:)
do l=intra_lmin,intra_lmax
   noise(l,:)    = intra_d2/ng_r(:)
   sig_auto(l,:) = sqrt(2./(2.*l+1.)/(intra_fsky))*(TH_auto(l,:)+intra_d2/ng_r(:))
   sig_cross(l,:)= sqrt(2./(2.*l+1.)/(intra_fsky))*(TH_cross(l,:))
end do
if(chatty)  then
   write(0,*) '...done!'
end if

deallocate(ng_r)
!--------------------------------------------------------------


!writing mock and dataset files--------------------------------
if(chatty)  then
   write(0,*) 'Writing output files'
end if
open(42,file='mock_output/'//trim(mock_root)//'.dataset')
write(42,'(A,I5)')'num_z_bins =',intra_numbins
do i=1,intra_numbins+1
   write(indbin,'(I2)') i
   write(42,'(A,1F10.2)')'zbin('//trim(adjustl(indbin))//') =',intra_zbins(i)
end do
write(42,'(A,I5)')'num_z_window =',intra_num_z_temp
write(42,'(A,1F10.2)')'max_z =',intra_zmax
write(42,'(A,1F10.2)')'kmax =',intra_kmax
write(42,'(A,I5)')'ell_max =',intra_lmax
write(42,'(A,I5)')'ell_min =',intra_lmin
write(42,'(A,I5)')'num_z_p =',num_z_win
write(42,'(A,1F10.2)')'fsky =',intra_fsky
write(42,'(A,A)')'WL_covmat =','%DATASETDIR%WL_mock/'//trim(mock_root)//'_covmat.dat'
write(42,'(A,A)')'WL_noise =','%DATASETDIR%WL_mock/'//trim(mock_root)//'_noise.dat'
write(42,'(A,A)')'window_file =','%DATASETDIR%WL_mock/'//trim(mock_root)//'_window.dat'
close(42)

open(42,file='mock_output/WL_mock/'//trim(mock_root)//'_noise.dat')
do l=intra_lmin,intra_lmax
   write(42,10) dble(l),(noise(l,bin),bin=1,intra_numbins)
end do
close(42)

open(42,file='mock_output/WL_mock/'//trim(mock_root)//'_covmat.dat')
do l=intra_lmin,intra_lmax
   do i=1,intra_numbins   
      write(42,10) (TH_matrix(l,i,bin),bin=1,intra_numbins)
   end do
end do
close(42)


open(42,file='mock_file_forplot/'//trim(mock_root)//'_fracbins.dat')
do bin=1,intra_numbins
   write(42,*)bin,intra_frac(bin)
end do
close(42)

open(42,file='mock_file_forplot/'//trim(mock_root)//'_autofile.dat')
do l=intra_lmin,intra_lmax
   if (l.ge.intra_lmin) then
       write(42,10) dble(l),(TH_auto(l,bin)*(l*(l+1))/(2*pi),bin=1,intra_numbins),(sig_auto(l,bin)*(l*(l+1))/(2*pi),bin=1,intra_numbins)
   end if
end do
close(42)

open(42,file='mock_file_forplot/'//trim(mock_root)//'_crossfile.dat')
do l=intra_lmin,intra_lmax
   if (l.ge.intra_lmin) then
       write(42,10) dble(l),(TH_cross(l,bin)*(l*(l+1))/(2*pi),bin=1,intra_numcross),(sig_cross(l,bin)*(l*(l+1))/(2*pi),bin=1,intra_numcross)
   end if
end do
close(42)

!---------------------------------------------------------------

if (chatty) then
   write(0,*) 'FINISHED!'
   write(0,*) 'produced files in mock_output folder'
   write(0,*) trim(mock_root)//'.dataset'
   write(0,*) trim(mock_root)//'_covmat.dat'
   write(0,*) trim(mock_root)//'_noise.dat'
   write(0,*) trim(mock_root)//'_window.dat'
   write(0,*) 'also some other stuff in mock_file_forplot folder'
   write(0,*) trim(mock_root)//'_autofile.dat'
   write(0,*) trim(mock_root)//'_crossfile.dat'
   write(0,*) trim(mock_root)//'_window_total.dat'
end if






write(0,*)'BELLA PE TUTTI!'



10 format(500(1e12.6,1x))
deallocate(win_z,win,win_tot)
deallocate(intra_zbins,intra_frac)
end subroutine wl_mock_maker


end module WL_MockMaker
