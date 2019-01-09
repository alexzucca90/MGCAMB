module WL_MockMaker

!----------------------------------------------!
!MOCK ESTIMATOR FOR SHEAR SPECTRA (MESS) v2.2  !
!                  07/09/2015                  !
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
real(mcp), allocatable, dimension(:) :: win_z, win_tot
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
real(mcp), dimension(num_z_win) :: diff


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
   SP%zmax=5.d0
   SP%zmin=0.d0
   SP%kmax=20.d0
   SP%sig_z0=0.05


else if (mock_type.eq.2) then

   write(0,*) 'Using Code Comparison specifications'

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
   SP%zmax=2.5d0
   SP%zmin=0.d0
   SP%kmax=20.d0
   SP%sig_z0=0.05


else if (mock_type.eq.3) then
   write(0,*) 'Using LSST specifications'
   write(0,*) 'This will be done eventually'
   write(0,*) 'now it is just the same Euclid thing but with Santiago specs'



   !SPECIFICATIONS
   SP%num_z_temp=num_z_win
   SP%numbins = 6
   SP%numcross = (SP%numbins*(SP%numbins+1)/2)-SP%numbins
   SP%lmin = 1
   SP%lmax = 5000
   SP%fsky= 15000/(4*pi*(180/pi)**2.) !0.3636
   SP%d2=2*0.0484
   SP%ng=30.
   SP%z0=0.9
   SP%zmax=3.d0
   SP%zmin=0.5d0
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
   allocate(win_z(num_z_win),win(num_z_win,SP%numbins),win_tot(num_z_win),win_temp(num_z_win,SP%numbins))
   do nz=1,num_z_win
      win_z(nz)=SP%zmin+(SP%zmax-SP%zmin)/(num_z_win)*(nz-1)
      win_tot(nz)=(win_z(nz)**2*exp(-(win_z(nz)/(SP%z0/1.412))**1.5))
   end do
   open(42,file='mock_output/'//trim(mock_root)//'_window_total.dat')
   do nz=1,num_z_win
      write(42,10) win_z(nz),win_tot(nz)
   end do
   close(42)

   totwin=0
   do nz=2,num_z_win
      if ((win_z(nz-1).ge.SP%zmin).and.(win_z(nz).le.SP%zmax)) then
         totwin=totwin+0.5*(win_tot(nz-1)+win_tot(nz))*(win_z(nz)-win_z(nz-1))
      end if
   end do



   !HERE BIN COMPUTATION
   if (bin_type.eq.1) then       !equipopulated bins
      write(0,*) 'Computing equipopulated redshift bins'
      SP%zbins(1)=SP%zmin
      do ib=1,SP%numbins
         zinf=SP%zbins(ib)
         if (ib.lt.SP%numbins) then
            diff(:)=1000
            do nz2=2,num_z_win
               zsup=win_z(nz2)
               do nz=1,num_z_win
                  sig_z=(1+win_z(nz))*SP%sig_z0
                   if ((win_z(nz).ge.zinf).and.(win_z(nz).le.zsup)) then
                      win_temp(nz,ib)=win_tot(nz)
                   else
                      win_temp(nz,ib)=0.d0
                   end if
               end do
               integral=0
               do nz=2,num_z_win
                  integral=integral+0.5*(win_temp(nz-1,ib)+win_temp(nz,ib))*(win_z(nz)-win_z(nz-1))
               end do
               diff(nz2)=abs(integral-totwin/SP%numbins)
               if (diff(nz2).lt.diff(nz2-1)) then
                  SP%zbins(ib+1)=zsup
                  win(:,ib)=win_temp(:,ib)
               end if
            end do
            write(0,*) 'redshift limit',ib+1,'=',SP%zbins(ib+1)
         else
            zsup=SP%zmax
            SP%zbins(ib+1)=zsup
            do nz=1,num_z_win
               sig_z=(1+win_z(nz))*SP%sig_z0
               win(nz,ib)=win_tot(nz)*0.5*(erf((win_z(nz)-zinf)/(sqrt(2.)*sig_z))-erf((win_z(nz)-zsup)/(sqrt(2.)*sig_z)))
            end do
            write(0,*) 'redshift limit',ib+1,'=',SP%zbins(ib+1)
         end if
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
   else
      write(0,*) 'Why are you doing this?'
      stop
   end if


   !writing windows to file
   open(42,file='mock_output/'//trim(mock_root)//'_window.dat')
   do nz=1,num_z_win
      write(42,10) win_z(nz),(win(nz,ib),ib=1,SP%numbins)
   end do
   close(42)

   intra_frac=0
   do ib=1,SP%numbins
      do nz=2,num_z_win
         intra_frac(ib)=intra_frac(ib)+0.5*(win(nz-1,ib)+win(nz,ib))*(win_z(nz)-win_z(nz-1))
      end do
      intra_frac(ib)=intra_frac(ib)/totwin
      write(0,*)'fraction of galaxies in bin',ib,'=',intra_frac(ib)
   end do

   intra_zbins=SP%zbins

   deallocate(win_temp)



10 format(500(1e12.6,1x))

end subroutine GetSpecs

subroutine wl_mock_maker(mock_root, mock_type, sys_root, sys_dir, TH_auto, TH_cross, THlmax, THbin, THncr, systematics)

implicit none
character(100) :: mock_root, sys_root, sys_dir
integer :: THlmax, THbin, THncr, mock_type
logical :: systematics
logical, parameter :: chatty=.true.
real(mcp), dimension(THlmax,THbin) :: TH_auto,sig_auto
real(mcp), dimension(THlmax,THncr) :: TH_cross,sig_cross
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
do l=1,intra_lmax!SISTEMA GLI ELLE
   sig_auto(l,:)=sqrt(2./(2.*l+1.)/(intra_fsky))*(TH_auto(l,:)+intra_d2/ng_r(:))
   sig_cross(l,:)=sqrt(2./(2.*l+1.)/(intra_fsky))*(TH_cross(l,:))
end do
if(chatty)  then
   write(0,*) '...done!'
end if

deallocate(ng_r)
!--------------------------------------------------------------


!Adding systematics -------------------------------------------
if(systematics) then

if (chatty) then
   write(0,*) 'Reading ',trim(sys_root),' systematics files'
end if

  do i=1,intra_numbins
     if (i.lt.10) then
        write(bchari,'(I1)') i
     else
        write(bchari,'(I2)') i
     endif

     do j=1,intra_numbins
        if (j.lt.10) then
           write(bcharj,'(I1)') j
        else
           write(bcharj,'(I2)') j
        endif

        if (i.le.j) then
           !reading additive systematics
           open(42,file=trim(sys_dir)//'/powspec-'//trim(sys_root)//'-B'//trim(adjustl(bchari))//'B'//trim(adjustl(bcharj))//'.dat')
           do l=1,intra_lmax
              read(42,*) fake, add_sys(l,i,j)
           end do
           close(42)

           !reading multiplicative systematics
           open(43,file=trim(sys_dir)//'/modelpar-'//trim(adjustl(sys_root))//'-B'//trim(adjustl(bchari))//'B'//trim(adjustl(bcharj))//'.dat')
           if (i.eq.j) then
              read(43,*) fake
              read(43,*) fake, fake, mul_sys(i)
              mul_fac(i,j) = (1+mul_sys(i))**2.
           else
              read(43,*) fake
              read(43,*) fake, fake, mul_sys(i), fake, fake, fake, fake, fake, fake, mul_sys(j)
              mul_fac(i,j) = (1+mul_sys(i)+mul_sys(j)+mul_sys(i)*mul_sys(j))
           end if
           close(43)
        end if
     end do
  end do

  !adding systematics to the mock spectra
  aa=0
  do i=1,intra_numbins
     do j=1,intra_numbins
        if (i.eq.j) then
           TH_auto(:,i) = mul_fac(i,j)*(TH_auto(:,i)+add_sys(:,i,j))
        else if (i.lt.j) then
           aa=aa+1
           TH_cross(:,aa) = mul_fac(i,j)*(TH_cross(:,aa)+add_sys(:,i,j))
        end if
     end do
  end do

if(chatty)  then
   write(0,*) 'including systematics in the power spectra'
end if


end if
if(chatty) then
   write(0,*)  '...done!'
end if
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
if (systematics) then
write(42,'(A)')'systematics = T'
else
write(42,'(A)')'systematics = F'
end if
write(42,'(A,A)')'WL_auto_file =','%DATASETDIR%WL_mock/'//trim(mock_root)//'_autofile.dat'
write(42,'(A,A)')'WL_cross_file =','%DATASETDIR%WL_mock/'//trim(mock_root)//'_crossfile.dat'
write(42,'(A,A)')'window_file =','%DATASETDIR%WL_mock/'//trim(mock_root)//'_window.dat'
close(42)

open(42,file='mock_output/'//trim(mock_root)//'_autofile.dat')
do l=1,intra_lmax
   if (l.ge.intra_lmin) then
       write(42,10) dble(l),(TH_auto(l,bin)*(l*(l+1))/(2*pi),bin=1,intra_numbins),(sig_auto(l,bin)*(l*(l+1))/(2*pi),bin=1,intra_numbins)
   end if
end do
close(42)

open(42,file='mock_output/'//trim(mock_root)//'_crossfile.dat')
do l=1,intra_lmax
   if (l.ge.intra_lmin) then
       write(42,10) dble(l),(TH_cross(l,bin)*(l*(l+1))/(2*pi),bin=1,intra_numcross),(sig_cross(l,bin)*(l*(l+1))/(2*pi),bin=1,intra_numcross)
   end if
end do
close(42)

!---------------------------------------------------------------

if (chatty) then
   write(0,*) 'FINISHED!'
   write(0,*) 'produced files in output folder'
   write(0,*) trim(mock_root)//'.dataset'
   write(0,*) trim(mock_root)//'_autofile.dat'
   write(0,*) trim(mock_root)//'_crossfile.dat'
   write(0,*) trim(mock_root)//'_window.dat'
end if






write(0,*)'BELLA PE TUTTI!'



10 format(500(1e12.6,1x))
deallocate(win_z,win,win_tot)
deallocate(intra_zbins,intra_frac)
end subroutine wl_mock_maker


end module WL_MockMaker
