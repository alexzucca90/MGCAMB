program wind
implicit none

integer :: i
integer,parameter :: nlines=100
real*8, dimension(nlines) :: w1, w2


open(42,file='W^k_z-CFHTLenS-1.dat')
do i=1,nlines
   read(42,*) w1(i)
end do
close(42)

open(42,file='W^k_z-CFHTLenS-2.dat')
do i=1,nlines
   read(42,*) w2(i)
end do
close(42)

open(42,file='stewin.dat')
do i=1,nlines
   write(42,*) (i-1)*0.01, w1(i), w2(i)
end do
close(42)


end program wind
