program wind
implicit none

integer :: i
integer,parameter :: nlines=30
real*8, dimension(nlines) :: ell, w1, w2, w12


open(42,file='C^k_l-CFHTLenS-11.dat')
do i=1,nlines
   read(42,*) ell(i), w1(i)
end do
close(42)
open(42,file='C^k_l-CFHTLenS-22.dat')
do i=1,nlines
   read(42,*) ell(i), w2(i)
end do
close(42)

open(42,file='C^k_l-CFHTLenS-12.dat')
do i=1,nlines
   read(42,*) ell(i), w12(i)
end do
close(42)

open(42,file='stecls.dat')
do i=1,nlines
   write(42,10) ell(i), w1(i)*ell(i)*(ell(i)+1)/(2*3.14159), w2(i)*ell(i)*(ell(i)+1)/(2*3.14159),w12(i)*ell(i)*(ell(i)+1)/(2*3.14159)
end do
close(42)

10 format(500(1e12.6,1x))
end program wind
