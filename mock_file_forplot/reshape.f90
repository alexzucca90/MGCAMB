program shuffle

implicit none

integer :: i
integer, parameter :: numlines = 2991
real*8, dimension(numlines) :: ell, auto11, auto22, cross12
real*8 :: fake


open(42,file='comparison_autofile.dat')
do i=1,numlines
   read(42,10) ell(i),auto11(i), fake, auto22(i)
end do
close(42)

open(42,file='comparison_crossfile.dat')
do i=1,numlines
   read(42,10) ell(i),fake, cross12(i)
end do
close(42)

open(42,file='comparison_cls_matteo.dat')
do i=1,numlines
   write(42,10) ell(i),auto11(i),auto22(i),cross12(i)
end do
close(42)


10 format(500(1e12.6,1x))
end 
