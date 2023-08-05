program test
use mcf_tipos
real(kind=dp), dimension(1000) :: gr
integer :: ig
real(kind=dp) :: a,b

do i=1, 1000
 gr(i) = 0.0_dp
end do

a=324.5
b=5.4
ig = idint(a/b)
print*, ig
print*, gr(ig) + 2
end program test
