module initial2
public :: Pb,Ub,Vb
private::erf_egonkortua
contains
!Presioa boundary------------------------------
function Pb(n)
use mcf_tipos
integer, intent(in) :: n                          !Mesh size
real(kind=dp), allocatable, dimension(:,:) :: Pb
allocate(Pb(n+2,n+2))                             !Domeinu kanpoko beste 2 zutabe eta errenkada
Pb=0.0_dp
!----------------------------------------------
end function Pb

!Abiaduran IC--------------------------------------------------------------

function Vb(n)
use mcf_tipos
integer, intent(in):: n
real(kind=dp), allocatable, dimension(:,:):: Vb

allocate(Vb(n+1,n+2))
Vb=0.0_dp

end function Vb
!----------------------------------------------------------------------------------------

!Abiadura U IC----------------------------------------------------------------------------------

function Ub(n)
use mcf_tipos
integer, intent(in):: n
real(kind=dp), allocatable, dimension(:,:):: Ub
integer:: nub, mub
integer:: i, j

allocate(Ub(n+2,n+1))

nub=size(Ub,1)                                                            !Errenkada kopurua.
mub=size(Ub,2)                                                            !Zutabe kopurua.

!Domeinu kanpoko balioak------------------------------------------------

do j=1, mub
 Ub(1,j)=erf_egonkortua(0.0_dp)
end do

do i=2,nub
 do j=1,mub
  Ub(i,j)=0.0_dp
 end do
end do
!------------------------------------------------------------------------

end function Ub
!-----------------------------------------------------------------------------------------------

!ERF FUNTZIOA EGONKORTUKO DUGU!
function erf_egonkortua(eps)  result(emaitza)
use mcf_tipos
real(kind=dp),intent(in)::eps
real(kind=dp),parameter:: xa=0.0_dp,xb=6.0_dp
real(kind=dp)::x,x1,emaitza,tarte,y1,y2
integer::i,n
n=1000
tarte=abs(xa-xb)/(n-1)
x1=xa
y2=0.0_dp
do
x=x1 + tarte
y1=erf(x)
 if (abs(y1-y2) == eps) then
    exit
 end if
 y2=y1
 x1=x
end do
emaitza=y2

end function erf_egonkortua


end module initial2
