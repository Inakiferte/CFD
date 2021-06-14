module mcf_Vt2
use mcf_tipos
!Vt-ren matrizea aldiunero
public :: Vt
!-----------------------------------------------------------------
contains
!Vt matrizea------------------------------------------------------
function Vt(P,U,V,tau,Re,deltaX)
use mcf_tipos
use mcf_F_G

real(kind=dp), dimension(:,:), intent(in) :: P                     !Aurreko aldiuneko presioaren matrizea
real(kind=dp), dimension(:,:), intent(in) :: U                     !Aurreko aldiuneko abiaduraren matrizea
real(kind=dp), dimension(:,:), intent(in) :: V                     !Aurreko aldiuneko abiaduraren matrizea

real(kind=dp), intent(in) :: tau, Re, deltaX                                   !Gauden pasuoa

real(kind=dp), dimension(size(P,1), size(P,2)) :: Ptp              !Presio matrize partziala, iteratzen joateko
real(kind=dp), dimension(size(V,1), size(V,2)) :: Vt               !Gure aldiuneko abiduraren matrizea

real(kind=dp) :: deltaY

integer :: j,k,nVt,mVt,i                                             !Iterazioak egiteko/dimentsioen balioak finkatzeko

nVt=size(Vt,1)                                                     ! n x m kasuko n-ren balioa
mVt=size(Vt,2)                                                     ! n x m kasuko m-ren balioa

deltaY=deltaX
Ptp=P                                                              !Partziala definitzen dugu, hau iteratuko dugu eta

!Domeinu barruko balioak-----------------------------------------
do j=2, nVt-1
 do k=2, mVt-1

  Vt(j,k)= G(U,V,j,k,tau,Re,deltaX) - (tau/deltaY) * (Ptp(j,k) - Ptp(j+1,k))

 end do
end do
!----------------------------------------------------------------

!Domeinu kanpoko balioak----------------------------------------
!OUTFLOW BOUNDARY
do i=1, nVt
 Vt(i,1)=-Vt(i,2)
 Vt(i,mVt)=-Vt(i,mVt-1)
end do
!----------------
!OUTFLOW BOUNDARY
do i=2, mVt-1
 Vt(nVt,i)=0.0_dp
end do
!----------------
!INFLOW BOUNDARY
do i=2, mVt-1
 Vt(1,i)=Vt(2,i)
end do
!----------------
!----------------------------------------------------------------
end function Vt
end module mcf_Vt2

