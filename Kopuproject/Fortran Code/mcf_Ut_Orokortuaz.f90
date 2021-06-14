module mcf_Ut2
use mcf_tipos
!Ut-ren matrizea aldiunero
public :: Ut
!-----------------------------------------------------------------
contains

!Ut matrizea------------------------------------------------------
function Ut(P,U,V,tau,Re,deltaX)
use mcf_tipos
use mcf_F_G

real(kind=dp), dimension(:,:), intent(in) :: P                     !Aurreko aldiuneko presioaren matrizea
real(kind=dp), dimension(:,:), intent(in) :: U                     !Aurreko aldiuneko abiaduraren matrizea
real(kind=dp), dimension(:,:), intent(in) :: V                     !Aurreko aldiuneko abiaduraren matrizea

real(kind=dp), intent(in) :: tau,Re,deltaX                         !Gauden pausoa

real(kind=dp), dimension(size(P,1), size(P,2)) :: Ptp              !Presio matrize partziala, iteratzen joateko
real(kind=dp), dimension(size(U,1), size(U,2)) :: Ut               !Gure aldiuneko abiduraren matrizea

integer :: j,k,nUt,mUt,i                                           !Iterazioak egiteko/dimentsioen balioak finkatzeko

nUt=size(Ut,1)                                                     ! n x m kasuko n-ren balioa
mUt=size(Ut,2)                                                     ! n x m kasuko m-ren balioa

Ptp=P                                                              !Partziala definitzen dugu, hau iteratuko dugu eta

!Domeinu barruko balioak-----------------------------------------
do j=2, nUt-1
 do k=2, mUt-1

  Ut(j,k)= F(U,V,j,k,tau,Re,deltaX) - (tau/deltaX) * (Ptp(j,k+1) - Ptp(j,k))

 end do
end do
!----------------------------------------------------------------

!Domeinu kanpoko balioak----------------------------------------

!OUTFLOW BOUNDARY
do i=2, nUt
 Ut(i,1)=0.0_dp         !Pareta
 Ut(i,mUt)=0.0_dp       !Pareta
end do

!----------------
!OUTFLOW BOUNDARY
do i=2, mUt-1
 Ut(nUt,i)= -Ut(nUt-1,i)
end do
!----------------
!U=KTE-----------
do i=1, mUt
 Ut(1,i)=0.9999999_dp   !Errekaren abiadura!
end do
!----------------

end function Ut

end module mcf_Ut2

