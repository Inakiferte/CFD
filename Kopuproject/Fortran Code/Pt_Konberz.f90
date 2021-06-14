module mcf_Pt2


public :: Pt

contains!=======================
        
function Pt(P,U,V,tau,Re,deltaX)
use mcf_tipos
use mcf_F_G

!MATRIZEAK--------------------------------------------------------------------------
real(kind=dp), dimension(:,:), intent(in) :: P,U,V                                   ! Aurreko aldiuneko matrizeak (in)
real(kind=dp), dimension(size(P,1), size(P,2)) :: Pp                                 ! Formula aplikaturiko matrizea
real(kind=dp), dimension(size(P,1), size(P,2)) :: P_last                             ! Aurreko iterazioko matrizea
real(kind=dp), dimension(size(P,1), size(P,2)) :: Pt                                ! Konbergitzen duen matrizea (out)
!-----------------------------------------------------------------------------------

!BEKTOREAK--------------------------------------------------------------------------
real(kind=dp), allocatable, dimension(:) :: yt                                       ! Pp-rekin elkarturiko elementuak
real(kind=dp), allocatable, dimension(:) :: y_last                                   ! P_last-ekin elkarturiko elementuak
!-----------------------------------------------------------------------------------

!GAINONTZEKO (IN)- AK---------------------------------------------------------------
real(kind=dp), intent(in) :: tau, deltaX,Re
!-----------------------------------------------------------------------------------

!REAL-AK----------------------------------------------------------------------------
real(kind=dp) :: deltaY
real(kind=dp), parameter :: eps=1.0E-7_dp                                           ! Konbergentzia parametroa
!-----------------------------------------------------------------------------------

!INTEGER-AK-------------------------------------------------------------------------
integer :: j,k                                                                      ! Matrizeen elemtuak hartzeko
integer :: npt,mpt,ny                                                               ! Matrize eta bektoreen dimentsioak finkatzeko
integer :: i_t, i_last                                                              ! Bektoreentzako iterazio parametroa
!-----------------------------------------------------------------------------------

!POISSON-SOLVER---------------------------------------------------------------------
npt=size(Pp,1)                                                                      ! Dimentsioak n x m moduan
mpt=size(Pp,2)                                                                      ! Dimentsioak n x m moduan
ny=npt*mpt
deltaY=deltaX                                                                       ! Simetria inposatu
allocate(yt(ny), y_last(ny))                                                        ! Dimentsioa finkatu elementu guztiak hartzeko
P_last=P                                                                            ! Aurreko aldiuneko matrizea finkatu
!ITERAZIOAK HASI KONBERGENTZIA LORTU ARTE-------------------------------------------
do 
!Presioaren matrizearen kalkuloa---------------------------------------------------
do j=2, npt-1
 do k=2, mpt-1

  Pp(j,k)= 0.25_dp*(P_last(j-1,k) + P_last(j+1,k) + P_last(j,k-1) + P_last(j,k+1))         &
         & - (1.0_dp/deltaX)*(F(U,V,j,k,tau,Re,deltaX) - F(U,V,j,k-1,tau,Re,deltaX))                &
         & - (1.0_dp/deltaY)*(G(U,V,j-1,k,tau,Re,deltaX) - G(U,V,j,k,tau,Re,deltaX))

 end do
end do
!Domeinu kanpoko balioak----------------------------------------------------------
do j=2,npt-1
 Pp(1,j)=0.0_dp
 Pp(npt,j)=0.0_dp
end do

do i=1, npt
 Pp(i,1)=0.0_dp
 Pp(i,npt)=0.0_dp
end do
!----------------------------------------------------------------------------------


!Bi matrizeak bektore batean gorde---------------------------------------------------
i_t=0
do j=1, npt
 do k=1, npt
  i_t=i_t+1
  yt(i_t)=Pp(j,k)
 end do
end do

i_last=0
do j=1, npt
 do k=1, npt
  i_last=i_last+1
  y_last(i_last)=P_last(j,k)
 end do
end do
!Konbergentzia---------------------------------------------------------------------
 
if (sum(abs(yt-y_last))/ny<eps) then
        Pt=Pp                                                                      ! Konbergitzen duen matrizea (out)
        print*, "Konbergitu du"
        exit
else
        P_last=Pp

end if
!----------------------------------------------------------------------------------
end do
end function Pt
end module mcf_Pt2

