!Bertsio honetan ohartu gara F eta G ez ditugula zertan eta bi ataletan banatu behar, Fa, Fb, Ga eta Gb hain zuzen ere.
!g95 WARNING g95
module mcf_F_G
public :: F,G
contains

!F funtzioa--------------------------------------------------------------------
function F(U,V,j,k,tau,Re,deltaX)
use mcf_tipos
integer, intent(in) :: j                               !Errenkaden kokapena
integer, intent(in) :: k                               !Zutabeen kokapena
real(kind=dp), dimension(:,:), intent(in) :: U         !Aurreko iterazioko U
real(kind=dp), dimension(:,:), intent(in) :: V         !Aurreko iterazioko V
real(kind=dp) :: F                                     !F-k balio bat j,k bakoitzerako
real(kind=dp), intent(in) :: tau, Re, deltaX
real(kind=dp) :: deltaY
deltaY=deltaX                                          !Simetria inposatu

F=U(j,k) + (tau/(Re*(deltaX**2)))*(U(j,k-1) - 2.0_dp*U(j,k) + U(j,k+1))          &
         & + (tau/(Re*(deltaY**2)))*(U(j-1,k)-2.0_dp*U(j,k)+U(j+1,k))                     &
         & - (tau/deltaX)*((0.5_dp*(U(j,k+1)+U(j,k)))**2 -( 0.5_dp*(U(j,k-1)+U(j,k)))**2) &
         & - (tau/deltaY)*((0.5_dp*(U(j,k)+U(j-1,k))*0.5_dp*(V(j-1,k)+V(j-1,k+1))))           &
         & + (tau/deltaY)*((0.5_dp*(U(j,k)+U(j+1,k))*0.5_dp*(V(j,k)+V(j,k+1))))

end function F
!G funtzioa------------------------------------------------
function G(U,V,j,k,tau,Re,deltaX)
use mcf_tipos
real(kind=dp), dimension(:,:), intent(in) :: U         !Aurreko iterazioko U
real(kind=dp), dimension(:,:), intent(in) :: V         !Aurreko iterazioko V
integer,intent(in)::j,k
real(kind=dp) :: G                                     !G j eta k bakoitzerako dagokion balioa izango da!
real(kind=dp), intent(in) :: tau, Re, deltaX
real(kind=dp) :: deltaY
deltaY=deltaX                                          !Simetria inposatu

G=V(j,k) + (tau/(Re*deltaX**2))*(V(j,k-1) - 2.0_dp*V(j,k) + V(j,k+1))                                                &
         & + (tau/(Re*deltaY**2))*(V(j-1,k)-2.0_dp*V(j,k)+V(j+1,k))                                                          &
         & - (tau/deltaX) * 0.25_dp*((U(j,k) + U(j+1,k))*(V(j,k) + V(j,k+1))-((U(j,k-1)+U(j+1,k-1))*(V(j,k-1)+V(j,k))))      &
         & - (tau/deltaY) * ((1.0_dp/2.0_dp*(V(j,k) + V(j-1,k)))**2.0_dp - (1.0_dp/2.0_dp*(V(j,k) + V(j+1,k)))**2.0_dp)

end function G


end module mcf_F_G

