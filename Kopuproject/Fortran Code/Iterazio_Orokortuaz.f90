program batpausua
!Programa honetako ideia pausu guztiak jarraitzea da eta aldiune baterako U,V eta P ploteatzea.
!JARRAITU BEHARREKO PAUSUAK:
!1)------------------------------------------------
!
!  READ MESH
!
!---------------------------------------------------
!
!2)------------------------------------------------
!
! SET IC
! IMPLEMENTATION OF BC
!
!--------------------------------------------------
!
!3)-------------------------------------------------
!
! COMPUTE F AND G
! pausu hau ez dugu egiten, ez dugu behar
!
!---------------------------------------------------
!
!4)-------------------------------------------------
!
! SOLVE POISSON ---> RELAXATION METHOD
!
!---------------------------------------------------
!
!5)-------------------------------------------------
!
! SOLVE U
! SOLVE V
!
!---------------------------------------------------
!
!6)-------------------------------------------------
!
! UPDATE BOUNDARY CONDITIONS
!
!---------------------------------------------------
use mcf_tipos                                                         !dp zehaztasunarekin egingo dugu lan
use initial2                                                          !Hasiera P, V eta U finkatzen ditu
use mcf_F_G                                                           !F(U,V,j,k) eta P(U,V,j,k) ditugu
use mcf_Pt2                                                           !Poisson solverra erlaxazio metodoarekin
use mcf_Ut2                                                           !U matrizea aldiunero
use mcf_Vt2                                                           !V matrizea aldiunero
integer, parameter ::           i=10000                              !Iterazio kopurua
integer, parameter ::           M=60                                 !Sentzu fisikodun zatiko mesh size
real(kind=dp),parameter ::      L=1.0_dp                              !Karratuaren aldeen luzera!
real(kind=dp),parameter ::      T=1.0_dp
real(kind=dp) :: tau                                                  !Pausuak
real(kind=dp) :: deltaX
real(kind=dp) :: deltaY
real(kind=dp), parameter:: Re=1.0_dp
real(kind=dp), allocatable, dimension(:,:) :: P_H, U_H, V_H
real(kind=dp), allocatable, dimension(:,:) :: P                       !Presioaren matrizea
real(kind=dp), allocatable, dimension(:,:) :: U                       !Abiaduraren matrizea
real(kind=dp), allocatable, dimension(:,:) :: V                       !Abiaduraren matrizea
integer :: nP,mP,nU,mU,nV,mV                                          !Erabiliko ditugun n x m dimentsioak, programa osoan berdinak
integer :: w,j,k,jj
real(kind=dp) :: prueb                                                !Probak egiteko

!PLOTEATZEKO--------------------------------------------------!
real(kind=dp), allocatable, dimension(:,:):: Uplot, Vplot
real(kind=dp):: x, y                                                  !Ploteatzeko puntuak
!-------------------------------------------------------------!

!--------------FINKATU REAL-AK------------------------------
deltaX=L/real(M)
deltaY=deltaX
tau=0.1_dp*Re*(deltaX)**2

!-----------------------------------------------------------

!--------------EGONKORTASUN BALDINTZA ZIURTATU-------------
if (tau/(Re*deltaX**2)<=0.24_dp) then
        print*, "Egonkortasun baldintzak betetzen ditu   "
else
 print*, "Ez ditu egonkortasun baldintzak betetzen"
        STOP
end if
!-----------------------------------------------------------

!1)MESH SIZE----------------------------------------------------------
print*, "Mesh size: ", M
!---------------------------------------------------------------------

!2)SET IC AND IMPLEMENTATION OF BC ON IC------------------------------
nP=size(Pb(M),1)
mP=size(Pb(M),2)
nU=size(Ub(M),1)
mU=size(Ub(M),2)
nV=size(Vb(M),1)
mV=size(Vb(M),2)
allocate(P_H(nP,mP))
allocate(U_H(nU,mU))
allocate(V_H(nV,mV))

P_H=Pb(M)
U_H=Ub(M)
V_H=Vb(M)
!----------------------------------------------------------------------

!-----------------------------------Ploteatzeko-----------------------

allocate(Uplot(M,M),Vplot(M,M))

!--------------------------------------------------------------------


!ITERAZIOA HASI-------------------------------------------------------
allocate(P(nP,mP))                                                     !Dimentsioak P-rentzat
allocate(U(nU,mU))                                                     !Dimentsioak U-rentzat
allocate(V(nV,mV))                                                     !Dimentsioak V-rentzat

do w=1, i
!4) SOLVE POISSON---->RELAXATION METHOD--------------------------------
P= Pt(P_H,U_H,V_H,tau,Re,deltaX)
!----------------------------------------------------------------------
!5) SOLVE U AND V  +  6) UPDATE BOUNDARY CONDITIONS  ------------------
U=Ut(P,U_H,V_H,tau,Re,deltaX)

V=Vt(P,U_H,V_H,tau,Re,deltaX)

!----------------------------------------------------------------------

!EGUNERATU MATRIZEAK---------------------------------------------------
P_H=P
U_H=U
V_H=V
!----------------------------------------------------------------------

print*, w, "Iterazioan goaz", i, "tik"


!FORT-ak EGITEKO DO-A JARRI HASIERAN ETA AMAIERA ETA AZPIKO ALGORTIMOA ERABILI======
!Uplot=0.0_dp
!Vplot=0.0_dp



!do j=1,M
 !do k=1,M
 !Uplot(j,k)=(U(j+1,k)+U(j+1,k+1))/2.0_dp
 !end do
!end do


!do j=1,M
 !do k=1,M
 !Vplot(j,k)=(V(j,k+1)+V(j+1,k+1))/2.0_dp
 !end do
!end do

!y=L
!do j=1,M
 !x=0
 !do k=1, M
 !write(unit=100+w,fmt="(100f12.6)") x, y, Uplot(j,k), Vplot(j,k)
 !x = x + deltaX
 !end do
!y = y - deltaY
!end do

end do

!====================================================================================

!ITERAZIOAK AMATITU----------------------------------------------------

!---------------------------------------------------------------------------!

!------------------Plot------------------------------------------------

Uplot=0.0_dp
Vplot=0.0_dp

do j=1,M
 do k=1,M
 Uplot(j,k)=(U(j+1,k)+U(j+1,k+1))/2.0_dp
 end do
end do


do j=1,M
 do k=1,M
 Vplot(j,k)=(V(j,k+1)+V(j+1,k+1))/2.0_dp
 end do
end do


open(unit=102, file="AbiK.dat", action="write", status="replace")

y=L
do j=1,M
 x=0
 do k=1, M
 write(unit=102,fmt="(100f12.6)") x, y, Uplot(j,k), Vplot(j,k)
 x = x + deltaX
 end do
y = y - deltaY
end do
close(unit=102)

open(unit=103, file="Pk.dat", action="write", status="replace")

y=L
do j=1,M
 x=0
 do k=1, M
 write(unit=103,fmt="(100f12.6)") x, y, P(j,k)
 x = x + deltaX
 end do
y = y - deltaY
end do
close(unit=103)

!--------------------------------------------------------------------



end program batpausua

