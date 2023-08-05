!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module containing physical constants
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module constants
implicit none
!Unit conversion
REAL*8, parameter :: J2eV=6.24150907e18 ! Joules to eV
REAL*8, parameter :: ev2cm1=8065.54429d0

!Mathematical constants
REAL*8,PARAMETER :: pi = 3.14159265358979324d0

!Physical constants
REAL*8, parameter :: kb=8.6173303d0*0.00001d0 ! Boltzman constant in eV/K
end module constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module containing physical parameters of the system
! also parameters from the calculation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module physical
use constants, only : J2eV
implicit none
INTEGER*8, parameter :: nbndmin=1          !minimum band index
INTEGER*8, parameter :: nbndmax=81        !maximum band index
INTEGER*8, parameter :: nelec=65   !electrons per spin componen
INTEGER*8, parameter :: nk=KKK*KKK      !number of k-points
INTEGER*8, parameter :: nq=QQQ*QQQ      !number of q-points
INTEGER*8, parameter :: nnu=21     !number of modes included in the gkk file
INTEGER*8, parameter :: nu_offset=OFFSET !to save memory 
REAL*8, parameter :: Suc=2.73e-19       ! Supercell are in m^2 
REAL*8, parameter :: Ef=4.67d0      ! Fermi level in eV
REAL*8, parameter :: Ex=3.1d0 ! energy of phothons in eV
REAL*8, parameter :: F = 38.d0*J2eV ! Fluence in Joules/m2
end module physical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module in files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module files
CHARACTER (len=7) :: gkkfile='gkk.dat'
CHARACTER (len=11) :: stretchfile='stretch.dat'
CHARACTER (len=8) :: tempsfile='temps.in'
CHARACTER (len=6) :: etafile='eta.in'
end module files


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM sum_eph
use constants
use physical
use files
IMPLICIT NONE
CHARACTER (len=10) :: filemode='eta.in'
INTEGER nbnd,ntemps,neta
INTEGER q,k,iFM,ifric,imass,ihot
REAL*8, DIMENSION(:,:,:,:,:,:) , ALLOCATABLE :: gkk
REAL*8, DIMENSION(:,:,:,:) , ALLOCATABLE :: stretch
REAL*8, DIMENSION(:) , ALLOCATABLE :: Te,Tl,EfT,etaim
REAL*8, DIMENSION(:,:) ,ALLOCATABLE :: massqnu
COMPLEX*16, DIMENSION(:,:,:,:) ,ALLOCATABLE :: sigmaFM
complex*16 :: eta

print*,"**************************************************"
print*,"In this program we sum to e-ph matrix elements"
print*,"To evaluate the first and second order correction"
print*,"to the CO stretch mode frequency and linewidth"
print*,"**************************************************"
print*, " "
print*,"We start allocating e-ph matrix elements"
print*, " "
!number of bands considered
nbnd=nbndmax-nbndmin+1	
! total number of k-an q-points
!nk=nk*nk
!nq=nq*nq

print*,"number of bands = ", nbnd
print*,"number of k-points = ", nk
print*,"number of q-points = ",  nq
print*,"number of modes printed =", nnu

ALLOCATE(gkk(nq,nk,nnu,nbnd,nbnd,4))
ALLOCATE(stretch(nk,nbnd,nbnd,4))

!initialize the gkk matrix
gkk = 0.d0
!initialize stretch matrix
stretch = 0.d0

print*," "
print*, "Fermi Level =", Ef, "eV"


print*," "
print*,"**************************************************"
! call subroutine read gkk and stretch files
CALL read_eph(nbnd,gkk,stretch)

print*,"**************************************************"
print*," "

!Read list of electronic and lattice temperatures
!Number of pairs of temperatures in file
print*, "Temperatures read from ", tempsfile
open(11,file=tempsfile, status = 'old')
read(11,*) ntemps
ALLOCATE(Te(ntemps),Tl(ntemps))
ALLOCATE(EfT(ntemps))
print*,"electronic and lattice temperatures that will be computed (in K)"
do k = 1, ntemps
read(11,*) Te(k),Tl(k)
print*, Te(k),Tl(k)
print*, " "
end do
eta=(0.d0,0.030d0)
print*, "value of eta = ", eta,"eV"

close(11)


print*," "
print*,"**************************************************"
!call subroutine that computes the chemical pontential as function of Te
CALL EF_T(nbnd,gkk(1,:,nnu,1,:,1),Te,ntemps,EfT)
print*,"**************************************************"
print*," "

print*," "
print*,"**************************************************"
!call subroutine that computes the firs order correction (bare phonon self-energy)
!CALL first_order(nbnd,gkk,Te,ntemps,eta)
CALL first_order_chem(nbnd,gkk,Te,ntemps,eta,EfT)
! In this routine we compute the dynamical susceptibility
CALL susceptibility(nbnd,gkk,Te,ntemps,eta,EfT)


print*," "
print*,"**************************************************"
!call subroutine that computes the second order 
! electron mediated phonon-phonon coupling (EMPP)
do k=1,nnu
write(filemode,'("mode",I2.2,".dat")') k
print*, "Computing mode:",k
CALL second_order_chem(nbnd,gkk,stretch,Te,Tl,ntemps,eta,k,k,filemode,EfT)
print*," "
end do
print*,"**************************************************"
print*," "



end program sum_eph

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE READ e-ph file
!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE read_eph(nbnd,gkk,stretch)
use physical, only : nk,nq,nnu,nu_offset
use files, only : gkkfile, stretchfile
IMPLICIT NONE
integer ppkq,ppkq2
INTEGER nbnd
INTEGER q,k,j,nu,m,n
REAL*8, DIMENSION(nq,nk,nnu,nbnd,nbnd,4):: gkk
REAL*8, DIMENSION(nk,nbnd,nbnd,4):: stretch
REAL*8 enk, enkq,omega,g
integer t1,t2,time

t1=time()
print*,"Reading the input file: ", gkkfile
open (12, file = gkkfile, status = 'old')


ppkq=nbnd*nbnd*nnu  !number points per k-q pair
print*,"number of elements per (k-q) pair = ", ppkq


do q = 1, nq
	print*,"q = ",q," of ",nq
	do k = 1, nk
		do j = 1, ppkq
			!read(12,*) m,n,nu,enk,enkq,omega,g
			read(12,*) n,m,nu,enk,enkq,omega,g
                        nu = nu-nu_offset
			! everythin in eV
			omega = omega*0.001d0
			g = g*0.001d0
			gkk(q,k,nu,m,n,1)= enk
			gkk(q,k,nu,m,n,2)= enkq
			gkk(q,k,nu,m,n,3)= omega
			gkk(q,k,nu,m,n,4)= g
		end do
	end do
end do
close(12)


print*,"Reading the input file: ", stretchfile
open (12, file = stretchfile, status = 'old')

ppkq2=nbnd*nbnd  !number points per k-q pair (nnu=1 in this case)
print*,"number of elements per (k-q) pair = ", ppkq2

! Now there is no loop in q as only q=0 is needed.
do k = 1, nk
	do j = 1, ppkq2
		!read(12,*) m,n,nu,enk,enkq,omega,g
		read(12,*) n,m,nu,enk,enkq,omega,g
		nu = nu-nu_offset
		! everythin in eV
		omega = omega*0.001d0
		g = g*0.001d0
		stretch(k,m,n,1)= enk
		stretch(k,m,n,2)= enkq
		stretch(k,m,n,3)= omega
		stretch(k,m,n,4)= g
	end do
end do
close(12)

print*,"Finished reading!"
t2=time()
print*,"elapsed time = ", t2-t1, "sec"
end subroutine read_eph


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPUTE sum of ocupancies
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sum_occ(nbnd,enk,Te,ntemps,it,EfT,Ne)
use constants
use physical
IMPLICIT NONE
INTEGER nbnd,ntemps
INTEGER q,k,j,nu,m,n,it
REAL*8, DIMENSION(nk,nbnd):: enk
REAL*8 enkq,omega,g, EfT
REAL*8, DIMENSION(ntemps) :: Te
REAL*8 kbTe
REAL*8 fmuk,occup,Ne
integer t1, t2,time


!initial time
!t1=time()
!:print*, " "
!print*, "Let's compute how the Fermi level changes with T_e"

!Lets start computing the sum of occupancies at each T.
! we do it for q=0 as it should be equal for each q
! It does not depend also on nu 
q = 1
nu= nnu
m = 1
!Initialize the values: sum of occ and E_F(T_e)
Ne = 0

kbTe=Te(it)*kb
do n = 1, nbnd
	do k = 1, nk
		fmuk=occup(enk(k,n)-EfT,kbTe)
                Ne= Ne + fmuk
	end do
end do
Ne=Ne/nk

end subroutine sum_occ

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compute the number of states in a window of energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE states_window(nbnd,enk,EfT,Ew,ns)
use physical, only : nk
IMPLICIT NONE
INTEGER nbnd,n,k
REAL*8, DIMENSION(nk,nbnd):: enk
REAL*8 EfT,ener,Ew,ns


ns = 0
if(Ew.lt.0) then
do n = 1, nbnd
	do k = 1, nk
                ener=enk(k,n)-EfT
                if ((ener.gt.Ew).and.(ener.lt.0)) then ! holes excited from [Ef-Ex,Ef]
	        ns= ns + 1.d0
                end if
        end do
end do
else ! Ex > 0 
do n = 1, nbnd
	do k = 1, nk
                ener=enk(k,n)-EfT
                if ((ener.gt.0).and.(ener.lt.Ew)) then ! holes excited from [Ef-Ex,Ef]
	        ns= ns + 1.d0
                end if
        end do
end do
end if

Ns=2*Ns/nk

end subroutine states_window


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPUTE E_F AS FUNCTION OF T_E
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EF_T(nbnd,enk,Te,ntemps,EfT)
use constants
use physical
IMPLICIT NONE
INTEGER nbnd,ntemps
INTEGER q,k,j,nu,m,n,it
!REAL*8, DIMENSION(nq,nk,nnu,nbnd,nbnd,4):: gkk
REAL*8, DIMENSION(nk,nbnd):: enk
REAL*8 enkq,omega,g
REAL*8, DIMENSION(ntemps) :: Te, EfT,NeT
REAL*8 kbTe
REAL*8 DeltaEf, docc, doccini
REAL*8  fmuk,fmukq,occup,fac1,Ne
REAL*8 Ef1, Ef2, docc1, docc2, Efnew, doccnew !Bolzano variables
integer t1, t2,time


!initial time
t1=time()
print*, " "
print*, "Let's compute how the Fermi level changes with T_e"


!Lets start computing the sum of occupancies at each T.
! we do it for q=0 as it should be equal for each q
! It does not depend also on nu 
q = 1
nu= nnu
m = 1
!Initialize the values: sum of occ and E_F(T_e)
Ne = 0
EfT = Ef
DeltaEf=0.01
!First we compute the occupancy at T=0
print*,"Computing the sum of occupancies at T=1 K..."
it=1
CALL sum_occ(nbnd,enk,Te,ntemps,it,Te(it),NeT(it))

print*,"The sum of occupancies is",NeT(it)
print*,"number of electrons/spin", nelec
print*,""
open (13, file = 'muT.dat', status = 'unknown', access = 'append')
!Loop in electronic temperatures
print*, "Te                    Re[Pi1]               Im[Pi1]"
do it = 1, ntemps
	!kbTe=Te(it)*kb
        !EfT(it)=EfT(it-1) !use result of previous Te to converge faster
        print*,"Electronic Temperature =",Te(it)
	!print*,"calculating for electronic temperature",Te(it),"K"
        CALL sum_occ(nbnd,enk,Te,ntemps,it,EfT(it),NeT(it))
        print*, "The initial sum of occupancies is and Ef="
        print*, NeT(it),EfT(it)
        !FINITE DIFFERENCES METHOD
        !docc= NeT(it)-nelec
        !doccini = docc
        !IF (docc.lt.0.d0) then
        !  DeltaEf = abs(DeltaEf)
        !ELSE
        !  DeltaEf = -abs(DeltaEf)
        !END IF
        !DO WHILE (docc*doccini.GT.0.d0)
        !EfT(it) = EfT(it) + DeltaEf
        !CALL sum_occ(nbnd,gkk,Te,ntemps,it,EfT(it),NeT(it))
        !docc= NeT(it) - nelec
        !END DO
        !print*,"Converged value for nocc, E_f(T_e) = ", NeT(it),EfT(it) 
        
        ! Bolzano mehod
        docc1= NeT(it)-nelec
        Ef1=EfT(it)
        IF (docc1.lt.0.d0) then
          Ef2 = Ef1 + 2.0d0
        ELSE
          Ef2 = Ef1 - 2.0d0
        END IF
        ! First Bolzano move
        CALL sum_occ(nbnd,enk,Te,ntemps,it,Ef2,NeT(it))
        docc2=NeT(it)-nelec
        IF(docc1*docc2.GT.0.d0) then
        print*,"Bolzano Faliled.."
        stop
        ELSE
        Efnew = (Ef1+Ef2)/2.d0
        END IF
         CALL sum_occ(nbnd,enk,Te,ntemps,it,EFnew,NeT(it))
         doccnew= NeT(it) - nelec
        j=0
        DO WHILE (abs(Ef1-Ef2).GT.1.d-6)

         IF(doccnew*docc1.le.0.d0) THEN
          Ef2=Efnew
          docc2=doccnew
          Efnew=(Ef1+Ef2)/2.d0
         ELSE
          Ef1=Efnew
          docc1=doccnew
          Efnew=(Ef1+Ef2)/2.d0
         END IF
         CALL sum_occ(nbnd,enk,Te,ntemps,it,EFnew,NeT(it))
         doccnew= NeT(it) - nelec
         j=j+1
        END DO
        EfT(it) = Efnew
       
        print*,"Converged value for nocc, E_f(T_e)" 
        write(13,*) Te(it),EfT(it)
        print*, NeT(it),EfT(it)
        print*,"Convergence achieved in", j,"iterations"
        print*,""
end do
close(13)

t2=time( )
print*,"elapsed time = ", t2-t1, "sec"

end subroutine EF_T


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE FIRST ORDER INCLUDING CHEMICAL POTENTIAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE first_order_chem(nbnd,gkk,Te,ntemps,eta,EfT)
use constants
use physical
IMPLICIT NONE
INTEGER nbnd,ntemps
INTEGER q,k,j,nu,m,n,it
REAL*8, DIMENSION(nq,nk,nnu,nbnd,nbnd,4):: gkk
REAL*8 enk, enkq,omega,g
REAL*8, DIMENSION(ntemps) :: Te,EfT
REAL*8 kbTe
REAL*8  fmuk,fmukq,occup,fac1
COMPLEX*16:: eta, fac2,pi_ep,pi_ep0,fac20
integer t1, t2,time


!initial time
t1=time()
print*, " "
print*, "Let's compute the first order"
print*, "here we include the chemical potential dependence with T mu(T)"


!In this routine we only use the stretch mode
! and only gamma point!!
q = 1
nu= nnu
!Loop in electronic temperatures
print*, "Te                    Re[Pi1]               Im[Pi1]"
do it = 1, ntemps
	kbTe=Te(it)*kb
	!print*,"calculating for electronic temperature",Te(it),"K"
	pi_ep=(0.d0,0.d0)
	pi_ep0=(0.d0,0.d0)
	do m = 1, nbnd
		do n = 1, nbnd
			if(m.ne.n) then
			do k = 1, nk
				enk= gkk(q,k,nu,m,n,1)
				enkq= gkk(q,k,nu,m,n,2)	
				omega= gkk(q,k,nu,m,n,3)
				g= gkk(q,k,nu,m,n,4)	
				fmuk=occup(enk-EfT(it),kbTe)
				fmukq=occup(enkq-EfT(it),kbTe)	
				fac1=fmuk-fmukq
				fac2 = omega + enk - enkq + eta
				fac20=	enk - enkq + eta
				pi_ep = pi_ep + g*g*fac1/fac2	
				pi_ep0 = pi_ep0 +g*g*fac1/fac20		
			end do
		end if
		end do
	end do
	open (13, file = 'results_first_order_muT.dat', status = 'unknown', access = 'append')
	pi_ep0=pi_ep-pi_ep0
	pi_ep0=pi_ep0/nk*eV2cm1
	pi_ep=pi_ep/nk*eV2cm1
	write(*,'(4f12.4)') Te(it), real(pi_ep0),real(pi_ep), imag(pi_ep)
	write(13,'(4f12.4)') Te(it), real(pi_ep0),real(pi_ep),imag(pi_ep)
	close(13)
end do
t2=time()
print*,"elapsed time = ", t2-t1, "sec"

end subroutine first_order_chem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE that computes susceptibility
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE susceptibility(nbnd,gkk,Te,ntemps,eta,EfT)
use constants
use physical
IMPLICIT NONE
INTEGER nbnd,ntemps
INTEGER q,k,j,nu,m,n,it
REAL*8, DIMENSION(nq,nk,nnu,nbnd,nbnd,4):: gkk
REAL*8 enk, enkq,omega,g
REAL*8, DIMENSION(ntemps) :: Te,EfT
REAL*8 kbTe
REAL*8  fmuk,fmukq,occup,fac1
COMPLEX*16:: eta, fac2,pi_ep,pi_ep0,fac20
integer t1, t2,time


!initial time
t1=time()
print*, " "
print*, "Let's compute the susceptibility"
print*, "the chemical potential dependence with T mu(T) is considered"


!In this routine we only use the stretch mode
! and only gamma point!!
q = 1
nu= nnu
!Loop in electronic temperatures
print*, "Te                    Re[Pi1]               Im[Pi1]"
do it = 1, ntemps
	kbTe=Te(it)*kb
	!print*,"calculating for electronic temperature",Te(it),"K"
	pi_ep=(0.d0,0.d0)
	pi_ep0=(0.d0,0.d0)
        do q = 1, nq
	do m = 1, nbnd
		do n = 1, nbnd
			if(m.ne.n) then
			do k = 1, nk
				enk= gkk(q,k,nu,m,n,1)
				enkq= gkk(q,k,nu,m,n,2)	
				omega= gkk(q,k,nu,m,n,3)
				!g= gkk(q,k,nu,m,n,4)	
				fmuk=occup(enk-EfT(it),kbTe)
				fmukq=occup(enkq-EfT(it),kbTe)	
				fac1=fmuk-fmukq
				fac2 = omega + enk - enkq + eta
				fac20= omega + enk - enkq + eta
				pi_ep = pi_ep + fac1/fac2	
				if(q.eq.1)then ! compute the q=0 part only
                                pi_ep0 = pi_ep0 +fac1/fac20		
                                end if
			end do
		end if
		end do
	end do
        end do
	open (13, file = 'susceptibility.dat', status = 'unknown', access = 'append')
	!pi_ep0=pi_ep-pi_ep0
	pi_ep0=pi_ep0/nk*eV2cm1
	pi_ep=pi_ep/nk*eV2cm1
	write(*,'(5E12.4)') Te(it), real(pi_ep0),imag(pi_ep0),real(pi_ep), imag(pi_ep)
	write(13,'(5E12.4)') Te(it), real(pi_ep0),imag(pi_ep0),real(pi_ep),imag(pi_ep)
	close(13)
end do
t2=time()
print*,"elapsed time = ", t2-t1, "sec"

end subroutine susceptibility




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE SECOND ORDER Including chemical potential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE second_order_chem(nbnd,gkk,stretch,Te,Tl,ntemps,eta,& 
mmin,mmax,filemode,EfT)
use constants
use physical
IMPLICIT NONE
CHARACTER (len=10) :: filemode
INTEGER mmin, mmax
INTEGER nbnd,ntemps
INTEGER q,k,j,nu,m,n,it
INTEGER qp,nup,s, sp
REAL*8, DIMENSION(nq,nk,nnu,nbnd,nbnd,4):: gkk
REAL*8, DIMENSION(nk,nbnd,nbnd,4):: stretch
REAL*8 enk, enkq,omega,g
REAL*8 emk, emkq,omegap,gp,deltaE, enkqEf,omegap_s,omegap_sps
REAL*8, DIMENSION(ntemps) :: Te, Tl,EfT
REAL*8 kbTe,kbTl
REAL*8  occup,fac1
REAL*8 bose,occ1, occ3sp,fac3
COMPLEX*16 :: eta, fac2,pi_ep2
COMPLEX*16 :: occ4sp,fac4,frac1,frac2
integer t1, t2,time


!initial time
t1=time()
print*, " "
print*, "Let's compute the second order"


!We compute the effect on the stretch mode
! and only gamma point!!
nu= nnu
!Loop in electronic temperatures
print*, "Te              Tl               Re[Pi2]            Im[Pi2]"
do it = 1, ntemps
	kbTe=Te(it)*kb
	kbTl=Tl(it)*kb
	!print*,"calculating for electronic temperature",Te(it),"K"
	pi_ep2=(0.d0,0.d0)
	do m = 1, nbnd-1 !n and m interchanged with respect to other routines
		do k = 1, nk
			!emk= gkk(1,k,nu,m,m,1)
			!emkq= gkk(1,k,nu,m,m,2)	
			!omega= gkk(1,k,nu,m,m,3)
			!g= gkk(1,k,nu,m,m,4)		
			emk   = stretch(k,m,m,1)
			emkq  = stretch(k,m,m,2)	
			omega = stretch(k,m,m,3)
			g     = stretch(k,m,m,4)	
                        occ1 = occup(emk-EfT(it),kbTe)
			frac2=(0.d0,0.d0)
			do n = 1, nbnd-1
				do qp=1, nq
					do nup=mmin,mmax
						enk= gkk(qp,k,nup,n,m,1)
						enkq= gkk(qp,k,nup,n,m,2)	
						omegap= gkk(qp,k,nup,n,m,3)
						gp = gkk(qp,k,nup,n,m,4)
						deltaE=emk-enkq
						enkqEf=enkq-EfT(it)
						frac1=(0.d0,0.d0)
						do sp=-1,1,2
							occ3sp=occup(sp*enkqEf,kbTe)
							occ4sp=sp*deltaE + omega + eta
							do s = -1, 1,2
							omegap_s=omegap*s
							omegap_sps=omegap_s*sp
							fac3=occ1-occup((enkqEf-omegap_sps),kbTe)
							fac3=fac3*s*(bose(omegap_s,kbTl)+occ3sp)
							fac4=deltaE+omegap_sps
							fac4=fac4*omega*(occ4sp+omegap_s)
							frac1=frac1+fac3/fac4
							end do !end s loop
						end do !end sp loop
					frac2=frac2+frac1*gp*gp	
					end do !end nup loop
				end do  !end qp loop
			end do  ! end n loop
			pi_ep2 = pi_ep2 - frac2*g*g
		end do	!end k loop
	end do	! end m loop
	pi_ep2=pi_ep2/nk/nq*ev2cm1
	open (14, file = filemode, status = 'unknown', access = 'append')
	write(*,'(4f12.4)') Te(it),Tl(it), real(pi_ep2), imag(pi_ep2)
	write(14,'(4f12.4)') Te(it),Tl(it), real(pi_ep2),imag(pi_ep2)
	close(14)
end do	!end of Te, Tl loop
t2=time()
print*,"elapsed time = ", t2-t1, "sec"

end subroutine second_order_chem


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fermi-Dirac occupation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Function occup(ener,kbTe)  
implicit none  
REAL*8 occup,x,kbTe,ener  
x=ener/kbTe
if (x.lt.-600.d0) then
	occup = 1.d0
else if (x.gt.600.d0) then
	occup = 0.d0
else
	occup = 1.d0/(exp(x)+1.d0)
end if
RETURN  
end function occup



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! BOSE occupation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Function bose(ener,kbTl)  
implicit none  
REAL*8 bose,x,kbTl,ener  
x=ener/kbTl
if (x.lt.-600) then
	bose = 1.d0
else if (x.gt.600) then
	bose = 0.d0
else
	bose = 1.d0/(exp(x)-1.d0)
end if
RETURN  
end function  bose


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Derivative of Fermi occupation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Function doccup(ener,kbTe)  
implicit none  
REAL*8 doccup,x,kbTe,ener  
x=ener/kbTe
if (x.lt.-600.d0) then
	doccup = 0.d0
else if (x.gt.600.d0) then
	doccup = 0.d0
else
	doccup = -1.d0/kbTe*exp(x)/((exp(x)+1.d0))**2.d0
end if
RETURN  
end function doccup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Lorentzian function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Function Lorentz(etar, x)
use constants, only : pi
implicit none
REAL*8 Lorentz, etar, x
if (x.lt.-600.d0) then
	Lorentz = 0.d0
else if (x.gt.600.d0) then
	Lorentz = 0.d0
else
Lorentz = etar/pi/(x*x+etar*etar)
end if
RETURN  
end function Lorentz
