
!----------------------------------------------------
!This program initialize all the parameters of our  -
!system. This program plots a 3D homogeneous cube,  -
!which is the initial state of the system.          -
!                                                   -
!It also adds the initial velocities, which would   -
!be the initial moementum since we assume           -
!that m = 1.                                        -
!                                                   -
!With the velocities, it calculates the histogram   -
!of how much atoms are with each velocities in a    -
!Delta range. It also calculates de histogram       -
!for the radial position, the so called radial      -
!distribution function.                             -
!                                                   -
!The total kinetic energy per particle and the      -
!total kinetic energy of the system is also         -
!computed.                                          -
!                                                   -
!COMPILE THE HOLE PROGRAM:                          -
!gfortran tipos.f90 Initialize-3DRan.f90            -
!.a.out                                             -
!                                                   -
!COMPILE(cube):                                     -
!gfortran tipos.f90 Initialize-3DRan.f90            -
!.a.out                                             -
!xmakemol -f fort.3                                 -
!                                                   -
!Positiion file:                                    -
!fort.4                                             -
!                                                   -
!Velocity file:                                     -
!fort.9                                             -
!                                                   -
!PLOT(histogram): Velocities                        - 
!gfortran tipos.f90 Initialize-3DRan.f90            -
!.a.out                                             -
!gnuplot                                            -
!plot 'fort.7' u 1:2 w l                            -
!                                                   -
!PLOT(histogram): Position                          -                       
!gfortran tipos.f90 Initialize-3DRan.f90            -
!.a.out                                             -
!gnuplot                                            -
!plot 'fort.8' u 1:2 w l                            -
!----------------------------------------------------

program initial
use mcf_tipos
implicit none
!Fixed values-------------------------
integer, parameter :: Nx = 7                   !For this format of creating the positions use an odd number of N
integer, parameter :: Ny = 7                   !For this format of creating the positions use an odd number of N
integer, parameter :: Nz = 7                   !For this format of creating the positions use an odd number of N
integer, parameter :: N = Nx * Ny * Nz
real(kind=dp), parameter :: L = 1.38_dp
character*2, parameter :: ELEM = 'N'
!-------------------------------------

!Other reals--------------------------
real(kind=dp) :: V                              !Volume
real(kind=dp) :: r                              !For the random number
!-------------------------------------

!Other integers-----------------------
integer :: i,jx,jy,jz,w,k,j
!integer, parameter :: iseed = 1234             !For generation of a random number
!-------------------------------------

!Position matrix(vector)--------------
real(kind=dp), allocatable, dimension(:) :: X
real(kind=dp), allocatable, dimension(:) :: Y
real(kind=dp), allocatable, dimension(:) :: Z
!-------------------------------------

!Position vector, hole system---------
real(kind=dp), allocatable, dimension(:) :: Xhl
real(kind=dp), allocatable, dimension(:) :: Yhl
real(kind=dp), allocatable, dimension(:) :: Zhl
!-------------------------------------

!Velocity vectors---------------------
real(kind=dp), allocatable, dimension(:) :: Vx
real(kind=dp), allocatable, dimension(:) :: Vy
real(kind=dp), allocatable, dimension(:) :: Vz
real(kind=dp), allocatable, dimension(:) :: Vsq !Square root of the velocity
!-------------------------------------

!Kinetic energy-----------------------
real(kind=dp) :: Ek                             !Total kinetic energy
real(kind=dp), allocatable, dimension(:) :: T   !Kinetic energy per particle
!-------------------------------------

!Velocity Histogram variables---------
real(kind=dp) :: Delta
real(kind=dp) :: va
integer, parameter :: nhis = 1000               !Times we divide the delta
real(kind=dp), allocatable, dimension(:) :: gv !Function to count histogram times
!-------------------------------------

!Radial function (histogram of positons)
real(kind=dp) :: Deltar
real(kind=dp) :: r_dis                         !Radial distance
real(kind=dp), allocatable, dimension(:) :: gr !Function to count histogram times
real(kind=dp) :: vol
real(kind=dp), parameter :: pi = acos(-1.0) 
!--------------------------------------



!Potential energy variables-----------
real(kind=dp) :: Epot, xx, yy, zz, rv, r_six
real(kind=dp), parameter :: d_zero = L
real(kind=dp), parameter :: sigma = d_zero * 0.5_dp ** (1 / 6)
real(kind=dp), parameter :: eps = 0.01_dp
!-------------------------------------


allocate(X(Nx))
allocate(Y(Ny))
allocate(Z(Nz))
allocate(Xhl(N))
allocate(Yhl(N))
allocate(Zhl(N))
allocate(Vx(N))
allocate(Vy(N))
allocate(Vz(N))
allocate(Vsq(N))
allocate(gv(nhis))
allocate(T(N))
allocate(gr(nhis))

!Define de volume---------------------
V = N * L
!-------------------------------------



!Generate the positions---------------

jx = ( Nx + 1) / 2
jy = ( Ny + 1) / 2
jz = ( Nz + 1) / 2


X(1) = 0.0_dp

do i=2, jx
 X(i) = X(i - 1) + L
 X(jx + i - 1) = - X(i)
end do

Y(1) = 0.0_dp

do i=2, jy
 Y(i) = Y(i - 1) + L
 Y(jy + i - 1) = - Y(i)
end do

Z(1) = 0.0_dp

do i=2, jz
 Z(i) = Z(i - 1) + L
 Z(jz + i - 1) = - Z(i)
end do

!-------------------------------------

!Generate velocities(randomly from -0.5 to 0.5)------------------
call random_seed()

do i=1, N
 call random_number(r)
 Vx(i) = r - 0.5_dp
end do

do i=1, N
 call random_number(r)
 Vy(i) = r - 0.5_dp
end do

do i=1, N
 call random_number(r)
 Vz(i) = r - 0.5_dp
end do

!---------------------------------------------------------------

!Write in a file (fort.3)-------------

write(3, '(i4,/)') N
do i=1, Nx
 do w=1, Ny
  do k=1, Nz
  write(3,'(a2,3F20.10)') ELEM, X(i), Y(w), Z(k)
  end do
 end do
end do
print*, "File for visualazing the cube in xmakemol has been generated: fort.3"
print*, "--------------------------------------------------------------------"
print*, " "
!-------------------------------------

!-----------------------------------------------------------------------------------
!                                                                                  -
!Now we are going to generate the hole position vectors, not Nx dimension but that -
!N dimension ones. This step can be done in just one do but for convenience we     -
!willl do it in two. Is not in any way more efficient to do it in two steps        -
!                                                                                  -
!-----------------------------------------------------------------------------------


!Readable position file---------------
do i=1, Nx
 do w=1, Ny
  do k=1, Nz
  write(4,'(3F20.10)') X(i), Y(w), Z(k)
  end do
 end do
end do
close(unit=4)
print*, "File with all the positions of the system has been generated: fort.4"
print*, "--------------------------------------------------------------------"
print*, " "
!-------------------------------------

!Readeable velocity file--------------
do i=1, N
   write(9,'(3F20.10)') Vx(i), Vy(i), Vz(i)
end do
close(unit=9)
print*, "File with all the velocities of the system has been generated: fort.9"
print*, "--------------------------------------------------------------------"
print*, " "
!-------------------------------------


!Read the hole positions--------------
open(unit=4)
do i=1, N
  read(4,'(3F20.10)') Xhl(i), Yhl(i), Zhl(i)
end do
close(unit=4)
!-------------------------------------

!--------------------------------------------------------------------------------
!                                                                               -
!Now we will generate the file for the plot of the  distribution of velocities. -
!An histogram for the velocities.                                               -
!                                                                               -  
!--------------------------------------------------------------------------------

!Compute the sqaure root of the velocities--------------------------------------

do i=1, N
 Vsq(i) = sqrt( Vx(i) ** 2 + Vy(i) ** 2 + Vz(i) ** 2 )
end do

!Define the delta where we will analyse how much Vsq-s are there-----------------

Delta = (sqrt(3.0_dp) / 2.0_dp) / nhis                                                 !sqrt(3)/2 is the maximun value of Vsq. So we
                                                                                       !divide it in nhis small lenght parts
!Generate gr---------------------------------------------------------------------

do i=1, N
 do j=1, nhis
  if (Vsq(i) >= ((j-1)*Delta) .and. Vsq(i) <= (j*Delta)) then
     gv(j) = gv(j) + 1.0_dp
  end if
 end do
end do

!Write the solution in a fort.7 file---------------------------------------------

do j=1, nhis
 write(7,'(2F20.10)') Delta * j /2.0, gv(j)
end do
close(unit=7)
print*, "File with the histogram of velocities has been generated: fort.7"
print*, "--------------------------------------------------------------------"
print*, " "
!---------------------------------------------------------------------------------
!                                                                                -
!Now we will compute the kinetic energy per particle and the total kinetic       -
!energy.                                                                         - 
!                                                                                -
!---------------------------------------------------------------------------------

!Kinetic energy per particle, assume m=1------------------------------------------
do i=1, N
 T(i) = 0.5_dp * Vsq(i) ** 2
end do

!Total kinetic energy-------------------------------------------------------------
Ek = sum(T)
print*, "Total kinetic energy of the system: ", Ek
print*, "--------------------------------------------------------------------"
print*, " "
!---------------------------------------------------------------------------------
!                                                                                -
!Now we will compute the potential energy of the system                          -
!energy.                                                                         - 
!                                                                                -
!Lennard-Jones potential: V(r) = 4*eps*((sigma/r)**12-(sigma/r)**6)              -
!                                                                                -
!---------------------------------------------------------------------------------

Epot = 0.0_dp

do i=1, N
 do j=i+1, N-1
  xx = Xhl(i) - Xhl(j)
  yy = Yhl(i) - Yhl(j)
  zz = Zhl(i) - Zhl(j)
  rv = sqrt(xx ** 2 + yy ** 2 + zz ** 2)
  r_six = (sigma / rv ) ** 6
  Epot = Epot + 4 * eps * r_six * (r_six - 1.0)
 end do
end do
print*, "Total potential (Lennard-Jones) energy of the system: ", Epot
print*, "--------------------------------------------------------------------"
print*, " "
!--------------------------------------------------------------------------------
!                                                                               -
!Now we will generate the radial function of the position (e.g. histogram       -
!of the positon).                                                               -
!Definition: https://en.wikipedia.org/wiki/Radial_distribution_function         -
!#Relations_involving_g(r)                                                      -
!                                                                               -
!--------------------------------------------------------------------------------

Deltar =sqrt( maxval(Xhl) ** 2 + maxval(Yhl) ** 2 + maxval(Zhl) ** 2) / nhis

do i=1, N
 do j=1, nhis
  r_dis = sqrt(Xhl(i) ** 2 + Yhl(i) ** 2 +Zhl(i) ** 2)
  if (r_dis >= ((j-1)*Deltar) .and. r_dis <= (j*Deltar)) then
     gr(j) = gr(j) + 1.0_dp
  end if
 end do
end do

do i=1, nhis
 vol = 4.0 / 3.0 * pi * Deltar
 gr(i) = gr(i) / (N * vol)
end do

do j=1, nhis
 write(8,'(2F20.10)') Deltar * j /2.0, gr(j)
end do
close(unit=8)
print*, "File with the histogram of positions has been generated: fort.8"
print*, "--------------------------------------------------------------------"
print*, " "

deallocate(X)
deallocate(Y)
deallocate(Z)
deallocate(Vx)
deallocate(Vy)
deallocate(Vz)
deallocate(Xhl)
deallocate(Yhl)
deallocate(Zhl)
deallocate(Vsq)
deallocate(gv)
deallocate(T)
deallocate(gr)
end program initial

