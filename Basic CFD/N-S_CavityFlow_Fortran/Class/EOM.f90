!==============================================================
!                                                             =
! <^> WARNING <^>:                                            =
! BEFORE USING THIS PROGRAM COMPILE Initialize-Final.f90      =
! WITH THE SAME VALUES OF Nx, Ny and Nz                       =
!                                                             =
!In this program we will start solving the EOM.               = 
!The Initialize-Final.f90 gives us all the initial info.      =
!The intial positions are coded in fort.4                     =
!The initial velocities are coded in fort.9                   =
!                                                             =
!Energies at each time step will be written in: fort.10       =
!Positions at each time step will be written in: fort.11      =
!                                                             =
!COMPILE:                                                     =
!gfortran tipos.90 EOM.f90                                    =
!./a.out                                                      =
!==============================================================

program EOM
use mcf_tipos
!Fixed values---------------------------
!This fixe values must be the same of  -
!Initialize-Final.f90                  -

integer, parameter :: Nx = 7                   !For this format of creating the positions use an odd number of N
integer, parameter :: Ny = 7                   !For this format of creating the positions use an odd number of N
integer, parameter :: Nz = 7                   !For this format of creating the positions use an odd number of N
integer, parameter :: N = Nx * Ny * Nz
real(kind=dp), parameter :: L = 1.38_dp
real(kind=dp), parameter :: pi = acos(-1.0)
character*2, parameter :: ELEM = 'N'           !Element of the particles
!--------------------------------------

!Some integers-------------------------
integer :: i,j,k
!--------------------------------------

!Position of the particles. Dim : N----
real(kind=dp), allocatable, dimension(:) :: Xhl
real(kind=dp), allocatable, dimension(:) :: Yhl
real(kind=dp), allocatable, dimension(:) :: Zhl
!--------------------------------------

!Velocities of the particles. Dim: N---
real(kind=dp), allocatable, dimension(:) :: Vx
real(kind=dp), allocatable, dimension(:) :: Vy
real(kind=dp), allocatable, dimension(:) :: Vz
!--------------------------------------

!Reals and integers for the equation solving
real(kind=dp), allocatable, dimension(:) :: X_zero !Bootstrap vector
real(kind=dp), allocatable, dimension(:) :: Y_zero !Bootstrap vector
real(kind=dp), allocatable, dimension(:) :: Z_zero !Bootstrap vector
real(kind=dp), allocatable, dimension(:) :: XX     !Time evolution position vector
real(kind=dp), allocatable, dimension(:) :: YY     !Time evolution position vector
real(kind=dp), allocatable, dimension(:) :: ZZ     !Time evolution position vector
real(kind=dp), allocatable, dimension(:) :: Fx     !Force in x direction
real(kind=dp), allocatable, dimension(:) :: Fy     !Force in y direction
real(kind=dp), allocatable, dimension(:) :: Fz     !Force in z direction
real(kind=dp) :: force                             !Potential force
real(kind=dp) :: r, r_six, xval, yval, zval        !Radial values for the force
real(kind=dp) :: Epot , Ek, Etot                   !Energy values
real(kind=dp), parameter :: dT = 0.00025_dp        !Time discretization
real(kind=dp), parameter :: r_min = 1.38_dp        !L-J potential minimun
real(kind=dp), parameter :: sigma = r_min * 0.5_dp **(1 / 6) !L-J potential sigma
integer, parameter :: n_step = 100000                !Time steps for the Eq. solver. Please even.
integer, parameter :: eps = 100.0_dp
!-------------------------------------------

!Radial distribution stuff------------------
real(kind=dp), allocatable, dimension(:) :: gr
integer, parameter :: n_his = 1000
real(kind=dp), parameter :: Edge = L * real(Nx)
real(kind=dp), parameter :: Density = real(N) / (edge ** 3) 
real(kind=dp) :: Deltagr 
real(kind=dp) :: vb
integer :: ig
!-------------------------------------------

!Define deltagr for radia distribution------

Deltagr = (edge * sqrt(3.0_dp) / 2.0_dp ) / real(n_his)

!-------------------------------------------


!Allocate things: Insert dim to vectors-----
allocate(Xhl(N))
allocate(Yhl(N))
allocate(Zhl(N))
allocate(Vx(N))
allocate(Vy(N))
allocate(Vz(N))
allocate(Fx(N))
allocate(Fy(N))
allocate(Fz(N))
allocate(X_zero(N))
allocate(XX(N))
allocate(Y_zero(N))
allocate(YY(N))
allocate(Z_zero(N))
allocate(ZZ(N))
allocate(gr(n_his))
!-------------------------------------------



!Start reading the inital positions and the velocities
open(unit=4)
do i=1, N
  read(4,'(3F20.10)') Xhl(i), Yhl(i), Zhl(i)
end do
close(unit=4)
open(unit=9)
do i=1, N
  read(9,'(3F20.10)') Vx(i), Vy(i), Vz(i)
end do
close(unit=9)
print*, "Initial positions and velocities have been readed succesfully"
print*, "============================================================="
print*, " "
!-----------------------------------------------------

!Fill the Bootstrap vector----------------------------
do i=1, N
 X_zero(i) = Xhl(i) - Vx(i) * dT
 Y_zero(i) = Yhl(i) - Vy(i) * dT
 Z_zero(i) = Zhl(i) - Vz(i) * dT
end do
!-----------------------------------------------------

!Radial distribution fucntion initialization----------
do i=1, n_his
 gr(i) = 0.0_dp
end do

!-----------------------------------------------------

!Solve the EOM----------------------------------------
print*, "Solving equation of motion. It might take a few minutes..."
print*, " "
do k=1, n_step
 if (k == n_step / 2) then
    print*, "Solving equation of motion. We are half way..."
    print*, "==============================================="
    print*, " "
 end if
 !Compute the forces for every i atom.
 Epot = 0.0_dp
 do i=1, N
  Fx(i) = 0.0_dp
  Fy(i) = 0.0_dp
  Fz(i) = 0.0_dp
 end do
 do i=1, N-1
  do j=i+1, N
  xval = Xhl(i) - Xhl(j)
  yval = Yhl(i) - Yhl(j)
  zval = Zhl(i) - Zhl(j)
  xval = xval - edge * anint(xval / edge)   ! PBC's
  yval = yval - edge * anint(yval / edge)   ! PBC's
  zval = zval - edge * anint(zval / edge)   ! PBC's
  r = sqrt(xval ** 2 + yval ** 2 + zval ** 2)
  r_six = (sigma / r )**6
  force = (48.0_dp * eps / (r ** 2)) * r_six * ( r_six - 0.5)
  Fx(i) = Fx(i) + xval * force
  Fx(j) = Fx(j) - xval * force
  Fy(i) = Fy(i) + yval * force
  Fy(j) = Fy(j) - yval * force
  Fz(i) = Fz(i) + zval * force
  Fz(j) = Fz(j) - zval * force
  Epot = Epot + r_six * (r_six - 1.0_dp) * 4.0_dp * eps
  end do
 end do
 !Forces done
 !Now we solve the equation of motion. Verlel Algorith
 Ekin = 0.0_dp
 do i=1, N
  XX(i) = 2 * Xhl(i) - X_zero(i) + Fx(i) * (dT**2)
  YY(i) = 2 * Yhl(i) - Y_zero(i) + Fy(i) * (dT**2)
  ZZ(i) = 2 * Zhl(i) - Z_zero(i) + Fz(i) * (dT**2)
  !Kinetic energy stuff
  Vx(i) = (XX(i) - X_zero(i)) / (2.0_dp * dT)
  Vy(i) = (YY(i) - Y_zero(i)) / (2.0_dp * dT)
  Vz(i) = (ZZ(i) - Z_zero(i)) / (2.0_dp * dT)
  Ek = Ek + (Vx(i) ** 2 + Vy(i) ** 2 + Vz(i) ** 2)
 end do
 !EOM solved for each time step
 !Rewrite the values for the next k step
 X_zero = Xhl
 Xhl = XX
 Y_zero = Yhl
 Yhl = YY
 Z_zero = Zhl
 Zhl = ZZ
 !Radial distrubution function stuff
 do i=1, N-1
  do j=1, N
  xval = Xhl(i) - Xhl(j)
  yval = Yhl(i) - Yhl(j)
  zval = Zhl(i) - Zhl(j)
  xval = xval - edge * anint(xval / edge)   ! PBC's
  yval = yval - edge * anint(yval / edge)   ! PBC's
  zval = zval - edge * anint(zval / edge)   ! PBC's
  r = sqrt(xval ** 2 + yval ** 2 + zval ** 2)
  if (r <= (edge * sqrt(3.0_dp) / 2.0_dp)) then
   ig = idint(r / Deltagr)
   gr(ig) = gr(ig) + 2.0_dp

  end if
  end do
 end do
 !Radial distrubition done
 !Energy stuff
 Ek = Ek / 2.0_dp
 Ek = Ek / N
 Epot = Epot / N
 Etot = Epot + Ek
 !Write as Jesus Ugalde wants
 !Energies at each time step(k*dT)
 open(unit=10)
 write(10,'(a6,4F20.10)')"Step:", k * dT, Ek, Epot, Etot
 !We write all the position in each every k * dT time step
 open(unit=11)
 write(11,'(i5,/)') N
 do i=1, N
  write(11, '(a2, 3(3x,F20.10))') ELEM, XX(i), YY(i), ZZ(i)
 end do
end do
!Normalization of g(r), average over trajectories and number of atoms
open(unit=12)
do i=1, n_his
 r = Deltagr * real(i) + 0.5_dp 
 vb = (4.0_dp / 3.0_dp) * pi * ((i+1) ** 3 - i ** 3) * Deltagr ** 3
 gr(i) = gr(i) / (real(n_step * N) * vb * Density)
 write(12,'(2d20.10)') r, gr(i)
end do
print*, "Updated energies at each time step in: fort.10"
print*, "=============================================="
print*, " "
print*, "Updated X, Y and Z positions at each time step in: fort.11"
print*, "=========================================================="
print*, " "
print*, "Total simulation time steps: ", n_step
print*, "======================================"
print*, " "
print*, "Total simulation time: ", n_step*dT, "s"
print*, "======================================="
print*, " "
print*, "Values of the g(r) function in: ", "fort.12"
print*, "======================================="
print*, " "



close(unit=10)
close(unit=11)
close(unit=12)
deallocate(Xhl)
deallocate(Yhl)
deallocate(Zhl)
deallocate(Vx)
deallocate(Vy)
deallocate(Vz)
deallocate(Fx)
deallocate(Fy)
deallocate(Fz)
deallocate(X_zero)
deallocate(XX)
deallocate(Y_zero)
deallocate(YY)
deallocate(Z_zero)
deallocate(ZZ)
deallocate(gr)
end program EOM
