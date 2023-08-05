!----------------------------------------------------
!This program plots a 2D homogeneous square         -
!gfortran tipos.f90 Initialize-2D.f90               -
!.a.out                                             -
!xmakemol -f fort.2                                 -
!----------------------------------------------------
program initial
use mcf_tipos

!Fixed values-------------------------
integer, parameter :: Nx = 5          !For this format of creating the positions use an odd number of N
integer, parameter :: Ny = 5          !For this format of creating the positions use an odd number of N
integer, parameter :: N = Nx * Ny
real(kind=dp), parameter :: L = 1.5_dp
character*2, parameter :: ELEM = 'N'
!-------------------------------------

!Other reals--------------------------
real(kind=dp) :: V
!-------------------------------------

!Other integers-----------------------
integer :: i,jx,jy,z
!-------------------------------------

!Position matrix(vector)--------------
real(kind=dp),allocatable,  dimension(:) :: X
real(kind=dp),allocatable,  dimension(:) :: Y
!-------------------------------------

allocate(X(Nx))
allocate(Y(Ny))


!Define de volume---------------------
V = N * L
!-------------------------------------



!Generate the positions---------------

jx = ( Nx + 1) / 2
jy = ( Ny + 1) / 2
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

!-------------------------------------


!Write in a file (fort.2)-------------

write(2, '(i4,/)') N
do i=1, Nx
 do z=1, Ny
  write(2,'(a2,3F20.10)') ELEM, X(i), Y(z), 0.0_dp
 end do
end do

!--------------------------------------

deallocate(X)
deallocate(Y)
end program initial

