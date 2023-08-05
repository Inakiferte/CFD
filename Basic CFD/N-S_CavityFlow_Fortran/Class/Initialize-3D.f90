
!----------------------------------------------------
!This program plots a 3D homogeneous cube           -
!gfortran tipos.f90 Initialize-3D.f90               -
!.a.out                                             -
!xmakemol -f fort.3                                 -
!----------------------------------------------------

program initial
use mcf_tipos

!Fixed values-------------------------
integer, parameter :: Nx = 5          !For this format of creating the positions use an odd number of N
integer, parameter :: Ny = 5          !For this format of creating the positions use an odd number of N
integer, parameter :: Nz = 5
integer, parameter :: N = Nx*Ny*Nz
real(kind=dp), parameter :: L = 1.5_dp
character*2, parameter :: ELEM = 'N'
!-------------------------------------

!Other reals--------------------------
real(kind=dp) :: V
!-------------------------------------

!Other integers-----------------------
integer :: i,jx,jy,jz,w,k
!-------------------------------------

!Position matrix(vector)--------------
real(kind=dp),allocatable,  dimension(:) :: X
real(kind=dp),allocatable,  dimension(:) :: Y
real(kind=dp),allocatable,  dimension(:) :: Z
!-------------------------------------

allocate(X(Nx))
allocate(Y(Ny))
allocate(Z(Nz))



!Define de volume---------------------
V = N*L
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


!Write in a file (fort.3)-------------

write(3, '(i4,/)') N
do, i=1, Nx
 do w=1, Ny
  do k=1, Nz
  write(3,'(a2,3F20.10)') ELEM, X(i), Y(w), Z(k)
  end do
 end do
end do

!--------------------------------------

deallocate(X)
deallocate(Y)
deallocate(Z)
end program initial

