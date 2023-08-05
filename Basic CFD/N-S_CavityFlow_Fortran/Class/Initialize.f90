!----------------------------------------------------
!This program plots a 1D array                      -
!gfortran tipos.f90 Initialize.f90                  -
!.a.out                                             -
!xmakemol -f fort.1                                 -
!----------------------------------------------------


program initial
use mcf_tipos

!Fixed values-------------------------
integer, parameter :: N = 5          !For this format of creating the positions use an odd number of N
real(kind=dp), parameter :: L = 1.5_dp
character*2, parameter :: ELEM = 'N'
!-------------------------------------

!Other reals--------------------------
real(kind=dp) :: V
!-------------------------------------

!Other integers-----------------------
integer :: i,j
!-------------------------------------

!Position matrix(vector)--------------
real(kind=dp),allocatable,  dimension(:) :: X
!-------------------------------------

allocate(X(N))

!Define de volume---------------------
V = N * L
!-------------------------------------



!Generate the positions---------------

j = ( N + 1) / 2
X(1) = 0.0_dp

do i=2, j
 X(i) = X(i - 1) + L
 X(j + i - 1) = - X(i)
end do
!-------------------------------------


!Write in a file (fort.1)-------------

write(1, '(i4,/)') N
do, i=1, N
 write(1,'(a2,3F20.10)') ELEM, X(i), 0.0_dp, 0.0_dp
end do

!--------------------------------------

deallocate(X)
end program initial


