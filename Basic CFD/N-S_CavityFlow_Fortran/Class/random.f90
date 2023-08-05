! random.f90
program main
    implicit none
    real :: r(5)

    call random_seed()
    call random_number(r)

    print '(5(f8.6, " "))', r
end program main
