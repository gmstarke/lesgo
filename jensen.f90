module jensen_m
use minimized_class
use types, only : rprec

type, extends(Minimized) :: jensen_t
    integer :: N
    real(rprec) :: D, alpha, Ui, Pref, Pfarm
    real(rprec), dimension(:), allocatable :: s, k, Ctp0, e, u
contains
    procedure, public :: eval
    procedure, public :: run
end type jensen_t

interface jensen_t
    module procedure :: constructor
end interface jensen_t

contains

function constructor(Ui, s, Ctp0, k, D, alpha, Pref, i_e) result(this)
    use types, only : rprec
    implicit none

    type(jensen_t) :: this
    real(rprec), dimension(:), intent(in) :: s, Ctp0, k, i_e
    real(rprec), intent(in) :: Ui, D, alpha, Pref

    this%N = size(s)
    allocate(this%s(this%N))
    allocate(this%k(this%N))
    allocate(this%Ctp0(this%N))
    allocate(this%e(this%N))
    allocate(this%u(this%N))

    this%D = D
    this%alpha = alpha
    this%Ui = Ui
    this%s = s
    this%k = k
    this%Ctp0 = Ctp0
    this%Pref = Pref
    this%e = i_e
end function constructor

subroutine eval(this, x, f, g)
    use types, only : rprec
    implicit none
    class(jensen_t), intent(inout) :: this
    real(rprec), dimension(:), intent(in) :: x
    real(rprec), intent(inout) :: f
    real(rprec), dimension(:), intent(inout) :: g
    real(rprec), dimension(:), allocatable :: Ctp
    type(jensen_t) :: jti
    integer :: i
    
    f = this%run(x)
    
    jti = jensen_t(this%Ui, this%s, this%Ctp0, this%k, this%D, this%alpha,     &
                   this%Pref, this%e)
    allocate(Ctp(this%N))
    do i = 1, this%N
        Ctp = x
        Ctp(i) = Ctp(i) + sqrt(epsilon(f))
        g(i) = jti%run(Ctp)
    end do
    
    g = (g - f) / sqrt(epsilon(f))

end subroutine eval

function run(this, Ctp) result(f)
    use types, only : rprec
    implicit none
    class(jensen_t), intent(inout) :: this
    real(rprec), dimension(:), intent(in) :: Ctp
    real(rprec) :: f
    real(rprec), dimension(:), allocatable :: a, du, P
    integer :: i, j
    
    allocate(a(this%N))
    allocate(du(this%N))
    allocate(P(this%N))
    a = Ctp / (4._rprec + Ctp)
    du = 0._rprec
    
    do i = 1, this%N
        do j = 1, i-1
            du(i) = du(i) + ( 2 * this%Ui * a(j) / (1._rprec + 2._rprec *      &
                this%k(j) * (this%s(i) - this%s(j)) / this%D) **2 )**2
        end do
    end do
    
    du = sqrt(du)
    this%u = this%Ui - du
    P = Ctp * (this%u + this%e)**3
    this%Pfarm = sum(P)
    f = (1._rprec - this%Pfarm / this%Pref)**2 + this%alpha * sum((Ctp - this%Ctp0)**2)

end function run 

end module jensen_m