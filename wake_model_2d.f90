!
!
!!  Copyright (C) 2016  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   Wake Model Base Class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*******************************************************************************
module wake_model_base_m
!*******************************************************************************
use types, only : rprec
use util, only : logistic, softplus, gaussian
use messages
implicit none

private
public wake_model_base

type wake_model_base
    real(rprec), dimension(:), allocatable :: k    ! wake expansion coefficient
    real(rprec), dimension(:), allocatable :: x    ! streamwise coordinate
    real(rprec), dimension(:), allocatable :: y    ! spanwise coordinate
    real(rprec), dimension(:), allocatable :: z    ! vertical coordinate
    real(rprec), dimension(:), allocatable :: yu_inf   ! freestream flow in y
    real(rprec), dimension(:,:), allocatable :: s      ! turbine location
    real(rprec), dimension(:,:), allocatable :: G      ! Gaussian forcing function (turbine, space)
    integer, dimension(:), allocatable :: Gstart, Gstop
    real(rprec), dimension(:,:), allocatable :: d      ! dimensionless wake diameter (turbine, space)
    real(rprec), dimension(:,:), allocatable :: dp     ! d/dx of d (turbine, space)
    real(rprec), dimension(:,:), allocatable :: w      ! wake expansion function (turbine, space)
    real(rprec), dimension(:,:), allocatable :: fp     ! forcing prefactor f = fp*Ctp*u_d  (turbine, space)
    real(rprec), dimension(:,:,:), allocatable :: disk ! matrix holding the disk faces
    real(rprec), dimension(:,:,:,:), allocatable :: diskw! matrix holding the disk faces
!    real(rprec), dimension(:,:,:,:), allocatable ::Wf  ! wake expansion function (turbine, space)
    real(rprec) :: U_infty = 0                         ! inlet velocity
    real(rprec) :: Delta   = 0                         ! Gaussian forcing width
    real(rprec) :: Dia     = 0                         ! rotor diameter
    real(rprec) :: dx      = 0                         ! rotor diameter
    real(rprec) :: dy      = 0                         ! delta in spanwise direction
    real(rprec) :: dz      = 0                         ! delta in vertical direction
    integer     :: Nx      = 0                         ! Number of streamwise points
    integer     :: Ny      = 0                         ! Number of spanwise points
    integer     :: dNy     = 32                        ! Number of spanwise cells in turbine disk
    integer     :: dNz     = 32                        ! Number of vertical cells in turbine disk
    integer     :: N       = 0                         ! Number of turbines
    integer     :: nfree   = 0                         ! Number of un-waked turb
    real(rprec) :: sigma_o = 0  
    integer, dimension(:), allocatable       :: wake_num
    integer, dimension(:), allocatable       :: free_turbines
    integer, dimension(:,:), allocatable     :: Cx
    integer, dimension(:,:), allocatable     :: px 
    integer, dimension(:,:), allocatable     :: ymin   ! Wake boundary vector L
    integer, dimension(:,:), allocatable     :: ymax   ! Wake boundary vector U
    logical     :: isDimensionless = .false.
    real(rprec) :: LENGTH=0, VELOCITY=0, TIME=0, FORCE=0
contains
    procedure, public :: initialize_val
    procedure, public :: print
    procedure, public :: makeDimensionless
    procedure, public :: makeDimensional
    procedure, public :: computeGaussians
    procedure, public :: computeWakeExpansionFunctions
    procedure, public :: compute2Dwakes
    procedure, public :: compute_wake_disks
    procedure, public :: woketurbines
end type wake_model_base

interface wake_model_base
    module procedure :: constructor_val
end interface wake_model_base

contains

!*******************************************************************************
function constructor_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ny) result(this)
!*******************************************************************************
! Constructor for wake model with values given
implicit none
type(wake_model_base)                   :: this
real(rprec), intent(in)                 :: i_U_infty, i_Delta, i_Dia
real(rprec), dimension(:), intent(in)   :: i_k
real(rprec), dimension(:,:), intent(in) :: i_s
integer, intent(in)                     :: i_Nx, i_Ny

call this%initialize_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ny)

end function constructor_val

!*******************************************************************************
subroutine initialize_val(this, i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ny)
!*******************************************************************************
use param, only : coord, comm, ierr
implicit none
class(wake_model_base), intent(inout)   :: this
real(rprec), intent(in)                 :: i_U_infty, i_Delta, i_Dia
real(rprec), dimension(:), intent(in)   :: i_k
real(rprec), dimension(:,:), intent(in) :: i_s
integer, intent(in)                     :: i_Nx, i_Ny
integer                                 :: i, j

! Allocate based on number of turbines
this%N = size(i_s(:,1))
if ( size(i_k) /= this%N ) then
    call error('wake_model_base.initialize','s and k must be the same size')
end if
allocate( this%s(this%N, 3)   )
allocate( this%k(this%N)      )

! Allocate based on number of gridpoints
this%Nx = i_Nx
this%Ny = i_Ny
allocate( this%x(this%Nx) )
allocate( this%y(this%Ny) )
allocate( this%z(this%dNz) )
allocate( this%yu_inf(this%Ny) )
allocate( this%wake_num(this%N) )
allocate( this%ymin(this%N, this%Nx) )
allocate( this%ymax(this%N, this%Nx) )
allocate( this%G(this%N,  this%Nx) )
allocate( this%Gstart(this%N), this%Gstop(this%N) )
allocate( this%d(this%N,  this%Nx) )
allocate( this%dp(this%N, this%Nx) )
allocate( this%w(this%N,  this%Nx) )
allocate( this%fp(this%N, this%Nx) )
allocate( this%disk(this%N, this%dNy, this%dNz) )
allocate( this%diskw(this%N, this%N, this%dNy, this%dNz) )
allocate( this%Cx(this%N, this%Nx) )
allocate( this%px(this%N, this%Nx) )
!allocate( this%Wf(this%N,  this%Nx, this%Ny, this%dNz) )

! Assign input arguments
this%s          = i_s
this%U_infty    = i_U_infty
this%Delta      = i_Delta
this%k          = i_k
this%Dia        = i_Dia

this%sigma_o   = this%Dia / 4

! Normalization constants
this%VELOCITY = i_U_infty
this%LENGTH   = i_Dia
this%TIME     = this%LENGTH / this%VELOCITY
this%FORCE    = this%VELOCITY / this%TIME

! Calculate other variables
this%dx = ( minval(this%s(:,1)) + maxval(this%s(:,1)) ) / this%Nx
this%x(1) = this%dx/2
this%dy = ( minval(this%s(:,2)) + maxval(this%s(:,2)) ) / this%Ny
this%y(1) = this%dy/2
this%dz = this%Dia / this%dNz
this%z(1) = this%dz + (this%s(1,3)-(this%Dia/2))

do i = 2, this%Nx
    this%x(i) = this%x(i-1) + this%dx
end do
do i = 2, this%Ny
    this%y(i) = this%y(i-1) + this%dy
end do
do i = 2, this%dNz
    this%z(i) = this%z(i-1) + this%dz
end do

call this%computeGaussians
call this%computeWakeExpansionFunctions
call this%compute2Dwakes
call this%compute_wake_disks
call this%woketurbines

end subroutine initialize_val

!*******************************************************************************
subroutine computeGaussians(this)
!*******************************************************************************
implicit none
class(wake_model_base), intent(inout) :: this
integer :: i
integer, dimension(1) :: temp_int

! Gaussian functions
this%G = 0._rprec
do i = 1, this%N
    temp_int = minloc(abs(this%x - this%s(i,1) + 4*this%Delta))
    this%Gstart(i) = temp_int(1)
    temp_int = minloc(abs(this%x -  this%s(i,1) - 4*this%Delta))
    this%Gstop(i) = temp_int(1)
    this%G(i,this%Gstart(i):this%Gstop(i)) =                                  &
        gaussian(this%x(this%Gstart(i):this%Gstop(i)), this%s(i,1), this%Delta)
    this%G(i,:) = this%G(i,:) / sum(this%G(i,:)) / this%dx
end do

end subroutine computeGaussians

!*******************************************************************************
subroutine computeWakeExpansionFunctions(this)
!*******************************************************************************
implicit none
class(wake_model_base), intent(inout) :: this
integer                               :: i, j, k, m

do i = 1, this%N
    this%d(i,:)  = 1.0 + this%k(i) * softplus(2.0 *                            &
        (this%s(i,1))/this%Dia, 2.0*this%x/this%Dia)
    this%dp(i,:) = 2.0 * this%k(i) * logistic(2.0 *                            &
        (this%s(i,1))/this%Dia, 2.0*this%x/this%Dia)/this%Dia
!    this%px(i,:) = 2 * (1 + (this%Dia / (max((this%x-this%s(i,1)),0.002))))
!    this%Cx(i,:) = (this%px(i,:)/(2 * gamma(2/this%px(i,:))))                  &
!         *(((this%Dia**2)/(8*this%sigma_o**2))**(2/this%px(i,:)))
!    do m = 1, this%Nx
!        do j = 1, this%Ny
!            do k = 1, this%dNz
!                rn = sqrt((this%y(j)-this%s(i,2))**2 + (this%z(k)-this%s(i,3))**2)
!                this%Wf(i,m,k,j)  = this%Cx(i,m) * exp(-(this%Dia**2)/(8*this%sigma_o**2)*     &
!                    (((2*rn)/(this%Dia * this%d(m,:)) )**this%px(i,m)))
!            enddo
!        enddo
!    enddo
    this%w(i,:)  = 2.0 * this%U_infty * this%dp(i,:) / this%d(i,:)
    this%fp(i,:) =  (1/2)*this%U_infty * this%G(i,:)
end do

end subroutine computeWakeExpansionFunctions

!*******************************************************************************
subroutine compute2Dwakes(this)
!*******************************************************************************
implicit none
class(wake_model_base), intent(inout)  :: this
integer                                :: i, j, dplus
real(rprec)                            :: y_val
real(rprec), dimension(:), allocatable :: dval,dmid
allocate( dval(this%Nx) )
allocate( dmid(this%Nx) )

this%ymin = 1
this%ymax = 1

do i = 1, this%N
    dmid = this%d(i,:) * this%Dia
    dval = this%s(i,2) - (dmid / 2)
    do j = 1, this%Nx
        dplus = 1
        y_val = this%y(dplus)
        do while (y_val < dval(j) .and. dplus < this%Ny)
           dplus = dplus + 1
           this%ymin(i,j) = dplus
           y_val = this%y(dplus)
        end do
    end do

    dval = this%s(i,2) + (dmid / 2)
    do j = 1, this%Nx
        dplus = 1
        y_val = this%y(dplus)
        do while (y_val < dval(j) .and. dplus < this%Ny)
           dplus = dplus + 1
           this%ymax(i,j) = dplus
           y_val = this%y(dplus)
        end do
    end do
!    this%wakeb(i, this%s_int(i,1), this%s_int(i,2)-d_points:this%s_int(i,2)+d_points) = 1
!    this%bounds(i, this%s_int(i,1),1) = this%s_int(i,2)-d_points
!    this%bounds(i, this%s_int(i,1),2) = this%s_int(i,2)+d_points
!    dval = floor((anint(this%d(i,:) * this%Dia - this%Dia)) / (2*this%dy) )
!    do j = this%s_int(i,1)+1, this%Nx
!        d_plus = d_points + dval(j)
!        this%wakeb(i, j, this%s_int(i,2)-d_plus:this%s_int(i,2)+d_plus) = 1
!        this%bounds(i, j, 1) = this%s_int(i,2)-d_plus
end do

end subroutine compute2Dwakes

!*******************************************************************************
subroutine compute_wake_disks(this)
!*******************************************************************************
    use functions, only: linear_interp
implicit none
class(wake_model_base), intent(inout)  :: this
integer                                :: i, j, k, m
real(rprec)                            :: cell_size, hcs, c_z, c_y, center, dist
real(rprec)                            :: ci, pix, rn, dc

cell_size = (this%Dia) / this%dNy
hcs = cell_size/2
c_z = hcs
c_y = hcs
center = (this%Dia) / 2

do i = 1, this%dNz
   c_z = (i-1)*cell_size + hcs
   do j = 1, this%dNy
       c_y = (j-1) * cell_size + hcs
       dist = sqrt((c_y-center)**2 + (c_z-center)**2)
       if (dist < (this%Dia/2)) then
          this%disk(1,j,i) = 1
       else
          this%disk(1,j,i) = 0
       endif
   enddo 
   c_y = hcs
enddo

c_y = 0
c_z = 0

do i =1, this%N
   this%disk(i,:,:) = this%disk(1,:,:)
!   pi = 2 * (1 + (this%Dia / 0.002))
!   ci = (pi/(2 * gamma(2/pi)))*(((this%Dia**2)/(8*this%sigma_o**2))**(2/pi))
  
   do m = 1, this%N
       pix = 2 * (1 + (this%Dia / max((this%s(m,1)-this%s(i,1)),0.02)))
       ci = (pix/(2 * gamma(2/pix)))*(((this%Dia**2)/(8*this%sigma_o**2))**(2/pix))
       c_y = this%s(m,2) - center + hcs
       dc = linear_interp(this%x,this%d(i,:),this%s(m,1))
       do j = 1, this%dNy
           c_z = this%s(m,3) - center + hcs
           do k = 1, this%dNz
               rn = sqrt((c_y-this%s(i,2))**2 + (c_z-this%s(i,3))**2)
               this%diskw(i,m,j,k) = ci * exp(-(this%Dia**2)/(8*this%sigma_o**2)*     &
                    ((((2*rn)/(this%Dia * dc)) )**pix))
               c_z = c_z + cell_size
          enddo
          c_y = c_y + cell_size
       enddo
   enddo
end do

end subroutine compute_wake_disks
!*********************************************************************************

! Find turbines not waked by other turbines
subroutine woketurbines(this)
    implicit none
    class(wake_model_base), intent(inout)  :: this
    integer                                :: i, j, loc(1), counter
    real(rprec)                            :: upper, lower, yup_i(1), yup, &
                                               ydown_i(1), ydown
!    integer, dimension(:), allocatable     :: wake_num
    this%wake_num = 0

    this%nfree = 0
    do i = 1,this%N
       counter = 0
       loc = minloc(abs(this%x - this%s(i,1)))
       do j = 1, this%N
!          print *, this%s(i,1), this%s(i,2)
          if (this%s(i,1) > this%s(j,1)) then
            upper = this%s(i,2) + this%Dia / 2
            lower = this%s(i,2) - this%Dia / 2
            yup_i = this%y(this%ymax(j,loc)); yup = yup_i(1)
            ydown_i = this%y(this%ymin(j,loc)); ydown =  ydown_i(1)
!            print *, lower, upper, ydown, yup
            if (ydown < upper .and. upper < yup) then
                 counter = counter + 1
            else if (ydown < lower .and. lower < yup) then
               counter = counter + 1
            else
            end if

          else
          end if
       end do
    if (counter == 0) then
       this%nfree = this%nfree + 1
!       print *, counter, this%nfree
       this%wake_num(this%nfree) = i
    end if
    end do

if ( allocated(this%free_turbines) ) deallocate(this%free_turbines)
allocate( this%free_turbines(this%nfree ) )
j = 1
do i = 1, this%N
    if (this%wake_num(i) > 0) then
        this%free_turbines(j) = this%wake_num(i)
        j = j + 1
    end if
end do
!write (*,*) 'free_turbines: ', this%free_turbines
!write(*,*) 'Number of free: ', this%nfree 

end subroutine woketurbines

!*******************************************************************************
subroutine makeDimensionless(this)
!*******************************************************************************
implicit none
class(wake_model_base), intent(inout)  :: this

if (.not.this%isDimensionless) then
    this%isDimensionless = .true.
    this%s       = this%s / this%LENGTH
    this%U_infty = this%U_infty / this%VELOCITY
    this%Delta   = this%Delta / this%LENGTH
    this%Dia     = this%Dia / this%LENGTH
    this%dx      = this%dx / this%LENGTH
    this%dy      = this%dy / this%LENGTH
    this%x       = this%x / this%LENGTH
    this%y       = this%y / this%LENGTH
    this%yu_inf  = this%yu_inf / this%VELOCITY
    this%G       = this%G * this%LENGTH         ! G has units 1/length
    this%dp      = this%dp * this%LENGTH        ! dp has units 1/length
    this%w       = this%w * this%TIME           ! w has units 1/time
    this%fp      = this%fp / this%FORCE
end if
end subroutine makeDimensionless

!*******************************************************************************
subroutine makeDimensional(this)
!*******************************************************************************
implicit none
class(wake_model_base), intent(inout)  :: this

if (this%isDimensionless) then
    this%isDimensionless = .false.
    this%s       = this%s * this%LENGTH
    this%U_infty = this%U_infty * this%VELOCITY
    this%Delta   = this%Delta * this%LENGTH
    this%Dia     = this%Dia * this%LENGTH
    this%dx      = this%dx * this%LENGTH
    this%dy      = this%dy * this%LENGTH
    this%x       = this%x * this%LENGTH
    this%y       = this%y * this%LENGTH
    this%yu_inf  = this%yu_inf * this%VELOCITY
    this%G       = this%G / this%LENGTH         ! G has units 1/length
    this%dp      = this%dp / this%LENGTH        ! dp has units 1/length
    this%w       = this%w / this%TIME           ! w has units 1/time
    this%fp      = this%fp * this%FORCE
end if
end subroutine makeDimensional

! Writes object to file
! subroutine write_to_file(this, fstring)
!     use open_file_fid_mod
!     use param, only : CHAR_BUFF_LENGTH
!     implicit none
!     class(wake_model_base), intent(in) :: this
!     character(CHAR_BUFF_LENGTH), intent(in)  :: fstring
!     integer :: i, fid
!
!     !  Open vel.out (lun_default in io) for final output
!     fid = open_file_fid(fstring, 'rewind', 'unformatted')
!     write(fid) this%N, this%Nx, this%dx, this%Dia, this%Delta,               &
!                this%U_infty, this%isDimensionless
!     write(fid) this%LENGTH, this%VELOCITY, this%TIME, this%FORCE
!     write(fid) this%s
!     write(fid) this%k
!     write(fid) this%x
!     write(fid) this%du
!     write(fid) this%u
!     write(fid) this%uhat
!     write(fid) this%Phat
!     write(fid) this%Ctp
!     close(fid)
!
! end subroutine write_to_file


!*******************************************************************************
subroutine print(this)
!*******************************************************************************
! Prints all variables of the class to standard output
implicit none
class(wake_model_base), intent(in) :: this
integer :: i

write(*,*) ' U_infty         = ', this%U_infty
write(*,*) ' Delta           = ', this%Delta
write(*,*) ' Dia             = ', this%Dia
write(*,*) ' Nx              = ', this%Nx
write(*,*) ' x               = ', this%x
write(*,*) ' Ny              = ', this%Ny
write(*,*) ' dNy             = ', this%dNy
write(*,*) ' dNz             = ', this%dNz
write(*,*) ' y               = ', this%y
write(*,*) ' z               = ', this%z
write(*,*) ' yu_inf          = ', this%yu_inf
write(*,*) ' dx              = ', this%dx
write(*,*) ' dy              = ', this%dy
write(*,*) ' isDimensionless = ', this%isDimensionless
write(*,*) ' LENGTH          = ', this%LENGTH
write(*,*) ' VELOCITY        = ', this%VELOCITY
write(*,*) ' TIME            = ', this%TIME
write(*,*) ' FORCE           = ', this%FORCE
do i = 1, this%N
    write(*,*) ' Wake', i,':'
    write(*,*) '  s       = ', this%s(i,:)
    write(*,*) '  k       = ', this%k(i)
    write(*,*) '  G       = ', this%G(i,:)
    write(*,*) '  d       = ', this%d(i,:)
    write(*,*) '  dp      = ', this%dp(i,:)
    write(*,*) '  w       = ', this%w(i,:)
end do
end subroutine print

end module wake_model_base_m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   Wake Model Class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module wake_model_class
use types, only : rprec
use util,  only : logistic, softplus, gaussian
use param, only : pi, coord, nproc
use wake_model_base_m
use messages
implicit none

private
public WakeModel

type, extends(wake_model_base) :: WakeModel
    real(rprec), dimension(:,:), allocatable   :: du   ! velocity deficit (turbine, space)
    real(rprec), dimension(:,:), allocatable   :: u    ! superimposed velocity (space)
    real(rprec), dimension(:), allocatable     :: uhat ! estimated local turbine velocity (turbine)
    real(rprec), dimension(:), allocatable     :: Phat ! estimated turbine power (turbine)
    real(rprec), dimension(:), allocatable     :: Ctp  ! local thrust coefficient (turbine)
contains
    procedure, public  :: initialize_val
    procedure, private :: initialize_file
    procedure, public  :: print
    procedure, public  :: write_to_file
    procedure, public  :: makeDimensionless
    procedure, public  :: makeDimensional
    procedure, public  :: advance
    procedure, private :: rhs

end type WakeModel

interface WakeModel
    module procedure :: constructor_val
    module procedure :: constructor_file
end interface WakeModel

contains

! Constructor for wake model with values given
function constructor_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ny) result(this)
    implicit none
    type(WakeModel)                         :: this
    real(rprec), intent(in)                 :: i_U_infty, i_Delta, i_Dia
    real(rprec), dimension(:), intent(in)   :: i_k
    real(rprec), dimension(:,:), intent(in) :: i_s
    integer, intent(in)                     :: i_Nx, i_Ny

    call this%initialize_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ny)
end function constructor_val

! Constructor for wake model that reads from file
function constructor_file(fstring) result(this)
    use param, only : CHAR_BUFF_LENGTH
    implicit none

    type(WakeModel)          :: this
    character(*), intent(in) :: fstring

    call this%initialize_file(fstring)

end function constructor_file

subroutine initialize_val(this, i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ny)
    use wake_model_base_m
    implicit none
    class(WakeModel), intent(inout)         :: this
    real(rprec), intent(in)                 :: i_U_infty, i_Delta, i_Dia
    real(rprec), dimension(:), intent(in)   :: i_k
    real(rprec), dimension(:,:), intent(in) :: i_s
    integer, intent(in)                     :: i_Nx, i_Ny

    ! Call base class initializer
    call this%wake_model_base%initialize_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ny)

    allocate( this%du(this%N, this%Nx) )
    allocate( this%u(this%Nx, this%Ny) )
    allocate( this%uhat(this%N) )
    allocate( this%Phat(this%N) )
    allocate( this%Ctp(this%N)  )

    this%du(:,:) = 0.d0
    this%u(:,:)    = this%U_infty
    this%uhat(:) = this%U_infty
    this%Phat(:) = 0.d0
    this%Ctp(:)  = 0.d0

end subroutine initialize_val

subroutine initialize_file(this, fstring)
    use param, only : CHAR_BUFF_LENGTH
    implicit none

    class(WakeModel), intent(inout) :: this
    character(*), intent(in)        :: fstring
    integer :: fid

    !  Open vel.out (lun_default in io) for final output
    open(newunit=fid , file=fstring, status='unknown', form='unformatted',     &
        action='read', position='rewind')
    read(fid) this%N, this%Nx, this%Ny, this%dNy, this%dNz, this%dx, this%dy, this%Dia, this%Delta,&
              this%U_infty, this%sigma_o, this%isDimensionless
    read(fid) this%LENGTH, this%VELOCITY, this%TIME, this%FORCE

    allocate( this%s(this%N, 3) )
    allocate( this%k(this%N)    )
    allocate( this%uhat(this%N) )
    allocate( this%wake_num(this%N) )
    allocate( this%Phat(this%N) )
    allocate( this%Ctp(this%N)  )
    allocate( this%x(this%Nx) )
    allocate( this%y(this%Ny) )
    allocate( this%z(this%dNz) )
    allocate( this%yu_inf(this%Ny) )
    allocate( this%ymin(this%N, this%Nx) )
    allocate( this%ymax(this%N, this%Nx) )
    allocate( this%u(this%Nx, this%Ny) )
    allocate( this%G(this%N,  this%Nx) )
    allocate( this%Gstart(this%N), this%Gstop(this%N) )
    allocate( this%d(this%N,  this%Nx) )
    allocate( this%dp(this%N, this%Nx) )
    allocate( this%w(this%N,  this%Nx) )
    allocate( this%fp(this%N, this%Nx) )
    allocate( this%du(this%N, this%Nx) )
    allocate( this%disk(this%N, this%dNy, this%dNz) )
    allocate( this%diskw(this%N, this%N, this%dNy, this%dNz) )
!    allocate( this%Wf(this%N, this%Nx, this%Ny, this%dNz) )

    read(fid) this%s
    read(fid) this%k
    read(fid) this%x
    read(fid) this%y
    read(fid) this%z
    read(fid) this%yu_inf
    read(fid) this%du
    read(fid) this%u
    read(fid) this%disk
    read(fid) this%uhat
    read(fid) this%Phat
    read(fid) this%Ctp
    close(fid)

    call this%computeGaussians
    call this%computeWakeExpansionFunctions
    call this%compute2Dwakes
    call this%compute_wake_disks
    call this%woketurbines

end subroutine initialize_file

subroutine makeDimensionless(this)
    implicit none
    class(WakeModel), intent(inout)  :: this

    if (.not.this%isDimensionless) then
        call this%wake_model_base%makeDimensionless
        this%du      = this%du / this%VELOCITY
        this%u       = this%u / this%VELOCITY
        this%uhat    = this%uhat / this%VELOCITY
        this%Phat    = this%uhat / this%VELOCITY**3
    end if
end subroutine makeDimensionless

subroutine makeDimensional(this)
    implicit none
    class(WakeModel), intent(inout)  :: this

    if (this%isDimensionless) then
        call this%wake_model_base%makeDimensional
        this%du      = this%du * this%VELOCITY
        this%u       = this%u * this%VELOCITY
        this%uhat    = this%uhat * this%VELOCITY
        this%Phat    = this%uhat * this%VELOCITY**3
    end if
end subroutine makeDimensional

! Writes object to file
subroutine write_to_file(this, fstring)
    use param, only : CHAR_BUFF_LENGTH
    implicit none
    class(WakeModel), intent(in) :: this
    character(CHAR_BUFF_LENGTH), intent(in)  :: fstring
    integer :: fid

    !  Open vel.out (lun_default in io) for final output
    open(newunit=fid , file=fstring, status='unknown', form='unformatted',     &
        action='write', position='rewind')
    write(fid) this%N, this%Nx, this%Ny, this%dNy, this%dNz, this%dx, this%dy, &
               this%Dia, this%Delta, this%U_infty, this%sigma_o, this%isDimensionless
    write(fid) this%LENGTH, this%VELOCITY, this%TIME, this%FORCE
    write(fid) this%s
    write(fid) this%k
    write(fid) this%x
    write(fid) this%y
    write(fid) this%z
    write(fid) this%yu_inf
    write(fid) this%du
    write(fid) this%u
    write(fid) this%disk
    write(fid) this%uhat
    write(fid) this%Phat
    write(fid) this%Ctp
    close(fid)

end subroutine write_to_file

! Prints all variables of the class to standard output
subroutine print(this)
    implicit none
    class(WakeModel), intent(in) :: this
    integer :: i

    write(*,*) ' U_infty         = ', this%U_infty
    write(*,*) ' Delta           = ', this%Delta
    write(*,*) ' Dia             = ', this%Dia
    write(*,*) ' Nx              = ', this%Nx
    write(*,*) ' Ny              = ', this%Ny
    write(*,*) ' dNy             = ', this%dNy
    write(*,*) ' dNz             = ', this%dNz
    write(*,*) ' x               = ', this%x
    write(*,*) ' y               = ', this%y
    write(*,*) ' z               = ', this%z
    write(*,*) ' yu_inf          = ', this%yu_inf
    write(*,*) ' dx              = ', this%dx
    write(*,*) ' dy              = ', this%dy
    write(*,*) ' isDimensionless = ', this%isDimensionless
    write(*,*) ' LENGTH          = ', this%LENGTH
    write(*,*) ' VELOCITY        = ', this%VELOCITY
    write(*,*) ' TIME            = ', this%TIME
    write(*,*) ' FORCE           = ', this%FORCE
    write(*,*) ' u               = ', this%u
    do i = 1, this%N
        write(*,*) ' Wake', i,':'
        write(*,*) '  uhat    = ', this%uhat(i)
        write(*,*) '  s       = ', this%s(i,:)
        write(*,*) '  k       = ', this%k(i)
        write(*,*) '  G       = ', this%G(i,:)
        write(*,*) '  d       = ', this%d(i,:)
        write(*,*) '  dp      = ', this%dp(i,:)
        write(*,*) '  w       = ', this%w(i,:)
    end do
end subroutine print

subroutine advance(this, Ctp, dt)
    use util, only : ddx_upwind1
    use functions, only: linear_interp
    implicit none
    class(WakeModel), intent(inout)          :: this
    real(rprec), intent(in)                  :: dt
    real(rprec), dimension(:), intent(in)    :: Ctp
    integer                                  :: i, j
    real(rprec)                              :: diff, a1, a2, du_loop
    real(rprec), dimension(:,:,:), allocatable :: new_ui
    real(rprec), dimension(:,:,:), allocatable :: du_sum, u_field

    allocate( new_ui(this%N,this%dNy,this%dNz) )

    if (size(Ctp) /= this%N) then
        call error('WakeModel.advance','Ctp must be size N')
    end if

    ! Compute new wake deficit and superimpose wakes
    do i = 1, this%N
        this%du(i,:) = this%du(i,:) +  dt * this%rhs(this%du(i,:),             &
            this%fp(i,:) * this%Ctp(i) * this%uhat(i), i)
!        if (i == 1 .and. coord == 0) then
!            write(*,*) 'dt: ', dt
!            write(*,*) 'rhs: ', this%rhs(this%du(i,:),             &
!               this%fp(i,:) * this%Ctp(i) / (4.0 + this%Ctp(i)), i)
!            write(*,*) 'inputs (forcing) = ', this%fp(i,:) * this%Ctp(i) / (4.0 + this%Ctp(i))
!        endif
    end do
    allocate( du_sum( this%N, this%dNy, this%dNz) )
    allocate( u_field( this%N, this%dNy, this%dNz) )
    du_sum = 0._rprec
    do i = 1, this%N
       do j = 1, this%N
           du_loop = linear_interp(this%x,this%du(j,:),this%s(i,1))
           du_sum(i,:,:) = du_sum(i,:,:) + du_loop * this%diskw(j,i,:,:)
!            du_superimposed(j,this%ymin(i,j):this%ymax(i,j)) =                 &
!                du_superimposed(j,this%ymin(i,j):this%ymax(i,j))               &
!                +this%du(i,j)**2
        end do
        new_ui(i,:,:) = linear_interp(this%y,this%yu_inf,this%s(i,2))
    end do

    ! Find the velocity field
    u_field = new_ui - du_sum

    ! Find estimated velocities
    this%uhat(:) = 0
    do i = 1, this%N
        a1 = sum(sum(this%disk(i,:,:),dim=1))
        a2 = sum(sum(u_field(i,:,:)*this%disk(i,:,:),dim=1))
        this%uhat(i) = a2 / a1
!        do j = this%ymin(i,1), this%ymax(i,1)
!            this%uhat(i) = this%uhat(i)                                        &
!                + sum(this%G(i,this%Gstart(i):this%Gstop(i))                   &
!                * this%u(this%Gstart(i):this%Gstop(i),j) * this%dx)
!        end do
!        diff = this%ymax(i,1) - this%ymin(i,1)
!        this%uhat(i) = this%uhat(i) / (diff + 1)
        this%Ctp(i) = Ctp(i)
        this%Phat(i) = this%Ctp(i) * this%uhat(i)**3
    end do

end subroutine advance

! Evaluates RHS of wake equation
function rhs(this, du, f, i) result(ddudt)
    use util, only : ddx_upwind1
    implicit none
    class(WakeModel), intent(in)           :: this
    real(rprec), dimension(:), intent(in)  :: f, du
    integer, intent(in)                    :: i
    real(rprec), dimension(:), allocatable :: ddudt, ddudx

    allocate(ddudt(this%Nx))
    allocate(ddudx(this%Nx))

    ddudx = ddx_upwind1(du, this%dx)
    ddudt = -this%U_infty * ddudx - this%w(i,:) * du + f
end function rhs

end module wake_model_class

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!                   Wake Model Adjoint Class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!module wake_model_adjoint_class
!use types, only : rprec
!use util,  only : logistic, softplus, gaussian
!use wake_model_base_m
!use messages
!implicit none
!
!private
!public WakeModelAdjoint
!
!type, extends(wake_model_base) :: WakeModelAdjoint
!   real(rprec), dimension(:,:), allocatable :: dustar ! adjoint velocity deficit (turbine, space)
!contains
!   procedure, public  :: initialize_val
!!     procedure, private :: initialize_file
!!     procedure, public  :: print
!!     procedure, public  :: write_to_file
!   procedure, public  :: makeDimensionless
!   procedure, public  :: makeDimensional
!   procedure, public  :: retract
!   procedure, private :: rhs
!
!end type WakeModelAdjoint
!
!interface WakeModelAdjoint
!   module procedure :: constructor_val
!!     module procedure :: constructor_file
!end interface WakeModelAdjoint
!
!contains
!
!! Constructor for wake model with values given
!function constructor_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ny) result(this)
!   implicit none
!   type(WakeModelAdjoint)                  :: this
!   real(rprec), intent(in)                 :: i_U_infty, i_Delta, i_Dia
!   real(rprec), dimension(:), intent(in)   :: i_k
!   real(rprec), dimension(:,:), intent(in) :: i_s
!   integer, intent(in)                     :: i_Nx, i_Ny
!
!   call this%initialize_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ny)
!end function constructor_val
!
!! Constructor for wake model that reads from file
!! function constructor_file(fstring) result(this)
!!     use open_file_fid_mod
!!     use param, only : CHAR_BUFF_LENGTH
!!     implicit none
!!
!!     type(WakeModel)          :: this
!!     character(*), intent(in) :: fstring
!!
!!     call this%initialize_file(fstring)
!!
!! end function constructor_file
!
!subroutine initialize_val(this, i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ny)
!   use wake_model_base_m
!   implicit none
!   class(WakeModelAdjoint), intent(inout)  :: this
!   real(rprec), intent(in)                 :: i_U_infty, i_Delta, i_Dia
!   real(rprec), dimension(:), intent(in)   :: i_k
!   real(rprec), dimension(:,:), intent(in) :: i_s
!   integer, intent(in)                     :: i_Nx, i_Ny
!
!   ! Call base class initializer
!   call this%wake_model_base%initialize_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ny)
!
!   allocate( this%dustar(this%N, this%Nx) )
!
!   this%dustar(:,:) = 0.d0
!
!end subroutine initialize_val
!!
!! subroutine initialize_file(this, fstring)
!!     use open_file_fid_mod
!!     use param, only : CHAR_BUFF_LENGTH
!!     implicit none
!!
!!     class(WakeModel), intent(inout) :: this
!!     character(*), intent(in)        :: fstring
!!     integer :: i, fid
!!
!!     !  Open vel.out (lun_default in io) for final output
!!     fid = open_file_fid(fstring, 'rewind', 'unformatted')
!!     read(fid) this%N, this%Nx, this%dx, this%Dia, this%Delta,                           &
!!               this%U_infty, this%isDimensionless
!!     read(fid) this%LENGTH, this%VELOCITY, this%TIME, this%FORCE
!!
!!     allocate( this%s(this%N)    )
!!     allocate( this%k(this%N)    )
!!     allocate( this%uhat(this%N) )
!!     allocate( this%Phat(this%N) )
!!     allocate( this%Ctp(this%N)  )
!!     allocate( this%x(this%Nx) )
!!     allocate( this%u(this%Nx) )
!!     allocate( this%G(this%N,  this%Nx) )
!!     allocate( this%d(this%N,  this%Nx) )
!!     allocate( this%dp(this%N, this%Nx) )
!!     allocate( this%w(this%N,  this%Nx) )
!!     allocate( this%fp(this%N, this%Nx) )
!!     allocate( this%du(this%N, this%Nx) )
!!
!!     read(fid) this%s
!!     read(fid) this%k
!!     read(fid) this%x
!!     read(fid) this%du
!!     read(fid) this%u
!!     read(fid) this%uhat
!!     read(fid) this%Phat
!!     read(fid) this%Ctp
!!     close(fid)
!!
!!     do i = 1, this%N
!!         this%G(i,:) = gaussian(this%x, this%s(i), this%Delta)
!!     end do
!!     call this%computeWakeExpansionFunctions
!!
!! end subroutine initialize_file
!
!subroutine makeDimensionless(this)
!   implicit none
!   class(WakeModelAdjoint), intent(inout)  :: this
!
!   if (.not.this%isDimensionless) then
!       call this%wake_model_base%makeDimensionless
!       this%dustar = this%dustar / this%VELOCITY
!   end if
!end subroutine makeDimensionless
!
!subroutine makeDimensional(this)
!   implicit none
!   class(WakeModelAdjoint), intent(inout)  :: this
!
!   if (this%isDimensionless) then
!       call this%wake_model_base%makeDimensional
!       this%dustar = this%dustar * this%VELOCITY
!   end if
!end subroutine makeDimensional
!
!! Writes object to file
!! subroutine write_to_file(this, fstring)
!!     use open_file_fid_mod
!!     use param, only : CHAR_BUFF_LENGTH
!!     implicit none
!!     class(WakeModel), intent(in) :: this
!!     character(CHAR_BUFF_LENGTH), intent(in)  :: fstring
!!     integer :: i, fid
!!
!!     !  Open vel.out (lun_default in io) for final output
!!     fid = open_file_fid(fstring, 'rewind', 'unformatted')
!!     write(fid) this%N, this%Nx, this%dx, this%Dia, this%Delta,                           &
!!                this%U_infty, this%isDimensionless
!!     write(fid) this%LENGTH, this%VELOCITY, this%TIME, this%FORCE
!!     write(fid) this%s
!!     write(fid) this%k
!!     write(fid) this%x
!!     write(fid) this%du
!!     write(fid) this%u
!!     write(fid) this%uhat
!!     write(fid) this%Phat
!!     write(fid) this%Ctp
!!     close(fid)
!!
!! end subroutine write_to_file
!
!! Prints all variables of the class to standard output
!! subroutine print(this)
!!     implicit none
!!     class(WakeModel), intent(in) :: this
!!     integer :: i
!!
!!     write(*,*) ' U_infty         = ', this%U_infty
!!     write(*,*) ' Delta           = ', this%Delta
!!     write(*,*) ' Dia             = ', this%Dia
!!     write(*,*) ' Nx              = ', this%Nx
!!     write(*,*) ' x               = ', this%x
!!     write(*,*) ' dx              = ', this%dx
!!     write(*,*) ' isDimensionless = ', this%isDimensionless
!!     write(*,*) ' LENGTH          = ', this%LENGTH
!!     write(*,*) ' VELOCITY        = ', this%VELOCITY
!!     write(*,*) ' TIME            = ', this%TIME
!!     write(*,*) ' FORCE           = ', this%FORCE
!!     write(*,*) ' u               = ', this%u
!!     do i = 1, this%N
!!         write(*,*) ' Wake', i,':'
!!         write(*,*) '  uhat    = ', this%uhat(i)
!!         write(*,*) '  s       = ', this%s(i)
!!         write(*,*) '  k       = ', this%k(i)
!!         write(*,*) '  G       = ', this%G(i,:)
!!         write(*,*) '  d       = ', this%d(i,:)
!!         write(*,*) '  dp      = ', this%dp(i,:)
!!         write(*,*) '  w       = ', this%w(i,:)
!!     end do
!! end subroutine print
!
!subroutine retract(this, fstar, dt, g)
!   implicit none
!   class(WakeModelAdjoint), intent(inout)   :: this
!   real(rprec), dimension(:,:), intent(in)  :: fstar
!   real(rprec), intent(in)                  :: dt
!   real(rprec), dimension(:), intent(inout) :: g
!   !real(rprec), dimension(:), allocatable   :: k1star, k2star, k3star, k4star
!   !real(rprec), dimension(:), allocatable   :: k1star_rhs, k2star_rhs, k3star_rhs, k4star_rhs
!   !real(rprec)                              :: gk1s, gk2s, gk3s, gk4s
!   integer                                  :: i
!
!   do i = 1, this%N
!       this%dustar(i,:) = this%dustar(i,:) - this%rhs(this%dustar(i,:), fstar(i,:), i) * dt
!       g(i) = sum(this%G(i,:) * this%dustar(i,:) / this%d(i,:) / this%d(i,:)) * this%dx
!   end do
!
!end subroutine retract
!
!! Evaluates RHS of wake equation using 3rd-order biased downwind differencing
!function rhs(this, du, f, i) result(ddudt)
!   use util, only : ddx_downwind1
!   implicit none
!   class(WakeModelAdjoint), intent(in)    :: this
!   real(rprec), dimension(:), intent(in)  :: f, du
!   real(rprec), dimension(:), allocatable :: ddudt, ddudx
!   integer, intent(in)                    :: i
!
!   allocate(ddudt(this%Nx))
!   allocate(ddudx(this%Nx))
!
!   ddudx = ddx_downwind1(du, this%dx)
!   ddudt = -this%U_infty * ddudx + this%w(i,:) * du - f
!end function rhs
!
!end module wake_model_adjoint_class

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                   Wake Model Estimator Class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module wake_model_estimator_class
use types, only : rprec
use messages
use wake_model_class
use param, only : pi, coord, nproc
use functions, only: linear_interp
#ifdef PPMPI
use param, only : MPI_RPREC, status, comm, ierr
use mpi
#endif
implicit none

private
public WakeModelEstimator

type :: WakeModelEstimator
    type(WakeModel), dimension(:), allocatable :: ensemble
    type(WakeModel)                            :: wm
    real(rprec)                                :: sigma_du, sigma_k, sigma_Phat
    integer                                    :: Ne, Ne_tot, Nm, Ns
    real(rprec), dimension(:,:), allocatable   :: A, Aprime, Ahat, Ahatprime, E, D, Dprime
    real(rprec), dimension(:), allocatable     :: Abar, Ahatbar
    real(rprec), dimension(:,:), allocatable   :: MM_array, SM_array, ME_array,&
         SE_array, DE_array
    real(rprec)                                :: tau_U_infty = 300
contains
    procedure, private :: initialize_val
    procedure, private :: initialize_file
    procedure, public  :: write_to_file
    procedure, public  :: generateInitialEnsemble
    procedure, public  :: advance
    procedure, private :: advanceEnsemble
    procedure, private :: sortFreeTurbines
end type WakeModelEstimator

interface WakeModelEstimator
    module procedure :: constructor_val
    module procedure :: constructor_file
end interface WakeModelEstimator

contains

! Constructor for wake model
function constructor_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ny, i_Ne,                &
                         i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau) result(this)
    implicit none
    type(WakeModelEstimator)                :: this
    real(rprec), intent(in)                 :: i_U_infty, i_Delta, i_Dia
    real(rprec), dimension(:), intent(in)   :: i_k
    real(rprec), dimension(:,:), intent(in) :: i_s
    integer, intent(in)                     :: i_Nx, i_Ny
    integer, intent(in)                     :: i_Ne
    real(rprec), intent(in)                 :: i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau


    call this%initialize_val(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ny, i_Ne,            &
                             i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau)


end function constructor_val

! Constructor for wake model that reads from file
function constructor_file(fpath, i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau) result(this)
    use param, only : CHAR_BUFF_LENGTH
    implicit none

    type(WakeModelEstimator) :: this
    character(*), intent(in) :: fpath
    real(rprec), intent(in)  :: i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau

    call this%initialize_file(fpath, i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau)

end function constructor_file

subroutine initialize_val(this, i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ny, i_Ne,         &
                          i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau)
    implicit none
    class(WakeModelEstimator), intent(inout)    :: this
    real(rprec), intent(in)                     :: i_sigma_du, i_U_infty, i_Delta, i_Dia
    real(rprec), dimension(:), intent(in)       :: i_k
    real(rprec), dimension(:,:), intent(in)     :: i_s
    integer, intent(in)                         :: i_Nx, i_Ny
    integer, intent(in)                         :: i_Ne
    real(rprec), intent(in)                     :: i_sigma_k, i_sigma_Phat, i_tau
    integer                                     :: i

    ! Set std deviations for noise
    this%sigma_du   = i_sigma_du
    this%sigma_k    = i_sigma_k
    this%sigma_Phat = i_sigma_Phat

    ! Filter time for U_infty
    this%tau_U_infty = i_tau
    
    ! Create ensemble
    this%Ne = ceiling( 1._rprec * i_Ne / nproc)
    this%Ne_tot = this%Ne * nproc

    if (coord == 0) write(*,*) "Resizing ensemble to have ", this%Ne_tot,      &
        "members"
!    if (coord == 0) write(*,*) "Resizing ensemble to have ", this%Ne_tot,      &
    
    ! Create wake model estim
    this%wm = WakeModel(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx, i_Ny)

    ! Create ensemble members
    allocate( this%ensemble(this%Ne_tot) )
    do i = 1, this%Ne
        this%ensemble(i+coord*this%Ne) = WakeModel(i_s, i_U_infty, i_Delta, i_k, i_Dia, i_Nx,&
             i_Ny)
    end do

    ! Allocate filter matrices
    this%Nm = this%wm%N                              ! Number of measurements
    this%Ns = this%wm%N * this%wm%Nx + this%wm%N     ! Number of states
    allocate( this%Abar(this%Ns) )
    allocate( this%Ahatbar(this%Nm) )
    allocate( this%A(this%Ns, this%Ne_tot) )
    allocate( this%Aprime(this%Ns, this%Ne_tot) )
    allocate( this%Ahat(this%Nm, this%Ne_tot) )
    allocate( this%Ahatprime(this%Nm, this%Ne_tot) )
    allocate( this%E(this%Nm, this%Ne_tot) )
    allocate( this%D(this%Nm, this%Ne_tot) )
    allocate( this%Dprime(this%Nm, this%Ne_tot) )
    allocate( this%MM_array(this%Nm, this%Nm) )
    allocate( this%SM_array(this%Ns, this%Nm) )
    allocate( this%ME_array(this%Nm, this%Ne_tot) )
    allocate( this%SE_array(this%Ns, this%Ne_tot) )
    allocate( this%DE_array(this%Nm, this%Ne_tot) )

end subroutine initialize_val

subroutine initialize_file(this, fpath, i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau)
    use param, only : CHAR_BUFF_LENGTH
    use string_util, only : string_splice
    implicit none

    class(WakeModelEstimator)   :: this
    character(*), intent(in)    :: fpath
    character(CHAR_BUFF_LENGTH) :: fstring
    integer                     :: i, fid
    real(rprec), intent(in)     :: i_sigma_du, i_sigma_k, i_sigma_Phat, i_tau

    ! Set std deviations for noise
    this%sigma_du   = i_sigma_du
    this%sigma_k    = i_sigma_k
    this%sigma_Phat = i_sigma_Phat

    ! Filter time for U_infty
    this%tau_U_infty = i_tau

    fstring = fpath // '/wm_est.dat'
    open(newunit=fid , file=fstring, status='unknown', form='unformatted',     &
        action='read', position='rewind')
    read(fid) this%Ne, this%Ne_tot, this%Nm, this%Ns

    allocate( this%Abar(this%Ns) )
    allocate( this%Ahatbar(this%Nm) )
    allocate( this%A(this%Ns, this%Ne_tot) )
    allocate( this%Aprime(this%Ns, this%Ne_tot) )
    allocate( this%Ahat(this%Nm, this%Ne_tot) )
    allocate( this%Ahatprime(this%Nm, this%Ne_tot) )
    allocate( this%E(this%Nm, this%Ne_tot) )
    allocate( this%D(this%Nm, this%Ne_tot) )
    allocate( this%Dprime(this%Nm, this%Ne_tot) )
    allocate( this%MM_array(this%Nm, this%Nm) )
    allocate( this%SM_array(this%Ns, this%Nm) )
    allocate( this%ME_array(this%Nm, this%Ne_tot) )
    allocate( this%SE_array(this%Ns, this%Ne_tot) )
    allocate( this%DE_array(this%Nm, this%Ne_tot) )

    read(fid) this%A
    read(fid) this%Aprime
    read(fid) this%Ahat
    read(fid) this%Ahatprime
    read(fid) this%E
    read(fid) this%D
    read(fid) this%Dprime
    read(fid) this%Abar
    read(fid) this%Ahatbar
    close(fid)

    fstring = fpath // '/wm.dat'
    this%wm = WakeModel(fstring)

    allocate( this%ensemble(this%Ne_tot) )
    do i = 1, this%Ne
        call string_splice( fstring, fpath // '/ensemble_', i+coord*this%Ne, '.dat' )
        this%ensemble(i+coord*this%Ne) = WakeModel(fstring)
    end do

end subroutine initialize_file

subroutine write_to_file(this, fpath)
    use param, only : CHAR_BUFF_LENGTH
    use string_util, only : string_splice
    implicit none

    class(WakeModelEstimator)   :: this
    character(*), intent(in)    :: fpath
    character(CHAR_BUFF_LENGTH) :: fstring
    integer                     :: i, fid

    call system('mkdir -vp ' // fpath)
    fstring = fpath // '/wm_est.dat'
  if (coord == 0) then
    open(newunit=fid , file=fstring, status='unknown', form='unformatted',     &
        action='write', position='rewind')
    write(fid) this%Ne, this%Ne_tot, this%Nm, this%Ns
    write(fid) this%A
    write(fid) this%Aprime
    write(fid) this%Ahat
    write(fid) this%Ahatprime
    write(fid) this%E
    write(fid) this%D
    write(fid) this%Dprime
    write(fid) this%Abar
    write(fid) this%Ahatbar
    close(fid)

    fstring = fpath // '/wm.dat'
    call this%wm%write_to_file(fstring)
   else
   endif

    do i = 1, this%Ne
        call string_splice( fstring, fpath // '/ensemble_', i+coord*this%Ne, '.dat' )
        call this%ensemble(i+coord*this%Ne)%write_to_file(fstring)
    end do

end subroutine write_to_file

subroutine generateInitialEnsemble(this, Ctp, Pm)
    use util, only : random_normal
    use functions, only: linear_interp
    implicit none
    class(WakeModelEstimator), intent(inout)    :: this
    real(rprec), dimension(:)                   :: Ctp, Pm
    real(rprec), allocatable, dimension(:)      ::  yu_turb, ylocs_sort
    real(rprec)                                 :: dt
    real(rprec), parameter                      :: cfl = 0.5
    real(rprec)                                 :: FTT
    integer                                     :: ii, i, j, N, Nx, jstart, jend


    if (size(Ctp) /= this%ensemble(1+coord*this%Ne)%N) then
        call error('WakeModelEstimator.generateInitialEnsemble','Ctp must be of size N')
    end if

    if (coord == 0) write(*,*) 'Generating initial wake model filter ensemble...'

    ! Initialize random number generator
    call init_random_seed

    ! Calculate safe dt.
    dt = cfl * this%wm%dx / this%wm%U_infty
    FTT = this%wm%x(this%wm%Nx) / this%wm%U_Infty

    ! set integer
    N = this%wm%N
    Nx = this%wm%Nx
    allocate( yu_turb(this%wm%nfree) )
    allocate( ylocs_sort(this%wm%nfree) )

    ! initial u_inf(y) setting
    do j = 1, size(this%wm%free_turbines)
       i = this%wm%free_turbines(j)  
       yu_turb(j) = (Pm(i) / Ctp(i))**(1._rprec/3._rprec)      &
           * (4._rprec + Ctp(i)) / 4._rprec
    end do
    
        call this%sortFreeTurbines(yu_turb, ylocs_sort)

    this%wm%yu_inf = linear_interp(ylocs_sort,yu_turb,this%wm%y)
    this%wm%U_infty = sum(this%wm%yu_inf) / size(this%wm%yu_inf)
    

    ! Do 1 FFT of k's
    do i = 1, this%Ne
        this%ensemble(i+coord*this%Ne)%yu_inf(:) = this%wm%yu_inf
        this%ensemble(i+coord*this%Ne)%U_infty = this%wm%U_infty
        do ii = 1, floor(FTT / dt)
            do j = 1, N
                this%ensemble(i+coord*this%Ne)%k(j) = max(this%ensemble(i+coord*this%Ne)%k(j)                  &
                   + sqrt(dt) * this%sigma_k * random_normal(), 0._rprec)
           end do
        end do
        call this%ensemble(i+coord*this%Ne)%computeWakeExpansionFunctions
        call this%ensemble(i+coord*this%Ne)%compute2Dwakes
        call this%ensemble(i+coord*this%Ne)%wokeTurbines
    end do

    ! Do at least 1 FTT of simulations to get good ensemble statistics
    do ii = 1, floor(FTT / dt)
        ! Advance step for objects
        ! always safeguard against negative k's, du's, and omega's
        do i = 1, this%Ne
            do j = 1, N
                this%ensemble(i+coord*this%Ne)%du(j,:) = max(this%ensemble(i+coord*this%Ne)%du(j,:)            &
                    + sqrt(dt) * this%sigma_du * random_normal()                   &
                    * this%ensemble(i+coord*this%Ne)%G(j,:) * sqrt(2*pi)                         &
                    * this%ensemble(i+coord*this%Ne)%Delta, 0._rprec)
           end do
            call this%ensemble(i+coord*this%Ne)%advance(Ctp, dt)
        end do
        call this%wm%advance(Ctp, dt)
    end do

    ! Place ensemble into a matrix with each member in a column
    this%A = 0._rprec
    this%Ahat = 0._rprec
    this%Abar = 0._rprec
    this%Ahatbar = 0._rprec
    do i = 1, this%Ne
        do j = 1, N
            jstart = (j-1)*Nx+1
            jend = j*Nx
            this%A(jstart:jend,i+coord*this%Ne) = this%ensemble(i+coord*this%Ne)%du(j,:)
        end do
        this%A((N*Nx+1):,i+coord*this%Ne) = this%ensemble(i+coord*this%Ne)%k(:)
        this%Ahat(:,i+coord*this%Ne) = this%ensemble(i+coord*this%Ne)%Phat
    end do

#ifdef PPMPI
    call mpi_allreduce(this%A, this%SE_array, this%Ns*this%Ne_tot, MPI_RPREC,      &
        MPI_SUM, comm, ierr)
    call mpi_allreduce(this%Ahat, this%ME_array, this%Nm*this%Ne_tot, MPI_RPREC,   &
        MPI_SUM, comm, ierr)
    this%A = this%SE_array
    this%Ahat = this%ME_array
#endif
  if (coord == 0) then
    do i = 1, this%Ne_tot
        this%Abar = this%Abar + this%A(:,i) / this%Ne_tot
        this%Ahatbar = this%Ahatbar + this%Ahat(:,i) / this%Ne_tot
    end do

    do i = 1, this%Ne_tot
        this%Aprime(:,i) = this%A(:,i) - this%Abar
        this%Ahatprime(:,i) = this%Ahat(:,i) - this%Ahatbar
    end do

!  if (coord == 0)

     write(*,*) "Done generating initial ensemble"
  endif
end subroutine generateInitialEnsemble

subroutine advance(this, dt, Pm, Ctp)
    use util, only : random_normal, inverse
    use functions, only: linear_interp
#ifdef PPIFORT
    use BLAS95
#endif
    implicit none
    class(WakeModelEstimator), intent(inout)    :: this
    real(rprec), intent(in)                     :: dt
    real(rprec), dimension(:), intent(in)       :: Pm, Ctp
    integer                                     :: i, j, N, Nx, jstart, jend, m
    real(rprec)                                 :: Uinftyi, alpha
    real(rprec), dimension(:), allocatable      :: yu_turb, yu_temp, ylocs_sort


    N = this%wm%N
    Nx = this%wm%Nx

    allocate( yu_turb(this%wm%nfree) )
    allocate( ylocs_sort(this%wm%nfree) )
    allocate( yu_temp(this%wm%Ny) )
    
    if (size(Pm) /= N ) then
        call error('WakeModelEstimator.advance','Pm must be of size N')
    end if
    if (size(Ctp) /= N ) then
        call error('WakeModelEstimator.advance','Ctp must be of size N')
    end if

   if (coord == 0) then
    ! Calculate noisy measurements
    this%E = 0._rprec
    do i = 1, N
        do j = 1, this%Ne_tot
            this%E(i, j) = this%sigma_Phat * random_normal()
            this%D(i, j) = Pm(i) + this%E(i, j)
        end do
    end do
    this%Dprime = this%D - this%Ahat

    ! Update Anew = A + A'*Ahat'^T * (Ahat'*Ahat'^T + E*E^T)^-1 * D'
    ! Since the dimension is small, we don't bother doing the SVD. If the matrix becomes
    ! singular, then this should be considered as in section 4.3.2 of Everson(2003)
!#ifdef PPIFORT
!    call gemm(this%E, this%E, this%MM_array, 'n', 't', 1._rprec, 0._rprec)
!    call gemm(this%Ahatprime, this%Ahatprime, this%MM_array, 'n', 't',         &
!        1._rprec, 1._rprec)
!    this%MM_array = inverse(this%MM_array)
!    call gemm(this%MM_array, this%Dprime, this%ME_array, 'n', 'n', 1._rprec,   &
!        0._rprec)
!    call gemm(this%Aprime, this%Ahatprime, this%SM_array, 'n', 't', 1._rprec,  &
!        0._rprec)
!    call gemm(this%SM_array, this%ME_array, this%A, 'n', 'n', 1._rprec,        &
!        1._rprec)
!#else
    this%A = this%A + matmul( matmul(this%Aprime, transpose(this%Ahatprime)),  &
        matmul(inverse(matmul(this%Ahatprime, transpose(this%Ahatprime)) +     &
        matmul(this%E, transpose(this%E))), this%Dprime))
!#endif
    else 
       this%A = 0._rprec
    endif

#ifdef PPMPI
    call mpi_allreduce(this%A, this%SE_array, this%Ns*this%Ne_tot, MPI_RPREC,      &
        MPI_SUM, comm, ierr)
    this%A = this%SE_array
#endif
    
 ! Compute mean
    this%Abar = 0._rprec
    do i = 1, this%Ne_tot
        this%Abar = this%Abar + this%A(:,i) / this%Ne_tot;
    end do

 ! Filter U_infty
    alpha = dt / (this%tau_U_infty + dt)
    Uinftyi = 0._rprec
    do j = 1, size(this%wm%free_turbines)
        i = this%wm%free_turbines(j)
!        Uinftyi = Uinftyi + (Pm(i) / this%wm%Ctp(i))**(1._rprec/3._rprec)      &
!            * (4._rprec + this%wm%Ctp(i)) / 4._rprec / size(this%wm%free_turbines)
        yu_turb(j) = (Pm(i) / this%wm%Ctp(i))**(1._rprec/3._rprec)      &
            * (4._rprec + this%wm%Ctp(i)) / 4._rprec
    end do

   call this%sortFreeTurbines(yu_turb, ylocs_sort)

    yu_temp = linear_interp(ylocs_sort,yu_turb,this%wm%y)
    this%wm%yu_inf = alpha * yu_temp + (1 - alpha) * this%wm%yu_inf
    this%wm%U_infty = sum(this%wm%yu_inf) / size(this%wm%yu_inf)
!    this%wm%U_infty = alpha * Uinftyi + (1 - alpha) * this%wm%U_infty
    this%wm%VELOCITY = this%wm%U_infty
    this%wm%TIME  = this%wm%LENGTH / this%wm%VELOCITY
    this%wm%FORCE = this%wm%VELOCITY / this%wm%TIME

    ! Fill into objects
    do i = 1, this%Ne
        do j = 1, N
            jstart = (j-1)*Nx+1
            jend   = j*Nx
            do m = 1, Nx
                this%ensemble(i+coord*this%Ne)%du(j,m) = max(this%A((m-1)+jstart,i+coord*this%Ne),0._rprec)
            enddo 
!            this%ensemble(i+coord*nproc)%du(j,:)  = max(this%A(jstart:jend,i+coord*nproc), 0._rprec)
            this%ensemble(i+coord*this%Ne)%k(j) = max(this%A(Nx*N+j,i+coord*this%Ne), 0._rprec)
        end do
        this%ensemble(i+coord*this%Ne)%U_infty  = this%wm%U_infty
        this%ensemble(i+coord*this%Ne)%yu_inf = this%wm%yu_inf
        this%ensemble(i+coord*this%Ne)%VELOCITY = this%wm%VELOCITY
        this%ensemble(i+coord*this%Ne)%TIME  = this%wm%TIME
        this%ensemble(i+coord*this%Ne)%FORCE = this%wm%FORCE
    end do
    do j = 1, N
        jstart = (j-1)*Nx+1
        jend = j*Nx
        this%wm%du(j,:) = max(this%Abar(jstart:jend), 0._rprec)
    end do
    this%wm%k(:) = max(this%Abar((N*Nx+1):), 0._rprec)

    ! Advance ensemble and mean estimate
    call this%advanceEnsemble(Ctp, dt)

    ! Place ensemble into a matrix with each member in a column
    this%A = 0._rprec
    this%Ahat = 0._rprec
    this%Abar = 0._rprec
    this%Ahatbar = 0._rprec
    do i = 1, this%Ne
        do j = 1, N
            jstart = (j-1)*Nx+1
            jend = j*Nx
            this%A(jstart:jend,i+coord*this%Ne) = this%ensemble(i+coord*this%Ne)%du(j,:)
        end do
        this%A((N*Nx+1):,i+coord*this%Ne) = this%ensemble(i+coord*this%Ne)%k(:)
        this%Ahat(:,i+coord*this%Ne) = this%ensemble(i+coord*this%Ne)%Phat
    end do
#ifdef PPMPI
    call mpi_allreduce(this%A, this%SE_array, this%Ns*this%Ne_tot, MPI_RPREC,      &
        MPI_SUM, comm, ierr)
    call mpi_allreduce(this%Ahat, this%ME_array, this%Nm*this%Ne_tot, MPI_RPREC,   &
        MPI_SUM, comm, ierr)
    this%A = this%SE_array
    this%Ahat = this%ME_array
#endif

    do i = 1, this%Ne_tot
        this%Abar = this%Abar + this%A(:,i) / this%Ne_tot
        this%Ahatbar = this%Ahatbar + this%Ahat(:,i) / this%Ne_tot
    end do

    do i = 1, this%Ne_tot
        this%Aprime(:,i) = this%A(:,i) - this%Abar
        this%Ahatprime(:,i) = this%Ahat(:,i) - this%Ahatbar
    end do

end subroutine

subroutine advanceEnsemble(this, Ctp, dt)
    use param, only : pi
    use util, only  : random_normal
    implicit none
    class(WakeModelEstimator), intent(inout) :: this
    real(rprec), dimension(:), intent(in)    :: Ctp
    real(rprec), intent(in)                  :: dt
    integer                                  :: i, j, N
    ! Advance step for objects
    N = this%wm%N
    do i = 1, this%Ne
        do j = 1,N
            this%ensemble(i+coord*this%Ne)%k(j) = max(this%ensemble(i+coord*this%Ne)%k(j)                      &
                + sqrt(dt) * this%sigma_k * random_normal(), 0._rprec)
            this%ensemble(i+coord*this%Ne)%du(j,:) = max(this%ensemble(i+coord*this%Ne)%du(j,:)                &
                + sqrt(dt) * this%sigma_du * random_normal()                       &
                * this%ensemble(i+coord*this%Ne)%G(j,:) * sqrt(2*pi) * this%ensemble(i+coord*this%Ne)%Delta,   &
    0._rprec)
        end do
        call this%ensemble(i+coord*this%Ne)%computeWakeExpansionFunctions
        call this%ensemble(i+coord*this%Ne)%compute2Dwakes
        call this%ensemble(i+coord*this%Ne)%wokeTurbines
        call this%ensemble(i+coord*this%Ne)%advance(Ctp, dt)
    end do
    call this%wm%computeWakeExpansionFunctions
    call this%wm%compute2Dwakes
    call this%wm%woketurbines
    call this%wm%compute_wake_disks
    call this%wm%advance(Ctp, dt)

end subroutine advanceEnsemble

subroutine sortFreeTurbines(this, yu_turb, ylocs_sort)
    use param, only : pi
    use util, only  : random_normal
    implicit none
    class(WakeModelEstimator), intent(inout) :: this
    real(rprec), dimension(:), intent(inout) :: yu_turb, ylocs_sort
    real(rprec), dimension(:)                :: buffer(2)
    real(rprec), dimension(:,:), allocatable :: sort_array
    integer                                  :: i, j

    ! Sort the free turbines in order from smallest y-location to largest
    allocate( sort_array(this%wm%nfree,2) )
    sort_array(:,2) = yu_turb
    sort_array(:,1) = this%wm%s(this%wm%free_turbines,2)
   
    do i = 1, this%wm%nfree
        j = minloc(sort_array(i:this%wm%nfree, 1), dim=1) + i - 1

        buffer(:) = sort_array(i,:)
        sort_array(i,:) = sort_array(j,:)
        sort_array(j,:) = buffer(:)
    enddo

    yu_turb(:) = sort_array(:,2)
    ylocs_sort(:) = sort_array(:,1)

end subroutine sortFreeTurbines

end module wake_model_estimator_class
