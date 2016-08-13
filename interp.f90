!!
!!  Copyright (C) 2009-2016  Johns Hopkins University
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

!*******************************************************************************
module interp
!*******************************************************************************
use messages
use types, only : rprec
implicit none
private

public linear_interp, bilinear_interp

character (*), parameter :: mod_name = 'interp'

interface linear_interp
    module procedure :: linear_interp_ss
    module procedure :: linear_interp_sa
    module procedure :: linear_interp_aa
end interface linear_interp

interface bilinear_interp
    module procedure :: bilinear_interp_ss
    module procedure :: bilinear_interp_sa
    module procedure :: bilinear_interp_aa
end interface bilinear_interp

contains

!**********************************************************************
real(rprec) function bilinear_interp_ss(u11,u21,u12,u22,dx,dy,xdiff,ydiff)
!**********************************************************************
!
!  This function performs linear interpolation 
!  
!  Inputs:
!  u11          - lower bound value in x direction for lower y
!  u21          - upper bound value in x direction for lower y
!  u12          - lower bound value in x direction for upper y
!  u22          - upper bound value in x direction for upper y
!  dx           - length delta for the grid in x direction
!  dy           - length delta for the grid in y direction
!  xdiff        - distance from the point of interest to the u11 node in x direction
!  xdiff        - distance from the point of interest to the u11 node in y direction
!
use types, only : rprec
implicit none

real(rprec), intent(in) :: u11, u12, u21, u22, dx, dy, xdiff, ydiff
real(rprec) :: v1, v2

v1 = linear_interp(u11, u21, dx, xdiff)
v2 = linear_interp(u12, u22, dx, xdiff)

bilinear_interp_ss = linear_interp(v1,v2,dy,ydiff)

return
end function bilinear_interp_ss

!**********************************************************************
function bilinear_interp_sa_nocheck(x, y, v, xq, yq) result(vq)
!**********************************************************************
!
!  This function performs the linear interpolation fo bilinear_interp_sa
!  and bilinear_interp_aa without checking bounds of the input arrays.
!  It should not be called directly.
!
implicit none
real(rprec), dimension(:), intent(in) :: x, y
real(rprec), dimension(:,:), intent(in) :: v
real(rprec), intent(in) :: xq, yq
real(rprec) :: vq
integer     :: i, j, Nx, Ny

Nx = size(x)
Ny = size(x)
i = binary_search(x, xq)
if (i == 0) then
    vq = linear_interp(y, v(1,:), yq)
else if (i == Nx) then
    vq = linear_interp(y, v(Nx,:), yq)
else
    j = binary_search(y, yq)
    if (j == 0) then
        vq = linear_interp(v(i,1), v(i+1,1), x(i+1)-x(i), xq - x(i))
    else if (j == Ny) then
        vq = linear_interp(v(i,Ny), v(i+1,Ny), x(i+1)-x(i), xq - x(i))    
    else
        vq = bilinear_interp_ss( v(i,j), v(i+1,j), v(i,j+1), v(i+1,j+1), &
             x(i+1) - x(i), y(j+1) - y(j), xq - x(i), yq - y(j) )
    end if
end if
    
end function bilinear_interp_sa_nocheck

!**********************************************************************
function bilinear_interp_sa(x, y, v, xq, yq) result(vq)
!**********************************************************************
!
!  This function performs linear interpolation from a set of points 
!  defined on a grid (x,y) with values v to a query point (xq, yq)
!  
!  Inputs:
!  x            - array of sample points
!  v            - array of values at sample points
!  xq           - query point
!
implicit none
real(rprec), dimension(:), intent(in) :: x, y
real(rprec), dimension(:,:), intent(in) :: v
real(rprec), intent(in) :: xq, yq
real(rprec) :: vq
character(*), parameter :: func_name = mod_name // '.bilinear_interp_sa'

if ( size(v,1) /= size(x) .or. size(v,2) /= size(y)) then
    call error(func_name, 'Array v must have a size of [size(x), size(y)]')
end if

vq = bilinear_interp_sa_nocheck(x, y, v, xq, yq)
    
end function bilinear_interp_sa

!**********************************************************************
function bilinear_interp_aa(x, y, v, xq, yq) result(vq)
!**********************************************************************
!
!  This function performs linear interpolation from a set of points 
!  defined on a grid (x,y) with values v to an array of query points
!  (xq, yq)
!  
!  Inputs:
!  x            - array of sample points
!  v            - array of values at sample points
!  xq           - array of query points
!
implicit none
real(rprec), dimension(:), intent(in) :: x, y
real(rprec), dimension(:,:), intent(in) :: v
real(rprec), dimension(:), intent(in) :: xq, yq
real(rprec), dimension(:), allocatable :: vq
integer :: i, N
character(*), parameter :: func_name = mod_name // '.bilinear_interp_sa'

if ( size(v,1) /= size(x) .or. size(v,2) /= size(y)) then
    call error(func_name, 'Array v must have a size of [size(x), size(y)]')
end if
if ( size(xq) /= size(yq) ) then
    call error(func_name, 'Arrays xq and yq must be the same size')
end if

N = size(xq)
allocate(vq(N))

do i = 1, N
    vq(i) = bilinear_interp_sa_nocheck(x, y, v, xq(i), yq(i))
end do
    
end function bilinear_interp_aa

!**********************************************************************
real(rprec) function linear_interp_ss(u1,u2,dx,xdiff)
!**********************************************************************
!
!  This function performs linear interpolation 
!  
!  Inputs:
!  u1           - lower bound value in the increasing index direction
!  u2           - upper bound value in the increasing index direction
!  dx           - length delta for the grid in the correct direction
!  xdiff        - distance from the point of interest to the u1 node
!
use types, only : rprec
implicit none

real(rprec), intent(in) :: u1, u2, dx, xdiff

linear_interp_ss = u1 + (xdiff) * (u2 - u1) / dx

return
end function linear_interp_ss

!**********************************************************************
function linear_interp_sa_nocheck(x, v, xq) result(vq)
!**********************************************************************
!
!  This function performs the linear interpolation fo linear_interp_sa
!  and linear_interp_aa without checking bounds of the input arrays.
!  It should not be called directly.
!
implicit none
real(rprec), dimension(:), intent(in) :: x, v
real(rprec), intent(in) :: xq
real(rprec) :: vq
integer :: i, N

N = size(v)
i = binary_search(x, xq)
if (i == 0) then
    vq = v(1)
else if (i == N) then
    vq = v(N)
else
    vq = linear_interp_ss(v(i), v(i+1), x(i+1)-x(i), (xq - x(i)))
end if
    
end function linear_interp_sa_nocheck

!**********************************************************************
function linear_interp_sa(x, v, xq) result(vq)
!**********************************************************************
!
!  This function performs linear interpolation from a set of points x
!  with values v to a query point xq
!  
!  Inputs:
!  x            - array of sample points
!  v            - array of values at sample points
!  xq           - query point
!
implicit none
real(rprec), dimension(:), intent(in) :: x, v
real(rprec), intent(in) :: xq
real(rprec) :: vq
character(*), parameter :: func_name = mod_name // '.linear_interp_sa'

if ( size(v) /= size(x) ) then
    call error(func_name, 'Arrays x and v must be the same size')
end if

vq = linear_interp_sa_nocheck(x, v, xq)
    
end function linear_interp_sa

!**********************************************************************
function linear_interp_aa(x, v, xq) result(vq)
!**********************************************************************
!
!  This function performs linear interpolation from a set of points x
!  with values v to an array of query points xq
!  
!  Inputs:
!  x            - array of sample points
!  v            - array of values at sample points
!  xq           - array of query points
!
implicit none
real(rprec), dimension(:), intent(in) :: x, v
real(rprec), dimension(:), intent(in) :: xq
real(rprec), dimension(:), allocatable :: vq
integer :: i, N
character(*), parameter :: func_name = mod_name // '.linear_interp_aa'

! Check array sizes
if ( size(v) /= size(x) ) then
    call error(func_name, 'Arrays x and v must be the same size')
end if

! Allocate output 
N = size(xq)
allocate(vq(N))

! For each element of the array perform interpolation
do i = 1, N
    vq(i) = linear_interp_sa_nocheck(x, v, xq(i))
end do
    
end function linear_interp_aa

!**********************************************************************
function binary_search(arr,val) result(low)
!**********************************************************************
!
!  This function performs a binary search on a sorted array. Given the 
!  provided value, adjacent low and high indices of the array are found
!  such that the provided value is bracketed. Guaranteed log2(N) search.
!  
!  Inputs:
!  arr          - sorted array of values to search
!  val          - value to be bracketed
!
!  Output:
!  low          - lower index of the array bracket 
!                 0 if val < arr(1), N if val < arr(N))
!
implicit none

real(rprec), dimension(:) :: arr
real(rprec) :: val
integer :: low, mid, high, N

! Size of array
N = size(arr)

! Check if value is outside bounds
if ( val < arr(1) ) then
    low = 0
    return
end if
if ( val > arr(N) ) then
    low = N
    return
end if

! Otherwise perform bisection
low = 1
high = N
do while (high - low > 1)
    mid = (low + high) / 2
    if ( arr(mid) > val ) then
        high = mid
    elseif ( arr(mid) < val ) then
        low = mid
    else
        low = mid
        return
    endif
end do

return
end function binary_search

end module interp