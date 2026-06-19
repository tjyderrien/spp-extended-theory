!    ==================================================================
!    interspectus -- coeffs.f90 : expansion coefficients for Taylor
!    ==================================================================
!
!    Copyright (C) 2009, Hans-Christian Weissker,
!
!    This program is free software; you can redistribute it and/or
!    modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation; either version 2 of
!    the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place - Suite 330,
!    Boston, MA 02111-1307, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! calc. the linear and the quadratic coefficients of the taylor expansion
! of the curves at some point X

subroutine coeffs (n_files,n_min, n_max, &
     f, & ! array of input functions: f(n_files,n_max_points_X) 
     q, & ! q vectors of n_files input files: q(n_files) FOR NOW: 3
     lincoef,quadcoef)
  
  implicit none
  integer n_files ! # of input functions (= # of files)
  integer n_min ! # points in smallest dataset (i.e., to be treated)
  integer n_max ! # points in longest dataset (for alloc)
  integer i,j
  double precision   q(n_files)   ! q1,q2,q3
  double precision   f(n_files,n_max)  ! f1(n_max),f2(n_max),f3(n_max)
  double precision lincoef(n_min), quadcoef(n_min)

  quadcoef=0.d0
  lincoef =0.d0
! version for quadratic procedure -- requires 3 input curves
  if(n_files.eq.3) then          
     write(*,*) '3 input files => quad. interpolation'
     do i=1,n_min  
! orig  -- this did produce the funny behavior, deviating from the inputs for the highest x  :       lincoef(i)  = ( f(3,i) - f(1,i) ) / ( q(3) - q(1) )
        lincoef(i)  = ( f(3,i) - f(1,i) ) / ( q(3) - q(1) )
        quadcoef(i) = ( f(3,i) + f(1,i) -  2*f(2,i) ) / ( (q(3)-q(2) )**2 )
!        write(*,*) i,quadcoef(i)
     end do
  end if
! version for linear procedure -- two input curves
  if(n_files.eq.2) then          
     do i=1,n_min  
        lincoef(i)  = ( f(2,i) - f(1,i) ) / ( q(2) - q(1) )
     end do
  end if
  
  return
end subroutine coeffs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine interpol_curves(n_files,n_min_points_X, &
     n_max_points_X, & !!! this is new!
     f_in, & 
     lincoef,quadcoef,  &
     q_in, &
     q, f_out)
  
  implicit none
  integer n_files, n_pts, n_max_points_X,n_min_points_X
  
  double precision  q_in(n_files)
  double precision  f_in(n_files,n_max_points_X)  
  double precision  q
  double precision  lincoef(n_min_points_X)
  double precision  quadcoef(n_min_points_X)
!! wrong  double precision  f_out(n_pts)
  double precision  f_out(n_max_points_X)
  
  integer i

!     write(*,*) 'n_files = ',n_files,'n_pts = ',n_pts
!     write(*,*) 'q_in = ',q_in
!     write(*,*) 'q = ',  q
  
!!! version for quadratic procedure -- requires 3 input curves
  if(n_files.eq.3) then          
     write(*,*) 'quadratic interpolation between curves'
     do i=1,n_min_points_X
        ! OK     write(*,*) i, lincoef(i) ,quadcoef(i)
        !      write(*,*) i, f_in(3,i) !, f_in(2,i)
        
        f_out(i) = f_in(2,i) &
             + lincoef(i)  * ( q - q_in(2) )  &
             + quadcoef(i) * ( ( q - q_in(2) )**2 ) / 2.d0

     enddo                  !   i => n_min_points_X      
  end if
  
  
!!! version for linear procedure -- 2 input curves
!!! extrapolation starts from the first input curve
  if(n_files.eq.2) then          
!     write(*,*) 'linear interpolation between curves'
     do i=1,n_min_points_X
        
        f_out(i) = f_in(1,i) &
             + lincoef(i)  * ( q - q_in(1) )  
!        write(*,*)  i , f_out(i)
     enddo                  ! n         
  end if
     
     

end subroutine interpol_curves


