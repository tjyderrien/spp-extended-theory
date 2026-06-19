!
!    ==================================================================
!    interspectus -- sum_rule.f90 : transformation of the respective sum
!                                   rule to the general form;
!       + + + This file must be edited according to your needs + + +
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! Function to take care of whatever happens in the integral which in the
! case of S(q,\omega) corresponds to f(omega) = S(q,\omega)*\omega
!
! => make sure the sum_rule_func_inv is consistent with that!!!

double precision function sum_rule_func(omega,f,q)
implicit none
double precision omega,f,q
!----------------------------------------------------------------------
!!!!!!   change here the sum-rule function             
sum_rule_func = f*omega ! don't forget to change also sum_rule_func_inv
!----------------------------------------------------------------------
return
end function sum_rule_func


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------------------------------------------------
! INVERSE of the Function to take care of whatever happens in the 
! integral which in the
! case of S(q,\omega) corresponds to f(omega) = S(q,\omega)*\omega

double precision function sum_rule_func_inv(omega,f,q)
implicit none
double precision omega,f,q

!----------------------------------------------------------------------
!!!!!!   change here the inverse of the sum-rule function             
sum_rule_func_inv = f /(omega+.0000000001)  
! don't forget to change also sum_rule_func
return
end function sum_rule_func_inv
