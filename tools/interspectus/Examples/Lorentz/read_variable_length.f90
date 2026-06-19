!
!    ==================================================================
!    interspectus -- read_variable_length.f90 : check input data files
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
!------------------------------------------------------------------------


subroutine read_variable_lenth(file_name,n_lines)

! Careful, as it is now, it reads the file
! 2 4 
! 3 
! 5 6
! as
! 2 4
! 3 5 etc. --- can I put some format that requires him to 
! read exactly two reals per line but without specifying 
! the format of those?

  implicit none
  
  integer :: IO, i, integer, n_lines
  real :: X,Y 
  character*200 file_name
  
  open(33,file=file_name)
  
  i = 0
  Read_Loop: DO
     READ( 33, *, IOSTAT = IO ) X,Y
!     write(*,*) IO
!     IF (IO < 0) THEN
     if (IO == -1) then
        close(33)
        EXIT Read_Loop
     END IF
     
     if (IO > 0 .OR. IO < -1) then
        close(33)
        STOP 'ERROR: Something seems to be wrong with the input file'
     END IF
     
     i = i + 1

  END DO Read_Loop
  
  n_lines = i
!  write(*,*) 'The file contains',i, 'lines of data'

  return
  
end subroutine read_variable_lenth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine check_order_omega(file_name,n_lines)
  implicit none
  
  integer :: IO, i, integer, n_lines
  DOUBLE PRECISION :: x(n_lines)  
  character*200 file_name
  
  open(33,file=file_name)

  read(33,*) x(1)
  do i = 2, n_lines
     read(33,*) x(i)     
     if(x(i).lt.x(i-1)) then
        write(*,*) 'ERROR: in input file ', trim(file_name)
        write(*,*) '.......  line' , i-1,' ===> ', i
        stop       ' .......  order of omegas in input files must be monotonously increasing'
     end if
  end do
  write(*,*) 'Order of omega in ', trim(file_name), ' is OK.'
  write(*,*) ''
  return
end subroutine check_order_omega
