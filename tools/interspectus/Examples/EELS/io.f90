!
!    ==================================================================
!    interspectus -- io.f90 : read input data + write function 
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


!   reads the input for file_in for the spectra 
subroutine read_file(file_in,dim_in,w,fkt)
  
  implicit none
  integer  dim_in, i
  double precision w(dim_in)
  double precision fkt(dim_in)
  character*80 file_in
  
  w = 0.d0
  fkt = 0.d0
  
  open( 20,file = file_in)  
  
  ! read the data  
  do i = 1,dim_in

     read(20,*) w(i),fkt(i)
     !  write(*,*) w(i),fkt(i)
  enddo
  !  print*,'Hola 2'
  close(20)
  
end subroutine read_file


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write an integer as a character with length as necessary
! without having to use a format
function char_write(i)
implicit none
integer j,i_pos,i
character*20  ctemp, cctemp,char_write

cctemp = '                    '
     write(ctemp, '(i10)' )  i
     i_pos = 1
     do j=1,20
        if (ctemp(j:j)/=' ') then
           cctemp(i_pos:i_pos) = ctemp(j:j) 
           i_pos = i_pos + 1         
        end if
     end do
     char_write = cctemp

return
end function char_write
