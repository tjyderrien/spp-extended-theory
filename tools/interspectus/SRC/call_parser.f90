!    ==================================================================
!    interspectus -- call_parser.f90 : interface to use parser of input file
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

 subroutine read_input(inputfile,interpol_type,n_inputs,q_flag_single,q_out_single,q_flag_series,q_out_array,stepsize_integral,stepsize_result,filename_compact,debug)
 
  implicit none
  character (len=*) :: inputfile
  integer           :: nl
  
  integer, external :: read_int, read_string, read_float, read_input_file
  integer, external :: read_float_array

  !!!! hansi
  integer interpol_type,n_inputs,i
  real q_out_array(3),q_out_single
  real r_inputs, stepsize_integral,stepsize_result
  logical q_flag_single,q_flag_series,filename_compact
  logical debug

  character*20, external ::  char_write
  !!! end hansi

  ! parse input file (or 'stdin')
  write(6,*) 'reading parameters from file: ', inputfile
  if(read_input_file(inputfile,10) /= 1) stop 'error in inputfile'

! first test of Ralf: exciton
  if(read_string('exciton',7)>0) then
    write(6,*) "exciton"
  endif

! the general output file form
  filename_compact = .TRUE. 
  if(read_string('outfile_sortable',16)>0) then
     filename_compact = .FALSE.
  endif

! output of files for debugging 
  debug = .FALSE. 
  if(read_string('debug',5)>0) then
     debug = .TRUE.
  endif

!!! which interpolation type to use
  interpol_type=77
  call read_int('interpol_type',13,interpol_type) 
  select case (interpol_type)
  case(1)
     write(*,*) 'linear interpolation, interpol_type = 1'
  case(2)
     write(*,*) 'polynomial interpolation, interpol_type = 2'
  case(3)
     write(*,*) 'cubic spline interpolation, interpol_type = 3'
  case(4)
     write(*,*) 'akima spline interpolation, interpol_type 4'
  case default
     STOP 'STOP -- invalid interpol_type' 
  end select
  

  ! number of input files
  n_inputs=0
  if(read_int('n_inputs',8,n_inputs)>0) then
     write(*,*) 'n_inputs = ',char_write(n_inputs)
  else
     write(*,*) 'ERROR: parameter n_inputs is not set properly'
  end if
  
!!! single q value for output
  q_flag_single = .FALSE.
  q_out_single=0
  if(read_float('q_out_single',12,q_out_single)>0) then
!     write(*,*) 'q_out_single = ',q_out_single
     q_flag_single = .TRUE.
     !  else
     !     write(*,*) 'ERROR: parameter q_out_single is not set properly'
  end if
  
!!! series of q values for output
  q_flag_series = .FALSE.
  q_out_array = 0
  if (read_float_array('q_out_series',12,q_out_array,3) == 1) then
!     write(*,*) 'use q_out_series = ',(q_out_array(i),i=1,3)
     q_flag_series = .TRUE.
     !  else 
     !     write(*,*) 'ERROR: parameter q_out_series is not set properly'
  end if
  
!!! test to make sure only single q value OR series is entered:
  if(q_flag_single.AND.q_flag_series) then 
     STOP 'ERROR: Specify only single q value OR series of q values for output'
  else if(.NOT.q_flag_single.AND..NOT.q_flag_series) then
     STOP 'ERROR: q for output not found; Specify single q value or series of q values for output'
  end if
  
!!!  stepsize along F(q,w), i.e., the NEW X (after inversion)
  stepsize_integral=0
  if(read_float('stepsize_integral',17,stepsize_integral)>0) then
     write(*,*) 'stepsize_integral = ',stepsize_integral
  else
     stepsize_integral=0.001
     write(*,*) 'stepsize_integral not in input file, using default =', stepsize_integral
  end if
  
!!!  stepsize along F(q,w), i.e., the NEW X (after inversion)
  stepsize_result=0
  if(read_float('stepsize_result',17,stepsize_result)>0) then
     write(*,*) 'stepsize_result = ',stepsize_result
  else
     stepsize_result=0.001
     write(*,*) 'stepsize_result not in input file, using default = ',stepsize_result
  end if
  
  
  
end subroutine read_input


