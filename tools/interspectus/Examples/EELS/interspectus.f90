!----------------------------------------------------------------------
!
!    ==================================================================
!    interspectus -- INTERpolate SPECtra Using Sum rules
!    ==================================================================
!
!    -- interpolates between positive-definite spectra that obey
!    a sum rule and depend on one parameter (e.q., momentum transfer)
!    
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
!-----------------------------------------------------------------------

program interspectus

  implicit none

  integer i, j, ifile, count
  
! for read-in files
  integer n_inputs     ! # of read-in files

! for interpolation type
  integer interpol_type

! dimension of the 2 or 3 read-in curves
  integer, dimension(:), ALLOCATABLE :: dim_in   ! number of points for each
  ! dataset seperately
  integer max_dim_in ! this is used for dimension of input valule array

! names of the 2 or 3 read-in curves
  character*200 , dimension(:), ALLOCATABLE :: file_in ! names of input files

! values of the read-in curves
  double precision, dimension(:,:), ALLOCATABLE :: w_in ! (n_dataset,values)
  double precision, dimension(:,:), ALLOCATABLE :: y_in ! (n_dataset,values)

! now make an array as well:  double precision q1, q2, q3
  double precision, dimension(:), ALLOCATABLE :: q_in
  double precision q_out_min, q_out_max

! array for "distribution-function-like" sum on input grid
  double precision, dimension(:,:), ALLOCATABLE ::   dist_in
  double precision, external ::  sum_rule_func ! for integrand 

! array for INVERSED "distribution-function-like" sum on input grid
  double precision, dimension(:,:), ALLOCATABLE ::  X_inv_int_in
  double precision, dimension(:,:), ALLOCATABLE ::  Y_inv_int_in

  double precision, dimension(:), ALLOCATABLE :: X_test
  double precision, dimension(:), ALLOCATABLE :: Y_test

  double precision, dimension(:), ALLOCATABLE :: X_interp_test
  double precision, dimension(:), ALLOCATABLE :: Y_interp_test

  integer n_max_points_X ! length of longest set after int along new X
  integer n_min_points_X ! length of shortest set after int along new X
  double precision max_X ! sort of dummy for comparison of length of sets
  double precision min_X ! sort of dummy for comparison of length of sets

! array for interpolated INVERSED sum on NEW grid
! (i.e., same grid for all curves -- then we can interpolate)
  !!! MUST BE REAL (not DOUBLE) for passing from c
  REAL stepsize_s ! this is the absolute stepsize along the NEW X
  double precision, dimension(:,:), ALLOCATABLE ::  X_inv_interp 
  double precision, dimension(:,:), ALLOCATABLE ::  Y_inv_interp

! for interpolation between curves
  double precision, dimension(:), ALLOCATABLE :: lincoef
  double precision, dimension(:), ALLOCATABLE :: quadcoef

  double precision, dimension(:,:), ALLOCATABLE :: f_out_inv
  double precision , dimension(:), ALLOCATABLE :: q_out
  integer nq_out, iq

! for returned fct., i.e., newly inversed back, but still on the
! irreagular grid x_int_out resulting from the function values
  double precision, dimension(:,:), ALLOCATABLE :: x_int_out
  double precision, dimension(:,:), ALLOCATABLE :: y_int_out
  double precision max_out , min_out
  integer imax_out, imin_out
! must become real for passing  double precision  stepsize_result

! for derivative of the returned distribution function
  double precision, dimension(:), ALLOCATABLE :: x_int_out_pass
  double precision, dimension(:), ALLOCATABLE :: y_int_out_pass

  double precision, dimension(:), ALLOCATABLE :: deriv_x
  double precision, dimension(:), ALLOCATABLE :: deriv_y
  
! for returning
  double precision, dimension(:,:), ALLOCATABLE ::   f_out

! inverse function of sum_rule_func
  double precision, external ::  sum_rule_func_inv

! declaration for gsl interpolation instead of Taylor

  logical TAYLOR
  double precision q_step
  integer i_point
  double precision, dimension(:), ALLOCATABLE :: Aq ! (n_inputs)
  double precision, dimension(:), ALLOCATABLE :: Xq ! (n_inputs)
  double precision, dimension(:), ALLOCATABLE :: Aq_interp ! (nq_out)
  double precision, dimension(:), ALLOCATABLE :: Xq_interp ! (nq_out)
  
!!! for writing to the files
!!  character*40 appendix,rm_cmd,tmp_file,outfile,outfile_trimmed
  logical filename_compact
  character*80  :: outfile = ' ', name_base, long, name_part
  character*20  ctemp, cctemp
  integer i_pos, length_name_base
  logical debug
  character*200 debug_file 
  character*20, external ::  char_write


!!! For the finalized version -- Einsortieren !
  LOGICAL q_flag_single,q_flag_series
  REAL q_out_single,q_out_array(3) ! must be REAL (not DOUBLE) for passing
  REAL stepsize_result_s
  !  REAL stepsize,stepsize_result
  DOUBLE PRECISION stepsize_result, stepsize
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TAYLOR = .FALSE.! False-> uses gsl; if  TAYLOR=.TRUE. my own Taylor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
! call input file for parameters
  call read_input("input.pars",interpol_type,n_inputs,q_flag_single,q_out_single,q_flag_series,q_out_array,stepsize_s,stepsize_result_s,filename_compact,debug)
  
! type conversion into DOUBLE -- is this necessary?
  stepsize = stepsize_s
  stepsize_result = stepsize_result_s

!!! the single q input does not yet work. (The problem is that for
!!! passing to the gsl, we need a step size for the function to be
!!! returned
  if (q_flag_single) then
     nq_out = 1   
     q_out_min = q_out_single
     q_out_max = q_out_single
  end if
  
  if (q_flag_series) then
     q_out_min = q_out_array(1)
     q_out_max = q_out_array(2)
     nq_out = int(q_out_array(3))   
  end if
  
! calculate the parameter values for the output
  ALLOCATE( q_out (nq_out) )
  if(q_flag_series) then
     do i=1,nq_out
        q_out(i) = q_out_min + (q_out_max - q_out_min) / (nq_out-1) * (i-1)
     end do
  end if
!!! the single q input does not yet work. (The problem is that for
!!! passing to the gsl, we need a step size for the function to be
!!! returned
!  if(q_flag_single) then
!        q_out(1) = q_out_single 
!  end if


  ALLOCATE(  q_in(n_inputs)  )
  ALLOCATE(  dim_in(n_inputs)  )
  ALLOCATE(  file_in(n_inputs)  )
  
  write(*,*) ''
! read input-file information from the file "input.files"
  write(*,*) 'reading file names from file: input.files'
  open(25,file='input.files')
  read(25,*) name_base
  do i = 1, n_inputs
     read(25,*) q_in(i)  ! q value for file
     read(25,*) file_in(i) ! file name

     write(*,*) 'file_in(',trim(char_write(i)),') = ',trim(file_in(i))          
     write(*,*) 'q_in(',trim(char_write(i)),') = ',q_in(i)
     write(*,*) ''
  end do
  close(25)
  
! determine the number of lines in the spectra input files
do ifile = 1, n_inputs
   call  read_variable_lenth(file_in(ifile),dim_in(ifile)) ! get dim_in(i)
   write(*,*) 'Input file  ',trim(file_in(ifile)), '  ..........'
   write(*,*) '.........  for parameter q = ',q_in(ifile)
   write(*,*) '.........  contains',dim_in(ifile) , 'line of data'
! and check if the x values are in order
   call check_order_omega(file_in(ifile),dim_in(ifile)) 
end do

! determine max dimension of largest input file
  max_dim_in = 0
  do i = 1, n_inputs
     if(max_dim_in.lt.dim_in(i)) then
        max_dim_in = dim_in(i) 
     end if
  end do
  if(max_dim_in.eq.0) STOP 'Problem max_dim_in' 

  ALLOCATE(w_in(n_inputs,max_dim_in))
  ALLOCATE(y_in(n_inputs,max_dim_in))

! read the data from the input files 
do ifile = 1, n_inputs
   call  read_file(file_in(ifile),dim_in(ifile),w_in(ifile,:),y_in(ifile,:))
!   write(*,*) file_in(ifile),'dim_in(',ifile,')=',dim_in(ifile)
end do ! ifile


! debug: echos the input files
debug_file = ''
if(debug) then
   do ifile=1,n_inputs
      write(debug_file,'(a,i5.5)')'debug_input_w_y_',ifile
      open(50,file=debug_file)
      do i=1,    max_dim_in
         write(50,*) w_in(ifile,i),y_in(ifile,i)
      end do
      close(50)
   end do
end if ! debug

! calculate distribution function dist_in (for distribution_in)
! dist_in(ifile,:) lives on the same grid w_in as the input
! NOTE 
! I do not use any sophisticated integration as there is no more information
! in the data anyway. We integrate naively and only then interpolate
! to get the functions all on the same grid in the old y direction, that is, 
! in the new x direction

ALLOCATE( dist_in(n_inputs,max_dim_in) )

dist_in = 0.d0
do ifile = 1, n_inputs
   dist_in(ifile,1) = 0 !sum_rule_func(w_in(ifile,1),y_in(ifile,1),q_in(ifile))
   do i = 2,dim_in(ifile)
      dist_in(ifile,i) = dist_in(ifile,i-1) + sum_rule_func( w_in(ifile,i),y_in(ifile,i),q_in(ifile) ) *  (w_in(ifile,i) - w_in(ifile,i-1))   &
           - ( sum_rule_func( w_in(ifile,i),y_in(ifile,i),q_in(ifile)) - sum_rule_func( w_in(ifile,i-1),y_in(ifile,i-1),q_in(ifile)) ) *  (w_in(ifile,i) - w_in(ifile,i-1)) * 0.5d0
   end do ! i
end do ! ifile

! debug:  write out the integrated \int f(x_in,y_in) dx, where f(.,.) is the
! function which creates the sum-rule integral
debug_file = ''
if(debug) then
   do ifile=1,n_inputs
      write(debug_file,'(a,i5.5)')'debug_integrated_w_dist_in_',ifile
      open(50,file=debug_file)
      do i=1,    max_dim_in
         write(50,*) w_in(ifile,i),dist_in(ifile,i)
      end do
      close(50)
   end do
end if ! debug


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Exchange x and y (in the sense of inverting the integral);
! the new functions live now in the "inverted" space -> X,Y

ALLOCATE(X_inv_int_in(n_inputs,max_dim_in))
ALLOCATE(Y_inv_int_in(n_inputs,max_dim_in))

Y_inv_int_in = w_in
X_inv_int_in = dist_in


! debug:  write out inverted integrals of the inputs
debug_file = ''
if(debug) then
   do ifile=1,n_inputs
      write(debug_file,'(a,i5.5)')'debug_inverted_integrals_X_inv_int_in_Y_inv_int_in_',ifile
      open(50,file=debug_file)
      do i=1,    max_dim_in
         write(50,*) X_inv_int_in(ifile,i),Y_inv_int_in(ifile,i)
      end do
      close(50)
   end do
end if ! debug



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! determine the number of points :
! 1) for the allocation of the arrays: Longest set / STEPSIZE
! 2) for interpolation between the interpolated
!    curves -- that is: min length in NEW X dir / STEPSIZE
!    NB: We discard the rest of the curves which are longer

max_X = 0.d0 ! for longest set -- to have sufficient allocated size
min_X = 1000000.d0 ! for shortest set -- to have the region where [...X...]
                   ! where all curves have values

do ifile=1,n_inputs

   if(max_X.lt.X_inv_int_in(ifile,dim_in(ifile))) then
      max_X = X_inv_int_in(ifile,dim_in(ifile))
      ! I insert +1 in order to always be on the safe side (???)
      n_max_points_X = int(max_X / STEPSIZE)  +1
   end if

   if(min_X.gt.X_inv_int_in(ifile,dim_in(ifile))) then
      min_X = X_inv_int_in(ifile,dim_in(ifile))
      n_min_points_X = int(min_X / STEPSIZE) +1 ! - 1  
   end if

end do

write(*,*) "n_max_points_X = ", trim(char_write(n_max_points_X))  , ";   max_X = ", max_X
write(*,*) "n_min_points_X = ", trim(char_write(n_min_points_X))  , ";   min_X = ", min_X


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate the inverted functions to have the same (inverse-) x grid,
! corresponding to the y grid in the original \int f(x_in,y_in) dx

ALLOCATE(X_inv_interp(n_inputs,n_max_points_X))
ALLOCATE(Y_inv_interp(n_inputs,n_max_points_X))
X_inv_interp = 0.d0
Y_inv_interp = 0.d0
! introduce arrays to pass to C routine:
! input only !!!

ALLOCATE(Y_test(max_dim_in))
ALLOCATE(X_test(max_dim_in))

Y_test = 0.d0
X_test = 0.d0
! output
ALLOCATE(X_interp_test(n_max_points_X))! introduces a dummy array for I don't
ALLOCATE(Y_interp_test(n_max_points_X))! want to (cannot???) call X(n,:) to C 
                                       ! (francesco advised against it)

!!!! We have to use now an absolute step size!
do i = 1,n_inputs

   X_interp_test = 0.d0       ! this is necessary, for 
   Y_interp_test = 0.d0       ! the lengths of Y_inv_int_in(i,:) are different
   
   Y_test = Y_inv_int_in(i,:)
   X_test = X_inv_int_in(i,:)
   
   call wrapper_int(Y_test, X_test, X_interp_test, Y_interp_test, dim_in(i), stepsize ) 

   
! Test
!   do j=1,max_X/stepsize +1 !                n_max_points_X
!      write(70+i,*)  X_interp_test(j),Y_interp_test(j)
!   end do
   
! recover the variable we want from the dummy array
   X_inv_interp(i,:)= X_interp_test
   Y_inv_interp(i,:)= Y_interp_test

end do ! i = 1,n_inputs


! debug: echos the inverted integrals belonging to the inputs
! after the interpolation that gives them the same (inverse-) x grid,
! corresponding to the y grid in the original \int f(x_in,y_in) dx
if(debug) then
debug_file = ''
   do ifile=1,n_inputs
      write(debug_file,'(a,i5.5)')'debug_X_inv_interp_Y_inv_interp_',ifile

      open(50,file=debug_file)
      do i=1,max_X/stepsize +1 ! ==  n_max_points_X
         write(50,*) X_inv_interp(ifile,i), Y_inv_interp(ifile,i)
      end do
      close(50)

   end do ! ifile
end if ! debug



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This is only if we want to use our Taylor instead of the gsl; this 
! might be useful if we want to EXtrapolate 
if(TAYLOR.eq..TRUE.) then
!!!
! calculate lin and quadr. coefficients for dependence on q:
ALLOCATE(lincoef(n_min_points_X))
ALLOCATE(quadcoef(n_min_points_X))

call coeffs(n_inputs,n_min_points_X,n_max_points_X,Y_inv_interp(:,:),&
     q_in,lincoef,quadcoef)

! test coeffs:
! 2007 looks reasonable, but need to check the sign of the quad. coefficients
! sign should now be OK
!do i = 1, n_min_points_X 
!!   write(85,*) X_inv_interp(1,i),lincoef(i),quadcoef(i)
!!   write(85,*) X_inv_interp(1,10),lincoef(10),quadcoef(10)
!end do
! following line to check a single point
!   write(85,*) X_inv_interp(1,10),(Y_inv_interp(i,10),i=1,3)

end if ! TAYLOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate new function for any q given explicitly
! for now, start from the middle data set

min_out = 10000000 ! just a large number which hopefully always > length sets

ALLOCATE( f_out_inv (nq_out, n_min_points_X) )

ALLOCATE( y_int_out (nq_out, n_min_points_X) )
ALLOCATE( x_int_out (nq_out, n_min_points_X) )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! use gsl interpolation instead of Taylor:
if(TAYLOR.eq..FALSE.) then

write(*,*) ''
write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
write(*,*) '++++++++++++++++ gsl interpolation +++++++++++++'
write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
write(*,*) ''

! the -.00001 at the end is to get the last point ... otherwise 
! it falls out of the interpolation 
q_step = (q_out_max - q_out_min) / float(nq_out - 1) - .000001

write(*,*) 'q_step =',q_step

ALLOCATE( Aq (n_inputs) )
ALLOCATE( Xq (n_inputs) )
ALLOCATE( Aq_interp (nq_out) )
ALLOCATE( Xq_interp (nq_out) )

Aq_interp = 0.d0
Xq_interp = 0.d0

do i_point = 2, n_min_points_X 
   
   Aq = Y_inv_interp(:,i_point)
   Xq =  q_in
   
!write(*,*) 'Aq=',Aq,'i_point=',i_point
!write(*,*) 'Xq=',Xq,'i_point=',i_point
   
! choose interpolation type
   select case (interpol_type)
   case(1) ! write(*,*) 'linear interpolation, interpol_type = 1'
      call wrapper_int_fq_lin(Aq, Xq, Xq_interp, Aq_interp, n_inputs, q_step ) 
   case(2) ! write(*,*) 'polynomial interpolation, interpol_type = 2'
      call wrapper_int_fq_pol(Aq, Xq, Xq_interp, Aq_interp, n_inputs, q_step ) 
   case(3) ! write(*,*) 'cubic spline interpolation, interpol_type = 3'
      call wrapper_int_fq_cspline(Aq, Xq, Xq_interp, Aq_interp, n_inputs,q_step) 
   case(4) ! write(*,*) 'akima spline interpolation, interpol_type 4'
      call wrapper_int_fq_akima(Aq, Xq, Xq_interp, Aq_interp, n_inputs, q_step )
    case default
      STOP 'STOP -- invalid interpolation type' 
   end select


   f_out_inv(:,i_point) =   Aq_interp
   
end do ! i_point = 1,n_min_points_X

write(*,*) 'nq_out =',nq_out

! Test interpolated inverse integrals (those which will become the outputs)
if(debug) then
debug_file = ''
   do j=1,nq_out
      write(debug_file,'(a,i5.5)')'debug_X_inv_out_Y_inv_out_',j
      open(50,file=debug_file)
      do i=1, n_min_points_X 
         write(50,*) X_inv_interp(1,i),f_out_inv(j,i), j
      end do ! i
      close(50)
   end do ! j
end if ! debug


end if ! switch TAYLOR/gsl
! end gsl interpolation instead of Taylor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do iq = 1, nq_out
   
   write(*,*) 'q_out( ',trim(char_write(iq)),' ) = ', q_out(iq)
   
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the following is the original Taylor-like interpolation that works
! I keep it in case we want to put later extrapolation
if(TAYLOR.eq..TRUE.) then

write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
write(*,*) '++++++++++++++++TAYLOR interpolation ++++++++++++++'
write(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++'

call  interpol_curves(n_inputs,n_min_points_X, &
           n_max_points_X, & !!! this is new!
           Y_inv_interp(:,:), &
           lincoef,quadcoef, &        
           q_in, &
           q_out(iq),&
           f_out_inv(iq,:))
end if ! TAYLOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! exchange again X and Y, to get the function to interpolate in 
!!! the original space -- only then do we interpolate / differenciate

!!! the X values are now in POINTS, the "units" come back in when
!!! we re-invert using the inverted X mesh

!!! here the X inv array as fkt of i
!!! This was wrong: X_inv_interp is the X of the interpolated input curves,
!!! Thus they should all be the same and we use simply X_inv_interp(1,.)
y_int_out(iq,:) = X_inv_interp(1,[1:n_min_points_X]) 
x_int_out(iq,:) = f_out_inv(iq,[1:n_min_points_X])

!!! to obtain the max length of the output data sets:
!!! use i_max_out = max_out / stepsize_result
if (max_out.lt.x_int_out(iq,n_min_points_X)) then
   max_out = x_int_out(iq,n_min_points_X)
end if

!!! to obtain the min length of the output data sets:
!!! i_min_out = min_out / stepsize_result
if (min_out.gt.x_int_out(iq,n_min_points_X)) then
   min_out = x_int_out(iq,n_min_points_X)
end if

end do ! iq


! debug: interpolated returned integrals
if(debug) then
debug_file = ''
   do iq=1,nq_out
      write(debug_file,'(a,i5.5)')'debug_interpolated_integrals_out_',iq
      open(50,file=debug_file)
      do i=1, n_min_points_X 
         write(50,*)  x_int_out(iq,i),y_int_out(iq,i)
      end do ! i
      close(50)
   end do ! iq
end if ! debug


!!! Determine length (# of data points) of output
write(*,*) ''
write(*,*) 'Characteristics of ouput spectra:'
imax_out = max_out / stepsize_result

write(*,*) 'imax_out = ', imax_out
write(*,'(a,f12.6)') ' max_out = ', max_out
write(*,'(a,f12.6)') ' stepsize_result = ', stepsize_result

imin_out = min_out / stepsize_result -1 ! to be on the safe side

write(*,*) 'imin_out = ',imin_out 
write(*,'(a,f12.6)') ' min_out = ', min_out 
write(*,'(a,f12.6)') ' stepsize_result = ', stepsize_result 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! obtain the desired spectra by differentiating and applying 
!! sum_rule_func ^ -1

ALLOCATE(f_out(nq_out,imax_out))
f_out = 0.d0
! the +1 in the following allocation avoids a crash in the deriv. (in which
! the program first finishes the loop but then haengs up -- probably a problem
! of the connection between c and fortran
ALLOCATE( deriv_x(imax_out+1) ) ! we leave deriv_*() local 
ALLOCATE( deriv_y(imax_out+1) ) ! in their loop over iq

ALLOCATE( x_int_out_pass(n_min_points_X) )
ALLOCATE( y_int_out_pass(n_min_points_X) )

do iq = 1, nq_out

x_int_out_pass = 0.d0
y_int_out_pass = 0.d0

deriv_x = 0.d0
deriv_y = 0.d0
y_int_out_pass = y_int_out(iq,:)
x_int_out_pass = x_int_out(iq,:)

! 2007 Test OK
!   do i = 1,  n_min_points_X
!      write(50+iq,*) x_int_out_pass(i),y_int_out_pass(i)
!      write(60+iq,*) x_int_out(iq,i),y_int_out(iq,i)
!   end do

call wrapper_differentiate(y_int_out_pass,x_int_out_pass,&
     deriv_x,deriv_y, &
     n_min_points_X, & ! this is # of points really read-in
     stepsize_result)  !

! apply the inverse sum_rule_func which should include the scaling
! for the output values of q if necessary
do i = 2,imin_out ! "2" to suppress the border effect at zero
   f_out(iq,i) = sum_rule_func_inv(deriv_x(i),deriv_y(i),q_out(iq))
   ! the result is now in: deriv_x(i),  f_out(iq,i)
end do

end do ! iq

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finished! Now write the curves to files;
! Attention: we use here the value of deriv_x(:) of the last loop over iq 
! -- they should all be the same for each q

!! file name creation : output to files
write(*,*) ''
write(*,*) 'Output files:'
do iq = 1,nq_out

   write(name_part,'(f10.5)') q_out(iq)
   long = trim(name_base) // trim(name_part)
   
   do i = 1,80
      if (name_base(i:i) == ' ') then
         length_name_base = i-1
         exit
      end if
   end do
! filename_compact=.TRUE. is the standard way for file naming   
! filename_compact=.FALSE. when 'outfile_sortable' is given in the 
! input file -> results in a bit more clumsy filenames which can, however,
! more easily be sorted numerically (sort -n)
   if(filename_compact) then
      i_pos = 1
      do j=1,80
         if (long(j:j)/=' ') then
            outfile(i_pos:i_pos) = long(j:j) 
            i_pos = i_pos + 1         
         end if
      end do
   else 
      do j=1,80
         if (long(j:j)/=' ') then
            outfile(j:j) = long(j:j)
         else if (long(j:j)==' '.AND.j<length_name_base+10) then
            outfile(j:j) = '0'
         end if
      end do
   end if

   write(*,*) trim(outfile)
   open(88,file=outfile)
   
   do i = 10,imin_out
      write(88,*)  deriv_x(i),f_out(iq,i) 
   end do
   
   close(88)
   
end do ! iq over output values

write(*,*) ''
write(*,*) 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
write(*,*) 'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX***********************'
write(*,*) 'XXXXXXXXXX                                          XXXX'
write(*,*) 'XXXXXX          Interspectus ended OK           XXXXXXXX'
write(*,*) 'XXXX                                      XXXXXXXXXXXXXX'
write(*,*) 'XXX           ________________xxxxxxxxxXXXXXXXXXXXXXXXXX'


end program interspectus


