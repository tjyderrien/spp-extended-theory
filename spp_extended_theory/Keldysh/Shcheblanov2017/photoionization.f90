    module constants  
      implicit none 
      real*8, parameter :: sr2 = dsqrt(2.0D0)      ! square root from two(=2)
      real*8, parameter :: pi = 4.0D0*datan(1.0D0) ! pi number = 3.1415...
      real*8, parameter :: c  = 2.998D+008         ! Light speed in [m/sec]
      real*8, parameter :: h  = 0.659D-015         ! Planck constant (h_bare) in [eV*sec]
      real*8, parameter :: hJ = 1.0546D-034        ! Planck constant (h_bare) in [J*sec]
      real*8, parameter :: eps0 = 8.854D-012       ! Electric permittivity of vacuum in [F/m]
      real*8, parameter :: Qe = 1.602D-019         ! Elementary charge in [Q]
      real*8, parameter :: m0 = 0.911D-030         ! Free e- mass in [kg]
      real*8, parameter :: mv = 1.0d0              ! 1.126D0,    ! Valence band e- mass in [m0]
      real*8, parameter :: mc = 1.0d0              ! 0.543D0,   ! Conduction band e- mass in [m0]
      real*8, parameter :: mvb = mv*m0             ! Valence band e- mass in [kg]
      real*8, parameter :: mcb = mc*m0             ! Conduction band e- mass in [kg]
      real*8, parameter :: mr = m0*0.9d0           ! (mv*mc)/(mv+mc) ! Reduced e- mass in [kg]
      real*8, parameter :: eJ = 1.602D-019         ! 1 [eV] = 1.6D-019 [J] 
      real*8, parameter :: gp = 9.0D0              ! Gap energy (SiO2) in [eV]
      real*8, parameter :: gpJ = gp*eJ             ! Gap energy (SiO2) in [J]

    end module constants  

    Program Photoionization
      use constants
      implicit none
        
      real*8, parameter :: nm_to_m = 1.0d-9
      real*8, parameter :: Eo = 1.0d8             ! Field step = 1 MV/cm in [V/m]
	  
	  real*8 :: lambda,omega,n_ref
	  real*8 :: rpi(4),lpi,Epi,egp,sig,Intensity,Eb
	  integer :: i 
 
      lambda = 800 ! laser wavelength in [nm]  
      omega  = 2.0d0*pi*c/(lambda*nm_to_m) ! laser frequency in [1/sec]
      n_ref  = 1.5916 ! refractive index of SiO2
  
      open(111,file='SiO2_PRA_96_2017.dat',status='unknown')
	  !write(111,'(5E16.8)')'E[V/m]  Intensity[TW/cm^2]  Keldysh[1/(sec cm^3)] corr_Keldysh[1/(sec cm^3)] Parabolic[1/(sec cm^3)]'
      do i=0,5000  
        Eb = Eo*i/10.0d0 ! incident electric field tension [V/m]
        if(Eb < Eo*0.1d0)Eb = Eo*0.1d0 
        Intensity=Eb**2*eps0*n_ref*c/2.d0 ! intensity inside material [W/m^2]
        
        call Krate(omega,Eb,rpi(1),lpi,Epi,egp,sig)      ! Keldysh model, Kane 
        call Krate_sin(omega,Eb,rpi(2),lpi,Epi,egp,sig)  ! Keldysh modified, Kane
        call Prate(omega,Eb,rpi(3),lpi,Epi,egp)          ! Parabolic
		call Prate_sin(omega,Eb,rpi(4),lpi,Epi,egp,sig)          ! Parabolic
	    
        write(111,'(6E16.8)')Eb,Intensity/1d16,rpi(:)/1.0d6
      enddo
        
    endprogram Photoionization
	
	
	
    SUBROUTINE Krate(w,Efi,wpi,lpi,Epi,egp,sig) ! Kane Band shape (Keldysh)
	  use constants
      implicit none
      real*8, INTENT(IN) :: Efi,w
      real*8, INTENT(OUT) :: wpi,lpi,Epi,egp,sig
        
      real*8 :: wpil,hw,kld,sqkld,gm1,gm2,sqgm1,sqgm2,fac,st,eel1,eel2,kel1,kel2,z,zz,tmp,ksi,eta,x,F,Di
	  integer :: l,n
     
      !!!!!!!!End Constants!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      hw = h*w         ! Photon energy in [eV]
      
      kld = w*dsqrt(gpJ*mr)/(Efi*Qe) ! Keldysh parameter (gamma)
      sqkld = kld*kld
      gm1 = 1.0D0/dsqrt(1.0D0+sqkld)
      gm2 = kld*gm1
      sqgm1 = gm1*gm1
      sqgm2 = gm2*gm2
      fac = (2.0D0*w/(9.0D0*pi))*((mr*w/(hJ*gm2))**(3.0D0/2.0D0)) ! pre-factor in final expression
      
      !!!! E(gm1),E(gm2),K(gm1),K(gm2) !!!!
      st = 0.5D0*pi/1000.0D0 ! integration step
      eel1 = 0.0D0
      eel2 = 0.0D0
      kel1 = 0.0D0
      kel2 = 0.0D0
      do l = 1,999
        z = dsin(st*l)
        zz = z*z
        eel1 = eel1 + dsqrt(1.0D0-sqgm1*zz)
        eel2 = eel2 + dsqrt(1.0D0-sqgm2*zz)
        kel1 = kel1 + 1.0D0/dsqrt(1.0D0-sqgm1*zz)
        kel2 = kel2 + 1.0D0/dsqrt(1.0D0-sqgm2*zz)
      enddo
      tmp = st/2.0D0
      eel1 = tmp*(1.0D0 + 2.0D0*eel1 + dsqrt(1.0D0-sqgm1))
      eel2 = tmp*(1.0D0 + 2.0D0*eel2 + dsqrt(1.0D0-sqgm2))
      kel1 = tmp*(1.0D0 + 2.0D0*kel1 + 1.0D0/dsqrt(1.0D0-sqgm1))
      kel2 = tmp*(1.0D0 + 2.0D0*kel2 + 1.0D0/dsqrt(1.0D0-sqgm2))
      !!!! E(gm1),E(gm2),K(gm1),K(gm2) !!!!
      
      ksi = (kel2-eel2)/eel1
      eta = 2.0D0*kel1*eel1
      
      egp = 2.0D0*gp*eel1/(pi*gm2)     ! Eff. gap (any gamma) in [eV]
      x = egp/hw
      lpi = IDINT(1.0D0 + x)           ! MIN # of photons (any gamma)
      Epi = (lpi*hw-gp)*(mr/mcb)       ! Kinetic energy due to PI in [eV]
      
      wpi = 0.0d0
      sig = 0.0d0 
      do n = 0,30 ! nMAX=3
        ! Dawson's integral calculation (Trapezium method)
        z = dsqrt(pi*pi*(lpi-x+n)/eta)
        zz = z*z
        st = z/500.0D0 ! integration step
        tmp = st*st
        F = 0.0D0
        do l = 1,499 ! Di(z) - Dawson's integral
          F = F + dexp(l*l*tmp)
        enddo
        Di = dexp(-zz)*(st/2.0D0)*(1.0D0+2.0D0*F+dexp(zz))
        ! End Dawson's integral calculation
        wpil = 2.0d0*fac*dexp(-pi*lpi*ksi)*dsqrt(0.5D0*pi/kel1)*Di*dexp(-pi*n*ksi) !2 for spin
        !if(lpi+n .ne. 2) wpil=0.0d0
        wpi = wpi+wpil
        
        sig = sig + 2.0d0*wpil*(lpi+n)*hJ*w/(eps0*Efi*Efi) ! equation (39) in description file
      enddo

      
      END SUBROUTINE Krate
    
      
      
      SUBROUTINE Krate_sin(w,Efi,wpi,lpi,Epi,egp,sig) ! Kane Band shape (Keldysh)
      use constants
	  implicit none
      real*8, INTENT(IN) :: Efi,w
      real*8, INTENT(OUT) :: wpi,lpi,Epi,egp,sig

      integer :: l,n
      real*8 :: cns36,wpil,hw,kld,sqkld,gm1,gm2,sqgm1,sqgm2,fac,st,eel1,eel2,kel1,kel2,z,zz,tmp,ksi,eta,x,F,Di
      
      !!!!!!!!End Constants!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      hw = h*w         ! Photon energy in [eV]
      
      kld = w*dsqrt(gpJ*mr)/(Efi*Qe) ! Keldysh parameter (gamma)
      sqkld = kld*kld
      gm1 = 1.0D0/dsqrt(1.0D0+sqkld)
      gm2 = kld*gm1
      sqgm1 = gm1*gm1
      sqgm2 = gm2*gm2
      fac = (4.0D0*w/(9.0D0*pi))*((mr*w/(hJ*gm2))**(3.0D0/2.0D0)) ! pre-factor in final expression
      
      !!!! E(gm1),E(gm2),K(gm1),K(gm2) !!!!
      st = 0.5D0*pi/1000.0D0 ! integration step
      eel1 = 0.0D0
      eel2 = 0.0D0
      kel1 = 0.0D0
      kel2 = 0.0D0
      do l = 1,999
        z = dsin(st*l)
        zz = z*z
        eel1 = eel1 + dsqrt(1.0D0-sqgm1*zz)
        eel2 = eel2 + dsqrt(1.0D0-sqgm2*zz)
        kel1 = kel1 + 1.0D0/dsqrt(1.0D0-sqgm1*zz)
        kel2 = kel2 + 1.0D0/dsqrt(1.0D0-sqgm2*zz)
      enddo
      tmp = st/2.0D0
      eel1 = tmp*(1.0D0 + 2.0D0*eel1 + dsqrt(1.0D0-sqgm1))
      eel2 = tmp*(1.0D0 + 2.0D0*eel2 + dsqrt(1.0D0-sqgm2))
      kel1 = tmp*(1.0D0 + 2.0D0*kel1 + 1.0D0/dsqrt(1.0D0-sqgm1))
      kel2 = tmp*(1.0D0 + 2.0D0*kel2 + 1.0D0/dsqrt(1.0D0-sqgm2))
      cns36 = 0.5d0*dsqrt(2.0d0*gp*kel1/(pi*hw*gm2))*(0.5d0*pi-2.0d0*datan(kld))
      
      ksi = (kel2-eel2)/eel1
      eta = 2.0D0*kel1*eel1
      
      egp = 2.0D0*gp*eel1/(pi*gm2)     ! Eff. gap (any gamma) in [eV]
      x = egp/hw
      lpi = IDINT(1.0D0 + x)           ! MIN # of photons (any gamma)
      Epi = (lpi*hw-gp)*(mr/mcb)       ! Kinetic energy due to PI in [eV]
      
      wpi = 0.0d0
      sig = 0.0d0 
      do n = 0,30 ! nMAX=3
        ! Dawson's integral calculation (Trapezium method)
        z = dsqrt(pi*pi*(lpi-x+n)/eta)
        zz = z*z
        st = z/500.0D0 ! integration step
        tmp = st*st
        F = 0.0D0
        do l = 1,499 ! Di(z) - Dawson's integral
          F = F + dexp(l*l*tmp)*(dsin(0.5d0*pi*(lpi+n)+cns36*l*st))**2
        enddo
        Di = dexp(-zz)*(st/2.0D0)*((dsin(0.5d0*pi*(n+lpi)))**2+2.0D0*F+dexp(zz)*(dsin(0.5d0*pi*(lpi+n)+cns36*z))**2)
        ! End Dawson's integral calculation
         wpil = 2.0d0*fac*dexp(-pi*lpi*ksi)*dsqrt(0.5D0*pi/kel1)*Di*dexp(-pi*n*ksi) ! 2 from spin and 2 was lost in prefactor
         !if(lpi+n .ne. 2.0d0) wpil = 0.0d0
         wpi = wpi+wpil
         sig = sig + 2.0d0*wpil*(lpi+n)*hJ*w/(eps0*Efi*Efi) ! equation (39) in description file
        enddo
    !    
        Epi = kld 
      
    !  sig = 2.0d0*wpi*gpJ/(eps0*Efi*Efi)

      
      END SUBROUTINE Krate_sin
      
      SUBROUTINE Prate(w,Efi,wpi,lpi,Epi,egp) ! Parabolic Band shape
      use constants
	  implicit none

      real*8, INTENT(IN) :: Efi,w
      real*8, INTENT(OUT) :: wpi,lpi,Epi,egp
	  
	  integer :: l,n
      real*8 :: wpil,hw,kld,sqkld,fac,st,z,zz,tmp,ksi,eta,x,F,Di,ars,QF
 
      hw = h*w          ! Photon energy in [eV]
      
      kld = w*dsqrt(gpJ*mr)/(Efi*Qe)    ! Keldysh parameter (gamma)
      sqkld = kld*kld
      egp = gp*(1.0D0 + 0.25D0/sqkld)    ! Eff. gap (any gamma) in [eV]
      x = egp/hw
      lpi = IDINT(1.0D0 + x)             ! MIN # of photons (any gamma)
      Epi = (lpi*hw-gp)/2.0D0            ! Kinetic energy due to PI in [eV]
      
      fac = (w/(4.0D0*pi))*((mr*w/hJ)**(3.0D0/2.0D0))  ! pre-factor in final expression
      
      eta = 1.0D0/dsqrt(1.0D0+0.5D0/sqkld)
      ksi = sqkld/(1.0D0+4.0D0*sqkld)
      ars = DLOG(sr2*kld+dsqrt(1.0D0+2.0D0*sqkld)) ! Arsh(sqrt(2.0d0)*gamma)
      
      QF = 0.0D0
      do n = 0,30 ! nMAX = 3
        ! Dawson's integral calculation (Trapezium method)
        z = dsqrt(2.0D0*eta*(lpi-x+n))
        zz = z*z
        st = z/1000.0D0 ! integration step
        tmp=st*st
        F = 0.0D0
        do l = 1,999
          F = F + dexp(l*l*tmp)
        enddo
        F = (st*0.5D0)*(1.0D0+2.0D0*F+dexp(zz))
        Di = dexp(-zz)*F ! Di(z) - Dawson's integral
        ! End Dawson's integral calculation
        QF = QF + Di*dexp(-2.0D0*n*(ars-eta))
      enddo
      wpi = fac*dsqrt(1.0D0/eta)*QF*dexp(-2.0D0*lpi*(ars-eta)-4.0D0*x*eta*ksi)
      
      END SUBROUTINE Prate
 
    
      SUBROUTINE Prate_sin(w,Efi,wpi,lpi,Epi,egp,sig) ! Parabolic Band shape
       use constants
 	   implicit none
           
       real*8 :: a19,wpil,alfa, beta
 
       real*8, INTENT(IN) :: Efi,w
       real*8, INTENT(OUT) :: wpi,lpi,Epi,egp,sig
	   
	   integer :: l,n
       real*8 :: hw,kld,sqkld,fac,st,z,zz,tmp,eta,x,F,Di,ars,qf
  
       hw = h*w          ! Photon energy in [eV]
       
       kld = w*dsqrt(gpJ*mr)/(Efi*Qe)    ! Keldysh parameter (gamma)
       sqkld = kld*kld
       egp = gp*(1.0D0 + 0.25D0/sqkld)    ! Effective gap (any gamma) in [eV]
       x = egp/hw
       lpi = IDINT(1.0D0 + x)             ! MIN # of photons (any gamma)
       Epi = (lpi*hw-gp)/2.0D0            ! Kinetic energy due to PI in [eV]
       
       fac = w/(4.0D0*pi)*(mr*w/hJ)**1.5D0  ! pre-factor in final expression
       ars = DLOG(sr2*kld+dsqrt(1.0D0+2.0D0*sqkld)) ! sinh^-1(sr2*kld)
       
       beta = 2.0d0/dsqrt(1.0d0 + 0.5d0/sqkld)
       alfa = 2.0d0*ars-beta
       eta  = 2.0D0*beta*sqkld/(1.0D0+4.0d0*sqkld)
         
       a19 = dsqrt(2.0d0*gp/(hw*beta))*(2.0d0*sr2/beta) ! eq(19) for a in parabolic
    
       QF = 0.0D0
       wpi =0.0d0
       sig=0.0d0
       do n = 0,30 ! nMAX = 3
         ! Dawson's integral calculation (Trapezium method)
         z = dsqrt(beta*(lpi-x+n))
         zz = z*z
         st = z/1000.0D0 ! integration step
         tmp=st*st
         F = 0.0D0
         do l = 1,999
           F = F + dexp(l*l*tmp)*(dsin(0.5d0*pi*(lpi+n)+a19*l*st))**2
         enddo
         Di = dexp(-zz)*(st*0.5D0)*((dsin(0.5d0*pi*(lpi+n)))**2+2.0D0*F+dexp(zz)*(dsin(0.5d0*pi*(lpi+n)+a19*z))**2)
         ! End Dawson's integral calculation
 
         wpil = 2.0d0*fac*dsqrt(2.0D0/beta)*Di*dexp(-alfa*n)*dexp(-alfa*lpi-eta*x) ! 2 for spin
  !       if(lpi+n .ne. 2.0d0) wpil = 0.0d0
         wpi = wpi+wpil
         sig = sig + 2.0d0*wpil*(lpi+n)*hJ*w/(eps0*Efi*Efi)
       enddo
       
       END SUBROUTINE Prate_sin

     