C **************************************************************
c     program name: dist_DF_zhou14.f 
C
C     ==========================================================
C     Calculate the ion distributions and moments at a virtual 
C     spacecraft location with a dipolarization front moving 
C     towards the spacecraft. See Zhou et al (JGR, 2014) 
C     
C **************************************************************

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION Be(3),Beini(3),vstart(6),vend(6)
	DOUBLE PRECISION PSD(80),EFLUX(80),PSDT(50,50,50)
	CHARACTER *1 NNENG
	CHARACTER *2 K3s, K4s
	CHARACTER *3 SECs
	INTEGER qmsign,fb

	COMMON /AB/Blobe,bz,HTx,XBD,ALPHA,PHII,Xm,Bm,Wm,IE
	COMMON /basics/Re,qm,PI,rad,v0,E0,B00,vc,Bio,qmsign
	COMMON /red/rp,ep,dpp
	COMMON /calcu/fb
	COMMON /outp/joutp,K3s,K4s,iof
	COMMON /DF/Bz0,Vw,tMAX,Xb,X0,tSEC,Y0,HF

! Basic parameters
	qe=1.60217646d-19  ! elementary charge (coulombs)
	pm=1.67262158d-27  ! proton mass (kg) 
	qm=qe/pm
	vc=3.d8            ! speed of light (m/s)
	Re=6.3712d6        ! Earth Radius (m)
      PI=3.14159265359D0
	Bio=5.6d6          ! B FIELD STRENGTH AT IONOSPHERE (LATITUDE=68 degree)

! Input Key Parameters
!	OPEN(UNIT=2,FILE='KeyPara.txt')
!	READ(2,*) Bz,Vth,Vw
!	CLOSE(2)

! Current sheet equilibrium parameters	

	XBD=-10.d0        ! X Location where b.c. are determined (Re)
	Blobe=30.d0		  ! Magnetic field in the Lobe at X=XBD (nT)
	HTx=0.5d0	      ! Current sheet Half-Thickness at X=XBD (Re)
	Bzxbd=2.d0        ! Magnetic Field Bz component at X=XBD, z=0 (nT)
	bz=2.d0 	      ! Magnetic field Bz component at X=-inf, z=0 (nT)

	IE=4			  ! THE EXPONENTIAL DECAY FACTOR OF Bz(x,0)
      Xm=-20.d0         ! X location of the Bz secondary peak (Re)
	Bm=0.d0			  ! The Bz value of the secondary peak (nT)
	Wm=2.d0           ! The width of the Bz secondary peak (Re)

	ALPHA=(Bzxbd-bz-Wm*Wm*Bm/((XBD-Xm)**2+Wm*Wm))
	1                     *XBD**IE/(Blobe*HTx**IE)
	PHII=HTx**(IE-1)*ALPHA/((IE-1.d0)*XBD**(IE-1))-bz*XBD/(Blobe*HTx)
	1          -Wm*Bm*datan((XBD-Xm)/Wm)/HTx/Blobe	

!  SO THAT Bz(x,0) = HTx^IE*Blobe*ALPHA/X^IE+bz

	Vth=700.d0		  ! Proton thermal velocity (km/s)
	Vsh=Vth*Vth/(HTx*Re*qm*Blobe/1.d3)*1.d9 ! Proton bulk velocity (km/s)
	DEN=0.35d0		  ! Plasma density at X=XBD, z=0 (cm^-3)
	DENf=DEN/(DEXP(Vsh**2/Vth**2)*PI**1.5*Vth**3)

	Vth2=400.d0		  ! Proton thermal velocity of background population
	Vsh2=0.d0		  ! no bulk velocity for background population
	DEN2=0.05d0		  
	DEN2f=DEN2/(PI**1.5*Vth2**3)

! Dipolarization front parameters
	Bz0=10.d0    ! Magnetic field Bz enhancement associated with the DF (nT)
	Vw=200.d0   ! Earthward propagating speed of the DF (km/s) 
! ENDTEST
	Xb=-12.d0   ! The initial location of the DF (Re)
	X0=0.1d0	! The DF half-thickness (Re). Non-applicable in sharp DF.
	HF=1.d0		! THE DF half-width in Y (Re)
	Y0=0.d0		! THE DF center in Y (Re)

! Normalization & calculation factors 
      fb=-1        ! backward tracing
	v0=1.d5
	E0=1.d-3     !(V/m)
      B00=1.d-9     !(nT)
	rp=qm*E0/v0
      ep=v0*B00/E0
      dpp=v0/Re	
	qmsign=1     ! proton
	t0=0.d0

! Spacecraft location 

	xini=-10.d0
	write(*,*) 'Input the Y location of the virtual spacecraft (Re)'
	read(*,*) yini
	zini=0.0047d0  !30 km

	vstart(4)=xini
	vstart(5)=yini
	vstart(6)=zini

! The time when DF arrives (s)
	tMAX=(xini-Xb+HF-DSQRT(HF**2-(yini-Y0)**2))*Re/1.d3/Vw  

	CALL VPot(xini,zini,Aini,Beini) 
	Bti=DSQRT(Beini(1)**2+Beini(2)**2+Beini(3)**2)

! Calculating the ion azimuthal angular spectra
	N4=80
	BL44=2.d0*PI
	H4=BL44/N4

	K4min=1
	K4max=N4

	DO NENG=6,7
	
	IF (NENG.EQ.6) THEN
	 OPEN(UNIT=NENG,FILE='Spec_5k22k_XZ14_Z0_cirF_V200.txt') 
	ELSE
	 OPEN(UNIT=NENG,FILE='Spec_30k100k_XZ14_Z0_cirF_V200.txt')
	ENDIF

	IF (NENG.EQ.6) THEN
	  NVMIN=1000
	  NVMAX=1900
	ELSEIF (NENG.EQ.7) THEN
	  NVMIN=2450
	  NVMAX=4350
	ELSE
	  NVMIN=1700
	  NVMAX=4700
	ENDIF
	NVH=100

	NAMIN=-25
	NAMAX=25
	NAH=10

	DHT=DSIN(PI*(NAMAX+0.5d0*NAH)/180.d0)
	1    -DSIN(PI*(NAMIN-0.5d0*NAH)/180.d0)

	DO tSEC=tMAX-1.5d0,0.d0,-3.d0
	  DO K4=K4min,K4max
	    PSD(K4)=0.d0
	    EFLUX(K4)=0.d0
	    DO NVsst=NVMIN,NVMAX,NVH
	      Vsst=DBLE(NVsst)
	      DO Nangle=NAMIN,NAMAX,NAH
	        ANG=PI*DBLE(Nangle)/180.d0
		    DH=DSIN(PI*(Nangle+0.5d0*NAH)/180.d0)
	1           -DSIN(PI*(Nangle-0.5d0*NAH)/180.d0)
	        ROW=Vsst
	        FIS=H4*(K4-1.d0)
		    WX=-ROW*DCOS(FIS)*DCOS(ANG)
		    X3=-ROW*DSIN(FIS)*DCOS(ANG)
	        vstart(1)=WX*1.d3/v0
		    vstart(2)=X3*1.d3/v0
		    vstart(3)=ROW*DSIN(ANG)*1.d3/v0
		    Pini=vstart(2)*v0/1.d3+Aini*1.d-9*Re*1.d-3*qm

	        CALL rkdumb(t0,TSEC,vstart,vend)

	   	    IF (iof.EQ.1) THEN
	          PSD(K4)=0.d0
	        ELSE  

		      X2=vend(1)*v0/1.d3
		      X3=vend(2)*v0/1.d3
	          X4=vend(3)*v0/1.d3
		      xend=vend(4)
		      zend=vend(6)

	          CALL VPot(xend,zend,Aend,Be)  
	          Pend=X3+Aend*1.d-9*Re*1.d-3*qm

	          PSD(K4)=DENf*DEXP((-X2**2-X3**2-X4**2)/Vth**2)
	1                      *DEXP(2.d0*Vsh*Pend/Vth**2)
	1                  +DEN2f*DEXP((-X2**2-X3**2-X4**2)/Vth2**2)

              ENDIF

	        EfH=5.d4*PSD(K4)*Vsst**4*((Vsst+NVH/2)**2-(Vsst-NVH/2)**2)
	        Eflux(K4)=Eflux(K4)+EfH*DH/DHT
		  ENDDO
	    ENDDO
	    Eflux(K4)=Eflux(K4)/((NVMAX+NVH/2)**2-(NVMIN-NVH/2)**2)
	  ENDDO
	  write(NENG,1001) EFLUX
	ENDDO

	CLOSE(UNIT=NENG)
	ENDDO

! Calculating the plasma bulk velocity ahead of the front arrival
	OPEN(9,FILE='VPERP_XZ14_V200_Z0.txt')

	Kmax=40
	Vmax=4000.d0
	SK=DBLE(Kmax+1)/2.d0
	SKm=SK-1.d0
	dv3=(Vmax/SKm)**3.d0

	DO tSEC=tMAX,tMAX-30.d0,-1.d0
	  VxDEN=0.d0
	  VyDEN=0.d0
	  VzDEN=0.d0
	  DDEN=0.d0
	  PDEN=0.d0
	  DO K2=1,Kmax
	    DO K3=1,Kmax
	      DO K4=1,Kmax
		    vstart(1)=((K2-SK)/SKm)*Vmax*1.d3/v0
	        vstart(2)=((K3-SK)/SKm)*Vmax*1.d3/v0
	        vstart(3)=((K4-SK)/SKm)*Vmax*1.d3/v0			

	        CALL rkdumb(t0,TSEC,vstart,vend)

		    X2=vend(1)*v0/1.d3
		    X3=vend(2)*v0/1.d3
	        X4=vend(3)*v0/1.d3

		    xend=vend(4)
		    zend=vend(6)
			CALL VPot(xend,zend,Aend,Be)  
	        Pend=X3+Aend*1.d-9*Re*1.d-3*qm

	        PSDT(K2,K3,K4)=DENf*DEXP((-X2**2-X3**2-X4**2)/Vth**2)
	1                           *DEXP(2.d0*Vsh*Pend/Vth**2)
	1                      +DEN2f*DEXP((-X2**2-X3**2-X4**2)/Vth2**2)

	        VxDEN=VxDEN+vstart(1)*PSDT(K2,K3,K4)
	        VyDEN=VyDEN+vstart(2)*PSDT(K2,K3,K4)
	        VzDEN=VzDEN+vstart(3)*PSDT(K2,K3,K4)
	        DDEN=DDEN+PSDT(K2,K3,K4)
	      ENDDO
	    ENDDO
	  ENDDO
	  Vxx=VxDEN*v0/1.d3/DDEN
	  Vyy=VyDEN*v0/1.d3/DDEN
	  Vzz=VzDEN*v0/1.d3/DDEN
	  DDEN=DDEN*dv3

	  DO K2=1,Kmax
	    DO K3=1,Kmax
	      DO K4=1,Kmax
		    Vx0=((K2-SK)/SKm)*Vmax-Vxx
	        Vy0=((K3-SK)/SKm)*Vmax-Vyy
	        Vz0=((K4-SK)/SKm)*Vmax-Vzz
			PDEN=PDEN+pm*(Vx0*Vx0+Vy0*Vy0+Vz0*Vz0)*PSDT(K2,K3,K4)	        
	      ENDDO
	    ENDDO
        ENDDO
	  PDEN=PDEN*dv3*1.d12*1.d9/3.d0

	  Vpara=(Vxx*Beini(1)+Vyy*Beini(2)+Vzz*Beini(3))/Bti
	  Vperpx=Vxx-Vpara*Beini(1)/Bti
	  Vperpy=Vyy-Vpara*Beini(2)/Bti
	  Vperpz=Vzz-Vpara*Beini(3)/Bti
	  write (9,1002) tSEC, Vxx, Vyy, Vzz, DDEN, PDEN, 
	1                  Vperpx, Vperpy, Vperpz

	ENDDO

1001  FORMAT(1X,240(E12.5,1X))  
1002  FORMAT(1X,9(E12.5,1X))
	END
c
c
	SUBROUTINE Bfield(posi,t,Bf)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION posi(3),Be(3),Bf(3)

      COMMON /basics/Re,qm,PI,rad,v0,E0,B00,vc,Bio,qmsign
	COMMON /DF/Bz0,Vw,tMAX,Xb,X0,tSEC,Y0,HF

	X=posi(1)
	Y=posi(2)
	Z=posi(3)

	CALL VPOT(X,Z,A,Be)

	Bf(1)=Be(1)
	Bf(2)=Be(2)

	tt=t+TSEC
	Xr=X-Xb-Vw*tt*1.d3/Re

	R=DSQRT((Xr+HF)**2+(Y-Y0)**2)
	Bp=Bz0*(1.d0-DTANH((R-HF)/X0))/2.d0

	Bf(3)=Be(3)+Bp

	return
	END
c
c
	SUBROUTINE Efield(posi,t,Ef)	

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION posi(3),Ef(3)

      COMMON /basics/Re,qm,PI,rad,v0,E0,B00,vc,Bio,qmsign
	COMMON /DF/Bz0,Vw,tMAX,Xb,X0,tSEC,Y0,HF

	X=posi(1)
	Y=posi(2)
	Z=posi(3)

	Ef(1)=0.d0
	Ef(3)=0.d0

	tt=t+TSEC
	Xr=X-Xb-Vw*tt*1.d3/Re

	R=DSQRT((Xr+HF)**2+(Y-Y0)**2)
	Bp=Bz0*(1.d0-DTANH((R-HF)/X0))/2.d0
	Ef(2)=Bp*Vw/1.d3

	return
	END

c
C **************************************************************
C     SUBROUTINE to calculate the magnetic vector potential in 
C     the equilibrium state of the current sheet.
C
C	Check Pritchett and Coroniti (2010,JGR) paper for details
C
C
C	UPDATES: FOR DISTRIBUTIONS OF A (VECTOR POTENTIAL), SEE
C  C:\TCS_Equilibrium\DipolarizationFront\protonaurora\initial4.m
c	which includes a Bz field proportional to 1/x^IE.
C
C **************************************************************
	SUBROUTINE VPot(x,z,A,Be)

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION Be(3)
	COMMON /AB/Blobe,bz,HTx,XBD,ALPHA,PHII,Xm,Bm,Wm,IE

	F=DEXP(bz*x/(Blobe*HTx)-ALPHA*HTx**(IE-1)/((IE-1)*x**(IE-1))
	1       +Wm*Bm*datan((x-xm)/Wm)/HTx/Blobe+PHII)
	A=-Blobe*HTx*DLOG(DCOSH(F*z/HTx)/F)

	Be(1)=Blobe*F*DTANH(z*F/HTx)
	Be(2)=0.d0
	Be(3)=(HTx**IE*Blobe*ALPHA/x**IE+Wm*Wm*Bm/((x-Xm)**2+Wm*Wm)+bz)
	1      *(1.d0-z*F*DTANH(z*F/HTx)/HTx)

	RETURN
	END
c	
	SUBROUTINE rkdumb(t0,TSEC,vstart,vend)

	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	DOUBLE PRECISION vstart(6),vend(6),v(6),dv(6)
	DOUBLE PRECISION Bf(3),Ef(3),posi(3)
	DOUBLE PRECISION p(8),dp(8)
	INTEGER qmsign,fb,joutp
	CHARACTER *2 K3s, K4s

      COMMON/basics/Re,qm,PI,rad,v0,E0,B00,vc,Bio,qmsign
	COMMON/calcu/fb
	COMMON/outp/joutp,K3s,K4s,iof

	IF (joutp.EQ.1) THEN
	  OPEN(111,FILE='TRAJ_YZ_'//K3s//'_'//K4s//'.txt')	  
	ENDIF

	IF (joutp.EQ.2) THEN
	  OPEN(111,FILE='TRAJ_XY_'//K3s//'_'//K4s//'.txt')	  
	ENDIF

	IF (joutp.EQ.3) THEN
	  OPEN(111,FILE='TRAJ_TA_'//K3s//'_'//K4s//'.txt')	  
	ENDIF

	IF (joutp.EQ.5) THEN
	  OPEN(UNIT=111,FILE='P4_SST_traj.txt',STATUS='unknown')
	ENDIF

	iof=0

	t=t0
	DO I=1,6
	  v(I)=vstart(I)
	ENDDO
      DO I=1,3
	 posi(I)=v(I+3)
	ENDDO

! Calculate initial s !//zhou 2002.10.12 
      vv=v(1)*v(1)+v(2)*v(2)+v(3)*v(3)
	gam=dsqrt(1-(vv*v0*v0/(vc*vc)))
	s=t0*gam*vc

! Calculate initial p by the value of v !//zhou 2002.10.12
      p(1)=(1.d0/gam)*vc  ! p(1)=gamma*c
      DO 3 i=1,3
	 p(i+1)=v(i)
3     continue
      p(5)=vc*t0
	DO 4 i=6,8
	 p(i)=v(i-2)
4     continue

	mt=0
	DO while(dabs(t).LE.TSEC)

        CALL Bfield(posi,t,Bf)

	  Btot2=Bf(1)**2+Bf(2)**2+Bf(3)**2
	  IF (Btot2.GE.Bio**2) THEN
          iof=1
	    return
	  ENDIF

	  CALL Efield(posi,t,Ef)

	  IF (joutp.EQ.1.OR.joutp.EQ.2.OR.joutp.EQ.3.OR.joutp.EQ.5) THEN
	    IF (IMOD(mt,1000).EQ.0) THEN

	      vtot=dsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
		  write(111,20) t+TSEC,v(1),v(2),v(3),v(4),v(5),v(6)
	1                    ,Bf(1),Bf(2),Bf(3),Ef(1),Ef(2),Ef(3)
		ENDIF
	  ENDIF

	  CALL derivs(s, p, dp)
        B_size=dsqrt(Bf(1)*Bf(1)+Bf(2)*Bf(2)+Bf(3)*Bf(3))
        vv=p(2)*p(2)+p(3)*p(3)+p(4)*p(4)
	  gam=dsqrt(1-(vv*v0*v0/(vc*vc)))

	  IF (B_size.LT.10.d0) THEN
 	    h=fb*0.02d0*(qmsign+1.1d0)
	  ELSE
	    h=fb*0.261799*2/(20*qm*B00*B_size)*(qmsign+1.1d0)
	  ENDIF
	  hs=h*vc*gam  ! written by zhou. The step of s. //2002.10.12
c calculate position and velocity at time t+h	
	  CALL rk4(s, hs, dp, p)

        s=s+hs ! //zhou 2002.10.12
        t=t+h

!  change dp and p into dv and v !//zhou 2002.10.12
        DO 6 i=1,3
	    dv(i)=dp(i+1)
	    dv(i+3)=dp(i+5)
          v(i)=p(i+1)
	    v(i+3)=p(i+5)
6       continue

C position and magnetic field at time t+h
        DO 5 i=1,3
          posi(i)=v(i+3)
5        CONTINUE

	ENDDO

	DO I=1,6
	  vend(I)=v(I)
	ENDDO

      CALL Bfield(posi,t,Bf)
	CALL Efield(posi,t,Ef)

	IF (joutp.EQ.1.OR.joutp.EQ.2.OR.joutp.EQ.3.OR.joutp.EQ.5) THEN	  
	  vtot=dsqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
	  write(111,20) t+TSEC,v(1),v(2),v(3),v(4),v(5),v(6)
	1                ,Bf(1),Bf(2),Bf(3),Ef(1),Ef(2),Ef(3)
	  CLOSE(111)
	ENDIF

	joutp=0

20    format(13(2x,e13.6))
	return
	END
c
      SUBROUTINE rk4(s, hs, dp, y)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION dp(8),y(8),dym(8),dyt(8),yt(8),dydx(8)
      INTEGER qmsign
	COMMON /basics/Re,qm,PI,rad,v0,E0,B00,vc,Bio,qmsign

	hh=hs*0.5d0
      h6=hs/6.0d0
      xh=s+hh

c  first step
      DO 10 i=1,8
        dydx(i)=dp(i)
10    CONTINUE
  	DO 20 i=1,8
        yt(i)=y(i)+hh*dydx(i)
20    CONTINUE

c  second step    
      CALL derivs(xh, yt, dp)
      DO 30 i=1,8
        dyt(i)=dp(i)
30    CONTINUE
      DO 40 i=1,8
        yt(i)=y(i)+hh*dyt(i)
40    CONTINUE

c  Third step
      CALL derivs(xh, yt, dp)
      DO 50 i=1,8
        dym(i)=dp(i)
50    CONTINUE
      DO 60 i=1,8
       yt(i)=y(i)+hs*dym(i)
       dym(i)=dym(i)+dyt(i)
60    CONTINUE
c  Fourth step
      CALL derivs(s+hs, yt, dp)
      DO 70 i=1,8
       dyt(i)=dp(i)
70    CONTINUE									   
c  Accumulate increments
      DO 80 i=1,8
       y(i)=y(i)+h6*(dydx(i)+dyt(i)+2.0d0*dym(i))
80    CONTINUE

	spe1=vc*dsqrt(1.d0-vc*vc/(y(1)*y(1)))/v0
	spe2=dsqrt(y(2)*y(2)+y(3)*y(3)+y(4)*y(4))
      EE=spe1/spe2

      DO 90 i=1,3
	 y(i+1)=y(i+1)*EE
90    continue

      return
      end
!
      SUBROUTINE derivs(s, p, dp)

	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	DIMENSION v(6),dv(6),posi(3),Bf(3),Ef(3),Ec(3),Ecr(3),Ei(3)
	INTEGER fb,qmsign
      DIMENSION p(8),dp(8)
	COMMON /basics/Re,qm,PI,rad,v0,E0,B00,vc,Bio,qmsign
	COMMON /red/rp,ep,dpp
!
      DO 10 i=1,3
        posi(i)=p(i+5)
10    CONTINUE
      vv=p(2)*p(2)+p(3)*p(3)+p(4)*p(4) !//zhou 2002.10.12
	gam=dsqrt(1-vv*v0*v0/(vc*vc))

	t=s/(vc*gam)
c  Calculate mfd and efd at each time and position
      CALL Bfield(posi,t,Bf)
!  Calculate electric field at each position !//zhou 2002.10.10
      CALL Efield(posi,t,Ef)

      dp(1)=qmsign*rp*v0*v0*(p(2)*Ef(1)+p(3)*Ef(2)+p(4)*Ef(3))/(vc*vc*
     !gam)
      dp(2)=qmsign*rp*(Ef(1)+ep*(p(3)*Bf(3)-p(4)*Bf(2)))/(vc*gam)					       
      dp(3)=qmsign*rp*(Ef(2)+ep*(p(4)*Bf(1)-p(2)*Bf(3)))/(vc*gam)
      dp(4)=qmsign*rp*(Ef(3)+ep*(p(2)*Bf(2)-p(3)*Bf(1)))/(vc*gam)

c  Dr/Dt  dpp=v0/Re	   !用来改变正推和反推的条件
!      dp(5)=dpp/gam
      dp(6)=dpp*p(2)/(vc*gam)
      dp(7)=dpp*p(3)/(vc*gam)
      dp(8)=dpp*p(4)/(vc*gam)

	continue
      return
      end
c