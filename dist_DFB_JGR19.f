C **************************************************************
c     program name: dist.f 
C
C     ==========================================================
C     Calculate the ion distribution function at a certain location 
C     of the current sheet with a dipolarization front propagating 
C     towards the observational location. 
C     
C **************************************************************

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION Be(3),Beini(3),Bf(3),vstart(6),vend(6),posiini(3)
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
	COMMON /DF/Bz0,Vw,tMAX,Xb,X0,tSECrd,Y0,HF

! Basic parameters
	qe=1.60217646d-19  ! elementary charge (coulombs)
	pm=1.67262158d-27  ! proton mass (kg) 
	qm=qe/pm
	vc=3.d8            ! speed of light (m/s)
	Re=6.3712d6        ! Earth Radius (m)
      PI=3.14159265359D0
	Bio=5.6d6          ! B FIELD STRENGTH AT IONOSPHERE (LATITUDE=68 degree)

! Current sheet equilibrium parameters	

	XBD=-10.d0        ! X Location where b.c. are determined (Re)
	Blobe=50.d0		  ! Magnetic field in the Lobe at X=XBD (nT)
	HTx=2.d0	      ! Current sheet Half-Thickness at X=XBD (Re)
	Bzxbd=15.d0       ! Magnetic Field Bz component at X=XBD, z=0 (nT)
	bz=5.d0 	      ! Magnetic field Bz component at X=-inf, z=0 (nT)
	IE=4			  ! THE EXPONENTIAL DECAY FACTOR OF Bz(x,0)
      Xm=-20.d0         ! X location of the Bz secondary peak (Re)
	Bm=0.d0			  ! The Bz value of the secondary peak (nT)
	Wm=2.d0           ! The width of the Bz secondary peak (Re)

	ALPHA=(Bzxbd-bz-Wm*Wm*Bm/((XBD-Xm)**2+Wm*Wm))
	1                     *XBD**IE/(Blobe*HTx**IE)
	PHII=HTx**(IE-1)*ALPHA/((IE-1.d0)*XBD**(IE-1))-bz*XBD/(Blobe*HTx)
	1          -Wm*Bm*datan((XBD-Xm)/Wm)/HTx/Blobe	

!  SO THAT Bz(x,0) = HTx^IE*Blobe*ALPHA/X^IE+bz

	Vth=1400.d0		  ! Proton thermal velocity (km/s)
	Vsh=Vth*Vth/(HTx*Re*qm*Blobe/1.d3)*1.d9 ! Proton bulk velocity (km/s)
	DEN=0.3d0		  ! Plasma density at X=XBD, z=0 (cm^-3)
	DENf=DEN/(DEXP(Vsh**2/Vth**2)*PI**1.5*Vth**3)

	Vth2=400.d0		  ! Proton thermal velocity of background population
	Vsh2=0.d0		  ! no bulk velocity for background population
	DEN2=0.05d0		  
	DEN2f=DEN2/(PI**1.5*Vth2**3)

! Dipolarization front parameters
	Bz0=20.d0    ! Magnetic field Bz enhancement associated with the DF (nT)
	Vw=400.d0   ! Earthward propagating speed of the DF (km/s) 
! ENDTEST
	Xb=-20.d0   ! The initial location of the DF (Re, at y=0)
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

	xini=-11.5d0
	yini=0.d0
	zini=-1.d0

	vstart(4)=xini
	vstart(5)=yini
	vstart(6)=zini
	posiini(1)=xini
	posiini(2)=yini
	posiini(3)=zini

	tMAX=(xini-Xb)*Re/1.d3/Vw  ! DFB center

	CALL VPot(xini,zini,Aini,Beini) 
	Bti=DSQRT(Beini(1)**2+Beini(2)**2+Beini(3)**2)

	N3=42
	BL33=4000.d0
	H3=BL33/DBLE(N3)

	N4=80
	BL44=2.d0*PI
	H4=BL44/N4

	NAMIN=-30
	NAMAX=30
	NAH=1

	NTW=12

! Distributions in the perp-para (YZ) plane
	DO nSEC=148,148,3
	 write(SECs,'(i3)') nSEC
	 tSEC=DBLE(nSEC)

	 CALL Bfield(posiini,tSEC,Bf)
	 Bangle=datan2(Bf(1),Bf(3))
	 sinBa=dsin(Bangle)
	 cosBa=dcos(Bangle)

	 OPEN(UNIT=nSEC,FILE='PP_Zm1_Y0_v400_'//SECs//'s_t.txt')
	 DO K3=1,N3
	  DO K4=1,N4
	   PSD(K4)=0.d0
	   KN=0
	   DO TS=-DBLE(NTW)/2.d0,DBLE(NTW)/2.d0,0.1
	    tSECrd=tSEC+TS
	    KN=KN+1
	    ANG=0.d0
	    ROW=H3*K3
	    FIS=H4*(K4-1.d0)

	    WX=ROW*DCOS(FIS)*DCOS(ANG)
	    X3=ROW*DSIN(FIS)*DCOS(ANG)

	    vperp=ROW*DSIN(ANG)*1.d3/v0
	    vperp2=WX*1.d3/v0
	    vpara=X3*1.d3/v0

	    vstart(1)=vpara*sinBa+vperp*cosBa
	    vstart(2)=vperp2
	    vstart(3)=vpara*cosBa-vperp*sinBa

	    Pini=vstart(2)*v0/1.d3+Aini*1.d-9*Re*1.d-3*qm

!	    IF (IABS(K3-35).LE.0.AND.IABS(K4-1).LE.0) THEN
!	      joutp=2
!	      write(K3s,'(i2)') K3
!	      write(K4s,'(i2)') K4
!	    ENDIF

	    CALL rkdumb(t0,TSECrd,vstart,vend)

	   	IF (iof.EQ.1) THEN
	      PSDANG=0.d0
	    ELSE  
		  X2=vend(1)*v0/1.d3
		  X3=vend(2)*v0/1.d3
	      X4=vend(3)*v0/1.d3
		  xend=vend(4)
		  zend=vend(6)

	      CALL VPot(xend,zend,Aend,Be)  
	      Pend=X3+Aend*1.d-9*Re*1.d-3*qm
		  PSDANG=DENf*DEXP((-X2**2-X3**2-X4**2)/Vth**2)
	1                  *DEXP(2.d0*Vsh*Pend/Vth**2)
	1            +DEN2f*DEXP((-X2**2-X3**2-X4**2)/Vth2**2)
	    ENDIF
	    PSD(K4)=PSD(K4)+PSDANG
	   ENDDO
	   PSD(K4)=PSD(K4)/KN
	  ENDDO
	  write(nSEC,1001) PSD
	 ENDDO
	 CLOSE(UNIT=nSEC)

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
	COMMON /DF/Bz0,Vw,tMAX,Xb,X0,tSECrd,Y0,HF

	X=posi(1)
	Y=posi(2)
	Z=posi(3)


	CALL VPOT(X,Z,A,Be)

	Bf(1)=Be(1)
	Bf(2)=Be(2)

	tt=t+TSECrd
	Xr=X-Xb-Vw*tt*1.d3/Re

!  Case 1: Planar DF with Y-width of 2*HF
!	GY=dexp(-(Y-Y0)**2/HF**2)
!	Bp=GY*Bz0*(1.d0-DTANH(Xr/X0))/2.d0

!  Case 2: Circular DF with radius of HF
	R=DSQRT((Xr+HF)**2+(Y-Y0)**2)
	Bp=Bz0*(1.d0-DTANH((R-HF)/X0))/2.d0

!  Case 3: Half-Circular DF with radius of HF
!	IF (DABS(Y-Y0).GT.HF) THEN
!	  Bp=0.d0
!	ELSE
!	  Bz0l=Bz0*(1.d0-((Y-Y0)/HF)**2)
!	  Bp=Bz0l*(1.d0-DTANH((Xr+HF-DSQRT(HF**2-(Y-Y0)**2))/X0))/2.d0
!	ENDIF

! Case 4: Sharp (non-continuous) DF with gradual decreasing Bz (with Y-width of 2*HF as well)
!         Also include Bx perturbations

!	HDX=1.d0  ! The characteristic length of the Dipolarized flux bundle
!	HDZ=1.d0  ! The characteristic thickness of the Dipolarized flux bundle

!	IF (Xr.GE.0.d0) THEN
!	  Bp=0.d0
!	ELSE
!	  GY=dexp(-(Y-Y0)**2/HF**2)
!	  Bp=GY*Bz0*DEXP(-Z**2/HDZ**2)/((DCOSH(Xr/HDX))**2)
!	  Bpx=2.d0*GY*Bz0*Z*HDX*DEXP(-Z**2/HDZ**2)*DTANH(Xr/HDX)/HDZ**2
!	  Bf(1)=Be(1)+Bpx
!	ENDIF

! Case 5: Parabolic DF with characteristic half-width of HF
!	Bp=Bz0*(1.d0-DTANH((Xr+(Y-Y0)**2/HF)/X0))/2.d0

! Case 6: Circular case (Case 2) multiplied by a Gaussian function
!         to generate a DFB with sharp Bz gradient at the leading edge and much weaker gradient at the trailing edge
!	R=DSQRT((Xr+HF)**2+(Y-Y0)**2)
!	Gau=DEXP(-Xr*Xr/(HF*HF))
!	Bp=Bz0*Gau*(1.d0-DTANH((R-HF)/X0))/2.d0

	Bf(3)=Be(3)+Bp

	return
	END
c
c
	SUBROUTINE Efield(posi,t,Ef)	

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION posi(3),Ef(3)

      COMMON /basics/Re,qm,PI,rad,v0,E0,B00,vc,Bio,qmsign
	COMMON /DF/Bz0,Vw,tMAX,Xb,X0,tSECrd,Y0,HF

	X=posi(1)
	Y=posi(2)
	Z=posi(3)

	Ef(1)=0.d0
	Ef(3)=0.d0

	tt=t+TSECrd
	Xr=X-Xb-Vw*tt*1.d3/Re

!  Case 1: Planar DF with Y-width of 2*HF
!	GY=dexp(-(Y-Y0)**2/HF**2)
!	Bp=GY*Bz0*(1.d0-DTANH(Xr/X0))/2.d0

!  Case 2: Circular DF with radius of HF
	R=DSQRT((Xr+HF)**2+(Y-Y0)**2)
	Bp=Bz0*(1.d0-DTANH((R-HF)/X0))/2.d0

!  Case 3: Half-Circular DF with radius of HF
!	IF (DABS(Y-Y0).GT.HF) THEN
!	  Bp=0.d0
!	ELSE
!	  Bz0l=Bz0*(1.d0-((Y-Y0)/HF)**2)
!	  Bp=Bz0l*(1.d0-DTANH((Xr+HF-DSQRT(HF**2-(Y-Y0)**2))/X0))/2.d0
!	ENDIF

! Case 4: Sharp (non-continuous) DF with gradual decreasing Bz (with Y-width of 2*HF as well)
!         Also include Bx perturbations

!	HDX=1.d0  ! The characteristic length of the Dipolarized flux bundle
!	HDZ=1.d0  ! The characteristic thickness of the Dipolarized flux bundle
!
!	IF (Xr.GE.0.d0) THEN
!	  Bp=0.d0
!	ELSE
!	  GY=dexp(-(Y-Y0)**2/HF**2)
!	  Bp=GY*Bz0*DEXP(-Z**2/HDZ**2)/((DCOSH(Xr/HDX))**2)
!	ENDIF

! Case 5: Parabolic DF with characteristic half-width of HF
!	Bp=Bz0*(1.d0-DTANH((Xr+(Y-Y0)**2/HF)/X0))/2.d0

! Case 6: Circular case (Case 2) multiplied by a Gaussian function
!         to generate a DFB with sharp Bz gradient at the leading edge and much weaker gradient at the trailing edge
!	R=DSQRT((Xr+HF)**2+(Y-Y0)**2)
!	Gau=DEXP(-Xr*Xr/(HF*HF))
!	Bp=Bz0*Gau*(1.d0-DTANH((R-HF)/X0))/2.d0

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
	  OPEN(111,FILE='TRAJ_XY_'//K3s//'_'//K4s//'_vs0.txt')	  
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