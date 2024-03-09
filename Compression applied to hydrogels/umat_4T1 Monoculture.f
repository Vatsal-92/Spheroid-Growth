
!DEC$ ATTRIBUTES ALIAS:"umat"::UMAT
	  SUBROUTINE umat(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)


C LOCAL ARRAYS
C ----------------------------------------------------------------
C XKIRCH1 - KIRCHHOFF STRESS
C DFP - INCREMENT OF THE PERTURBED DEF. GRAD.
C DFGRD_PERT - PERTURBED DEF. GRAD.
C XKIRCH_PERT - PERTURBED KIRCHHOFF STRESS
C CMJ - (:,:,I,J) COMPONENTS OF MATERIAL JACOBIAN
C CMJVEC - ABOVE IN VECTOR FORM
C ILIST, JLIST - SET OF THE COMPONENTS TO BE PERTURBED UPON
C DUMSTRSS - DUMMY STRESS TENSOR
C XMMAT - ARRAY OF THE UNIT VECTORS DESCRIBING THE FIBRE DIRECTION
C ----------------------------------------------------------------
C
      DIMENSION XKIRCH1(3,3), DFP(3,3), DFGRD_PERT(3,3),XKIRCH_PERT(3,3)
     1 ,CMJ(3,3), CMJVEC(NTENS), ilist(6), jlist(6), DUMSTRSS(3,3) 
C     
C     ARRAYS FOR THE POLAR DECOMPOSITION     
      dimension spolF(3,3), spolC(3,3), spolCS(3,3), spolUINV(3,3),
     & U(3,3), R(3,3), RT(3,3), VAR(1), xR(3,3), DFGRD1x(3,3), 
     & xFloc(3,3), DFGRD1e(3,3), DFGRD0x(3,3), xF0inv(3,3), xF0invT(3,3),
     & xtemp(3,3), xS(3,3), Xparams(10), xSTRESS(3,3)
	 
      REAL ystart(1), x1time, x2time, rkeps, xh1, xhmin, dydx(1)
C      
C
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0,
     1 SIX=6.D0)
C
C ----------------------------------------------------------------
C UMAT FOR GROWTH WITH NEO-HOOKE MATERIAL
C
C 3D CONTINUUM ELEMENTS
C 
C LATEST VERSION 13/05/2019 - EOIN MCEVOY
C -Cleaned up old NH UMAT
C ----------------------------------------------------------------
C
C    ** PROPERTIES **
C
C PROPS(1) - C10
C PROPS(2) - D1
C PROPS(3) - switch for local or global coordinate system 
C            PROPS(8)=1 for local, PROPS(8)=2 for global
C
C    ** STATE DEPENDENT VARIABLES **   
C
C STATEV(1) - I1 
C STATEV(2) - I2
C STATEV(3) - I3
C STATEV(4) - J
C
C ----------------------------------------------------------------
C     MATERIAL PROPERTIES
C ----------------------------------------------------------------

        C10     = PROPS(1)
        D1      = PROPS(2)
        iswitch = int(PROPS(3))

        xtau = PROPS(4)   ! change timescale
	    p0   = PROPS(5)      ! pressure sensitivity
	    xb   = PROPS(6)

        XPI = 3.14159265359
		
!       Runge-Kutta params
        x1time = time(2)
        x2time = x1time + dtime
        rkeps = 1e-3
        xh1 = dtime/10.   
        xhmin = 0.
		
c-----------------------------------------------------------------------
c                Growth decomp
C--------------------------------------------------------------------  
				
C    IF *Orientation, then set F local equal to U
      IF (iswitch.EQ.1) then
      CALL spolar(dfgrd1, DFGRD1x)
      CALL spolar(dfgrd0, DFGRD0x) 
		 
      ELSE IF (iswitch.EQ.2) then
	     DFGRD1x = DFGRD1
	     DFGRD0x = DFGRD0
      ELSE
	     WRITE (6,*) '*** ERROR: Global or local coordinate system must be'
	     WRITE (6,*) 'specIFied in PROPS(3),=1 for local, =2 for global'
         CALL EXIT(1)
		 
      ENDIF     
	  
C-----------------------------------------------------------------
C    Growth
C-----------------------------------------------------------------
	  
      IF ((kstep.EQ.1).AND.(kinc.EQ.1)) THEN      
          STATEV(5) = 1.0d0   
          STATEV(6) = 0.0d0
      ENDIF
	  
! 	  Value from previous iteration
      ystart(1) = STATEV(5)
	  
        xDET=DFGRD0x(1, 1)*DFGRD0x(2, 2)*DFGRD0x(3, 3)
     1   -DFGRD0x(1, 2)*DFGRD0x(2, 1)*DFGRD0x(3, 3)
     2   +DFGRD0x(1, 2)*DFGRD0x(2, 3)*DFGRD0x(3, 1)
     3   +DFGRD0x(1, 3)*DFGRD0x(3, 2)*DFGRD0x(2, 1)
     4   -DFGRD0x(1, 3)*DFGRD0x(3, 1)*DFGRD0x(2, 2)
     5   -DFGRD0x(2, 3)*DFGRD0x(3, 2)*DFGRD0x(1, 1)
	  
      CALL KFINDINV(DFGRD0x, xF0inv, 3, ierrorflag)
      DO I=1,3
         DO j=1,3
	        xF0invT(I,J) = xF0inv(J,I)
			xSTRESS(I,J) = 0.0
         ENDDO
      ENDDO
	  
	  xSTRESS(1,1) = STRESS(1)
	  xSTRESS(2,2) = STRESS(2)
	  xSTRESS(3,3) = STRESS(3)
	  xSTRESS(1,2) = STRESS(4)
	  xSTRESS(2,1) = STRESS(4)
	  
! 	  Use stress from previous increment
      CALL KMTMS(3, 3, 3, xF0inv, 3, xSTRESS, 3, xtemp, 3)
      CALL KMTMS(3, 3, 3, xtemp, 3, xF0invT, 3, xS, 3)
	  
      xtrS = xDET*(xS(1,1)+xS(2,2)+xS(3,3))
	  ! write(6,*) xS(2,2)
	  ! write(6,*) xS(3,3)
	  ! write(6,*) xS(3,3)
	  
      Xparams(1) = xtrS
      Xparams(2) = xtau   ! timescale
	  Xparams(3) = p0     ! pressure sensitivity
	  Xparams(4) = xb     ! ratio apoptosis / proliferation

         CALL odeint(ystart,1,x1time,x2time,rkeps,xh1,xhmin,nok,
     1 nbad, Xparams)

      xg = ystart(1)

      xlamg = (xg)**(0.333)
      ! xlamg = (1+time(1))
	  
! 	  Elastic deformation gradient
      dfgrd1e = (1/xlamg)*DFGRD1x
	  
! 	  Store current value of g
      STATEV(5) = ystart(1)
      STATEV(6) = xtrS
	  
C-----------------------------------------------------------------
C  ZERO THE TANGENT MATRIX
C-----------------------------------------------------------------        
C
      DO I=1,NTENS
      DO J=1,NTENS
        DDSDDE(I,J) = ZERO
      END DO
      END DO
	  
C-----------------------------------------------------------------
C CALCULATE THE STRESS FOR THE NEO-HOOKEAN MODEL
C-----------------------------------------------------------------
C     
      call kstress_calc(dfgrd1e, C10, D1, XKIRCH1, XJ1, 
     1 XI1, XI2, XI3)
     
        STATEV(1) = XI1
        STATEV(2) = XI2
        STATEV(3) = XI3
        STATEV(4) = XJ1     
C             
C-----------------------------------------------------------------
C CONVERT KIRCHHOFF STRESS TO CAUCHY STRESS               
C-----------------------------------------------------------------      
      DUMSTRSS = XKIRCH1 / (XJ1)    

      call kmatrix2vector(DUMSTRSS, STRESS, nshr)                  
C
C*****************************************************************     
C          
C-----------------------------------------------------------------
C CALCULATE THE PERTURBATION OF THE KIRCHHOFF STRESS
C-----------------------------------------------------------------
C
      eps = 1.0e-08
C                
      ilist(1) = 1; ilist(2) = 2; ilist(3) = 3
      ilist(4) = 1; ilist(5) = 1; ilist(6) = 2
      jlist(1) = 1; jlist(2) = 2; jlist(3) = 3
      jlist(4) = 2; jlist(5) = 3; jlist(6) = 3
C      
C      
      Perturbation: DO iter = 1,NTENS
C      
        ii = ilist(iter)
        jj = jlist(iter)     
C     
          call kdelF(ii, jj, dfgrd1e, eps, DFP)
C          
C-----------------------------------------------------------------
C CREATE THE PERTURBATION OF THE DEFORMATION GRADIENT
C-----------------------------------------------------------------
C
          DFGRD_PERT = dfgrd1e + DFP        
C
C-----------------------------------------------------------------
C CALCULATE THE STRESS BASED ON THIS NEW DEFORMATION GRADIENT
C-----------------------------------------------------------------
C
          call kstress_calc(DFGRD_PERT, C10, D1, 
     1 XKIRCH_PERT, XJP, XI1, XI2, XI3)          
C
C-----------------------------------------------------------------
C DIFFERENCE BETWEEN THE PERTURBED(i,j) AND UNPERT. STRESS 
C-----------------------------------------------------------------      
C            
          do i = 1,3
          do j = 1,3        
            CMJ(i,j) = XKIRCH_PERT(i,j) - XKIRCH1(i,j)
          end do
          end do
C         
          CMJ = CMJ/(XJ1*eps)
C      
C      
C-----------------------------------------------------------------
C VECTORISE AND INSERT INTO THE DDSDDE MATRIX
C-----------------------------------------------------------------      
C
          call kmatrix2vector(CMJ, CMJVEC, NSHR)
          
          do insert = 1,NTENS
          
            DDSDDE(insert,iter) = CMJVEC(insert)
            
          end do
		  
C
       end do Perturbation  
	   
C*****************************************************************	  
C
      RETURN
      contains
C
C
C
C
C      
C-----------------------------------------------------------------------------
C           SUBROUTINES
C-----------------------------------------------------------------------------
C
C * KSTRESS_CALC -   Calculate the Kirchhoff stress based on the deformation 
C                   gradient and the elastic constants C10 and D1.
C
C * KDELF        -   Calculate the increment of the deformation gradient for
C                   a given perturbation in (i,j), with epsilon
C
C * KPRINTER     -   Print out a matrix of any size
C
C * KMTMS     -   Multiply two 2nd order tensors
C
C * KMATRIX2VECTOR     -   Convert a 3x3 matrix to a 6x1 vector
C                   
C * KDOTPROD           - Dot product of two vectors
C
C-----------------------------------------------------------------------------      
C      
C      
      subroutine kstress_calc(DGRAD, C10, D1, XKIRCH, 
     1 DET, XI1, XI2, XI3)
C     
      INCLUDE 'ABA_PARAM.INC'
C      
      intent(in) :: DGRAD, C10, D1
      intent(out):: XKIRCH, DET, XI1, XI2, XI3   
C      
      dimension DGRAD(3,3), BMAT(3,3), XKIRCH(3,3),
     1 B2MAT(3,3), U(3,3)
C     
C      
C JACOBIAN AND DISTORTION TENSOR
C
        DET=DGRAD(1, 1)*DGRAD(2, 2)*DGRAD(3, 3)
     1   -DGRAD(1, 2)*DGRAD(2, 1)*DGRAD(3, 3)
     2   +DGRAD(1, 2)*DGRAD(2, 3)*DGRAD(3, 1)
     3   +DGRAD(1, 3)*DGRAD(3, 2)*DGRAD(2, 1)
     4   -DGRAD(1, 3)*DGRAD(3, 1)*DGRAD(2, 2)
     5   -DGRAD(2, 3)*DGRAD(3, 2)*DGRAD(1, 1)
C     
C        
C     ZERO MATRICES
        DO I=1,3
        DO J=1,3
            BMAT(I,J) = ZERO          
        END DO
        END DO                               
C                    	  

C CALCULATE LEFT CAUCHY-GREEN DEFORMATION TENSOR
C        
        DO I=1,3
            DO J=1,3
                DO K = 1,3
                    BMAT(I,J) = BMAT(I,J) + DGRAD(I,K)*DGRAD(J,K)
                END DO
            END DO
        END DO
C   
C-----------------------------------------------------------------------
C     CALCULCATE THE INVARIANTS                
C-----------------------------------------------------------------------
C     
        XI1 = BMAT(1,1)+BMAT(2,2)+BMAT(3,3)
C  
C        KMTMS (M, N, L, A, KA, B, KB, C, KC)
        CALL KMTMS(3, 3, 3, BMAT, 3, BMAT, 3, B2MAT, 3)    
		
        TRB2 = B2MAT(1,1)+B2MAT(2,2)+B2MAT(3,3)       
        XI2 = 0.5*((XI1**2.0) - TRB2)     
		
        XI3=BMAT(1, 1)*BMAT(2, 2)*BMAT(3, 3)
     1   -BMAT(1, 2)*BMAT(2, 1)*BMAT(3, 3)
     2   +BMAT(1, 2)*BMAT(2, 3)*BMAT(3, 1)
     3   +BMAT(1, 3)*BMAT(3, 2)*BMAT(2, 1)
     4   -BMAT(1, 3)*BMAT(3, 1)*BMAT(2, 2)
     5   -BMAT(2, 3)*BMAT(3, 2)*BMAT(1, 1)
C                 
C--------------------------------------------------------------------
C     CALCULATE THE ISOTROPIC PORTION OF THE KIRCH STRESS
C--------------------------------------------------------------------
C
C     PART 1
      COEFF1 = TWO*C10/(DET**(TWO/THREE))
C
C     KIRCHHOFF STRESS PART 2
      TRBMAT= COEFF1*(BMAT(1,1)+BMAT(2,2)+BMAT(3,3))/THREE
C
C     KIRCHHOFF STRESS PART 1      
      DO I=1,3
        DO J=1,3
            XKIRCH(I,J) = BMAT(I,J) * COEFF1
        END DO
      END DO
C
C     SUBTRACT THE PART 2   
      DO I = 1,3
        XKIRCH(I,I) = XKIRCH(I,I) - TRBMAT
      END DO
C
C     FORM PART 3     
      COEFF3 = 2*(DET-ONE)*DET/D1
C
C     ADD TO THE PREVIOUS PARTS     
      DO I = 1,3
        XKIRCH(I,I) = XKIRCH(I,I) + COEFF3
      END DO          		
C        
      return
      end subroutine kstress_calc
 
C-----------------------------------------------------------------------------

      SUBROUTINE derivs(x, y, dydx, Xparams)

      INCLUDE 'ABA_PARAM.INC'

      intent(in) x, y, Xparams
      intent(out) dydx
      REAL x, y, dydx
      DIMENSION Xparams(10), y(1), dydx(1)
	 
C 	    x = time
C       y = val
	 
C 	  Growth Parameters (only needed in SUBROUTINE)
      xtau = Xparams(2)   ! change timescale
	  p0   = Xparams(3)      ! pressure sensitivity
	  xb   = Xparams(4) ! ratio of Apoptosis rate to Proliferation rate 
	  
      xtrS = Xparams(1) ! trace of stress
      ct   = y(1) ! number of cells
	  
      IF (xtrS.LE.0.0) THEN
		dgdt = (ct/xtau)*(1.0d0 + (xtrS/(3.0d0*p0)) - xb)
      ELSE
		dgdt = (ct/xtau)*(1.0d0 - xb)
      ENDIF
C 	  
      dydx(1) = dgdt
	  
      END SUBROUTINE derivs
	  
C-----------------------------------------------------------------------------

      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad
     1 ,Xparams)
	  
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      REAL eps,h1,hmin,x1,x2,ystart(nvar),TINY
      DIMENSION Xparams(10)
C       EXTERNAL derivs,rkqs
      Parameter (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.e-30)
C Runge-Kutta driver with adaptive stepsize control. Integrate the starting values ystart(1:nvar)
C from x1 to x2 with accuracy eps, storing intermediate results in the common block /path/.
C h1 should be set as a guessed first stepsize, hmin as the minimum allowed stepsize (can
C be zero). On output nok and nbad are the number of good and bad (but retried and
C fixed) steps taken, and ystart is replaced by values at the END of the integration interval.
C derivs is the user-supplied SUBROUTINE for calculating the right-hand side derivative, while
C rkqs is the name of the stepper routine to be used. /path/ contains its own information
C about how often an intermediate value is to be stored.
      INTEGER i,kmax,kount,nstp
      REAL dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),y(NMAX),
     1 yp(NMAX,KMAXX),yscal(NMAX)
C       COMMON /path/ kmax,kount,dxsav,xp,yp
C       User storage for intermediate results. Preset dxsav and kmax.
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
	  
      DO i=1,nvar
         y(i)=ystart(i)
      ENDDO
	  
      IF (kmax.gt.0) xsav=x-2.*dxsav    !Assures storage of first step.

      DO nstp=1,MAXSTP    !Take at most MAXSTP steps.
         CALL derivs(x,y,dydx,Xparams)    ! CALL of derivs function
		 
         DO i=1,nvar
C          Scaling used to monitor accuracy. This general-purpose choice can be modIFed IF need be.
           yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
         ENDDO
		 
         IF(kmax.gt.0)then
            IF(abs(x-xsav).gt.abs(dxsav)) then !Store intermediate results.
               IF(kount.lt.kmax-1)then
                  kount=kount+1
                  xp(kount)=x
				  
                  DO i=1,nvar
                     yp(i,kount)=y(i)
                  ENDDO
				  
                  xsav=x
               ENDIF
            ENDIF
         ENDIF
		 
         IF((x+h-x2)*(x+h-x1).gt.0.) h=x2-x    !IF stepsize can overshoot, decrease.

         CALL rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,Xparams)

         IF(hdid.eq.h)then
            nok=nok+1
         ELSE
            nbad=nbad+1
         ENDIF
		 
         IF((x-x2)*(x2-x1).ge.0.)then    !Are we DOne?
		 
            DO i=1,nvar
               ystart(i)=y(i)
            ENDDO
			
            IF(kmax.ne.0)then
               kount=kount+1    !Save final step.
               xp(kount)=x
			   
               DO i=1,nvar
                  yp(i,kount)=y(i)
               ENDDO
			   
            ENDIF
            RETURN    !Normal exit.
         ENDIF
		 
         IF(abs(hnext).lt.hmin) then
		    pause 'stepsize smaller than min odeint' 
		    CALL EXIT(1)
		 ENDIF
	     h=hnext
			
      ENDDO
	  
      pause 'too many steps in odeint' 
	  CALL EXIT(1)
	  
      RETURN
      END SUBROUTINE odeint
	  
c-----------------------------------------------------------------------------
   
      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,Xparams)
      INTEGER n,NMAX
      REAL eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n)
      REAL yout(n)
      DIMENSION Xparams(10)
C       EXTERNAL derivs
      Parameter (NMAX=50)    !Maximum number of equations.
C    USES derivs,rkck
C FIFth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy
C and adjust stepsize. Input are the depENDent variable vector y(1:n) and its derivative
C dydx(1:n) at the starting value of the indepENDent variable x. Also input are the stepsize
C to be attempted htry, the required accuracy eps, and the vector yscal(1:n) against
C which the error is scaled. On output, y and x are replaced by their new values, hdid is the
C stepsize that was actually accomplished, and hnext is the estimated next stepsize. derivs
C is the user-supplied SUBROUTINE that computes the right-hand side derivatives.
      INTEGER i
      REAL errmax,h,htemp,xnew,yerr(n),ytemp(n),SAFETY,PGROW,
     1 PSHRNK,ERRCON
      REAL ak2(n),ak3(n),ak4(n),ak5(n),ak6(n),
     1 A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51,
     1 B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,
     1 DC4,DC5,DC6
      Parameter (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.89e-4)
      Parameter (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40.,
     1 B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5,
     1 B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512.,
     1 B63=575./13824.,B64=44275./110592.,B65=253./4096.,
     1 C1=37./378.,C3=250./621.,C4=125./594.,C6=512./1771.,
     1 DC1=C1-2825./27648.,DC3=C3-18575./48384.,
     1 DC4=C4-13525./55296.,DC5=-277./14336.,DC6=C6-.25)
C       The value ERRCON equals (5/SAFETY)**(1/PGROW), see use below.

      h=htry    !Set stepsize to the initial trial value.	  
C 	  SUBROUTINE HAS BEEN INCLUDED - rkck not used
C 1     CALL rkck(y,dydx,n,x,h,ytemp,yerr,Xparams,derivs)    !Take a step.

c-----------------------------------------------------------------------------
1     DO i=1,n    !First step.
         ytemp(i)=y(i)+B21*h*dydx(i)
      ENDDO
	  
      CALL derivs(x+A2*h,ytemp,ak2, Xparams)    !Second step.

      DO i=1,n
         ytemp(i)=y(i)+h*(B31*dydx(i)+B32*ak2(i))
      ENDDO
	  
      CALL derivs(x+A3*h,ytemp,ak3, Xparams)    !Third step.

      DO i=1,n
         ytemp(i)=y(i)+h*(B41*dydx(i)+B42*ak2(i)+B43*ak3(i))
      ENDDO
	  
      CALL derivs(x+A4*h,ytemp,ak4, Xparams)    !Fourth step.

      DO i=1,n
         ytemp(i)=y(i)+h*(B51*dydx(i)+B52*ak2(i)+B53*ak3(i)+
     1 B54*ak4(i))
      ENDDO

      CALL derivs(x+A5*h,ytemp,ak5, Xparams)    !FIFth step.

      DO i=1,n
         ytemp(i)=y(i)+h*(B61*dydx(i)+B62*ak2(i)+B63*ak3(i)+
     1 B64*ak4(i)+B65*ak5(i))
      ENDDO

      CALL derivs(x+A6*h,ytemp,ak6, Xparams)    !Sixth step.

      DO i=1,n    !Accumulate increments with proper weights.
         yout(i)=y(i)+h*(C1*dydx(i)+C3*ak3(i)+C4*ak4(i)+
     1 C6*ak6(i))
      ENDDO

      DO i=1,n
      !   Estimate error as dfference between fourth and fIFth order methods.
         yerr(i)=h*(DC1*dydx(i)+DC3*ak3(i)+DC4*ak4(i)+DC5*ak5(i)
     1 +DC6*ak6(i))
      ENDDO
c-----------------------------------------------------------------------------

      errmax=0.    !   Evaluate accuracy.
	  
      DO i=1,n
         errmax=max(errmax,abs(yerr(i)/yscal(i)))
      ENDDO
	  
      errmax=errmax/eps    !Scale relative to required tolerance.
      IF(errmax.gt.1.)then    !Truncation error too large, reduce stepsize.
         htemp=SAFETY*h*(errmax**PSHRNK)
         h=sign(max(abs(htemp),0.1*abs(h)),h)     !No more than a factor of 10.
         xnew=x+h
         IF(xnew.eq.x) then
		    pause 'stepsize underflow in rkqs' 
			CALL EXIT(1)
	     ENDIF
         goto 1    !For another try.
      ELSE     !Step succeeded.  Compute size of next step.
         IF(errmax.gt.ERRCON)then
            hnext=SAFETY*h*(errmax**PGROW)
         ELSE    !No more than a factor of 5 increase.
            hnext=5.*h
         ENDIF
         hdid=h
         x=x+h
		 
         DO i=1,n
            y(i)=yout(i)
         ENDDO
		 
         RETURN
      ENDIF
	  
      END SUBROUTINE rkqs
	  
c----------------------------------------------------------------------------- 
C-----------------------------------------------------------------------------

      subroutine kdelF(m, n, DGRAD, eps, DF)
      
      INCLUDE 'ABA_PARAM.INC'
      
      intent (in) :: DGRAD, eps, m, n
      intent (out):: DF      
        
C Input: the index's i & j; The current deformation gradient (DGRAD). The perturbation
C        increment (eps)
C
C Output: The perturbed increment DF
C
      dimension dyad1(3,3), dyad2(3,3), DGRAD(3,3), DF(3,3), DFp1(3,3)
        
c Zero the dyad matrices
c
      do i = 1,3
        do j = 1,3
            dyad1(i,j) = zero
            dyad2(i,j) = zero
        end do
      end do
      
c Place the 1's in the correct location        
      dyad1(m,n) = 1.0;
c
      dyad2(n,m) = 1.0;
c      
c      KMTMS (M, N, L, A, KA, B, KB, C, KC)  
      call KMTMS(3, 3, 3, dyad1, 3, DGRAD, 3, DFp1, 3)
      DF = DFp1
      
      
      call KMTMS(3, 3, 3, dyad2, 3, DGRAD, 3, DFp1, 3)
      DF = DF + DFp1
           
      
      DF = 0.5*DF*eps            
      
      end subroutine kdelF
      
      
c-----------------------------------------------------------------------------            
       
      subroutine kprinter(tens, m, n)
      
      INCLUDE 'ABA_PARAM.INC'
      
      intent(in):: tens, m, n      
      
      dimension tens(m,n)
        
        write(6,*)
        do i = 1,m
        do j = 1,n
            write(6,'(e19.9)',advance='no'),tens(i,j)
        end do
        write(6,*)
        end do
        write(6,*)
      return
      end subroutine kprinter
      
c------------------------------------------------------------------------------      
      
      SUBROUTINE KMTMS (M, N, L, A, KA, B, KB, C, KC)
      
      INCLUDE 'ABA_PARAM.INC'
C      
      intent(in) :: M, N, L, A, KA, B, KB, KC
      intent(out):: C      
C      
C
C    PRODUCT OF REAL MATRICES
C
      DIMENSION A(KA,N), B(KB,L), C(KC,L)
      DOUBLE PRECISION W
C       
C
      DO 30 J = 1,L
         DO 20 I = 1,M
            W = 0.D0
            DO 10 K = 1,N
               W = W + A(I,K) * B(K,J)
   10       CONTINUE
            C(I,J) = W
   20    CONTINUE
   30 CONTINUE
      RETURN
      END SUBROUTINE


c-----------------------------------------------------------------------

      
      subroutine kmatrix2vector(XMAT, VEC, NSHR)
      
      INCLUDE 'ABA_PARAM.INC'
      
      intent(in) :: XMAT, NSHR
      intent(out):: VEC
      
      dimension xmat(3,3), vec(6)
  
        do i=1,3
            vec(i) = xmat(i,i);
        end do
               
        vec(4) = xmat(1,2);
        
        IF (NSHR==3) then
            vec(5) = xmat(1,3);
            vec(6) = xmat(2,3);
        END IF
      
      end subroutine kmatrix2vector
      
C-----------------------------------------------------------------------      
      
      subroutine kdotprod(A, B, dotp, n)
      
      INCLUDE 'ABA_PARAM.INC'
      
      intent(in) :: A, B, n
      intent(out):: dotp      
      
      dimension A(n), B(n)
      dotp = 0.0
      
      do i = 1,n
        dotp = dotp + A(i)*B(i)
      end do
      
      end subroutine kdotprod      
      
      
      
      END

C---------------------------------------------------------------------------
C     POLAR DECOMPOSITION
C--------------------------------------------------------------------------- 

      SUBROUTINE spolar(dfgrd, xFloc) 
C
      INCLUDE 'ABA_PARAM.INC'
	  
      DIMENSION spolF(3,3), spolC(3,3), spolCS(3,3), spolUINV(3,3),
     & U(3,3), R(3,3), RT(3,3), xFloc(3,3), dfgrd(3,3)
       

      spolF=dfgrd
      DO  J=1,3
        DO  I = 1,3
          spolC(I,J)= 0.0D0
          DO  K=1,3
            spolC(I,J)= spolC(I,J) + spolF(K,I)*spolF(K,J)
          ENDDO
        ENDDO
      ENDDO
      
      DO  J=1,3
        DO  I=1,3
          spolCS(I,J)= 0.0D0
          DO K=1,3
            spolCS(I,J)= spolCS(I,J) + spolC(I,K)*spolC(K,J)
          ENDDO
        ENDDO
      ENDDO
      
      spolC1= spolC(1,1) + spolC(2,2) + spolC(3,3)
      spolC2= 0.5D0 * (spolC1**2.0D0 - (spolCS(1,1)+spolCS(2,2)+
     +    spolCS(3,3)))
      spolC3= spolC(1,1) * (spolC(2,2)*spolC(3,3)-spolC(2,3)*spolC(3,2))
     1    +spolC(1,2) * (spolC(2,3)*spolC(3,1)-spolC(2,1)*spolC(3,3))
     2    +spolC(1,3) * (spolC(2,1)*spolC(3,2)-spolC(2,2)*spolC(3,1))
      spolU3= SQRT(spolC3)
      spolX1=2.0**5.0 /27.0*(2.0*spolC1**3.0-9.0*spolC1*spolC2+
     &27.0*spolC3)
      spolX2= 2.**10.0 /27.0*(4.0*spolC2**3.0-spolC1**2.0*spolC2**2.0 +
     &4.0*spolC1**3.0*spolC3-18.0*spolC1*spolC2*spolC3+27.0*spolC3**2.0)
     
      IF (spolX2.LT.0.) spolX2= 0.0
      F1= spolX1 + SQRT(spolX2)
      IFLAG2= 0
      IFLAG1= 0
      IF (F1.LT.0.0) IFLAG1= 1
      F2= spolX1 - SQRT(spolX2)
      IF (F2.LT.0.0) IFLAG2= 1
      IF (IFLAG2.EQ.1) F2= -F2
      IF (IFLAG1.EQ.1) F1= -F1
      spolX3= -2.0/3.0*spolC1 + F1**(1.0/3.0) + F2**(1.0/3.0)
      IF (IFLAG1.EQ.1) spolX3= -2.0/3.0*spolC1 + F2**(1.0/3.0) -
     +                     F1**(1.0/3.0)
      IF (IFLAG2.EQ.1) spolX3= -2.0/3.0*spolC1 + F1**(1.0/3.0) -
     +                     F2**(1.0/3.0)
      B= -2.0*spolC1
      IF (spolX3.EQ.B) THEN
      U1= SQRT(spolC1+2.0*SQRT(spolC2))
      ELSE 
      U1= 0.5 * (SQRT(2.0*spolC1+spolX3) + SQRT(2.0*spolC1 -
     +    spolX3+16.0*SQRT(spolC3)/SQRT(2.0*spolC1+spolX3)))
      ENDIF
      U2= SQRT(spolC2+2.0*spolU3*U1)
      B1= spolU3**2.0 * (spolU3+U1*spolC1) + 
     &U1**2.0 * (U1*spolC3+spolU3*spolC2)
      B2= U1 * (U1*U2-spolU3) / B1
      B3=-(U1*U2-spolU3) * (spolU3+U1*spolC1) / B1
      B4= (U2*spolU3*(spolU3+U1*spolC1) + 
     &U1**2.0 * (U2*spolC2+spolC3))/B1
     
      DO  J=1,3
        DO  I=1,3
          spolUINV(I,J)= B2*spolCS(I,J) + B3*spolC(I,J)
          IF (I.EQ.J)   spolUINV(I,J)= spolUINV(I,J) + B4
        ENDDO
      ENDDO
      
      DO  J=1,3
        DO  I=1,3
          R(I,J)=0.0
          DO  K=1,3
            R(I,J)= R(I,J) + spolF(I,K)*spolUINV(K,J)
          ENDDO
        ENDDO
      ENDDO
      DO  I=1,3
        DO  J=1,3
          RT(I,J) = R(J,I)
        ENDDO
      ENDDO
      
      U=0.0
      DO I=1,3
           DO J=1,3
               DO M=1,3
		          U(I,J)=U(I,J)+RT(I,M)*spolF(M,J)
               ENDDO
           ENDDO
       ENDDO  
	   
       DO I=1,3
          DO J=1,3
              xFloc(I,J) = 0.0
          END DO
       END DO	  
	  
C 	  F (local) = F (abaqus) * RT
      DO I=1,3
           DO J=1,3
               DO M=1,3
                  xFloc(I,J)=xFloc(I,J) + DFGRD(I,M)*RT(M,J)
               ENDDO
           ENDDO
       ENDDO  
C
      END SUBROUTINE
	  
C -----------------------------------------------------------------------------
C
      SUBROUTINE KFINDINV(xmatrix, xinverse, n, ierrorflag)
      INCLUDE 'ABA_PARAM.INC'
C      
      INTENT(in) :: xmatrix, n
      INTENT(out):: xinverse, ierrorflag
C      
C
      DIMENSION xmatrix(n,n), xinverse(n,n), augmatrix(n,2*n)
C	
	  LOGICAL :: flag = .true.
C	
C	
C	augment input matrix with an identity matrix
		DO i = 1, n
			DO j = 1, 2*n
				IF (j <= n ) THEN
					augmatrix(i,j) = xmatrix(i,j)
				ELSE IF ((i+n) == j) THEN
					augmatrix(i,j) = 1
				ELSE
					augmatrix(i,j) = 0
				ENDIF
			END DO
		END DO
C	
C	
		DO k =1, n-1
			IF (augmatrix(k,k) == 0) THEN
				Write(6,*),'IF statement', k
				flag = .false.
				DO i = k+1, n
					IF (augmatrix(i,k) /= 0) THEN
						DO j = 1,2*n
							augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
						END DO
						flag = .true.
						exit
					ENDIF
					IF (flag .eqv. .false.) THEN
						Write(6,*),"matrix is non - invertible 1"
						xinverse = 0
						ierrorflag = -1
						RETURN
					ENDIF
				END DO
			ENDIF
			DO j = k+1, n			
				xm = augmatrix(j,k)/augmatrix(k,k)
				DO i = k, 2*n
					augmatrix(j,i) = augmatrix(j,i) - xm*augmatrix(k,i)
				END DO
			END DO
		END DO
C
C	
C	  test for invertibility
	  DO i = 1, n
		  IF (augmatrix(i,i) == 0) THEN
			  Write(6,*), "matrix is non - invertible 2"
			  xinverse = 0
			  ierrorflag = -1
			  RETURN
		  ENDIF
	  END DO
C	
C	  make diagonal elements as 1
	  DO i = 1 , n
		  xm = augmatrix(i,i)
		  DO j = i , (2 * n)				
			   augmatrix(i,j) = (augmatrix(i,j) / xm)
		  END DO
	  END DO
C	
C	reduced right side half of augmented matrix to identity matrix
	  DO k = n-1, 1, -1
		  DO i =1, k
		  xm = augmatrix(i,k+1)
			  DO j = k, (2*n)
				augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * xm
		      END DO
		  END DO
	  END DO				
C	
C	  store answer
	  DO i =1, n
		  DO j = 1, n
			xinverse(i,j) = augmatrix(i,j+n)
		  END DO
	  END DO
	  ierrorflag = 0
      END SUBROUTINE KFINDINV   