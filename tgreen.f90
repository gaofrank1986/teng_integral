!  
!  G=-1/4/pi*1/r 
!
!
        SUBROUTINE TGRN(H,X,X0,Y,Y0,Z,Z0,XHF) 
! 
        IMPLICIT NONE
		INTEGER I

	    REAL*8,INTENT(IN)::  H,X,Y,Z,X0,Y0,Z0
		REAL*8,INTENT(OUT):: XHF(4)
	    REAL*8 DX,DY,DZ,DZ1
		REAL*8 RXY2,DZ02,DZ12,R02,R12,SR,SR1,PI,PI4
!
	    DATA PI/3.14159265358979/ 
!C 
!C H: water depth, negtive for infinity water depth 
!C
        PI4 = 4.0D0*PI

        DX=X-X0
        DY=Y-Y0
        DZ  =Z- Z0


		RXY2=DX*DX+DY*DY
	    DZ02=DZ*DZ
		DZ12=DZ1*DZ1
!
	   IF(H .GT. 0) THEN
        R02=RXY2+DZ02                      
        R12=RXY2+DZ12           
        SR =DSQRT(R02) 
        SR1=DSQRT(R12)   
! 
        XHF(1)= 1.D0/SR  
        XHF(2)=-DX/SR**3  
        XHF(3)=-DY/SR**3  
        XHF(4)=-DZ/SR**3     
! 
	   ELSE
        R02=RXY2+DZ02                      
        SR =DSQRT(R02) 
! 
        XHF(1)= 1.D0/SR  
        XHF(2)=-DX/SR**3  
        XHF(3)=-DY/SR**3  
        XHF(4)=-DZ/SR**3     
! 
	   ENDIF
!
       DO 120   I=1,  4
         XHF(I)=-XHF(I)/PI4 
120	   CONTINUE
!
        RETURN 
        END 
		           

      SUBROUTINE MGREEN(ROU,ZZ,ZP,H,GF) 
      IMPLICIT NONE
	  INTEGER I
	  REAL*8,INTENT(IN)::  ROU,ZZ,ZP,H
	  REAL*8,INTENT(OUT):: GF(6) ! 6?
      REAL*8 DZ1,DZ12,SR1,R12,PI,PI4,ROU1
	  DATA PI/3.14159265358979/ 
	  PI4 = 4.0D0*PI

      DZ1 = ZZ+ZP+2.0*H
	  DZ12 = DZ1*DZ1
	  R12 = ROU+DZ12
	  SR1 = DSQRT(R12)
	  ROU1 = DSQRT(ROU)

	  GF(1) = 1.D0/SR1 
	  GF(2) = -ROU1/SR1**3 
	  GF(3) = -DZ1/SR1**3

      DO 130   I=1,  3
         GF(I)=-GF(I)/PI4 
130	   CONTINUE

       RETURN 
       END 

      SUBROUTINE MGREEN1(ROU,ZZ,ZP,H,GF) 
      IMPLICIT NONE
	  INTEGER I
	  REAL*8,INTENT(IN)::  ROU,ZZ,ZP,H
	  REAL*8,INTENT(OUT):: GF(6) ! 6?
      REAL*8 DZ1,DZ12,SR1,R12,PI,PI4,ROU1,DZ,DZ02,R02,SR
	  DATA PI/3.14159265358979/ 
	  PI4 = 4.0D0*PI
      DZ = ZZ-ZP
      DZ1 = ZZ+ZP+2.0*H
      DZ02 = DZ*DZ
	  DZ12 = DZ1*DZ1
	  R02 = ROU+DZ02
	  R12 = ROU+DZ12
	  SR = DSQRT(R02)
	  SR1 = DSQRT(R12)
	  ROU1 = DSQRT(ROU)

	  GF(1) = 1.0D0/SR+1.D0/SR1 
	  GF(2) = -ROU1/SR**3-ROU1/SR1**3 
	  GF(3) = -DZ/SR**3-DZ1/SR1**3

      DO 230   I=1,  3
         GF(I)=-GF(I)/PI4 
230	   CONTINUE

       RETURN 
       END 


       SUBROUTINE GREEN1(DX,DY,Z,Z0,XFF) !,D,V,XF)
IMPLICIT REAL*8(A-H, O-Z)
! XF(2) :dG/dx;XF(3):dG/dy£»XF(4),dG/dz
! XF(5):dG/dx/dx0; XF(6):dG/dx/dy0;XF(7):dG/dx/dz0
! XF(8):dG/dy/dx0; XF(9):dG/dy/dy0;XF(10):dG/dy/dz0
! XF(11):dG/dz/dx0; XF(12):dG/dz/dy0;XF(13):dG/dz/dz0
REAL*8  XFF(13)  
PI  = 4.0*DATAN(1.0D0) 
PI4 = 4.*PI
DX2 = DX*DX
DY2 = DY*DY
RXY2= DX2 + DY2
RXY = DSQRT(RXY2)
DZ  = Z-Z0
!DZ1 = Z+Z0+2*D  
DZ02= DZ*DZ
!DZ12= DZ1*DZ1
R02=RXY2+DZ02                     
!R12=RXY2+DZ12          
R =DSQRT(R02)
!R1=DSQRT(R12)  

               XFf(:)=0.0d0
	
               
               !	XF(1) =(-1./R -1./R1)/PI4
!	XF(2) =(DX/R**3+ DX/R1**3)/PI4
!	XF(3) =(DY/R**3+ DY /R1**3)/PI4
!	XF(4) =(DZ/R**3+ DZ1/R1**3)/PI4   
!
!	XF(5) =((3.*DX2-R02)/R**5+(3.*DX2-R12)/R1**5)/PI4
!    	XF(6) =(3.*DX*DY*(1./R**5+1./R1**5))/PI4
              	XFf(7) =(3.*DX*(DZ/R**5))/PI4  !-DZ1/R1**5))/PI4

!               XF(8) =(3.*DX*DY*(1./R**5+1./R1**5))/PI4
!               XF(9)=((3.*DY2-R02)/R**5+(3.*DY2-R12) /R1**5)/PI4
               XFf(10)=(3.*DY*(DZ/R**5))/PI4 !-DZ1/R1**5))/PI4

!               XF(11)=(3.*DX*(DZ/R**5+DZ1/R1**5))/PI4
!               XF(12)=(3.*DY*(DZ/R**5+DZ1/R1**5))/PI4
               XFf(13)=((3.*DZ02-R02)/R**5)/PI4 !-(3.*DZ12-R12) /R1**5)/PI4
	RETURN
END 