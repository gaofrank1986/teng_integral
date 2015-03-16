!       MESHBD +MESHT
!


C ***************************************************************
C *                                                             *
C *  Generate nodal and element data on the body surface of     *
C *  an arbitrary body                                          *
C *                                                             *
C ***************************************************************
C 
        SUBROUTINE MESHBD(NCOR)

	  USE MVAR_MOD
        IMPLICIT   NONE  

!	  INTEGER NELEM,NELEMB,NELEMF,NNODE,NNODED,NNB,NNBD,NNF,NNTCH


	  INTEGER NCOR,IPL,IPOLAR(50)
	  INTEGER I,IE,M,INODE,NCNN,K,KK
	  REAL*8  A1,A2,R1,R2,Z1
        REAL*8  XOFSET(50),YOFSET(50),ZOFSET(50)
!
! --------------------------------------------
!
        DO 5 I=1, NCOR
        READ(2,*)  M, IPOLAR(I), XOFSET(I),YOFSET(I),ZOFSET(I)
5       CONTINUE
!
        
        DO 10 INODE=1, NNB
        READ(2,*) M, IPL, R1, A1, Z1
        IF(IPOLAR(IPL).EQ.0) THEN
         XYZ(1,INODE)=R1+XOFSET(IPL)
         XYZ(2,INODE)=A1+YOFSET(IPL)
        ELSE IF(IPOLAR(IPL).EQ.1) THEN
         XYZ(1,INODE)=R1*DCOS(A1*PI/180.0D0)+XOFSET(IPL)
         XYZ(2,INODE)=R1*DSIN(A1*PI/180.0D0)+YOFSET(IPL)
        ENDIF
	   XYZ(3,INODE)=Z1+ZOFSET(IPL)
10      CONTINUE
!       
        DO 11 INODE=1, NNBD
        READ(2,*) M,IPL,R2,A2,DXYZ(3,INODE)
        IF(IPOLAR(IPL).EQ.0) THEN
          DXYZ(1,INODE)= R2
          DXYZ(2,INODE)= A2
        ELSE IF(IPOLAR(IPL).NE.0) THEN
          DXYZ(1,INODE)=R2*DCOS( A2*PI/180.0D0 )
          DXYZ(2,INODE)=R2*DSIN( A2*PI/180.0D0 )
        ENDIF
11      CONTINUE
C
C      
        DO 20 I=1, NELEMB
	  IE=I
	  IETYPE(IE)=1
        READ(2,*) M, NCN(IE)    
        READ(2,*) (NCON(I,K), K=1, NCN(IE))
20    CONTINUE
!
	 DO 30 I=1,  NELEMB
	  IE=I
        READ(2,*)  M, NCNN
        READ(2,*) (NCOND(I,K), K=1, NCN(IE))
30     CONTINUE
!

	   WRITE(10,102)  

         DO INODE=1,NNODE
	   WRITE(10,114)  INODE,XYZ(1,INODE),XYZ(2,INODE),XYZ(3,INODE)
         ENDDO
        
	   WRITE(10,104)  
         DO INODE=1, NNODE
	   WRITE(10,114)  INODE,DXYZ(1,INODE),DXYZ(2,INODE),DXYZ(3,INODE)
         ENDDO

102	  FORMAT( '   INODE       X           Y        Z')
104	  FORMAT( '   INODE       nx          ny       nz')
114	  FORMAT(I6,4(3x,E11.4))

!
       RETURN
       END


C  
C *******************************************************
C *                                                     *
C *  Compute mesh position at T-time                    *
C *                                                     *
C *******************************************************
C 
        SUBROUTINE MESHT

	  USE MVAR_MOD
	  USE PVar_mod
        IMPLICIT   NONE  

	  INTEGER IE,J
C
C    Rotation 
C
        DO 100 IE=1, NELEM
	   DO 80 J=1, NCN(IE)
          TXYZE(1,J,IE)=(TXYZE(1,J,IE)-XC)*TRXYZ(1,1)+
	1			      (TXYZE(2,J,IE)-YC)*TRXYZ(1,2)+
     2                  (TXYZE(3,J,IE)-ZC)*TRXYZ(1,3)   
          TXYZE(2,J,IE)=(TXYZE(1,J,IE)-XC)*TRXYZ(2,1)+
	1	              (TXYZE(2,J,IE)-YC)*TRXYZ(2,2)+
     2                  (TXYZE(3,J,IE)-ZC)*TRXYZ(2,3)   
          TXYZE(3,J,IE)=(TXYZE(1,J,IE)-XC)*TRXYZ(3,1)+
	1	              (TXYZE(2,J,IE)-YC)*TRXYZ(3,2)+
     2                  (TXYZE(3,J,IE)-ZC)*TRXYZ(3,3)        
80      CONTINUE
C
100    CONTINUE
C
C  Translation
C
        DO 200 IE=1, NELEM

	  DO 180 J=1, NCN(IE)
         TXYZE(1,J,IE)=TXYZE(1,J,IE)+XC+DISP(1)
         TXYZE(2,J,IE)=TXYZE(2,J,IE)+YC+DISP(2)  
         TXYZE(3,J,IE)=TXYZE(3,J,IE)+ZC+DISP(3)  
180      CONTINUE
C      
200    CONTINUE
C
      RETURN
      END


!C         
!C ******************************************************
!C *  Attention: the order of rotations                 *
!C *                                                    *
!C ******************************************************
!C
        SUBROUTINE TRMAX
        USE MVar_mod
	  USE PVar_mod
	  IMPLICIT NONE

	  REAL*4 CSAX,SNAX,CSAY,SNAY,CSAZ,SNAZ
        REAL*4 TRX(3,3),TRY(3,3),TRZ(3,3),TRXY(3,3)
!C
        CSAX=COS(DISP(4))
        SNAX=SIN(DISP(4))
        CSAY=COS(DISP(5))
        SNAY=SIN(DISP(5))
        CSAZ=COS(DISP(6))
        SNAZ=SIN(DISP(6))

!                           
        TRXYZ(1,1)= CSAY*CSAZ
        TRXYZ(2,1)= CSAX*SNAZ+SNAX*SNAY*CSAZ
        TRXYZ(3,1)= SNAX*SNAZ-CSAX*SNAY*CSAZ

        TRXYZ(1,2)=-CSAY*SNAZ
        TRXYZ(2,2)= CSAX*CSAZ-SNAX*SNAY*SNAZ
        TRXYZ(3,2)= SNAX*CSAZ+CSAX*SNAY*SNAZ

        TRXYZ(1,3)= SNAY
        TRXYZ(2,3)=-SNAX*CSAY
        TRXYZ(3,3)= CSAX*CSAY
!
        RETURN
        END

