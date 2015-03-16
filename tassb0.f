C  TASSB0+ TINBOD0+ NSBD0
C *********************************************************
C *                                                       *
C * Calculate the element contribution and assembly the   *
C * coefficients of the corresponding system of equationn *
C *                                                       *
C *********************************************************
C
        SUBROUTINE TASSB0
C  
	  USE MVAR_MOD
	  USE PVAR_MOD
	  USE MFUNC_mod
        USE SEBSM_MOD
        IMPLICIT   NONE  
	  INTEGER  INODE,IELEM, J,JNODE,IND,IP,
	1	       I,II,IS,JNCON,KNCON,L
	  REAL*8  XP,YP,ZP,RSN(4,4)
        REAL*8  VALG(4,8),VALDG(4,8)
	  REAL*8  S_ANGLE     
	  REAL*8  DSIGN
	  REAL*8  IN_ANGLE(NNODE)
C
        DATA RSN /1.,  1.,  1.,  1., 
     1            1., -1.,  1., -1.,
     2            1.,  1., -1., -1., 
     3            1., -1., -1.,  1./ 
C
!	  Print *,' IN TASSB0 '
        IN_ANGLE(:)=0.0D0
!        DO   INODE=1,NNODE
!          READ(8,*) IN_ANGLE(INODE)
!	  END DO
        DO 50 INODE=1, NNODE 
        L=0
        DO 40 IELEM=1,  NELEM
        DO 30 J=1,      NCN(IELEM)
        IF(INODE.EQ.NCON(IELEM,J)) THEN
        L=L+1
        NODELE(INODE,L)=IELEM
        NODELJ(INODE,L)=J
        ENDIF
30      CONTINUE
40      CONTINUE
        NODNOE(INODE)=L
C                          
        NODQUA(INODE)=0
        IF( NSYS .GE. 2) THEN
          IF( DABS(XYZ(2,INODE)).LT.1.0E-06 ) THEN
          NODQUA(INODE)=2
          END IF
        END IF
C
        IF( NSYS .EQ. 4) THEN
          IF( DABS(XYZ(1,INODE)).LT.1.0E-06.AND.
     1        DABS(XYZ(2,INODE)).LT.1.0E-06) THEN
           NODQUA(INODE)=5
          ELSE IF( DABS(XYZ(1,INODE)).LT.1.0E-06 ) THEN
           NODQUA(INODE)=4
          ENDIF
        END IF
!
50      CONTINUE
!                       
!
!                
! ***************************************
!
        AMATA(:,:,:)=(0.0D0,0.0D0)
        CMATA(:,:,:)=(0.0D0,0.0D0)
!        LEFT(:,:,:)=(0.0D0,0.0D0)
!	  RIGHT(:,:,:)=(0.0D0,0.0D0)
        
!        DO  500   INODE=1,  NNODE  
           
           INODE=Node_Singular
           

!
	   XP=XYZ(1,INODE)
         YP=XYZ(2,INODE)
         ZP=XYZ(3,INODE) 
        
         line_sum=0.0d0
         Area0_sum=0.0d0
         Area1_sum=0.0d0
         Area2_sum=0.0d0
!	  
!
!        IF(INODE .EQ. 1) STOP
!     	  CALL SOLIDANGLE(INODE,S_ANGLE)     
!        S_ANGLE=IN_ANGLE(INODE)
	   WRITE(9,102)  INODE, XP, YP, ZP
	   WRITE(10,102)  INODE, XP, YP, ZP
      	 
      	 WRITE(17,*) ' ============================================='
      	 WRITE(17,102)  INODE, XP, YP, ZP
        

102	  FORMAT( /,' INODE=',I4,' X=',F10.4,' Y=',F10.4,' Z=',F10.4)
!
!
!
! ** BODY INTEGRAL
!
        DO  400   IELEM=1,  NELEM

        II=0   
C
C Using TSING if the source point is in the element or its mirror
C     elements about any symmetrical axis, otherwise using TINBOD
C
        IF(NODNOE(INODE) .NE. 0) THEN          
         DO 100 I=1, NODNOE(INODE)
         IF(IELEM .EQ. NODELE(INODE,I)) THEN
         II=II+1
         CALL TSING0(INODE,IELEM,NODQUA(INODE),XP,YP,ZP,VALG,VALDG)
         ENDIF
100      CONTINUE
        ENDIF
C
        IF (II .EQ. 0)   THEN 
         CALL TINBOD0(IELEM,XP,YP,ZP,VALG,VALDG)
        END IF                


!
!        DO   J=1,  NCN(IELEM) 
!	   WRITE(6, 214) J, VALDG(1,J)
!	   WRITE(13, 214) J, VALDG(1,J)
!	   ENDDO
214	  FORMAT(I3,' VALDG=',F12.6)
	  
!  --------------------------
!

!
! ---------------------------
!
!
400     CONTINUE
!
        Area_sum=Area0_sum+Area1_sum+Area2_sum

       write(17,*) 
       write(17,*) '线积分之和',line_sum
       write(17,*) '面积分F1之和',Area1_sum
       write(17,*) '面积分F2之和',Area2_sum
       write(17,*) '面积分F0之和',Area0_sum
       write(17,*) '面积分之和',Area_sum
       write(17,*) ' 积分之和',Area_sum+Line_sum
!
!
500     CONTINUE
!
         STOP
! =============================================
!
        IF( NSYS .EQ. 2) THEN
!
	   DO INODE=1, NNF
	    IF(NODQUA(INODE) .EQ. 2) THEN
	     AMATA(INODE,INODE,2)=1.0E20	 
		ENDIF
	   ENDDO
!
	  ELSE IF( NSYS .EQ. 4) THEN
!
	   DO INODE=1, NNF
	    IF(NODQUA(INODE) .EQ. 2) THEN
	     AMATA(INODE,INODE,2)=1.0E20
	     AMATA(INODE,INODE,4)=1.0E20	    
	    ELSE IF(NODQUA(INODE) .EQ. 4) THEN
	     AMATA(INODE,INODE,3)=1.0E20
	     AMATA(INODE,INODE,2)=1.0E20
		ELSE IF(NODQUA(INODE) .EQ. 5) THEN
	     AMATA(INODE,INODE,2)=1.0E20
	     AMATA(INODE,INODE,3)=1.0E20	    
	     AMATA(INODE,INODE,4)=1.0E20	    
		ENDIF
	   ENDDO
!
	  ENDIF
!
! =============================================
!
	 DO 600 IP=1, NSYS
        WRITE(11, *) '  IP=',IP
	  WRITE(11, *) '    INODE=',nNODE,'      AMATA' 
       DO 600   INODE=1,  NNODE
        WRITE(11, *) '  INODE=',INODE
       DO 600   IND= 1,  NNODE !nnf+1,nnf+1
        WRITE(11, 620) IND,XYZ(1,IND),XYZ(2,IND),
	1              XYZ(3,IND),AMATA(INODE,IND,IP)
 600   CONTINUE 
C
C
      ALLOCATE( A_SEBSM(NNODE,NNODE),B_SEBSM(NNODE))
       A_SEBSM(:,:)=AMATA(:,:,1)
       N_SEBSM=NNODE
         
!          call  timestamp
	  
	  DO IP=1, NSYS
          WRITE(6, *) '  IP=',IP,'    Before RLUDCMP'
	    CALL RLUDCMP(IP,AMATA,NNODE,NNODE,NSYS,INDX,DSIGN)  
	  ENDDO
!
c 610     FORMAT(10x,'NNODE=',I6,/,2X,'IND',10X,'AMATA(IND,IND,1)')
 620     FORMAT(1X,I4,2(2X,F13.6,2X,F20.12))                            
C      
C                    
      RETURN
      END




! ======================================================
C
C   Integration on an element without source in itself
C   nor its symmetrical ones
C
! ======================================================
C
        SUBROUTINE TINBOD0(IELEM,XP,YP,ZP,VALG,VALDG)
	  USE MVAR_MOD
        IMPLICIT   NONE 
	  
	  INTEGER  NSAMB,IS,ND,J,NP,IELEM,NCNE
	  REAL*8  XP,YP,ZP 
        REAL*8  VALG(4,8),VALDG(4,8)
C   
C            
	   NCNE=NCN(IELEM)
	                 
c	 Print *,' In TINBOD   PI=',PI             
        NSAMB=16 
        IF(NCNE.EQ.6) NSAMB=4
C
          VALG =(0.0D0,0.0D0)
          VALDG=(0.0D0,0.0D0)
C
        DO 100   IS=1,   NSYS  
          CALL NSBD0(IS,IELEM,NSAMB,NCNE,XP,YP,ZP,VALG,VALDG)
100     CONTINUE
C
        RETURN
        END           

C
! ======================================================
C   Integration on an element with source in itself
C   or its mirror ones about any symmetrical axis
C
! ======================================================
C
        SUBROUTINE TSING0(INODE,IELEM,NUMQUA,XP,YP,ZP,VALG,VALDG)
	  USE MVAR_MOD
	  USE MFUNC_mod
!
        IMPLICIT   NONE  
!
	  INTEGER I,J,IS,IELEM,INODE,NODNUM,ND,NP,NSAMB,NUMQUA
        REAL*8  XP,YP,ZP,XYZT(3,8),DXYZT(3,8)
        REAL*8 VALG(4,8),VALDG(4,8)
!

        DO 5     I=1,  NCN(IELEM)
        XYZT(1, I)  =  XYZ(1,  NCON(IELEM,I))  
        XYZT(2, I)  =  XYZ(2,  NCON(IELEM,I))  
        XYZT(3, I)  =  XYZ(3,  NCON(IELEM,I))
	    
        DXYZT(1, I) = DXYZ(1,  NCOND(IELEM,I))  
        DXYZT(2, I) = DXYZ(2,  NCOND(IELEM,I))  
        DXYZT(3, I) = DXYZ(3,  NCOND(IELEM,I))

        IF(INODE.EQ.NCON(IELEM,I)) NODNUM=I
5       CONTINUE
C
!        CALL TRIPOL(NODNUM,NCN(IELEM),XYZT,DXYZT)
c	  PRINT *,' AFTER TRIPOL'
C
          VALG = 0.0d0  !(0.0D0,0.0D0)
          VALDG= 0.0d0 !(0.0D0,0.0D0)
C
        NSAMB=16
        IF(NCN(IELEM).EQ.6)   NSAMB=4
C
        !write(*,*) 'NUMQUA=',NUMQUA,'nsys=',nsys
        !write(*,*) 'NODNUM=',nodnum
        !pause
        IF(NUMQUA.EQ.0)       THEN
         DO 100 IS=1,  NSYS
          IF( IS.NE.1 ) THEN   
            CALL NSBD0(IS,IELEM,NSAMB,NCN(IELEM),XP,YP,ZP,
     1          VALG,VALDG)
          ELSE IF(IS.EQ.1) THEN 
         !     CALL SGBD0(IS,IELEM,XP,YP,ZP,VALG,VALDG)
              CALL SGBD0_1(IS,INODE,IELEM,NODNUM,XP,YP,ZP,
     1     VALG,VALDG)
           ! write(*,*) 'After Subroutine SGBD0_1'
          END IF
100      CONTINUE
C
        ELSE IF(NUMQUA.EQ.2) THEN
         DO 200 IS=1,NSYS     
          IF(IS.EQ.3.OR.IS.EQ.4) THEN  
            CALL NSBD0(IS,IELEM,NSAMB,NCN(IELEM),XP,YP,ZP,VALG,VALDG)
          ELSE IF(IS.EQ.1.OR.IS.EQ.2) THEN  
              CALL SGBD0(IS,IELEM,XP,YP,ZP,VALG,VALDG)
          END IF
C
200      CONTINUE  
C
        ELSE IF(NUMQUA.EQ.4) THEN
         DO 300  IS=1,  NSYS
          IF(IS.EQ.2.OR.IS.EQ.3) THEN  
            CALL NSBD0(IS,IELEM,NSAMB,NCN(IELEM),XP,YP,ZP,VALG,VALDG)
          ELSE IF(IS.EQ.1.OR.IS.EQ.4) THEN   
              CALL SGBD0(IS,IELEM,XP,YP,ZP,VALG,VALDG)
          END IF
300      CONTINUE
C
        ELSE IF(NUMQUA.EQ.5) THEN
         DO 400 IS=1, NSYS  
              CALL SGBD0(IS,IELEM,XP,YP,ZP,VALG,VALDG)
400      CONTINUE
        ENDIF
C
        RETURN
        END
                           
C ======================================================
C
C Integration on an element without source point
C 
C ======================================================
C                      
        SUBROUTINE NSBD0(IS,IELEM,NSAMB,NCNE,XP,YP,ZP,VALG,VALDG)
	  USE MVAR_MOD
        IMPLICIT   NONE  
 
	  INTEGER IS,IELEM,N,NSAMB,NCNE,J,IP
        REAL*8  XP,YP,ZP,EX(4,4),EY(4,4)
	  REAL*8  X,X0,Y,Y0,Z,Z0       

        REAL*8  DPOX,DPOY,DPOZ,DGN
        REAL*8  VALG(4,8),VALDG(4,8),GXF(4)
C   
        DATA EX /  1.,  1., -1., -1.,  
     1             1.,  1., -1., -1.,
     2            -1., -1.,  1.,  1., 
     3            -1., -1.,  1.,  1./
C                                                  
        DATA EY /  1., -1., -1.,  1.,
     1            -1.,  1.,  1., -1.,
     2            -1.,  1.,  1., -1.,
     3             1., -1., -1.,  1./
C
c	 PRINT *,' IN  NSBD0'
C
        DO 100    N=1,   NSAMB     

	  X =SAMBXY(IELEM,N,1)
	  Y =SAMBXY(IELEM,N,2)
	  Z =SAMBXY(IELEM,N,3)
	  X0=EX(IS,1)*XP
	  Y0=EY(IS,1)*YP
	  Z0= ZP

	  CALL	TGRN (H,X,X0,Y,Y0,Z,Z0,GXF) 
!
!	  CALL	TGRN(H,DLTX,DLTY,DZ,DZ1,GXF) 
!                                                
        DGN=-(GXF(2)*DSAMB(IELEM,N,1)+
	1		GXF(3)*DSAMB(IELEM,N,2)+
     2        GXF(4)*DSAMB(IELEM,N,3))           
C
        DO 20    J=1,   NCNE
        VALG(IS,J) =VALG(IS,J) +GXF(1)*SAMB(IELEM,N,J)
        VALDG(IS,J)=VALDG(IS,J)+   DGN*SAMB(IELEM,N,J)
20      CONTINUE 
C
100     CONTINUE
C
        RETURN
        END
C
C
C Integration on an element with source point
C
        SUBROUTINE SGBD0(IS,IELEM,XP,YP,ZP,VALG,VALDG)
	  USE MVAR_MOD
	  USE TRVAR_MOD
        IMPLICIT   NONE  
!
	  INTEGER IS,IELEM,N,J,IP       
	  REAL*8  XP,YP,ZP,EX(4,4),EY(4,4)
        REAL*8  X,X0,Y,Y0,Z,Z0       
	  REAL*8  DPOX,DPOY,DPOZ,DGN
        REAL*8  VALG(4,8),VALDG(4,8),GXF(4)
C
        DATA EX  /  1., 1.,-1.,-1.,
     1              1., 1.,-1.,-1., 
     2             -1.,-1., 1., 1., 
     3             -1.,-1., 1., 1./
        DATA EY  /  1.,-1.,-1., 1.,
     1             -1., 1., 1.,-1.,
     2             -1., 1., 1.,-1.,
     3              1.,-1.,-1., 1./
C   
!	  write(6,201) XP,YP,ZP
201	  FORMAT(' In SGBD0    XP, YP, ZP=',3F12.4)

 
        DO 130 N=1, NOSAMP
	
	  X =XYNOD(1,N)
	  Y =XYNOD(2,N)
	  Z =XYNOD(3,N)
	  X0=EX(IS,1)*XP
	  Y0=EY(IS,1)*YP
	  Z0= ZP

	  CALL	TGRN (H,X,X0,Y,Y0,Z,Z0,GXF) 
C
!	  write(6,*) ' GXF(1)=',GXF(1)
!	  write(6,*) ' GXF(2)=',GXf(2)
!	  write(6,*) ' GXF(3)=',GXF(3)
!	  write(6,*) ' GXf(4)=',GXF(4)
!
        DGN=-(GXF(2)*DXYNOD(1,N)+GXF(3)*DXYNOD(2,N)+
     1        GXF(4)*DXYNOD(3,N))             

!	  write(6,*) ' DGN=',DGN


        DO 140 J=1, NCN(IELEM)
         VALG(IS,J) =VALG(IS,J) +GXF(1)*SAMNOD(N,J)
         VALDG(IS,J)=VALDG(IS,J)+   DGN*SAMNOD(N,J)

!	   write(6,*) '  J=',J,' SAM=',SAMNOD(N,J),' VALG=',VALG(1,J)

140     CONTINUE
C
!	 pause

130     CONTINUE
C

        RETURN
        END


