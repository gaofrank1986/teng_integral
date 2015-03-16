!
! *************************************************************
! *                                                           *
! *  The source point is in the mesh                          *
! *                                                           *
! *************************************************************
!
        SUBROUTINE SGBD0_1(IS,INODE,IELEM,NODN,XP,YP,ZP,VALG,VALDG) 
        USE MVAR_MOD
        USE TRVar_mod    
        USE  MFUNC_mod   
        
        IMPLICIT NONE
        
        INTEGER IS,IELEM,N,J,IP,NODN   
        integer loop1,loop2,i,NSAMB
        INTEGER INODD,LI,LJ,LK,INODE

	    REAL*8  EX(4,4),EY(4,4)
        REAL*8  XP,YP,ZP,X,Y,Z,r       
	    REAL*8  DPOX,DPOY,DPOZ,DGN
        REAL*8  VALG(4,8),VALDG(4,8),GXF(4)
        REAL*8  XIQSI(8),XIQET(8),XITSI(6),XITET(6),SI,ETA,DUMX
        real*8  Xiq(8),Wiq(8),Xit(7),EtaT(7),WIT(7)
        REAL*8  SF(8),DSF(2,8),DDSF(3,8),DET1,DET2,DET3,DET 
        REAL*8  SISM,ETASM,FITG,SI1,ETA1
        REAL*8  XJ(2,3),XJSM(2,3),XXJ(3,3)
        REAL*8  jk0(3),jk1c(3),jk1s(3),jk1(3),n0,n1c,n1s
        REAL*8  DX,DY,DZ,nx,ny,nz
        real*8  a1,a2,a3,b1,b2,b3,ast,bst,tot,totF,totf_1,totf_2
        real*8  n1,b30,a30,b31,a31,g31,f,f1,f2,gv,xff(13)    
        REAL*8  CSST,SNST,PLO,BETAs,GAMMA
        REAL*8  AS_3,AS_2

!                
	  DATA EX/  1.0d0,  1.0d0, -1.0d0, -1.0d0,    &
                1.0d0,  1.0d0, -1.0d0, -1.0d0,    &
               -1.0d0, -1.0d0,  1.0d0,  1.0d0,    &
               -1.0d0, -1.0d0,  1.0d0,  1.0d0/
!                                                  
	  DATA EY/  1.0d0, -1.0d0, -1.0d0,  1.0d0,    &
               -1.0d0,  1.0d0,  1.0d0, -1.0d0,    &
               -1.0d0,  1.0d0,  1.0d0, -1.0d0,    &
                1.0d0, -1.0d0, -1.0d0,  1.0d0/
!
	!DATA CI  /(0.0D0,1.0D0)/
!
! ** XIQSI AND XIQET: NODE COORDINATES FOR QUADRILATERAL ELEMENTS
      	DATA XIQSI/-1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0,-1.0d0,-1.0d0/
      	DATA XIQET/-1.0d0,-1.0d0,-1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0/    
!
! ** XITSI AND XITET: NODE COORDINATES FOR TRIANGULAR ELEMENTS
!
     	DATA XITSI/0.00d0, 1.00d0, 0.00d0, 0.50d0, 0.50d0, 0.00d0/
      	DATA XITET/0.00d0, 0.00d0, 1.00d0, 0.00d0, 0.50d0, 0.50d0/
      	


! ==================================================================
!
	DATA XIT/0.101286507323456D0,0.797426985353087D0,		&
			 0.101286507323456D0,0.470142064105115D0,		&
			 0.470142064105115D0,0.059715871789770D0,		&
		     0.333333333333333D0/

	DATA ETAT/0.101286507323456D0,0.101286507323456D0,		&
			  0.797426985353087D0,0.059715871789770D0,		&
			  0.470142064105115D0,0.470142064105115D0,		&
			  0.333333333333333D0/

	DATA WIT/0.062969590272414D0,0.062969590272414D0,		&
			 0.062969590272414D0,0.066197076394253D0,		&
			 0.066197076394253D0,0.066197076394253D0,		&
			 0.112500000000000D0/
!
!  -----------------------------------------------

	DATA XIQ/ 0.960289856497536D+00, 0.796666477413626D+00,	&
              0.525532409916329D+00, 0.183434642495650D+00,	&
			 -0.183434642495650D+00,-0.525532409916329D+00,	&
			 -0.796666477413626D+00,-0.960289856497536D+00/
  
  	DATA WIQ/ 0.101228536290376D+00, 0.222381034453374D+00,	&
			  0.313706645877887D+00, 0.362683783378362D+00,	&
			  0.362683783378362D+00, 0.313706645877887D+00,	&
			  0.222381034453374D+00, 0.101228536290376D+00/   	
!    	
!    ============================================    

                 
	 INODD=NCOND(IELEM,NODN)
!
! **** FOR QUADRILATERIAL ELEMENTS  ****************
!    
      write(10,*) 
      write(10,*) '单元',ielem
      write(17,*) '单元',ielem
  
       WRITE(10,*) ' Inside SGBD0_1'
       WRITE(10,*) ' IS,INODE,IELEM,NODN:',IS,INODE,IELEM,NODN
       WRITE(10,*) ' XP,YP,ZP:',XP,YP,ZP
!
      
       IF(NCN(IELEM).EQ.8)  THEN 
!               
        SI =XIQSI(NODN)
	    ETA=XIQET(NODN)
        CALL SPFUNC8_1(SI,ETA,SF,DSF,DDSF) 

       ELSE IF(NCN(IELEM).EQ.6)  THEN
!C
        SI =XITSI(NODN)
	    ETA=XITET(NODN)
        CALL SPFUNC6_1(SI,ETA,SF,DSF,DDSF) 

       ENDIF
!      	
! ----------------------------------------------  (C2)      	
!
      	DO  30  LI=1,2
      	DO  30  LJ=1,3
      	DUMX=0.0D0
      	DO   LK=1, NCN(IELEM)
         DUMX=DUMX+DSF(LI,LK)*XYZ(LJ,NCON(IELEM,LK))
        ENDDO
     	XJ(LI,LJ)=DUMX
30	   CONTINUE
!
       write(10,*) ' XJ(1,1)=',XJ(1,1), ' XJ(2,1)=',XJ(2,1)
       write(10,*) ' XJ(1,2)=',XJ(1,2), ' XJ(2,2)=',XJ(2,2)
       write(10,*) ' XJ(1,3)=',XJ(1,3), ' XJ(2,3)=',XJ(2,3)

!
! ** compute the determinant of the Jacobian matrix at (SI,ETA), DET
!
! **计算Jk0                         ------ (C12)

      	JK0(1)=XJ(1,2)*XJ(2,3)-XJ(1,3)*XJ(2,2)    !  J1
      	JK0(2)=XJ(1,3)*XJ(2,1)-XJ(1,1)*XJ(2,3)    !  J2
      	JK0(3)=XJ(1,1)*XJ(2,2)-XJ(1,2)*XJ(2,1)    !  J3
!  
       write(10,*) ' JK0=',JK0


!  --------------------------------------------  (C3)
!
      DO  35 LI=1,3
      	DO  35 LJ=1,3
      	DUMX=0.0D0
      	DO   LK=1,NCN(IELEM)
      	DUMX=DUMX+DDSF(LI,LK)*XYZ(LJ,NCON(IELEM,LK))   
        ENDDO
   	  XXJ(LI,LJ)=DUMX
35	 CONTINUE

       write(10,41) XXJ(1,1), XXJ(2,1), XXJ(3,1)
       write(10,42) XXJ(1,2), XXJ(2,2), XXJ(3,2)
       write(10,43) XXJ(1,3), XXJ(2,3), XXJ(3,3)

41     FORMAT('XXJ(1,1)=',E15.8, ' XXJ(2,1)=',E15.8, ' XXJ(3,1)=',E15.8)
42     FORMAT('XXJ(1,2)=',E15.8, ' XXJ(2,2)=',E15.8, ' XXJ(3,2)=',E15.8)
43     FORMAT('XXJ(1,3)=',E15.8, ' XXJ(2,3)=',E15.8, ' XXJ(3,3)=',E15.8)


! **计算  JK1=Jk1C*cos（theta）+Jk1S*sin（theta） ------ (C12)
     
       JK1C(1)=(XXJ(1,2)*XJ(2,3)+XJ(1,2)*XXJ(3,3)) -   &
               (XXJ(1,3)*XJ(2,2)+XJ(1,3)*XXJ(3,2))
       JK1C(2)=(XXJ(1,3)*XJ(2,1)+XJ(1,3)*XXJ(3,1)) -   &
               (XXJ(1,1)*XJ(2,3)+XJ(1,1)*XXJ(3,3))
       JK1C(3)=(XXJ(1,1)*XJ(2,2)+XJ(1,1)*XXJ(3,2)) -  &
               (XXJ(1,2)*XJ(2,1)+XJ(1,2)*XXJ(3,1))
              
       JK1S(1)=(XXJ(3,2)*XJ(2,3)+XJ(1,2)*XXJ(2,3)) -   &
               (XXJ(3,3)*XJ(2,2)+XJ(1,3)*XXJ(2,2))
       JK1S(2)=(XXJ(3,3)*XJ(2,1)+XJ(1,3)*XXJ(2,1)) -   &
               (XXJ(3,1)*XJ(2,3)+XJ(1,1)*XXJ(2,3))
       JK1S(3)=(XXJ(3,1)*XJ(2,2)+XJ(1,1)*XXJ(2,2)) -  &
               (XXJ(3,2)*XJ(2,1)+XJ(1,2)*XXJ(2,1))    


       write(10,*) ' JK1C=',JK1C
       write(10,*) ' JK1S=',JK1S

! **计算N0，N1=N1C*cos（theta）+N1S*sin（theta）   (C13)

      n0=0.0d0
      n1c=0.0d0
      n1s=0.0d0
       WRITE(10,*) '   LK,  SF(LK),  DSF(1,LK),   DSF(2,LK)' 

!      DO LK=1,NCN(IELEM)     
!        N0=SF(LK)+N0
!        N1C=DSF(1,LK)+N1C
!        N1S=DSF(2,LK)+N1S

!        N0=SF(NODN)
 !       N1C=DSF(1,NODN)
 !       N1S=DSF(2,NODN)

        N0=1.0d0
        N1C=0.0d0
        N1S=0.0d0


!       WRITE(10,*)   LK, SF(LK),  DSF(1,LK),   DSF(2,LK) 

!      END DO

      
       WRITE(10,*) '  N0=',N0
       WRITE(10,*) '  N1C=',N1C
       WRITE(10,*) '  N1S=',N1S
       

!  **计算单重积分

      	CALL CIRBOD_2(NODN,INODE,IELEM,SI,ETA,FITG,   &
                      JK0,JK1C,JK1S,N0,N1C,N1S,XJ,XXJ)
!                                  
       line_sum=line_sum+fitg

      write(17,*) '        线积分        ',FITG

!
! ======================================================================
!
! ** 1000循环计算双重积分中函数（F2/plo/plo+F1/plo）项积分

       tot=0.0d0    
       totf_1=0.0d0    
       totf_2=0.0d0    
       totF=0.0d0                     
       NSAMB=0

       write(10,*) ' Surface Integration'
 
       IF(NCN(IELEM).EQ.8)  THEN 

      DO 500 LOOP1=1, 8 
      DO 500 LOOP2=1, 8      
       NSAMB=NSAMB+1     
       SISM=XIQ(loop1)
       ETASM=XIQ(loop2)  
       PLO=DSQRT((SISM-SI)*(SISM-SI)+(ETASM-ETA)*(ETASM-ETA))
       CSST=(SISM-SI)/PLO 
       SNST=(ETASM-ETA)/PLO 
        
        CALL AREA_COEF(IELEM,SISM,ETASM,CSST,SNST,JK0,JK1C,JK1S,   &
                   N0,N1C,N1S,XJ,XXJ,X,Y,Z,DET,F1,F2)


        F=F2/PLO/PLO+F1/PLO                       ! (C19)

      
        TOT=TOT+F/PLO*WIQ(LOOP1)*WIQ(LOOP2)
        TOTf_2=TOTf_2+F2/PLO/PLO/PLO*WIQ(LOOP1)*WIQ(LOOP2)
        TOTf_1=TOTf_1+F1/PLO/PLO*WIQ(LOOP1)*WIQ(LOOP2)
         
   
        DX=X-XP
        DY=Y-YP       
        DZ=Z-ZP       
        r=DSQRT(DX**2+DY**2+DZ**2)
        GV=1.0/(r**3)
      
        TOTF=TOTF+GV*DET*WIQ(LOOP1)*WIQ(LOOP2)
                            

500    CONTINUE
!
! ======================================================================
!
! **** FOR TRIANGULAR ELEMENTS **********************
! 
      	ELSE IF(NCN(IELEM).EQ.6)  THEN
!
 

      DO 1000 NSAMB=1, 7

      
       SISM=XIT(NSAMB)
       ETASM=ETAT(NSAMB) 
       PLO=DSQRT((SISM-SI)*(SISM-SI)+(ETASM-ETA)*(ETASM-ETA))
       CSST=(SISM-SI)/PLO 
       SNST=(ETASM-ETA)/PLO 
      
       CALL SPFUNC6_1(SISM,ETASM,SF,DSF,DDSF)

       CALL AREA_COEF(IELEM,SISM,ETASM,CSST,SNST,JK0,JK1C,JK1S,   &
                   N0,N1C,N1S,XJ,XXJ,X,Y,Z,DET,F1,F2)

        F=F2/PLO/PLO+F1/PLO                       ! (C19)

        TOT=TOT+F/PLO*WIT(NSAMB)
        TOTf_2=TOTf_2+F2/PLO/PLO/PLO*WIT(NSAMB)
        TOTf_1=TOTf_1+F1/PLO/PLO*WIT(NSAMB)
 
        DX=X-XP
        DY=Y-YP       
        DZ=Z-ZP       
        r=DSQRT(DX**2+DY**2+DZ**2)
        GV=1.0/(r**3)
        TOTF=TOTF+GV*DET*WIT(NSAMB)

1000    CONTINUE
      
!C
      	ENDIF
!C

         write(17,*)
         write(17,*) '       面积分 F1/plo=    ',totf_1
         write(17,*) '       面积分 F2/plo/plo=',totf_2
         write(17,*) '       面积分 F=',totF
         write(17,*) '       面积分 F-F2/plo/plo+F1/plo=',totF-tot
         write(17,*)

        Area1_sum=Area1_sum-totf_1

        Area2_sum=Area2_sum-totF_2

        Area0_sum=Area0_sum+totF

        write(17,*) '    Area+Line=',totF-tot+FITG
        write(17,*)


      	RETURN
      	END
               
!
!C
!C **************************************************************
!C *                                                            *
!C *  Calculating the one-dimension integral part of the        *
!C *  singular integral                                         * 
!C *  The singularity is on the body surface                    *
!C *     by M. Guiggiani and A. Gigante's Method (Dec. 1990)    *
!C *                                                            *
!C **************************************************************
!C
      	SUBROUTINE CIRBOD_2(NODN,INODE,IELEM,SI,ETA,FITG,   &
                            JK0,JK1C,JK1S,N0,N1C,N1S,XJ,XXJ)
        USE MVAR_MOD
        IMPLICIT NONE 
!
       Integer li,lj,lk
       INTEGER NEWNUM,NODN,IN,IELEM,I,INODE
       INTEGER NEW6(6),NEW8(8)
!
       real*8 XIQSI(32),XIQET(32),Wiq1(17),WIQ2(25)
       REAL*8 XITSI(24),XITET(24),WIT1(9),WIT2(17)
       REAL*8 FITG,FITG1,FITG2
       REAL*8 SI,ETA,SISM,ETASM,SISM1,ETASM1
       REAL*8 CSST,SNST,PLO,BETAs,GAMMA
       REAL*8 A1,A2,A3,AST,B1,B2,B3,BST
       REAL*8 A30,B30,jk0(3),jk1c(3),jk1s(3),jk1(3),n0,n1,n1c,n1s
       REAL*8 A31,B31,AS_3,AS_2
       REAL*8 F1,F2,SUM1,SUM2,G31
       REAL*8 XJ(2,3),XXJ(3,3)
       REAL*8 det1,det2,det3,dumx,sf(8),dsf(2,8),ddsf(3,8)      	 
!    
! ** SAMPLING POINTS AND WEIGHTING FACTORS FOR QUADRILATERAL ELEMENT
!
!   8 nodes at each side
!       
         DATA NEW8/1, 5, 9, 13, 17, 21, 25, 29/
     	 DATA XIQSI/-1.,-.75,-.5,-.25,0.,.25,.5,.75,1.,   &
                   1.,1.,1.,1.,1.,1.,1.,1.,    &
                   .75,.5,.25,0.,-.25,-.5,-.75,-1.,   &
                   -1., -1.,-1.,-1.,-1.,-1.,-1./
         DATA XIQET/-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,  &
                     -.75,-.5,-.25,0.,.25,.5,.75,1.,    &
                     1.,1.,1.,1.,1.,1.,1.,1.,   &
                    .75,.5,.25,0.,-.25,-.5,-.75/ 
      
      
       DATA WIQ1/.0621774973,.122489332,.117207838,.109334473,   & 
                 .099914323,.08992675,.080115342,    &
         .070948527,.066568164,.070948527,.080115342,.08992675,  &
         .099914323,.109334473,.117207838,.122489332,.0621774973/
     

       DATA WIQ2/0.1224893316,0.2318238045,0.1992612228,0.1608752772,0.1262771379,   &
                 0.0986977799,0.0777974140,0.0621774973,0.0801877220,0.1093344729,   &
                 0.1172078379,0.1224893316,0.1243549945,0.1224893316,0.1172078379,   &
                 0.1093344729,0.0801877220,0.0621774973,0.0777974140,0.0986977799,   &
                 0.1262771379,0.1608752772,0.1992612228,0.2318238045,0.1224893316 /


!    
! ** SAMPLING POINTS AND WEIGHTING FACTORS FOR TRIANGULAR ELEMENT
!   8 nodes at each side
!
     	 DATA NEW6/1, 9, 17, 5, 13, 21/     
     	 DATA XITSI/0.0, 0.125, 0.25, 0.375, 0.5, 0.615, 0.75, 0.875,   & 
     	            1.0, 0.875, 0.75, 0.615, 0.5, 0.375, 0.25, 0.125,   &
     	            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  /
         DATA XITET/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,       &
                    0.125, 0.25, 0.375, 0.5, 0.615, 0.75, 0.875, 1.0,   &
                    0.875, 0.75, 0.615, 0.5, 0.375, 0.25, 0.125/
                     
	     DATA WIT1/0.0574805, 0.1212819, 0.1326451, 0.1405174,            &
                   0.1433476, 0.1405174, 0.1326451, 0.1212819, 0.0574805/

	     DATA WIT2/0.1212819, 0.2617994, 0.2810349, 0.2617994, 0.2163447,  &
	               0.1667366, 0.1255904, 0.0950628, 0.0822924, 0.0950628,  &
                   0.1255904, 0.1667366, 0.2163447, 0.2617994, 0.2810349,  &
                   0.2617994, 0.1212819/
!
       FITG=0.0D0
       FITG1=0.0d0
       FITG2=0.0d0
       SUM1=0.0d0
       SUM2=0.0d0          

        WRITE(10,*) ' Inside CIRBOD_2'

!  =====================================================================         
!C **** FOR QUADRILATERIAL ELEMENTS  **************************
!C       
	IF (NCN(IELEM) .EQ. 8) THEN 
       NEWNUM=NEW8(NODN) 
!C                          	     
!C
     IF(NODN==1 .OR. NODN==3 .OR. NODN==5 .OR. NODN==7)  THEN  
!C            
	 DO  200   i=1, 17
	 WRITE(10,*) '  I=1',I
       IN=MOD(NEWNUM+I+6, 32) + 1 
        SISM  = XIQSI(IN) 
     	ETASM = XIQET(IN)   
!C
	 PLO=DSQRT((SISM-SI)*(SISM-SI)+(ETASM-ETA)*(ETASM-ETA))
       CSST=(SISM-SI)/PLO 
       SNST=(ETASM-ETA)/PLO 
       
       JK1(1)=JK1C(1)*CSST+ JK1S(1)*SNST    !(C12)
       JK1(2)=JK1C(2)*CSST+ JK1S(2)*SNST
       JK1(3)=JK1C(3)*CSST+ JK1S(3)*SNST     
   
        WRITE(10,*) ' JK1=',JK1
   
       N1= N1C*CSST+ N1S*SNST               !(C13)
      
	 A1=XJ(1,1)*CSST+XJ(2,1)*SNST
	 A2=XJ(1,2)*CSST+XJ(2,2)*SNST
	 A3=XJ(1,3)*CSST+XJ(2,3)*SNST
      
       B1=XXJ(1,1)*CSST*CSST*0.5+XXJ(3,1)*CSST*SNST +  &
          XXJ(2,1)*SNST*SNST*0.5
       B2=XXJ(1,2)*CSST*CSST*0.5+XXJ(3,2)*CSST*SNST +  &
          XXJ(2,2)*SNST*SNST*0.5
       B3=XXJ(1,3)*CSST*CSST*0.5+XXJ(3,3)*CSST*SNST +  &
          XXJ(2,3)*SNST*SNST*0.5
	
	 AST=DSQRT(A1*A1 + A2*A2 + A3*A3)   
       BST=DSQRT(B1*B1 + B2*B2 + B3*B3)        !(C7)

       AS_3=1.0d0/AST**3
       AS_2=-3.0d0*(A1*B1+A2*B2+A3*B3)/AST**5

       B30=-JK0(3)
       A30=B30*N0
       
       G31=B1*JK0(1)+B2*JK0(2)+B3*JK0(3)+A1*JK1(1)+A2*JK1(2)+A3*JK1(3)
       G31=(A3/AST/AST)*G31
       B31=3*G31-JK1(3)
       A31=B31*N0+B30*N1
      !WRITE(*,*) 'jk0(3)=',jk0(3)

        
        F2=AS_3*JK0(3)
        F1=AS_2*JK0(3)+AS_3*JK1(3)     

!        F1=AS_2*A30+AS_3*A31     
!        F2=AS_3*A30

! --------------------------------------------
        GAMMA=-(A1*B1+A2*B2+A3*B3)/AST**4
        BETAs=1.0d0/AST
  
         WRITE(10,*) ' F1=',F1,' GAMMA=',GAMMA
         WRITE(10,*) ' PLO=',PLO,'  BETAs=',BETAs
         WRITE(10,*) ' DLOG(PLO/BETAs)=',DLOG(PLO/BETAs)
        FITG2=FITG2-F2*(GAMMA/BETAs**2+1.0/PLO)*WIQ1(I) 
        FITG1=FITG1+F1*DLOG(PLO/BETAs)*WIQ1(I) 

         WRITE(10,*) ' FITG1=',FITG1
         WRITE(10,*)
    
! --------------------------------------------
!
!**  线积分中sum1应为0， sum2应为2*Pi    
      
       SUM1=SUM1+F1*WIQ1(I)                      !(25)
       SUM2=SUM2+F2/BETAs*WIQ1(I)                !(26)
!       SUM1=SUM1+F1*DLOG(1.0/AST)*WIQ11(I)                      !(25)
!       SUM2=SUM2+F2*(-1.0*GAMMA/BETAs**2)*WIQ11(I)              !(26)
!
200	 CONTINUE
!
	   ELSE	IF(NODN==2 .OR. NODN==4 .OR. NODN==6 .OR. NODN==8)   THEN  
     
       DO  400   i=1, 25
	   WRITE(10,*) '  I=1',I
       IN=MOD(NEWNUM+I+2, 32) + 1 
       SISM  = XIQSI(IN) 
       ETASM = XIQET(IN)   
!C
	   PLO=DSQRT((SISM-SI)*(SISM-SI)+(ETASM-ETA)*(ETASM-ETA))
       CSST=(SISM-SI)/PLO 
       SNST=(ETASM-ETA)/PLO 
       
       JK1(1)=JK1C(1)*CSST+ JK1S(1)*SNST    !(C12)
       JK1(2)=JK1C(2)*CSST+ JK1S(2)*SNST
       JK1(3)=JK1C(3)*CSST+ JK1S(3)*SNST     
   
        WRITE(10,*) ' JK1=',JK1
   
       N1= N1C*CSST+ N1S*SNST               !(C13)
      
	 A1=XJ(1,1)*CSST+XJ(2,1)*SNST
	 A2=XJ(1,2)*CSST+XJ(2,2)*SNST
	 A3=XJ(1,3)*CSST+XJ(2,3)*SNST
      
       B1=XXJ(1,1)*CSST*CSST*0.5+XXJ(3,1)*CSST*SNST +  &
          XXJ(2,1)*SNST*SNST*0.5
       B2=XXJ(1,2)*CSST*CSST*0.5+XXJ(3,2)*CSST*SNST +  &
          XXJ(2,2)*SNST*SNST*0.5
       B3=XXJ(1,3)*CSST*CSST*0.5+XXJ(3,3)*CSST*SNST +  &
          XXJ(2,3)*SNST*SNST*0.5
	
	   AST=DSQRT(A1*A1 + A2*A2 + A3*A3)   
       BST=DSQRT(B1*B1 + B2*B2 + B3*B3)        !(C7)

       AS_3=1.0d0/AST**3
       AS_2=-3.0d0*(A1*B1+A2*B2+A3*B3)/AST**5

       B30=-JK0(3)
       A30=B30*N0
       
       G31=B1*JK0(1)+B2*JK0(2)+B3*JK0(3)+A1*JK1(1)+A2*JK1(2)+A3*JK1(3)
       G31=(A3/AST/AST)*G31
       B31=3*G31-JK1(3)
       A31=B31*N0+B30*N1
      !WRITE(*,*) 'jk0(3)=',jk0(3)

        
        F2=AS_3*JK0(3)
        F1=AS_2*JK0(3)+AS_3*JK1(3)     

!        F1=AS_2*A30+AS_3*A31     
!        F2=AS_3*A30

! --------------------------------------------
        GAMMA=-(A1*B1+A2*B2+A3*B3)/AST**4
        BETAs=1.0d0/AST
  
         WRITE(10,*) ' F1=',F1,' GAMMA=',GAMMA
         WRITE(10,*) ' PLO=',PLO,'  BETAs=',BETAs
         WRITE(10,*) ' DLOG(PLO/BETAs)=',DLOG(PLO/BETAs)
        FITG2=FITG2-F2*(GAMMA/BETAs**2+1.0/PLO)*WIQ2(I) 
        FITG1=FITG1+F1*DLOG(PLO/BETAs)*WIQ2(I) 

         WRITE(10,*) ' FITG1=',FITG1
         WRITE(10,*)
    
! --------------------------------------------
!
!**  线积分中sum1应为0， sum2应为2*Pi    
      
       SUM1=SUM1+F1*WIQ2(I)                      !(25)
       SUM2=SUM2+F2/BETAs*WIQ2(I)                !(26)
!       SUM1=SUM1+F1*DLOG(1.0/AST)*WIQ11(I)                      !(25)
!       SUM2=SUM2+F2*(-1.0*GAMMA/BETAs**2)*WIQ11(I)              !(26)
!
400	 CONTINUE 
       
       
       END IF
!      
      
!  =====================================================================         
! **** FOR TRIANGULAR ELEMENTS *****************************
!            
	ELSE  IF (NCN(IELEM) .EQ. 6) THEN            
! 
       NEWNUM=NEW6(NODN) 
                     	     
!
     IF(NODN==1 .OR. NODN==2 .OR. NODN==3) THEN  
!            
	 DO  600   i=1, 9
	 WRITE(10,*) '  I=1',I
       IN=MOD(NEWNUM+I+6, 24) + 1 
       SISM  = XITSI(IN) 
       ETASM = XITET(IN)   
!
	   PLO=DSQRT((SISM-SI)*(SISM-SI)+(ETASM-ETA)*(ETASM-ETA))
       CSST=(SISM-SI)/PLO 
       SNST=(ETASM-ETA)/PLO 
       
       JK1(1)=JK1C(1)*CSST+ JK1S(1)*SNST    !(C12)
       JK1(2)=JK1C(2)*CSST+ JK1S(2)*SNST
       JK1(3)=JK1C(3)*CSST+ JK1S(3)*SNST     
   
        WRITE(10,*) ' JK1=',JK1
   
       N1= N1C*CSST+ N1S*SNST               !(C13)
      
	   A1=XJ(1,1)*CSST+XJ(2,1)*SNST
	   A2=XJ(1,2)*CSST+XJ(2,2)*SNST
	   A3=XJ(1,3)*CSST+XJ(2,3)*SNST
      
       B1=XXJ(1,1)*CSST*CSST*0.5+XXJ(3,1)*CSST*SNST +  &
          XXJ(2,1)*SNST*SNST*0.5
       B2=XXJ(1,2)*CSST*CSST*0.5+XXJ(3,2)*CSST*SNST +  &
          XXJ(2,2)*SNST*SNST*0.5
       B3=XXJ(1,3)*CSST*CSST*0.5+XXJ(3,3)*CSST*SNST +  &
          XXJ(2,3)*SNST*SNST*0.5
	
	   AST=DSQRT(A1*A1 + A2*A2 + A3*A3)   
       BST=DSQRT(B1*B1 + B2*B2 + B3*B3)        !(C7)

       AS_3=1.0d0/AST**3
       AS_2=-3.0d0*(A1*B1+A2*B2+A3*B3)/AST**5

       B30=-JK0(3)
       A30=B30*N0
       
       G31=B1*JK0(1)+B2*JK0(2)+B3*JK0(3)+A1*JK1(1)+A2*JK1(2)+A3*JK1(3)
       G31=(A3/AST/AST)*G31
       B31=3*G31-JK1(3)
       A31=B31*N0+B30*N1
      !WRITE(*,*) 'jk0(3)=',jk0(3)

        
        F2=AS_3*JK0(3)
        F1=AS_2*JK0(3)+AS_3*JK1(3)     

!        F1=AS_2*A30+AS_3*A31     
!        F2=AS_3*A30

! --------------------------------------------
        GAMMA=-(A1*B1+A2*B2+A3*B3)/AST**4
        BETAs=1.0d0/AST
  
         WRITE(10,*) ' F1=',F1,' GAMMA=',GAMMA
         WRITE(10,*) ' PLO=',PLO,'  BETAs=',BETAs
         WRITE(10,*) ' DLOG(PLO/BETAs)=',DLOG(PLO/BETAs)
         FITG2=FITG2-F2*(GAMMA/BETAs**2+1.0/PLO)*WIT1(I)
         FITG1=FITG1+F1*DLOG(PLO/BETAs)*WIT1(I) 

         WRITE(10,*) ' FITG1=',FITG1
         WRITE(10,*)
    
! --------------------------------------------
!
!**  线积分中sum1应为0， sum2应为2*Pi    
      
       SUM1=SUM1+F1*WIT1(I)                      !(25)
       SUM2=SUM2+F2/BETAs*WIT1(I)                !(26)
!       SUM1=SUM1+F1*DLOG(1.0/AST)*WIT1(I)                      !(25)
!       SUM2=SUM2+F2*(-1.0*GAMMA/BETAs**2)*WIT1(I)              !(26)
!
600	 CONTINUE
!
	   ELSE	IF(NODN==4 .OR. NODN==5 .OR. NODN==6) THEN  
     
       DO  800   i=1, 17
	   WRITE(10,*) '  I=1',I
       IN=MOD(NEWNUM+I+2, 24) + 1 
       SISM  = XITSI(IN) 
       ETASM = XITET(IN)   
!C
	   PLO=DSQRT((SISM-SI)*(SISM-SI)+(ETASM-ETA)*(ETASM-ETA))
       CSST=(SISM-SI)/PLO 
       SNST=(ETASM-ETA)/PLO 
       
       JK1(1)=JK1C(1)*CSST+ JK1S(1)*SNST    !(C12)
       JK1(2)=JK1C(2)*CSST+ JK1S(2)*SNST
       JK1(3)=JK1C(3)*CSST+ JK1S(3)*SNST     
   
        WRITE(10,*) ' JK1=',JK1
   
       N1= N1C*CSST+ N1S*SNST               !(C13)
      
	 A1=XJ(1,1)*CSST+XJ(2,1)*SNST
	 A2=XJ(1,2)*CSST+XJ(2,2)*SNST
	 A3=XJ(1,3)*CSST+XJ(2,3)*SNST
      
       B1=XXJ(1,1)*CSST*CSST*0.5+XXJ(3,1)*CSST*SNST +  &
          XXJ(2,1)*SNST*SNST*0.5
       B2=XXJ(1,2)*CSST*CSST*0.5+XXJ(3,2)*CSST*SNST +  &
          XXJ(2,2)*SNST*SNST*0.5
       B3=XXJ(1,3)*CSST*CSST*0.5+XXJ(3,3)*CSST*SNST +  &
          XXJ(2,3)*SNST*SNST*0.5
	
	   AST=DSQRT(A1*A1 + A2*A2 + A3*A3)   
       BST=DSQRT(B1*B1 + B2*B2 + B3*B3)        !(C7)

       AS_3=1.0d0/AST**3
       AS_2=-3.0d0*(A1*B1+A2*B2+A3*B3)/AST**5

       B30=-JK0(3)
       A30=B30*N0
       
       G31=B1*JK0(1)+B2*JK0(2)+B3*JK0(3)+A1*JK1(1)+A2*JK1(2)+A3*JK1(3)
       G31=(A3/AST/AST)*G31
       B31=3*G31-JK1(3)
       A31=B31*N0+B30*N1
      !WRITE(*,*) 'jk0(3)=',jk0(3)

        
        F2=AS_3*JK0(3)
        F1=AS_2*JK0(3)+AS_3*JK1(3)     

!        F1=AS_2*A30+AS_3*A31     
!        F2=AS_3*A30

! --------------------------------------------
        GAMMA=-(A1*B1+A2*B2+A3*B3)/AST**4
        BETAs=1.0d0/AST
  
         WRITE(10,*) ' F1=',F1,' GAMMA=',GAMMA
         WRITE(10,*) ' PLO=',PLO,'  BETAs=',BETAs
         WRITE(10,*) ' DLOG(PLO/BETAs)=',DLOG(PLO/BETAs)
        FITG2=FITG2-F2*(GAMMA/BETAs**2+1.0/PLO)*WIT2(I) 
        FITG1=FITG1+F1*DLOG(PLO/BETAs)*WIT2(I)

         WRITE(10,*) ' FITG1=',FITG1
         WRITE(10,*)
    
! --------------------------------------------
!
!**  线积分中sum1应为0， sum2应为2*Pi    
      
       SUM1=SUM1+F1*WIT2(I)                      !(25)
       SUM2=SUM2+F2/BETAs*WIT2(I)                !(26)
!       SUM1=SUM1+F1*DLOG(1.0/AST)*WIT2(I)                      !(25)
!       SUM2=SUM2+F2*(-1.0*GAMMA/BETAs**2)*WIT2(I)              !(26)
!
800	 CONTINUE 
       
       
       END IF
!      
             
	END IF             
!
!  =====================================================================         

       FITG=FITG2+FITG1
       SUMF_1=SUMF_1+SUM1
       SUMF_2=SUMF_2+SUM2
       
       write(10,*) ' Inside CIRBOD_2   ftig=', fitg   
       WRITE(10,*) ' FITG1= ',FITG1,'   FITG2=',FITG2
!       WRITE(17,*) ' FITG1= ',FITG1,'   FITG2=',FITG2
       WRITE(10,*) ' SUM1= ',SUM1,'   SUM2=',SUM2




      	RETURN
      	END




!
! ================================================================================
!
!

        SUBROUTINE AREA_COEF(IELEM,SISM,ETASM,CSST,SNST,JK0,JK1C,JK1S,   &
                   N0,N1C,N1S,XJ,XXJ,X,Y,Z,DET,F1,F2)
!	Input:	IELEM--element src point in
!		SISM---src point ksi 
!		ETASM--src point eta
!		CSST --
!		SNST --
!		JK0  --
!		JK1C --
!		JK1S --
!		N0   --
!		N1S --
!		N1C --
!		XJ --
!		XXJ ---
!	OUTPUT:	X,Y,Z
!		DET
!		F1,F2


        USE MVAR_MOD
        USE  MFUNC_mod   

        IMPLICIT NONE 
!
       Integer li,lj,lk
       INTEGER IN,IELEM,I
!
       REAL*8 SI,ETA,SISM,ETASM
       REAL*8 CSST,SNST,BETAs,GAMMA
       REAL*8 jk0(3),jk1c(3),jk1s(3),jk1(3),n0,n1,n1c,n1s
!
        REAL*8  X,Y,Z       
	    REAL*8  DPOX,DPOY,DPOZ,DGN
        REAL*8  VALG(4,8),VALDG(4,8),GXF(4)
        REAL*8  XIQSI(8),XIQET(8),XITSI(6),XITET(6)
        REAL*8  SF(8),DSF(2,8),DDSF(3,8),DET1,DET2,DET3,DET,DUMX 
        REAL*8  XJ(2,3),XJp(2,3),XJSM(2,3),XXJ(3,3)
        REAL*8  DX,DY,DZ
        real*8  a1,a2,a3,b1,b2,b3,ast,bst,tot,totF,totf_1,totf_2
        real*8  b30,a30,b31,a31,g31,f,f1,f2,gv,xff(13)    
        REAL*8  AS_3,AS_2
!

! ==================================================================
!	File 10:	OUTPUT.TXT
!
       write(10,*) ' Inside AREA_COEF'
       WRITE(10,*) ' SISM,ETASM=',SISM,ETASM
       
       
	   IF (NCN(IELEM) .EQ. 8) THEN            
         CALL SPFUNC8_1(SISM,ETASM,SF,DSF,DDSF)
	   ELSE IF (NCN(IELEM) .EQ. 6) THEN            
        CALL SPFUNC6_1(SISM,ETASM,SF,DSF,DDSF)
       ENDIF
!

      X=0.0D0
      Y=0.0D0
      Z=0.0D0
      DO  LK=1,  NCN(IELEM)
        X=X+SF(LK)*XYZ(1,NCON(IELEM,LK))    
        Y=Y+SF(LK)*XYZ(2,NCON(IELEM,LK))
        Z=Z+SF(LK)*XYZ(3,NCON(IELEM,LK))
      END DO
      
!     
       WRITE(10,*) '   X,Y,Z=',X,Y,Z

!           
       JK1(1)=JK1C(1)*CSST+ JK1S(1)*SNST    !(C12)
       JK1(2)=JK1C(2)*CSST+ JK1S(2)*SNST
       JK1(3)=JK1C(3)*CSST+ JK1S(3)*SNST     
 
        write(10,*) ' JK1=',JK1
        
       N1= N1C*CSST+ N1S*SNST               !(C13)
!       
!   -----------------(C6)
!             
	 A1=XJ(1,1)*CSST+XJ(2,1)*SNST
	 A2=XJ(1,2)*CSST+XJ(2,2)*SNST
	 A3=XJ(1,3)*CSST+XJ(2,3)*SNST
       
       B1=XXJ(1,1)*CSST*CSST*0.5+XXJ(3,1)*CSST*SNST +   &
          XXJ(2,1)*SNST*SNST*0.5
       B2=XXJ(1,2)*CSST*CSST*0.5+XXJ(3,2)*CSST*SNST +   &
          XXJ(2,2)*SNST*SNST*0.5
       B3=XXJ(1,3)*CSST*CSST*0.5+XXJ(3,3)*CSST*SNST +   &
          XXJ(2,3)*SNST*SNST*0.5
	

	   AST=DSQRT(A1*A1 + A2*A2 + A3*A3)         !  (C7)  
       BST=DSQRT(B1*B1 + B2*B2 + B3*B3)

       AS_3=1.0d0/AST**3                        !  (C10)
       AS_2=-3.0d0*(A1*B1+A2*B2+A3*B3)/AST**5


       G31=B1*JK0(1)+B2*JK0(2)+B3*JK0(3)+A1*JK1(1)+A2*JK1(2)+A3*JK1(3)
       G31=(A3/AST/AST)*G31                    !  (C16)

       B30=-JK0(3)                             !  (C17)            
       B31=3*G31-JK1(3)

       A30=B30*N0                              !  (C18)
       A31=B31*N0+B30*N1

        F2=AS_3*JK0(3)
        F1=AS_2*JK0(3)+AS_3*JK1(3)     
      
!        F1=AS_2*A30+AS_3*A31     
!        F2=AS_3*A30
      
      
      !F=F/PI4
      
!
!  ---------------------------------------------------
!
!       
    
       	DO  230  LI=1,2
      	DO  230  LJ=1,3
      	DUMX=0.0D0
      	DO    LK=1,NCN(IELEM)
         DUMX=DUMX+DSF(LI,LK)*XYZ(LJ,NCON(IELEM,LK))
        enddo
  	     XJP(LI,LJ)=DUMX
230	   CONTINUE
!
       	DET1=XJP(1,2)*XJP(2,3)-XJP(1,3)*XJP(2,2)    !  J1
      	DET2=XJP(1,3)*XJP(2,1)-XJP(1,1)*XJP(2,3)    !  J2
      	DET3=XJP(1,1)*XJP(2,2)-XJP(1,2)*XJP(2,1)    !  J3
      	DET=DSQRT(DET1*DET1+DET2*DET2+DET3*DET3)
!      
                    
       END
