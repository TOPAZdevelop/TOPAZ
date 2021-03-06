      SUBROUTINE SUG_TTBUG(P1,ANS)
C  
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : u g -> t t~ u g  
C  
C Crossing   1 is u g -> t t~ u g  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NEXTERNAL,   NCOMB,     NCROSS         
      PARAMETER (NEXTERNAL=6, NCOMB= 64, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 UG_TTBUG
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL)
      LOGICAL GOODHEL(NCOMB,NCROSS)
      DATA GOODHEL/THEL*.FALSE./
      DATA NTRY/0/
      DATA (NHEL(IHEL,  1),IHEL=1,6) / -1, -1, -1, -1, -1, -1/
      DATA (NHEL(IHEL,  2),IHEL=1,6) / -1, -1, -1, -1, -1,  1/
      DATA (NHEL(IHEL,  3),IHEL=1,6) / -1, -1, -1, -1,  1, -1/
      DATA (NHEL(IHEL,  4),IHEL=1,6) / -1, -1, -1, -1,  1,  1/
      DATA (NHEL(IHEL,  5),IHEL=1,6) / -1, -1, -1,  1, -1, -1/
      DATA (NHEL(IHEL,  6),IHEL=1,6) / -1, -1, -1,  1, -1,  1/
      DATA (NHEL(IHEL,  7),IHEL=1,6) / -1, -1, -1,  1,  1, -1/
      DATA (NHEL(IHEL,  8),IHEL=1,6) / -1, -1, -1,  1,  1,  1/
      DATA (NHEL(IHEL,  9),IHEL=1,6) / -1, -1,  1, -1, -1, -1/
      DATA (NHEL(IHEL, 10),IHEL=1,6) / -1, -1,  1, -1, -1,  1/
      DATA (NHEL(IHEL, 11),IHEL=1,6) / -1, -1,  1, -1,  1, -1/
      DATA (NHEL(IHEL, 12),IHEL=1,6) / -1, -1,  1, -1,  1,  1/
      DATA (NHEL(IHEL, 13),IHEL=1,6) / -1, -1,  1,  1, -1, -1/
      DATA (NHEL(IHEL, 14),IHEL=1,6) / -1, -1,  1,  1, -1,  1/
      DATA (NHEL(IHEL, 15),IHEL=1,6) / -1, -1,  1,  1,  1, -1/
      DATA (NHEL(IHEL, 16),IHEL=1,6) / -1, -1,  1,  1,  1,  1/
      DATA (NHEL(IHEL, 17),IHEL=1,6) / -1,  1, -1, -1, -1, -1/
      DATA (NHEL(IHEL, 18),IHEL=1,6) / -1,  1, -1, -1, -1,  1/
      DATA (NHEL(IHEL, 19),IHEL=1,6) / -1,  1, -1, -1,  1, -1/
      DATA (NHEL(IHEL, 20),IHEL=1,6) / -1,  1, -1, -1,  1,  1/
      DATA (NHEL(IHEL, 21),IHEL=1,6) / -1,  1, -1,  1, -1, -1/
      DATA (NHEL(IHEL, 22),IHEL=1,6) / -1,  1, -1,  1, -1,  1/
      DATA (NHEL(IHEL, 23),IHEL=1,6) / -1,  1, -1,  1,  1, -1/
      DATA (NHEL(IHEL, 24),IHEL=1,6) / -1,  1, -1,  1,  1,  1/
      DATA (NHEL(IHEL, 25),IHEL=1,6) / -1,  1,  1, -1, -1, -1/
      DATA (NHEL(IHEL, 26),IHEL=1,6) / -1,  1,  1, -1, -1,  1/
      DATA (NHEL(IHEL, 27),IHEL=1,6) / -1,  1,  1, -1,  1, -1/
      DATA (NHEL(IHEL, 28),IHEL=1,6) / -1,  1,  1, -1,  1,  1/
      DATA (NHEL(IHEL, 29),IHEL=1,6) / -1,  1,  1,  1, -1, -1/
      DATA (NHEL(IHEL, 30),IHEL=1,6) / -1,  1,  1,  1, -1,  1/
      DATA (NHEL(IHEL, 31),IHEL=1,6) / -1,  1,  1,  1,  1, -1/
      DATA (NHEL(IHEL, 32),IHEL=1,6) / -1,  1,  1,  1,  1,  1/
      DATA (NHEL(IHEL, 33),IHEL=1,6) /  1, -1, -1, -1, -1, -1/
      DATA (NHEL(IHEL, 34),IHEL=1,6) /  1, -1, -1, -1, -1,  1/
      DATA (NHEL(IHEL, 35),IHEL=1,6) /  1, -1, -1, -1,  1, -1/
      DATA (NHEL(IHEL, 36),IHEL=1,6) /  1, -1, -1, -1,  1,  1/
      DATA (NHEL(IHEL, 37),IHEL=1,6) /  1, -1, -1,  1, -1, -1/
      DATA (NHEL(IHEL, 38),IHEL=1,6) /  1, -1, -1,  1, -1,  1/
      DATA (NHEL(IHEL, 39),IHEL=1,6) /  1, -1, -1,  1,  1, -1/
      DATA (NHEL(IHEL, 40),IHEL=1,6) /  1, -1, -1,  1,  1,  1/
      DATA (NHEL(IHEL, 41),IHEL=1,6) /  1, -1,  1, -1, -1, -1/
      DATA (NHEL(IHEL, 42),IHEL=1,6) /  1, -1,  1, -1, -1,  1/
      DATA (NHEL(IHEL, 43),IHEL=1,6) /  1, -1,  1, -1,  1, -1/
      DATA (NHEL(IHEL, 44),IHEL=1,6) /  1, -1,  1, -1,  1,  1/
      DATA (NHEL(IHEL, 45),IHEL=1,6) /  1, -1,  1,  1, -1, -1/
      DATA (NHEL(IHEL, 46),IHEL=1,6) /  1, -1,  1,  1, -1,  1/
      DATA (NHEL(IHEL, 47),IHEL=1,6) /  1, -1,  1,  1,  1, -1/
      DATA (NHEL(IHEL, 48),IHEL=1,6) /  1, -1,  1,  1,  1,  1/
      DATA (NHEL(IHEL, 49),IHEL=1,6) /  1,  1, -1, -1, -1, -1/
      DATA (NHEL(IHEL, 50),IHEL=1,6) /  1,  1, -1, -1, -1,  1/
      DATA (NHEL(IHEL, 51),IHEL=1,6) /  1,  1, -1, -1,  1, -1/
      DATA (NHEL(IHEL, 52),IHEL=1,6) /  1,  1, -1, -1,  1,  1/
      DATA (NHEL(IHEL, 53),IHEL=1,6) /  1,  1, -1,  1, -1, -1/
      DATA (NHEL(IHEL, 54),IHEL=1,6) /  1,  1, -1,  1, -1,  1/
      DATA (NHEL(IHEL, 55),IHEL=1,6) /  1,  1, -1,  1,  1, -1/
      DATA (NHEL(IHEL, 56),IHEL=1,6) /  1,  1, -1,  1,  1,  1/
      DATA (NHEL(IHEL, 57),IHEL=1,6) /  1,  1,  1, -1, -1, -1/
      DATA (NHEL(IHEL, 58),IHEL=1,6) /  1,  1,  1, -1, -1,  1/
      DATA (NHEL(IHEL, 59),IHEL=1,6) /  1,  1,  1, -1,  1, -1/
      DATA (NHEL(IHEL, 60),IHEL=1,6) /  1,  1,  1, -1,  1,  1/
      DATA (NHEL(IHEL, 61),IHEL=1,6) /  1,  1,  1,  1, -1, -1/
      DATA (NHEL(IHEL, 62),IHEL=1,6) /  1,  1,  1,  1, -1,  1/
      DATA (NHEL(IHEL, 63),IHEL=1,6) /  1,  1,  1,  1,  1, -1/
      DATA (NHEL(IHEL, 64),IHEL=1,6) /  1,  1,  1,  1,  1,  1/
      DATA (  IC(IHEL,  1),IHEL=1,6) /  1,  2,  3,  4,  5,  6/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  96/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
      ANS(IPROC) = 0D0
      DO IHEL=1,NCOMB
          IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
             T=UG_TTBUG(P ,NHEL(1,IHEL),JC(1))            
             ANS(IPROC)=ANS(IPROC)+T
              IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                  GOODHEL(IHEL,IPROC)=.TRUE.
C             WRITE(*,*) IHEL,T
              ENDIF
          ENDIF
      ENDDO
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION UG_TTBUG(P,NHEL,IC)
C  
C FUNCTION GENERATED BY MADGRAPH
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u g -> t t~ u g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN,    NEXTERNAL       
      PARAMETER (NGRAPHS=  38,NEIGEN= 20,NEXTERNAL=6)   
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  36, NCOLOR=  20) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(6,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      INCLUDE 'coupl.inc'
C  
C COLOR DATA
C  
      DATA Denom(1  )/           18/                                       
      DATA (CF(i,1  ),i=1  ,6  ) /    64,   56,   -7,  -16,   -8,   -8/    
      DATA (CF(i,1  ),i=7  ,12 ) /     2,    2,   -7,    1,   -6,   -6/    
      DATA (CF(i,1  ),i=13 ,18 ) /    -6,   12,  -16,    2,    2,   56/    
      DATA (CF(i,1  ),i=19 ,20 ) /    -7,   -7/                            
C               T[5,4]T[3,1,2,6]                                           
      DATA Denom(2  )/           18/                                       
      DATA (CF(i,2  ),i=1  ,6  ) /    56,   64,   -8,   -8,  -16,   -7/    
      DATA (CF(i,2  ),i=7  ,12 ) /     1,   10,    1,    2,   -7,   -7/    
      DATA (CF(i,2  ),i=13 ,18 ) /     2,    2,   -6,   -6,   12,   57/    
      DATA (CF(i,2  ),i=19 ,20 ) /    -6,   -6/                            
C               T[5,4]T[3,1,2,6]                                           
      DATA Denom(3  )/           18/                                       
      DATA (CF(i,3  ),i=1  ,6  ) /    -7,   -8,   64,    1,    2,   -7/    
      DATA (CF(i,3  ),i=7  ,12 ) /    10,    1,   -8,    2,   56,   -7/    
      DATA (CF(i,3  ),i=13 ,18 ) /   -16,    2,   -6,   12,   -6,   -6/    
      DATA (CF(i,3  ),i=19 ,20 ) /    57,   -6/                            
C               T[5,4]T[3,1,6,2]                                           
      DATA Denom(4  )/           18/                                       
      DATA (CF(i,4  ),i=1  ,6  ) /   -16,   -8,    1,   64,   56,    2/    
      DATA (CF(i,4  ),i=7  ,12 ) /    -8,    1,   10,   -7,   -7,   56/    
      DATA (CF(i,4  ),i=13 ,18 ) /     2,  -16,   12,   -6,   -6,   -6/    
      DATA (CF(i,4  ),i=19 ,20 ) /    -6,   57/                            
C               T[5,4,6]T[3,1,2]                                           
      DATA Denom(5  )/           18/                                       
      DATA (CF(i,5  ),i=1  ,6  ) /    -8,  -16,    2,   56,   64,    1/    
      DATA (CF(i,5  ),i=7  ,12 ) /    -7,   -7,    2,   -8,   -6,   57/    
      DATA (CF(i,5  ),i=13 ,18 ) /    -6,   -6,    2,    2,  -16,   -7/    
      DATA (CF(i,5  ),i=19 ,20 ) /    -7,   56/                            
C               T[5,4,6]T[3,1,2]                                           
      DATA Denom(6  )/           18/                                       
      DATA (CF(i,6  ),i=1  ,6  ) /    -8,   -7,   -7,    2,    1,   64/    
      DATA (CF(i,6  ),i=7  ,12 ) /     2,  -16,   56,   -8,   -6,   -6/    
      DATA (CF(i,6  ),i=13 ,18 ) /    57,   -6,   56,   -7,   -7,  -16/    
      DATA (CF(i,6  ),i=19 ,20 ) /     2,    2/                            
C               T[5,4,2]T[3,1,6]                                           
      DATA Denom(7  )/           18/                                       
      DATA (CF(i,7  ),i=1  ,6  ) /     2,    1,   10,   -8,   -7,    2/    
      DATA (CF(i,7  ),i=7  ,12 ) /    64,   -8,    1,   -7,    2,  -16/    
      DATA (CF(i,7  ),i=13 ,18 ) /    -7,   56,   -6,   57,   -6,   -6/    
      DATA (CF(i,7  ),i=19 ,20 ) /    12,   -6/                            
C               T[5,4,2,6]T[3,1]                                           
      DATA Denom(8  )/           18/                                       
      DATA (CF(i,8  ),i=1  ,6  ) /     2,   10,    1,    1,   -7,  -16/    
      DATA (CF(i,8  ),i=7  ,12 ) /    -8,   64,   -8,   56,    2,    2/    
      DATA (CF(i,8  ),i=13 ,18 ) /    -7,   -7,   -6,   -6,   57,   12/    
      DATA (CF(i,8  ),i=19 ,20 ) /    -6,   -6/                            
C               T[5,4,6,2]T[3,1]                                           
      DATA Denom(9  )/           18/                                       
      DATA (CF(i,9  ),i=1  ,6  ) /    -7,    1,   -8,   10,    2,   56/    
      DATA (CF(i,9  ),i=7  ,12 ) /     1,   -8,   64,  -16,  -16,    2/    
      DATA (CF(i,9  ),i=13 ,18 ) /    56,   -7,   57,   -6,   -6,   -6/    
      DATA (CF(i,9  ),i=19 ,20 ) /    -6,   12/                            
C               T[5,4,2]T[3,1,6]                                           
      DATA Denom(10 )/           18/                                       
      DATA (CF(i,10 ),i=1  ,6  ) /     1,    2,    2,   -7,   -8,   -8/    
      DATA (CF(i,10 ),i=7  ,12 ) /    -7,   56,  -16,   64,   12,   -6/    
      DATA (CF(i,10 ),i=13 ,18 ) /    -6,   -6,   -7,   -7,   56,    2/    
      DATA (CF(i,10 ),i=19 ,20 ) /     2,  -16/                            
C               T[3,1]T[5,4,6,2]                                           
      DATA Denom(11 )/           18/                                       
      DATA (CF(i,11 ),i=1  ,6  ) /    -6,   -7,   56,   -7,   -6,   -6/    
      DATA (CF(i,11 ),i=7  ,12 ) /     2,    2,  -16,   12,   64,   -8/    
      DATA (CF(i,11 ),i=13 ,18 ) /    -8,    1,   -7,    2,    2,   -7/    
      DATA (CF(i,11 ),i=19 ,20 ) /    56,  -16/                            
C               T[5,4]T[3,1,6,2]                                           
      DATA Denom(12 )/           18/                                       
      DATA (CF(i,12 ),i=1  ,6  ) /    -6,   -7,   -7,   56,   57,   -6/    
      DATA (CF(i,12 ),i=7  ,12 ) /   -16,    2,    2,   -6,   -8,   64/    
      DATA (CF(i,12 ),i=13 ,18 ) /     1,   -8,    2,   -7,   -7,    2/    
      DATA (CF(i,12 ),i=19 ,20 ) /   -16,   56/                            
C               T[5,4,6]T[3,1,2]                                           
      DATA Denom(13 )/           18/                                       
      DATA (CF(i,13 ),i=1  ,6  ) /    -6,    2,  -16,    2,   -6,   57/    
      DATA (CF(i,13 ),i=7  ,12 ) /    -7,   -7,   56,   -6,   -8,    1/    
      DATA (CF(i,13 ),i=13 ,18 ) /    64,   -8,   56,  -16,    2,   -7/    
      DATA (CF(i,13 ),i=19 ,20 ) /    -7,    2/                            
C               T[5,4,2]T[3,1,6]                                           
      DATA Denom(14 )/           18/                                       
      DATA (CF(i,14 ),i=1  ,6  ) /    12,    2,    2,  -16,   -6,   -6/    
      DATA (CF(i,14 ),i=7  ,12 ) /    56,   -7,   -7,   -6,    1,   -8/    
      DATA (CF(i,14 ),i=13 ,18 ) /    -8,   64,  -16,   56,   -7,    2/    
      DATA (CF(i,14 ),i=19 ,20 ) /     2,   -7/                            
C               T[3,1]T[5,4,2,6]                                           
      DATA Denom(15 )/           18/                                       
      DATA (CF(i,15 ),i=1  ,6  ) /   -16,   -6,   -6,   12,    2,   56/    
      DATA (CF(i,15 ),i=7  ,12 ) /    -6,   -6,   57,   -7,   -7,    2/    
      DATA (CF(i,15 ),i=13 ,18 ) /    56,  -16,   64,   -8,    1,   -8/    
      DATA (CF(i,15 ),i=19 ,20 ) /     1,   10/                            
C               T[3,1,6]T[5,4,2]                                           
      DATA Denom(16 )/           18/                                       
      DATA (CF(i,16 ),i=1  ,6  ) /     2,   -6,   12,   -6,    2,   -7/    
      DATA (CF(i,16 ),i=7  ,12 ) /    57,   -6,   -6,   -7,    2,   -7/    
      DATA (CF(i,16 ),i=13 ,18 ) /   -16,   56,   -8,   64,   -8,    1/    
      DATA (CF(i,16 ),i=19 ,20 ) /    10,    1/                            
C               T[3,1]T[5,4,2,6]                                           
      DATA Denom(17 )/           18/                                       
      DATA (CF(i,17 ),i=1  ,6  ) /     2,   12,   -6,   -6,  -16,   -7/    
      DATA (CF(i,17 ),i=7  ,12 ) /    -6,   57,   -6,   56,    2,   -7/    
      DATA (CF(i,17 ),i=13 ,18 ) /     2,   -7,    1,   -8,   64,   10/    
      DATA (CF(i,17 ),i=19 ,20 ) /     1,   -8/                            
C               T[3,1]T[5,4,6,2]                                           
      DATA Denom(18 )/           18/                                       
      DATA (CF(i,18 ),i=1  ,6  ) /    56,   57,   -6,   -6,   -7,  -16/    
      DATA (CF(i,18 ),i=7  ,12 ) /    -6,   12,   -6,    2,   -7,    2/    
      DATA (CF(i,18 ),i=13 ,18 ) /    -7,    2,   -8,    1,   10,   64/    
      DATA (CF(i,18 ),i=19 ,20 ) /    -8,    1/                            
C               T[5,4]T[3,1,2,6]                                           
      DATA Denom(19 )/           18/                                       
      DATA (CF(i,19 ),i=1  ,6  ) /    -7,   -6,   57,   -6,   -7,    2/    
      DATA (CF(i,19 ),i=7  ,12 ) /    12,   -6,   -6,    2,   56,  -16/    
      DATA (CF(i,19 ),i=13 ,18 ) /    -7,    2,    1,   10,    1,   -8/    
      DATA (CF(i,19 ),i=19 ,20 ) /    64,   -8/                            
C               T[5,4]T[3,1,6,2]                                           
      DATA Denom(20 )/           18/                                       
      DATA (CF(i,20 ),i=1  ,6  ) /    -7,   -6,   -6,   57,   56,    2/    
      DATA (CF(i,20 ),i=7  ,12 ) /    -6,   -6,   12,  -16,  -16,   56/    
      DATA (CF(i,20 ),i=13 ,18 ) /     2,   -7,   10,    1,   -8,    1/    
      DATA (CF(i,20 ),i=19 ,20 ) /    -8,   64/                            
C               T[3,1,2]T[5,4,6]                                           
C ----------
C BEGIN CODE
C ----------
      CALL IXXXXX(P(0,1   ),ZERO ,NHEL(1   ),+1*IC(1   ),W(1,1   ))        
      CALL VXXXXX(P(0,2   ),ZERO ,NHEL(2   ),-1*IC(2   ),W(1,2   ))        
      CALL OXXXXX(P(0,3   ),TMASS ,NHEL(3   ),+1*IC(3   ),W(1,3   ))       
      CALL IXXXXX(P(0,4   ),TMASS ,NHEL(4   ),-1*IC(4   ),W(1,4   ))       
      CALL OXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))        
      CALL VXXXXX(P(0,6   ),ZERO ,NHEL(6   ),+1*IC(6   ),W(1,6   ))        
      CALL FVOXXX(W(1,3   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,7   ))     
      CALL JIOXXX(W(1,4   ),W(1,7   ),GG ,ZERO    ,ZERO    ,W(1,8   ))     
      CALL FVOXXX(W(1,5   ),W(1,8   ),GG ,ZERO    ,ZERO    ,W(1,9   ))     
      CALL IOVXXX(W(1,1   ),W(1,9   ),W(1,6   ),GG ,AMP(1   ))             
      CALL JIOXXX(W(1,1   ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,10  ))     
      CALL JVVXXX(W(1,6   ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,11  ))      
      CALL FVOXXX(W(1,3   ),W(1,11  ),GG ,TMASS   ,TWIDTH  ,W(1,12  ))     
      CALL IOVXXX(W(1,4   ),W(1,12  ),W(1,10  ),GG ,AMP(2   ))             
      CALL FVOXXX(W(1,3   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,13  ))     
      CALL FVOXXX(W(1,13  ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,14  ))     
      CALL IOVXXX(W(1,4   ),W(1,14  ),W(1,10  ),GG ,AMP(3   ))             
      CALL FVIXXX(W(1,4   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,15  ))     
      CALL IOVXXX(W(1,15  ),W(1,7   ),W(1,10  ),GG ,AMP(4   ))             
      CALL FVIXXX(W(1,1   ),W(1,8   ),GG ,ZERO    ,ZERO    ,W(1,16  ))     
      CALL IOVXXX(W(1,16  ),W(1,5   ),W(1,6   ),GG ,AMP(5   ))             
      CALL FVOXXX(W(1,7   ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,17  ))     
      CALL IOVXXX(W(1,4   ),W(1,17  ),W(1,10  ),GG ,AMP(6   ))             
      CALL VVVXXX(W(1,6   ),W(1,10  ),W(1,8   ),G ,AMP(7   ))              
      CALL FVIXXX(W(1,4   ),W(1,2   ),GG ,TMASS   ,TWIDTH  ,W(1,18  ))     
      CALL FVIXXX(W(1,1   ),W(1,6   ),GG ,ZERO    ,ZERO    ,W(1,19  ))     
      CALL JIOXXX(W(1,18  ),W(1,3   ),GG ,ZERO    ,ZERO    ,W(1,20  ))     
      CALL IOVXXX(W(1,19  ),W(1,5   ),W(1,20  ),GG ,AMP(8   ))             
      CALL FVOXXX(W(1,3   ),W(1,10  ),GG ,TMASS   ,TWIDTH  ,W(1,21  ))     
      CALL IOVXXX(W(1,4   ),W(1,21  ),W(1,11  ),GG ,AMP(9   ))             
      CALL IOVXXX(W(1,18  ),W(1,13  ),W(1,10  ),GG ,AMP(10  ))             
      CALL IOVXXX(W(1,15  ),W(1,21  ),W(1,2   ),GG ,AMP(11  ))             
      CALL FVOXXX(W(1,5   ),W(1,6   ),GG ,ZERO    ,ZERO    ,W(1,22  ))     
      CALL IOVXXX(W(1,1   ),W(1,22  ),W(1,20  ),GG ,AMP(12  ))             
      CALL FVIXXX(W(1,18  ),W(1,6   ),GG ,TMASS   ,TWIDTH  ,W(1,23  ))     
      CALL IOVXXX(W(1,23  ),W(1,3   ),W(1,10  ),GG ,AMP(13  ))             
      CALL JVVXXX(W(1,6   ),W(1,10  ),G ,ZERO    ,ZERO    ,W(1,24  ))      
      CALL IOVXXX(W(1,18  ),W(1,3   ),W(1,24  ),GG ,AMP(14  ))             
      CALL JIOXXX(W(1,4   ),W(1,3   ),GG ,ZERO    ,ZERO    ,W(1,25  ))     
      CALL JVVXXX(W(1,25  ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,26  ))      
      CALL IOVXXX(W(1,19  ),W(1,5   ),W(1,26  ),GG ,AMP(15  ))             
      CALL VVVXXX(W(1,10  ),W(1,11  ),W(1,25  ),G ,AMP(16  ))              
      CALL JVVXXX(W(1,10  ),W(1,2   ),G ,ZERO    ,ZERO    ,W(1,27  ))      
      CALL IOVXXX(W(1,4   ),W(1,13  ),W(1,27  ),GG ,AMP(17  ))             
      CALL IOVXXX(W(1,15  ),W(1,3   ),W(1,27  ),GG ,AMP(18  ))             
      CALL IOVXXX(W(1,1   ),W(1,22  ),W(1,26  ),GG ,AMP(19  ))             
      CALL JVVXXX(W(1,6   ),W(1,25  ),G ,ZERO    ,ZERO    ,W(1,28  ))      
      CALL VVVXXX(W(1,10  ),W(1,2   ),W(1,28  ),G ,AMP(20  ))              
      CALL VVVXXX(W(1,24  ),W(1,2   ),W(1,25  ),G ,AMP(21  ))              
      CALL GGGGXX(W(1,10  ),W(1,2   ),W(1,25  ),W(1,6   ),G ,AMP(22  ))    
      CALL GGGGXX(W(1,25  ),W(1,10  ),W(1,2   ),W(1,6   ),G ,AMP(23  ))    
      CALL GGGGXX(W(1,2   ),W(1,25  ),W(1,10  ),W(1,6   ),G ,AMP(24  ))    
      CALL FVOXXX(W(1,5   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,29  ))     
      CALL IOVXXX(W(1,19  ),W(1,29  ),W(1,25  ),GG ,AMP(25  ))             
      CALL FVIXXX(W(1,1   ),W(1,25  ),GG ,ZERO    ,ZERO    ,W(1,30  ))     
      CALL IOVXXX(W(1,30  ),W(1,5   ),W(1,11  ),GG ,AMP(26  ))             
      CALL JIOXXX(W(1,1   ),W(1,29  ),GG ,ZERO    ,ZERO    ,W(1,31  ))     
      CALL IOVXXX(W(1,4   ),W(1,13  ),W(1,31  ),GG ,AMP(27  ))             
      CALL IOVXXX(W(1,15  ),W(1,3   ),W(1,31  ),GG ,AMP(28  ))             
      CALL IOVXXX(W(1,30  ),W(1,22  ),W(1,2   ),GG ,AMP(29  ))             
      CALL IOVXXX(W(1,1   ),W(1,29  ),W(1,28  ),GG ,AMP(30  ))             
      CALL FVOXXX(W(1,29  ),W(1,6   ),GG ,ZERO    ,ZERO    ,W(1,32  ))     
      CALL IOVXXX(W(1,1   ),W(1,32  ),W(1,25  ),GG ,AMP(31  ))             
      CALL FVOXXX(W(1,5   ),W(1,25  ),GG ,ZERO    ,ZERO    ,W(1,33  ))     
      CALL IOVXXX(W(1,19  ),W(1,33  ),W(1,2   ),GG ,AMP(32  ))             
      CALL IOVXXX(W(1,1   ),W(1,33  ),W(1,11  ),GG ,AMP(33  ))             
      CALL FVIXXX(W(1,1   ),W(1,2   ),GG ,ZERO    ,ZERO    ,W(1,34  ))     
      CALL JIOXXX(W(1,34  ),W(1,5   ),GG ,ZERO    ,ZERO    ,W(1,35  ))     
      CALL IOVXXX(W(1,4   ),W(1,13  ),W(1,35  ),GG ,AMP(34  ))             
      CALL IOVXXX(W(1,15  ),W(1,3   ),W(1,35  ),GG ,AMP(35  ))             
      CALL FVIXXX(W(1,34  ),W(1,25  ),GG ,ZERO    ,ZERO    ,W(1,36  ))     
      CALL IOVXXX(W(1,36  ),W(1,5   ),W(1,6   ),GG ,AMP(36  ))             
      CALL VVVXXX(W(1,6   ),W(1,35  ),W(1,25  ),G ,AMP(37  ))              
      CALL IOVXXX(W(1,34  ),W(1,33  ),W(1,6   ),GG ,AMP(38  ))             
      JAMP(   1) = -AMP(   1)-AMP(   7)-AMP(  15)+AMP(  21)
      JAMP(   2) = +AMP(   2)-AMP(   6)-AMP(  16)-AMP(  23)+AMP(  24)
      JAMP(   3) = -AMP(   2)-AMP(   3)+AMP(  16)+AMP(  17)-AMP(  22)
     &             +AMP(  23)
      JAMP(   4) = -AMP(   4)+AMP(  18)+AMP(  22)-AMP(  24)
      JAMP(   5) = -AMP(   5)+AMP(   7)-AMP(  19)-AMP(  21)
      JAMP(   6) = -AMP(   8)-AMP(  14)+AMP(  15)-AMP(  21)
      JAMP(   7) = +AMP(   9)-AMP(  11)+AMP(  16)-AMP(  18)-AMP(  22)
     &             +AMP(  23)
      JAMP(   8) = -AMP(   9)-AMP(  13)-AMP(  16)-AMP(  23)+AMP(  24)
      JAMP(   9) = -AMP(  10)-AMP(  17)+AMP(  22)-AMP(  24)
      JAMP(  10) = -AMP(  12)+AMP(  14)+AMP(  19)+AMP(  21)
      JAMP(  11) = -AMP(  20)-AMP(  34)
      JAMP(  12) = +AMP(  20)-AMP(  35)
      JAMP(  13) = +AMP(  20)-AMP(  27)+AMP(  30)
      JAMP(  14) = -AMP(  20)-AMP(  28)-AMP(  30)
      JAMP(  15) = -AMP(  25)
      JAMP(  16) = +AMP(  26)-AMP(  31)
      JAMP(  17) = -AMP(  26)-AMP(  29)
      JAMP(  18) = -AMP(  32)+AMP(  33)
      JAMP(  19) = -AMP(  33)-AMP(  37)-AMP(  38)
      JAMP(  20) = -AMP(  36)+AMP(  37)
      UG_TTBUG = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          UG_TTBUG =UG_TTBUG+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
