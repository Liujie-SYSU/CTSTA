module misc
contains
  function getroot(N,i)
     integer, dimension(:), intent(in) :: n
     integer, intent(in) :: i
     integer :: getroot
     integer :: tmp
     tmp=N(i)
     do while (N(tmp).ne.tmp)
       tmp=N(tmp)
     end do
     getroot=tmp
     return
   end function getroot
end module misc

module sort
contains
recursive subroutine QsortC(A,idx)
  integer*8, dimension(:), intent(inout) :: A
  integer, dimension(:), intent(inout) :: idx
  integer :: iq

  if(size(A) > 1) then
     call Partition(A,idx,iq)
     call QsortC(A(:iq-1),idx(:iq-1))
     call QsortC(A(iq:),idx(iq:))
  endif
end subroutine QsortC

subroutine Partition(A,idx, marker)
  integer*8, dimension(:) :: A
  integer, dimension(:) :: idx
  integer, intent(out) :: marker
  integer :: i, j
  integer*8 :: tempA
  integer :: tempidx
  integer*8 :: x      ! pivot point

  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) >= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) <= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        tempA = A(i)
        tempidx = idx(i)
        A(i) = A(j)
        idx(i) = idx(j)
        A(j) = tempA
        idx(j) = tempidx
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition
end module sort
!
!
       PROGRAM CTSTA10
!
! A program conducts statistical analysis of a rectangular volume of microtomographic data
!           using percolation theory,
!           and exports the size, shape and orientation of all clusters in the model.
!
! ********************************************************************************************************
!
! The modification from certsta1.f90 to certsta2.f90 is in analysing clusters in different groups. 
!      Certsta1.f90 is in groups of every 10 clusters. When the number of clusters is big, it runs slow.
!      Certsta2.f90 is in groups of variable clusters, 
!                   for the largest clusters, there may be 1 cluster or several cluster in a group,
!                   for smaller cluster, there may be hundreds thousands clusters in a group.
!
! The modification from certsta2.f90 to certsta3.f90 is to add the calculation of square radius of every cluster.
!      The square radius of clusters are exported in the file "*.clus".
!
! The modification from certsta3.f90 to certsta4.f90 is to read file 10 "CERTSTA.DAT" continually
!      until the end of the file. Once an input file name is read from CERTSTA.DAT, the program 
!      begins the calculation, then reads the next input file name, and so on. All allocatable arrays
!      are deallocated when a calculation is finished. 
!
! The modification from certsta4.f90 to certsta5.f90 is to output all percolating clusters in an extra file (.perc).
!      This is specially designed for tree branch samples which include lots of percolating clusters and they are
!      are necessary to be identified individually. 
!  
! The modification from certsta5.f90 to certsta6.f90 is simply to have the dimension of analysed volume
!      in the output file names.
!
! The modification from ctsta6.f90 to ctsta7.f90 is to have different input file: 
!           1) .bin1, 
!           2) .raw, 
!           3) .header and .binary
!!!!!! Be ware !
!      In the .raw data conducted from Paul's code and the .binary data exported from Avizo  
!               the target phase are generally labelled as 1. 
!      In this code, the target phase is necessary to be 0. 
!      Thus the code reverses all labels after reading .raw and .binary !!!!!
!
! The modification from ctsta7.f90 to ctsta8.f90 is to modify the calculation of orientation tensors.
! 1) each component of the orientation tensor is defined as the sqare root of arithmatic mean of 
!      the summation of dyadic product of the site vector a_i and the transposed site vector a_i^T, or
!      sqrt(sum(a_i*a_i^T)/n).
! 2) the output of eigenvalues is replaced by the output of finite dimensions. The finite dimensions
!      are defined as 4 times of the eigenvalues of the new orientation tensor.
! 3) in previous versions, for small cluster with any tmt(i,i)=0 (means there is only one voxel in a 
!      direction), I modified the tensor by force in adding 0.5 in the diagonal components. This way avoided
!      the unreasonable results of isotropy and elongation indices but caused some other problems, such as
!      finite dimensions calculated from the original .out2 file are not correct for these small clusters. 
!      Now the calculation of each orientation tensor is based on sqrt(sum(a_i*a_i^T)/n). For the situation 
!      of tmt(i,i)=0, a very small value (1.0e-4) is added in diagonal components. At the same time, the output
!      of finite dimensions for fdim(i)<=0.1 is modified to fdim(i)=1.0.
!
! The modification from CTSTA8.f90 to CTSTA9.f90 is the allocation of some arrays of calcuating orientation tensors.
!      For a very large cluster with more than 1 billion voxels, it needs several GB memory to load the array.
!      This new version avoids the extra memory demand simply by a "reading and calculating" scheme 
!                   particularly for the very large cluster. 
!
! The improvement from CTSTA9.f90 to CTSTA10.f90 is 1) using some subroutines to replace blocks in main program; 
!      2) the segment of calculating anisotropy of cluster is simplified -- in the previous version there were
!      two similar blocks doing the same thing for the largest cluster and the rest clusters, respectively, which
!      is not necessary; and 3) the output file .ptd is renamed as .psd to match 'pore-size-distribution' or 
!      'particle-size-distribution'.
!       
! There are 7 output files:
!      1) '.out1' includes general features of the model
!      2) '.out2' includes size and shape features of all clusters in the model
!      3) '.out3' includes orientation matrixs of all the clusters in the model
!      4) '.out4' includes normalized orientation matrix of clusters for stereonet projection
!      5) '.clus' -- clusters' file 
!      6) '.psd' -- statistics of pore-size-distribution or particle-size-distribution
!      7) '.perc' -- the file of percolating clusters
!
       use sort
       use misc
       implicit none
!
    INTERFACE
       SUBROUTINE LABELLING(MAXNUM,MAX_X,MAX_Y,MAX_Z,LABEL,LS,ICI,NCLUSTER,NON0CLU,KLASS)
           use sort
           use misc
           integer, dimension(:,:,:), allocatable, INTENT(INOUT) :: LABEL
           integer, dimension(:,:,:), allocatable, INTENT(OUT) :: LS
           INTEGER*8, DIMENSION(MAXNUM), INTENT(OUT) :: ICI
           integer :: N(MAXNUM),KLASS(0:MAXNUM)
           integer, dimension(:), allocatable :: rank, lglb
           integer :: info,neibors,lleft,lbelow,lback,root
           integer :: iper(6), NCLU
           integer :: error,i,j,k,m,IM,inix,iniy,iniz
           INTEGER :: MAX_X,MAX_Y,MAX_Z
           INTEGER :: ncluster,non0clu,loop,nfi,nc
           integer ::  NLARGE,SEED,SUBLARG,outsides0,samelab,outsides1
        END SUBROUTINE LABELLING
    END INTERFACE
!
    INTERFACE
       SUBROUTINE PSD(ISTEP,MAX_X,MAX_Y,MAX_Z,LABEL)
          INTEGER :: MAX_X,MAX_Y,MAX_Z
          INTEGER :: LEN_X,LEN_Y,LEN_Z
          INTEGER :: ISTEP,LLMAX,I,J,K
          INTEGER :: PORE_TH,MPORE_X,MPORE_Y,MPORE_Z,MAPERT
          INTEGER, DIMENSION(:), ALLOCATABLE :: PTX,PTY,PTZ
          INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: LABEL
!          character*80 labfl,outname
       END SUBROUTINE PSD
    END INTERFACE
!
    INTERFACE
	   SUBROUTINE JUDGE_PER(MAX_X,MAX_Y,MAX_Z,LABEL,nlarge,ixperc,iyperc,izperc,outname)
          INTEGER :: MAX_X,MAX_Y,MAX_Z
          integer, dimension(:,:,:), allocatable :: label
          integer, dimension(:), allocatable :: lglbx0, lglby0, lglbz0
          integer, dimension(:), allocatable :: lglbxm, lglbym, lglbzm
          integer, dimension(:), allocatable :: ixperc,  iyperc, izperc
          integer :: SIXPERC,SIYPERC,SIZPERC
          integer :: iper(6), NCLU
          integer :: i,j,k,NC,NLARGE
          character*80 outname
       END SUBROUTINE
    END INTERFACE
!
    INTERFACE
      SUBROUTINE ORIAN_TMT(MAXNUM,INGP,NTURN,NLARGE,ICI,LABEL,RS_SQ,MAX_X,MAX_Y,MAX_Z)
        INTEGER :: MAXNUM, INGP
        INTEGER :: MAX_X,MAX_Y,MAX_Z
        integer, dimension(:,:,:), allocatable :: label
        INTEGER*8, DIMENSION(MAXNUM) :: ICI
        integer :: iper(6), NCLU
        integer*8 , dimension(:), ALLOCATABLE :: NUMCL
        integer*8 , dimension(:), ALLOCATABLE :: XNUM,YNUM,ZNUM
        REAL :: XC,YC,ZC,ev(3),work(8),EANI_SUM,PROIND1,PROIND2,RR,fdim(3)
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: X,Y,Z
        REAL, DIMENSION(:), ALLOCATABLE :: AX,AY,AZ,RS_SQ
        INTEGER*4 :: INDEX1,INDEX2,NUM_ANI
        real :: elmd(6),EGNVLU(3),EANI(3),avssa,efssa
        REAL,DIMENSION(3,3) :: TMT0,TMT,V0,V,FMT
        integer :: error,i,j,k,m,IM,inix,iniy,iniz
        INTEGER :: ncluster,non0clu,loop
        real :: pbden, emu, FI_XL,TR_TMT,temp
        integer ::  NLARGE,SEED,SUBLARG,outsides0,samelab,outsides1
        integer*8 , dimension(:) , ALLOCATABLE ::  outsurf1
        INTEGER :: NTURN, NGP
      END SUBROUTINE
    END INTERFACE
!
! *****************************************************************************
! ***************************  MAIN PROGRAM  **********************************

       INTEGER, PARAMETER ::MAXNUM=30000000,INGP=10000000
       integer*1, dimension(:,:,:), allocatable :: labmat
       integer, dimension(:,:,:), allocatable :: label
       integer, dimension(:,:,:), allocatable :: ls
       INTEGER*8, DIMENSION(MAXNUM) :: ICI
       integer*8 :: FI_XL0,TV,outsurf0,rlen
       INTEGER*8 , dimension(:), ALLOCATABLE :: ACCUPORE
       integer :: info,neibors,lleft,lbelow,lback,root
       integer :: N(MAXNUM),KLASS(0:MAXNUM)
       integer, dimension(:), allocatable :: rank, lglb
       integer, dimension(:), allocatable :: lglbx0, lglby0, lglbz0
       integer, dimension(:), allocatable :: lglbxm, lglbym, lglbzm
       integer :: iper(6), NCLU
       integer*8 , dimension(:), ALLOCATABLE :: NUMCL
       integer*8 , dimension(:), ALLOCATABLE :: XNUM,YNUM,ZNUM
       REAL :: XC,YC,ZC,ev(3),work(8),EANI_SUM,PROIND1,PROIND2,RR,fdim(3)
       INTEGER, DIMENSION(:,:), ALLOCATABLE :: X,Y,Z
       REAL, DIMENSION(:), ALLOCATABLE :: AX,AY,AZ,RS_SQ
       INTEGER*4 :: INDEX1,INDEX2,NUM_ANI
       real :: elmd(6),EGNVLU(3),EANI(3),avssa,efssa
       REAL,DIMENSION(3,3) :: TMT0,TMT,V0,V,FMT
       CHARACTER*80 LABFL,PERFL,DISFL,ANIFL,OUTNAME
       CHARACTER*4 CHX0,CHXM,CHY0,CHYM,CHZ0,CHZM
       integer :: error,i,j,k,m,IM,inix,iniy,iniz
       INTEGER :: MAX_X0,MAX_Y0,MAX_Z0,VOLX0,VOLXM,VOLY0,VOLYM,VOLZ0,VOLZM
       integer*8 :: lnx, lny,lnz
       INTEGER :: MAX_X,MAX_Y,MAX_Z
       INTEGER :: smlvol,llmax,istep,ncluster,non0clu,loop,nfi,nc
       integer, dimension(:), allocatable :: ixperc,  iyperc, izperc
       integer :: SIXPERC,SIYPERC,SIZPERC
       real :: pbden, emu, FI_XL,TR_TMT,temp
       integer ::  NLARGE,SEED,SUBLARG,outsides0,samelab,outsides1
       integer*8 , dimension(:) , ALLOCATABLE ::  outsurf1
       INTEGER :: PORE_TH,MPORE_X,MPORE_Y,MPORE_Z,MAPERT, NTURN, NGP
       INTEGER, DIMENSION(:), ALLOCATABLE :: PTX,PTY,PTZ
       REAL, DIMENSION(:), ALLOCATABLE ::  MAINPORE
       LOGICAL :: fexist1,fexist2,fexist3,fexist4,fexist5
       character*15 :: char, char0
       integer*1 :: nbit
!
       OPEN(10,FILE='CTSTA.DAT',STATUS='OLD')
1      READ(10,*,end=555) LABFL
!
!************************* Block 1 *********************************************
! This block finds the right input file and prepares output files.
!
       inquire(file=trim(LABFL)//'.bin1',exist=fexist1)
       inquire(file=trim(LABFL)//'.header',exist=fexist2)
       inquire(file=trim(LABFL)//'.binary',exist=fexist3)
       inquire(file=trim(LABFL)//'.raw',exist=fexist4)
       if ((fexist1.and.fexist3).or.(fexist1.and.fexist4).or.(fexist3.and.fexist4)) then
         write(*,*) ' There are more than one data files available. Check and keep only one!'
         stop
       endif
!
       if (fexist1) then
! Case 1: input data is from previously created .bin1 format.
         open(1,FILE=trim(LABFL)//'.bin1',action='read',&
         &form='unformatted',iostat=error)
         WRITE(*,*) ' ***** CALCULATION OF FILE ',trim(LABFL),'.bin1' 
         read(1) max_x0
         read(1) max_y0
         read(1) max_z0
         allocate(labmat(max_x0,max_y0,max_z0))
         read(1) labmat
         close(1)
       elseif (fexist2.and.fexist3) then
! Case 2: input data is from Avizo binary .am and separated to .header (ASCII) and .binary 
         OPEN(1,FILE=trim(LABFL)//'.header',STATUS='OLD')
	       DO I=1,3
             READ(1,*)
	       ENDDO
         READ(1,'(a15,a15)') char,char0
         CLOSE(1)
         open(2,file='max',status='unknown')
         write(2,*) char0
         close(2)
         open(2,file='max',status='old')
         read(2,*) MAX_X0,MAX_Y0,MAX_Z0
         close(2,status='delete')
         WRITE(*,*) ' ***** CALCULATION OF FILE ',trim(LABFL),'.binary' 
         lnx=MAX_X0; lny=MAX_Y0; lnz=MAX_Z0
         rlen=lnx*lny*lnz
         allocate(labmat(max_x0,max_y0,max_z0))
         OPEN(1,FILE=trim(LABFL)//'.binary',STATUS='OLD',access='direct',form='Binary',recl=rlen)
         READ(1,rec=1) LABMAT
         CLOSE(1)
         labmat=1-labmat
       elseif (fexist4) then
! Case 3: input data is .raw data with 13-bit header defining model size and data type
         OPEN(1,FILE=trim(LABFL)//'.raw',STATUS='OLD',access='direct',form='Binary',recl=13)
         READ(1,rec=1) nbit,MAX_Z0,MAX_Y0,MAX_X0
         write(*,*) ' nbit,NX,NY,NZ=',nbit,MAX_X0,MAX_Y0,MAX_Z0
         CLOSE(1)
         if (nbit.ne.0.and.nbit.ne.1) stop ' Label-data is not integer*1, cannot read correctly!'
         allocate(labmat(MAX_X0,MAX_Y0,MAX_Z0))
         WRITE(*,*) ' ***** CALCULATION OF FILE ',trim(LABFL),'.raw' 
         lnx=MAX_X0; lny=MAX_Y0; lnz=MAX_Z0
         rlen=lnx*lny*lnz+13
         write(*,*) ' Length of record is:', rlen
         OPEN(1,FILE=trim(LABFL)//'.raw',STATUS='OLD',access='direct',form='Binary',recl=rlen)
         READ(1,rec=1) nbit,MAX_Z0,MAX_Y0,MAX_X0,LABMAT
         CLOSE(1)
         labmat=1-labmat
       else 
         write(*,*) ' Input data file not exist. Check CTSTA.DAT ...'
         stop
       end if
!
       WRITE(*,*) 'READ LABEL DATA PASS.'
	   WRITE(*,*) ' THE SIZE OF THE MODEL IS:',MAX_X0,MAX_Y0,MAX_Z0
       READ(10,*) VOLX0,VOLXM,VOLY0,VOLYM,VOLZ0,VOLZM
!
       if (volx0.lt.10) write(chx0,'(i1)') volx0
       if (volx0.ge.10.and.volx0.lt.100) write(chx0,'(i2)') volx0
       if (volx0.ge.100.and.volx0.lt.1000) write(chx0,'(i3)') volx0
       if (volx0.ge.1000.and.volx0.lt.10000) write(chx0,'(i4)') volx0
!
       if (volxm.lt.10) write(chxm,'(i1)') volxm
       if (volxm.ge.10.and.volxm.lt.100) write(chxm,'(i2)') volxm
       if (volxm.ge.100.and.volxm.lt.1000) write(chxm,'(i3)') volxm
       if (volxm.ge.1000.and.volxm.lt.10000) write(chxm,'(i4)') volxm
!
       if (voly0.lt.10) write(chy0,'(i1)') voly0
       if (voly0.ge.10.and.voly0.lt.100) write(chy0,'(i2)') voly0
       if (voly0.ge.100.and.voly0.lt.1000) write(chy0,'(i3)') voly0
       if (voly0.ge.1000.and.voly0.lt.10000) write(chy0,'(i4)') voly0
!
       if (volym.lt.10) write(chym,'(i1)') volym
       if (volym.ge.10.and.volym.lt.100) write(chym,'(i2)') volym
       if (volym.ge.100.and.volym.lt.1000) write(chym,'(i3)') volym
       if (volym.ge.1000.and.volym.lt.10000) write(chym,'(i4)') volym
!
       if (volz0.lt.10) write(chz0,'(i1)') volz0
       if (volz0.ge.10.and.volz0.lt.100) write(chz0,'(i2)') volz0
       if (volz0.ge.100.and.volz0.lt.1000) write(chz0,'(i3)') volz0
       if (volz0.ge.1000.and.volz0.lt.10000) write(chz0,'(i4)') volz0
!
       if (volzm.lt.10) write(chzm,'(i1)') volzm
       if (volzm.ge.10.and.volzm.lt.100) write(chzm,'(i2)') volzm
       if (volzm.ge.100.and.volzm.lt.1000) write(chzm,'(i3)') volzm
       if (volzm.ge.1000.and.volzm.lt.10000) write(chzm,'(i4)') volzm
!
       outname=trim(labfl)//'('//trim(chx0)//','//trim(chy0)//','//trim(chz0)//')('   &
               //trim(chxm)//','//trim(chym)//','//trim(chzm)//')'
       write(*,*) 'outname=',outname
!
       OPEN(2,FILE=TRIM(outname)//'.out1',STATUS='UNKNOWN')
	   WRITE(2,*) ' THE SIZE OF THE MODEL ',TRIM(LABFL),'  IS:'
       WRITE(2,*) MAX_X0,MAX_Y0,MAX_Z0
       WRITE(2,*) ' THE VOLUME FOR PROBING IS: X = (',VOLX0,',',VOLXM,')'
       WRITE(2,*) '                            Y = (',VOLY0,',',VOLYM,')'
       WRITE(2,*) '                            Z = (',VOLZ0,',',VOLZM,')'
       WRITE(*,*) ' THE VOLUME FOR PROBING IS: X = (',VOLX0,',',VOLXM,')'
       WRITE(*,*) '                            Y = (',VOLY0,',',VOLYM,')'
       WRITE(*,*) '                            Z = (',VOLZ0,',',VOLZM,')'
!
!******************** End of Block 1 *******************************************
!
!************************* Block 2 *********************************************
! This block simply calculates porosity and iso-surface.
!
       MAX_X=VOLXM-VOLX0+1
       MAX_Y=VOLYM-VOLY0+1
       MAX_Z=VOLZM-VOLZ0+1
!	   smlvol=MIN(MAX_X,MAX_Y,MAX_Z)
!	   LLMAX=MAX(MAX_X,MAX_Y,MAX_Z)
!
! INITIALIZATION 
       LNX=MAX_X; LNY=MAX_Y; LNZ=MAX_Z
       TV = LNX * LNY * LNZ
       WRITE(2,*) ' TOTAL VOXELS OF THE VOLUME = ',TV
       WRITE(*,*) ' TOTAL VOXELS OF THE VOLUME = ',TV
       allocate(label(0:MAX_X+1,0:MAX_Y+1,0:MAX_Z+1))
!
       LABEL(1:MAX_X, 1:MAX_Y, 1:MAX_Z)=LABMAT(VOLX0:VOLXM, VOLY0:VOLYM, VOLZ0:VOLZM)
       DEALLOCATE(LABMAT)
!
!C CALCULATING POROSITY OF THE VOLUME (VARIABLE 'FI_XL0')
! Calculating the IsoSurface Area
       write(*,*) ' Calculating the IsoSurface Area ......' 
       LABEL(0, 1:MAX_Y, 1:MAX_Z) = LABEL(1, 1:MAX_Y, 1:MAX_Z)
       LABEL(1:MAX_X, 0, 1:MAX_Z) = LABEL(1:MAX_X, 1, 1:MAX_Z)
       LABEL(1:MAX_X, 1:MAX_Y, 0) = LABEL(1:MAX_X, 1:MAX_Y, 1)
       FI_XL0=0
       OUTSURF0=0
       DO K=1,MAX_Z
         DO J=1,MAX_Y
           DO I=1,MAX_X
             IF (LABEL(I,J,K).EQ.0) FI_XL0=FI_XL0+1
		     OUTSIDES0 = ABS(LABEL(I,J,K)-LABEL(I-1,J,K))&
		             & + ABS(LABEL(I,J,K)-LABEL(I,J-1,K))&
		             & + ABS(LABEL(I,J,K)-LABEL(I,J,K-1))
		     OUTSURF0=OUTSURF0+OUTSIDES0
           ENDDO
         ENDDO
       ENDDO
!
!******************** End of Block 2 *******************************************
!
!************************* Block 3 *********************************************
! This block calculates pore size distribution (outputs the file .psd).
! ISTEP is the interval of slices to be counted
!
         READ(10,*) ISTEP
         write(*,*) ' Calculating and outputing the distribution of pore throat ......'
         OPEN(6,FILE=TRIM(outname)//'.psd',STATUS='UNKNOWN')
	     WRITE(6,*) ' THE SIZE OF THE MODEL ',TRIM(LABFL),'  IS:'
         WRITE(6,*) MAX_X0,MAX_Y0,MAX_Z0
         WRITE(6,*) ' THE VOLUME FOR PROBING IS: X = (',VOLX0,',',VOLXM,')'
         WRITE(6,*) '                            Y = (',VOLY0,',',VOLYM,')'
         WRITE(6,*) '                            Z = (',VOLZ0,',',VOLZM,')'
         WRITE(6,*) '  THE DISTRIBUTION OF PORE THROAT IN THIS VOLUME :'
         CALL PSD(ISTEP,MAX_X,MAX_Y,MAX_Z,LABEL)
!       
!******************** End of Block 3 *******************************************
!
!************************* Block 4 *********************************************
! This block performs labelling clusters, 
!            which is the most important part of the program.
!
         CALL LABELLING(MAXNUM,MAX_X,MAX_Y,MAX_Z,LABEL,LS,ICI,NCLUSTER,NON0CLU,KLASS)
!
!******************** End of Block 4 *******************************************
!
!************************* Block 5 *********************************************
! This block sorts all the clusters,
!            is the second important part of the program.
! NOTE: CALLING QSORTC IN A SUBROUTINE WITH ALLOCATABLE ARRAY(S) MAY CAUSE TROUBLE.
! TO ENSURE A SMOOTH PROCESSING, KEEP THIS PART IN MAIN PROGRAM!
!
        allocate(rank(ncluster),lglb(0:ncluster))
        lglb=0
        do i=1,ncluster
          rank(i) = i
        end do
        call QsortC(ici(:ncluster),rank(:ncluster))
!
!    lglb(some_cluster_index) tells you how many clusters are bigger than this cluster
        DO I=1,NON0CLU
          lglb(rank(i)) = i
        ENDDO
        READ(10,*) smlvol 
        NLARGE=0
        do i=1,non0clu
          IF (ICI(I).GE.smlvol) THEN
            NLARGE=I
          ENDIF
        enddo
! 
! THE LARGEST CLUSTER IS LABELLED AS CLUSTER 1;
! THE SECOND LARGEST CLUSTER IS  CLUSTER 2; 
! AND SO ON ......
!
           DO K=1,MAX_Z
             DO J=1,MAX_Y
               DO I=1,MAX_X
                 if (ls(i,j,k).ne.0) then
                   LABEL(i,j,k) = lglb(klass(ls(i,j,k)))
                 else
                   LABEL(i,j,k) = 0
                 end if
               ENDDO
             ENDDO
           ENDDO
           DEALLOCATE(RANK,LGLB) 
!        write(*,*) '  Sorting clusters and rewriting labels, finished. '
!******************** End of Block 5 *******************************************
!
!************************* Block 6 *********************************************
! Printing block: print to output files and screen for some major parameters.
!
         WRITE(2,*)
         WRITE(2,*) ' Voxels of target fabric in this model is:', FI_XL0
         WRITE(2,*) ' The volume percentage of the target fabric is:', REAL(FI_XL0)/REAL(TV)
         WRITE(2,*)
         IF (REAL(FI_XL0)/REAL(TV).GE.0.5) THEN
           WRITE(2,*) ' Porosity >= 50% !!! Check the target phase is labelled as 0 or 1.'
           WRITE(2,*) ' A easy way to reverse the label is to comment out labmat=1-labmat.'
         ENDIF 
!
         WRITE(2,*) ' Surface of all target fabric:', OUTSURF0
         WRITE(2,*) ' The specific surface area (SSA) is:', REAL(OUTSURF0)/REAL(TV)
         WRITE(2,*)
!
         WRITE(*,*) '       Total number of clusters in this model is:', NON0CLU
         WRITE(*,'(a17,i9,a51,i6)') '        There are', NLARGE,'  large clusters are &
                                     &considered and their voxels >=',smlvol
         WRITE(2,*) ' Total number of clusters in this model is:', NON0CLU
         WRITE(2,*) ' The voxel-number for large cluster is truncated to:', smlvol
         WRITE(2,*) ' Number of large clusters analysed (in .out2 & .out3) is:', NLARGE
         WRITE(2,*)
!
         ALLOCATE (accupore(non0clu), mainpore(non0clu))
         accupore=0
         mainpore=0.0
         accupore(1)=ici(1)
         mainpore(1)=real(accupore(1))/real(fi_xl0)
         do i=2,non0clu
           accupore(i)=accupore(i-1)+ici(i)
           mainpore(i)=real(accupore(i))/real(fi_xl0)
           if ((mainpore(i-1)-0.50).lt.0.0.and.(mainpore(i)-0.50).ge.0.0) then
             write(*,*) '  50% target fabric in clusters 1 to',i,', their voxels >=', ici(i)
             write(2,*) '  50% target fabric in clusters 1 to',i,', their voxels >=', ici(i)
           endif
           if ((mainpore(i-1)-0.60).lt.0.0.and.(mainpore(i)-0.60).ge.0.0) then
             write(*,*) '  60% target fabric in clusters 1 to',i,', their voxels >=', ici(i)
             write(2,*) '  60% target fabric in clusters 1 to',i,', their voxels >=', ici(i)
           endif
           if ((mainpore(i-1)-0.70).lt.0.0.and.(mainpore(i)-0.70).ge.0.0) then
             write(*,*) '  70% target fabric in clusters 1 to',i,', their voxels >=', ici(i)
             write(2,*) '  70% target fabric in clusters 1 to',i,', their voxels >=', ici(i)
           endif
           if ((mainpore(i-1)-0.80).lt.0.0.and.(mainpore(i)-0.80).ge.0.0) then
             write(*,*) '  80% target fabric in clusters 1 to',i,', their voxels >=', ici(i)
             write(2,*) '  80% target fabric in clusters 1 to',i,', their voxels >=', ici(i)
           endif
           if ((mainpore(i-1)-0.90).lt.0.0.and.(mainpore(i)-0.90).ge.0.0) then
             write(*,*) '  90% target fabric in clusters 1 to',i,', their voxels >=', ici(i)
             write(2,*) '  90% target fabric in clusters 1 to',i,', their voxels >=', ici(i)
           endif
           if ((mainpore(i-1)-0.95).lt.0.0.and.(mainpore(i)-0.95).ge.0.0) then
             write(*,*) '  95% target fabric in clusters 1 to',i,', their voxels >=', ici(i)
             write(2,*) '  95% target fabric in clusters 1 to',i,', their voxels >=', ici(i)
           endif
           if ((mainpore(i-1)-0.99).lt.0.0.and.(mainpore(i)-0.99).ge.0.0) then
             write(*,*) '  99% target fabric in clusters 1 to',i,', their voxels >=', ici(i)
             write(2,*) '  99% target fabric in clusters 1 to',i,', their voxels >=', ici(i)
           endif
         enddo
         close(5)
!
!******************** End of Block 6 *******************************************
!
!************************* Block 7 *********************************************
! Judging the percolation of the target volume
!
        CALL JUDGE_PER(MAX_X,MAX_Y,MAX_Z,LABEL,nlarge,ixperc,iyperc,izperc,outname)
!
!        write(*,*) '  Judging percolation, finished. '
!******************** End of Block 7 *******************************************
!
!************************* Block 8 *********************************************
!  CALCULATING THE ANISOTROPY OF LARGE CLUSTERS IN THE MODEL
!
         OPEN(3,FILE=TRIM(outname)//'.out2',STATUS='UNKNOWN')
         WRITE(3,'(A88)') '   Cluster    Voxels   Surface        Dim-1     Dim-2     Dim-3&
                          &   Iso-index  Elong-index' 
         OPEN(4,FILE=TRIM(outname)//'.out3',STATUS='UNKNOWN')
         WRITE(4,*) VOLX0,VOLXM,VOLY0,VOLYM,VOLZ0,VOLZM
         OPEN(5,FILE=TRIM(outname)//'.out4',STATUS='UNKNOWN')
         WRITE(5,*) VOLX0,VOLXM,VOLY0,VOLYM,VOLZ0,VOLZM
!
         write(*,*) '   Analysing all large clusters in groups of clusters ......'
         write(2,*)
         write(2,*) '   Analysing all large clusters in groups of clusters ......'
!
         LABEL(0,1:MAX_Y,1:MAX_Z) = LABEL(1,1:MAX_Y,1:MAX_Z)
         LABEL(1:MAX_X,0,1:MAX_Z) = LABEL(1:MAX_X,1,1:MAX_Z)
         LABEL(1:MAX_X,1:MAX_Y,0) = LABEL(1:MAX_X,1:MAX_Y,1)
         LABEL(MAX_X+1,1:MAX_Y,1:MAX_Z) = LABEL(MAX_X,1:MAX_Y,1:MAX_Z)
         LABEL(1:MAX_X,MAX_Y+1,1:MAX_Z) = LABEL(1:MAX_X,MAX_Y,1:MAX_Z)
         LABEL(1:MAX_X,1:MAX_Y,MAX_Z+1) = LABEL(1:MAX_X,1:MAX_Y,MAX_Z)
!
! Sometimes, the first cluster may be very big (~ 10^9 voxels), and much bigger than other  
!     clusters, then the arrays of X, Y, Z, AX, AY, AZ need tens GB of memory to load. 
!     In this case, we can only process this biggest cluster in one turn and leave other 
!     clusters to be processed later. In the following:
! NGP is the the number of clusters to be dealt with in one turn.
! ICI(i) is the number of voxels included in the i-th cluster.
! INGP (=10000000) is a large integer estimated to calculate NGP. It could be increased.
!
         ALLOCATE (RS_SQ(non0clu))
         NTURN=0
!
        CALL ORIAN_TMT(MAXNUM,INGP,NTURN,NLARGE,ICI,LABEL,RS_SQ,MAX_X,MAX_Y,MAX_Z)
!
        DEALLOCATE(LABEL)
        CLOSE(3)
        CLOSE(4)
        CLOSE(5)
!
!******************** End of Block 8 *******************************************
!
        OPEN(7,FILE=TRIM(outname)//'.clus',STATUS='UNKNOWN')
        WRITE(7,*) NON0CLU
        WRITE(7,*) '         NO.           SIZE           RS^2     ACCUPORE      PERCENT'
        DO I=1,non0clu
          WRITE(7,'(I13,I15,F15.2,I13,F13.5)') I,ICI(I),RS_SQ(I),ACCUPORE(I),MAINPORE(I)
        ENDDO
        DEALLOCATE(ACCUPORE, MAINPORE, RS_SQ)
        CLOSE(7)
!
         write(*,*)
         GOTO 1
!
555      WRITE(*,*) '  ALL DONE !'
         WRITE(2,*)
         WRITE(2,*) '  ALL DONE !'
         CLOSE(2)
!C
       STOP 
       END PROGRAM CTSTA10
!
!
! ******************************************************
       SUBROUTINE LABELLING(MAXNUM,MAX_X,MAX_Y,MAX_Z,LABEL,LS,ICI,NCLUSTER,NON0CLU,KLASS)
!
       use sort
       use misc
!
! Label clusters
! LABEL describes whether the point is porus or not
! LS is the cluster number
! When N(i,j,k)!=LS(i,j,k) the point was labeled as a new cluster before
!      it was realised it was connected to an old cluster.
!
       integer, dimension(:,:,:), allocatable, INTENT(INOUT) :: LABEL
       integer, dimension(:,:,:), allocatable, INTENT(OUT) :: LS
       INTEGER*8, DIMENSION(MAXNUM), INTENT(OUT) :: ICI
       integer :: N(MAXNUM),KLASS(0:MAXNUM)
       integer, dimension(:), allocatable :: rank, lglb
       integer :: info,neibors,lleft,lbelow,lback,root
       integer :: iper(6), NCLU
       integer :: error,i,j,k,m,IM,inix,iniy,iniz
       INTEGER :: MAX_X,MAX_Y,MAX_Z
       INTEGER :: ncluster,non0clu,loop,nfi,nc
       integer :: NLARGE,SEED,SUBLARG,outsides0,samelab,outsides1
!
       write(*,*) ' Labelling clusters ......'
       IF (ALLOCATED(LABEL)) THEN
         WRITE(*,*) ' IN SUB LABELLING LABEL() HAS BEEN ALLOCATED.'
         WRITE(*,*) ' LABEL() SIZE: ', SIZE(LABEL,1),SIZE(LABEL,2),SIZE(LABEL,3)
       ELSE
         allocate(label(0:MAX_X+1,0:MAX_Y+1,0:MAX_Z+1))
         WRITE(*,*) ' ALLOCATING LABEL() SIZE: ', SIZE(LABEL,1),SIZE(LABEL,2),SIZE(LABEL,3)
       ENDIF
       IF (ALLOCATED(LS)) THEN
         WRITE(*,*) ' IN SUB LABELLING LS() HAS BEEN ALLOCATED.'
         WRITE(*,*) ' LS() SIZE: ', SIZE(LS,1),SIZE(LS,2),SIZE(LS,3)
       ELSE
         allocate(ls(1:MAX_X,1:MAX_Y,1:MAX_Z))
         WRITE(*,*) ' ALLOCATING LS() SIZE: ', SIZE(LS,1),SIZE(LS,2),SIZE(LS,3)
       ENDIF
           LABEL(0,1:MAX_Y,1:MAX_Z) = 1
           LABEL(1:MAX_X,0,1:MAX_Z) = 1
           LABEL(1:MAX_X,1:MAX_Y,0) = 1
           LS = MAXNUM
           N = MAXNUM
           KLASS = 0
           ICI = 0
           NCLUSTER=0
!
           DO K=1,MAX_Z
             DO J=1,MAX_Y
               DO I=1,MAX_X
!
! KEY PROCEDURE FOR PERCOLATION ANALYSIS
! LS contains a cluster index for each point where label=0
! N(LS) contains either the cluster index LS, or a lower valued cluster
! index to which the cluster LS is linked.
! ICI contains a count of the number of points with a given cluster index.
                 NEIBORS=LABEL(I-1,J,K)+LABEL(I,J-1,K)+LABEL(I,J,K-1)
                 IF (LABEL(I,J,K).EQ.0) THEN
                   IF (NEIBORS.EQ.3) THEN
                     NCLUSTER=NCLUSTER+1
                     IF (NCLUSTER.GT.MAXNUM) THEN
                       WRITE(*,*) 'Number of cluster exceeds the MAXNUM !!!'
                       STOP
                     ENDIF 
                     LS(I,J,K)=NCLUSTER
                     N(LS(I,J,K))=NCLUSTER
                   ELSE
                     if (label(i,j,k-1).eq.0) then
                        LLEFT=getroot(N,LS(I,J,K-1))
                     else
                        lleft=MAXNUM
                     end if
                     if (label(i,j-1,k).eq.0) then
                       LBELOW=getroot(N,LS(I,J-1,K))
                     else
                       lbelow=MAXNUM
                     end if
                     if (label(i-1,j,k).eq.0) then
                       LBACK=getroot(N,LS(I-1,J,K))
                     else
                       lback=MAXNUM
                     end if
                     M=MIN(LBACK,LBELOW,LLEFT)
                     if (m.eq.MAXNUM) then
                       write(*,*) 'error identifying connected clusters'
		               write(*,*) ' i,j,k,label(i,j,k)=',i,j,k,label(i,j,k)
!                       stop
                     end if
                     LS(I,J,K)=M
! If the given point brings two clusters into contact with each other,
! update N to indicate the cluster linkage.
                     IF (LABEL(I-1,J,K).EQ.0) then
                       root=getroot(N,LS(i-1,j,k))
                       if (root.gt.M) then
                         N(root)=m
                       end if
                     end if
                     IF (LABEL(I,J-1,K).EQ.0) then
                       root=getroot(N,LS(i,j-1,k))
                       if (root.gt.M) then
                         N(root)=m
                       end if
                     ENDIF
                     IF (LABEL(I,J,K-1).EQ.0)then
                       root=getroot(N,LS(i,j,k-1))
                       if (root.gt.M) then
                         N(root)=m
                       end if
                     ENDIF
                   ENDIF
                   ICI(N(LS(I,J,K)))=ICI(N(LS(I,J,K)))+1
                 ENDIF
               ENDDO
             ENDDO
           ENDDO
!
!        write(*,*) '  Step 1 --- labelling clusters,  finished. '
! HOSHEN-KOPELMAN ALGORITHM FOR LABELS OF CLUSTERS
! I think KLASS(N(I))) ends up as the root label for a given label list
! Starting from the highest valued cluster, collapse the count of the
! number of points
!
           NON0CLU=0
           KLASS=0
           DO I=1,NCLUSTER
             IM=I
             DO WHILE (N(IM).NE.IM)
               ICI(N(IM))=ICI(N(IM))+ICI(IM)
               ICI(IM)=0
               IM=N(IM)
             END DO
             KLASS(I)=IM
             IF (ICI(I).NE.0) NON0CLU=NON0CLU+1
           ENDDO
!
        RETURN
        END SUBROUTINE LABELLING
!
! ***********************************************
       SUBROUTINE PSD(ISTEP,MAX_X,MAX_Y,MAX_Z,LABEL)
!
       INTEGER :: MAX_X,MAX_Y,MAX_Z
       INTEGER :: ISTEP,LLMAX,I,J,K
       INTEGER :: PORE_TH,MPORE_X,MPORE_Y,MPORE_Z,MAPERT
       INTEGER, DIMENSION(:), ALLOCATABLE :: PTX,PTY,PTZ
       INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: LABEL
!       character*80 labfl,outname

	   LLMAX=MAX(MAX_X,MAX_Y,MAX_Z)
!
       IF (ALLOCATED(LABEL)) THEN
         WRITE(*,*) ' IN SUB PSD LABEL() HAS BEEN ALLOCATED.'
         WRITE(*,*) ' LABEL() SIZE: ', SIZE(LABEL,1),SIZE(LABEL,2),SIZE(LABEL,3)
       ELSE
         allocate(label(0:LEN_X+1,0:LEN_Y+1,0:LEN_Z+1))
         WRITE(*,*) ' ALLOCATING LABEL() SIZE IN SUB PSD: ', SIZE(LABEL,1),SIZE(LABEL,2),SIZE(LABEL,3)
       ENDIF
       ALLOCATE(PTX(0:LLMAX),PTY(0:LLMAX),PTZ(0:LLMAX))
          PTX=0; PTY=0; PTZ=0
          MPORE_X=0; MPORE_Y=0; MPORE_Z=0
!
           DO K=1,MAX_Z,ISTEP
             DO J=1,MAX_Y
               PORE_TH=0
               DO I=1,MAX_X
                 IF (LABEL(I,J,K).EQ.0) THEN
                   PORE_TH=PORE_TH+1
                   PTX(PORE_TH)=PTX(PORE_TH)+1
                   PTX(PORE_TH-1)=PTX(PORE_TH-1)-1
                   PTX(0)=0
                 ELSE
                   PORE_TH=0
                 ENDIF
                 MPORE_X=MAX(MPORE_X,PORE_TH)
               ENDDO
               PTX(0)=0
             ENDDO
           ENDDO
!
           DO I=1,MAX_X,ISTEP
             DO K=1,MAX_Z
               PORE_TH=0
               DO J=1,MAX_Y
                 IF (LABEL(I,J,K).EQ.0) THEN
                   PORE_TH=PORE_TH+1
                   PTY(PORE_TH)=PTY(PORE_TH)+1
                   PTY(PORE_TH-1)=PTY(PORE_TH-1)-1
                   PTY(0)=0
                 ELSE
                   PORE_TH=0
                 ENDIF
                 MPORE_Y=MAX(MPORE_Y,PORE_TH)
               ENDDO
               PTY(0)=0
             ENDDO
           ENDDO
!
           DO J=1,MAX_Y,ISTEP
             DO I=1,MAX_X
               PORE_TH=0
               DO K=1,MAX_Z
                 IF (LABEL(I,J,K).EQ.0) THEN
                   PORE_TH=PORE_TH+1
                   PTZ(PORE_TH)=PTZ(PORE_TH)+1
                   PTZ(PORE_TH-1)=PTZ(PORE_TH-1)-1
                   PTZ(0)=0
                 ELSE
                   PORE_TH=0
                 ENDIF
                 MPORE_Z=MAX(MPORE_Z,PORE_TH)
               ENDDO
               PTZ(0)=0
             ENDDO
           ENDDO
!
         MAPERT=MAX(MPORE_X,MPORE_Y,MPORE_Z)
         WRITE(6,*) ' PORE_THROAT    X-DIRECTION    Y-DIRECTION    Z-DIRECTION'
         DO I=1,MAPERT
           WRITE(6,'(I12,3I15)') I,PTX(I),PTY(I),PTZ(I)
         ENDDO
         WRITE(6,*)
         WRITE(6,*) ' THE MAXIMAL APERTURES IN X-, Y- AND Z-DIRECTIONS ARE:'
         WRITE(6,'(I27,2I15)') MPORE_X,MPORE_Y,MPORE_Z
         DEALLOCATE(PTX,PTY,PTZ)
         CLOSE(6)
!       
        RETURN
        END SUBROUTINE PSD
!
! *********************************************** 
	    SUBROUTINE JUDGE_PER(MAX_X,MAX_Y,MAX_Z,LABEL,nlarge,ixperc,iyperc,izperc,outname)
        INTEGER :: MAX_X,MAX_Y,MAX_Z
        integer, dimension(:,:,:), allocatable :: label
        integer, dimension(:), allocatable :: lglbx0, lglby0, lglbz0
        integer, dimension(:), allocatable :: lglbxm, lglbym, lglbzm
        integer, dimension(:), allocatable :: ixperc,  iyperc, izperc
        integer :: SIXPERC,SIYPERC,SIZPERC
        integer :: iper(6), NCLU
        integer :: i,j,k,NC,NLARGE
        character*80 outname
!
        allocate(lglbx0(0:NLARGE+1),lglbxm(0:NLARGE+1))
        allocate(lglby0(0:NLARGE+1),lglbym(0:NLARGE+1))
        allocate(lglbz0(0:NLARGE+1),lglbzm(0:NLARGE+1))
        lglbx0=0
        lglbxm=0
        DO K=1,MAX_Z
          DO J=1,MAX_Y
            IF (LABEL(1,J,K).GT.NLARGE) LABEL(1,J,K)=NLARGE+1
            IF (LABEL(MAX_X,J,K).GT.NLARGE) LABEL(MAX_X,J,K)=NLARGE+1    
            IF (LABEL(1,J,K).NE.0) LGLBX0(LABEL(1,J,K))=LGLBX0(LABEL(1,J,K))+1
            IF (LABEL(MAX_X,J,K).NE.0) LGLBXM(LABEL(MAX_X,J,K))=LGLBXM(LABEL(MAX_X,J,K))+1
          ENDDO
        ENDDO
        lglby0=0
        lglbym=0
        DO K=1,MAX_Z
          DO I=1,MAX_X
            IF (LABEL(I,1,K).GT.NLARGE) LABEL(I,1,K)=NLARGE+1
            IF (LABEL(I,MAX_Y,K).GT.NLARGE) LABEL(I,MAX_Y,K)=NLARGE+1 
            IF (LABEL(I,1,K).NE.0) LGLBY0(LABEL(I,1,K))=LGLBY0(LABEL(I,1,K))+1
            IF (LABEL(I,MAX_Y,K).NE.0) LGLBYM(LABEL(I,MAX_Y,K))=LGLBYM(LABEL(I,MAX_Y,K))+1
          ENDDO
        ENDDO
        lglbz0=0
        lglbzm=0
        DO J=1,MAX_Y
          DO I=1,MAX_X
            IF (LABEL(I,J,1).GT.NLARGE) LABEL(I,J,1)=NLARGE+1
            IF (LABEL(I,J,MAX_Z).GT.NLARGE) LABEL(I,J,MAX_Z)=NLARGE+1
            IF (LABEL(I,J,1).NE.0) LGLBZ0(LABEL(I,J,1))=LGLBZ0(LABEL(I,J,1))+1
            IF (LABEL(I,J,MAX_Z).NE.0) LGLBZM(LABEL(I,J,MAX_Z))=LGLBZM(LABEL(I,J,MAX_Z))+1
          ENDDO
        ENDDO
!
        allocate(ixperc(nlarge),iyperc(nlarge),izperc(nlarge))
        IXPERC=0; IYPERC=0; IZPERC=0
        SIXPERC=0; SIYPERC=0; SIZPERC=0
        DO NC=1,NLARGE
          IF (LGLBX0(NC).gt.0 .and. LGLBXM(NC).gt.0) THEN
            IXPERC(NC)=1
            SIXPERC=SIXPERC+1
          ENDIF
        ENDDO
!
        DO NC=1,NLARGE
          IF (LGLBY0(NC).gt.0 .and. LGLBYM(NC).gt.0) THEN
            IYPERC(NC)=1
            SIYPERC=SIYPERC+1
          ENDIF
        ENDDO
!
        DO NC=1,NLARGE
          IF (LGLBZ0(NC).gt.0 .and. LGLBZM(NC).gt.0) THEN
            IZPERC(NC)=1
            SIZPERC=SIZPERC+1
          ENDIF
        ENDDO
!
        WRITE(2,*)
        WRITE(2,*) ' Percolation result in every direction:'
        WRITE(2,*) '    (=0---nonpercolating, >0---number of percolating clusters )'
        WRITE(2,'(A18,2i5)') '      X-direction:',SIXPERC
        WRITE(2,'(A18,2i5)') '      Y-direction:',SIYPERC
        WRITE(2,'(A18,2i5)') '      Z-direction:',SIZPERC
!
        OPEN(6,FILE=TRIM(outname)//'.perc',STATUS='UNKNOWN')
        WRITE(6,*) ' Percolation result in every direction:'
        WRITE(6,*) '    (=0---nonpercolating, >0---number of percolating clusters )'
        WRITE(6,'(A18,2i5)') '      X-direction:',SIXPERC
        WRITE(6,'(A18,2i5)') '      Y-direction:',SIYPERC
        WRITE(6,'(A18,2i5)') '      Z-direction:',SIZPERC
        WRITE(6,*)
!
        WRITE(6,*) ' Percolating clusters include:'
        WRITE(6,*) '    Cluster#      X-perc      Y-perc      Z-perc'
        DO NC=1,NLARGE
          IF ((IXPERC(NC)+IYPERC(NC)+IZPERC(NC)).NE.0) THEN
!            WRITE(2,*) NC, IXPERC(NC),IYPERC(NC),IZPERC(NC)
            WRITE(6,*) NC, IXPERC(NC),IYPERC(NC),IZPERC(NC)
          ENDIF
        ENDDO
!
        deallocate(lglbx0,lglbxm,lglby0,lglbym,lglbz0,lglbzm)
        DEALLOCATE(IXPERC,IYPERC,IZPERC)
        CLOSE(6)
!
         RETURN
         END SUBROUTINE JUDGE_PER
!
!  **********************************************
      SUBROUTINE ORIAN_TMT(MAXNUM,INGP,NTURN,NLARGE,ICI,LABEL,RS_SQ,MAX_X,MAX_Y,MAX_Z)
        INTEGER :: MAXNUM, INGP
        INTEGER :: MAX_X,MAX_Y,MAX_Z
        integer, dimension(:,:,:), allocatable :: label
        INTEGER*8, DIMENSION(MAXNUM) :: ICI
        integer :: iper(6), NCLU
        integer*8 , dimension(:), ALLOCATABLE :: NUMCL
        integer*8 , dimension(:), ALLOCATABLE :: XNUM,YNUM,ZNUM
        REAL :: XC,YC,ZC,ev(3),work(8),EANI_SUM,PROIND1,PROIND2,RR,fdim(3)
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: X,Y,Z
        REAL, DIMENSION(:), ALLOCATABLE :: AX,AY,AZ,RS_SQ
        INTEGER*4 :: INDEX1,INDEX2,NUM_ANI
        real :: elmd(6),EGNVLU(3),EANI(3),avssa,efssa
        REAL,DIMENSION(3,3) :: TMT0,TMT,V0,V,FMT
        integer :: error,i,j,k,m,IM,inix,iniy,iniz
        INTEGER :: ncluster,non0clu,loop
        real :: pbden, emu, FI_XL,TR_TMT,temp
        integer ::  NLARGE,SEED,SUBLARG,outsides0,samelab,outsides1
        integer*8 , dimension(:) , ALLOCATABLE ::  outsurf1
        INTEGER :: NTURN, NGP
!
        IF (ALLOCATED(LABEL)) THEN
          WRITE(*,*) ' IN SUB ORIAN_TMT LABEL() HAS BEEN ALLOCATED.'
          WRITE(*,*) ' LABEL() SIZE: ', SIZE(LABEL,1),SIZE(LABEL,2),SIZE(LABEL,3)
        ELSE
          allocate(label(0:MAX_X+1,0:MAX_Y+1,0:MAX_Z+1))
          WRITE(*,*) ' ALLOCATING LABEL() SIZE: ', SIZE(LABEL,1),SIZE(LABEL,2),SIZE(LABEL,3)
        ENDIF
!
         DO WHILE ((NTURN+1).LE.NLARGE)
           NGP=INT(INGP/ICI(NTURN+1))
           IF (NGP.EQ.0) NGP=1
           IF ((NTURN+NGP).GT.NLARGE) NGP=NLARGE-NTURN
           WRITE(*,*) '       From cluster No.',NTURN+1,'    to cluster No.',NTURN+NGP
           WRITE(2,*) '       From cluster No.',NTURN+1,'    to cluster No.',NTURN+NGP
           ALLOCATE (NUMCL(NTURN+1:NTURN+NGP),OUTSURF1(NTURN+1:NTURN+NGP))
           ALLOCATE (XNUM(NTURN+1:NTURN+NGP),YNUM(NTURN+1:NTURN+NGP),ZNUM(NTURN+1:NTURN+NGP))
!
            ALLOCATE (X(NTURN+1:NTURN+NGP,ICI(NTURN+1)))
            ALLOCATE (Y(NTURN+1:NTURN+NGP,ICI(NTURN+1)))
            ALLOCATE (Z(NTURN+1:NTURN+NGP,ICI(NTURN+1)))
            ALLOCATE (AX(ICI(NTURN+1)),AY(ICI(NTURN+1)),AZ(ICI(NTURN+1)))
            NUMCL=0
            XNUM=0;  YNUM=0;     ZNUM=0
            X=0;     Y=0;        Z=0
            OUTSURF1=0
!           
            DO K=1,MAX_Z
              DO J=1,MAX_Y
                DO I=1,MAX_X
                  NCLU=LABEL(I,J,K)
                  IF (NCLU.GE.(NTURN+1).AND.NCLU.LE.(NTURN+NGP)) THEN
                    NUMCL(NCLU)=NUMCL(NCLU)+1
                    XNUM(NCLU)=XNUM(NCLU)+I
                    YNUM(NCLU)=YNUM(NCLU)+J
                    ZNUM(NCLU)=ZNUM(NCLU)+K
                    X(NCLU,NUMCL(NCLU))=I
                    Y(NCLU,NUMCL(NCLU))=J
                    Z(NCLU,NUMCL(NCLU))=K
                    SAMELAB=(LABEL(I-1,J,K)+LABEL(I+1,J,K)+LABEL(I,J-1,K)+LABEL(I,J+1,K)&
                           &+LABEL(I,J,K-1)+LABEL(I,J,K+1))/LABEL(I,J,K)
		            OUTSIDES1=6-SAMELAB
		            OUTSURF1(NCLU)=OUTSURF1(NCLU)+OUTSIDES1
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
!           write(*,*) '  Assigning arrays of all clusters, finished. '
            DO I= NTURN+1, NTURN+NGP
              TMT=0.0
              EGNVLU=0.0
              FMT=0.0
              XC=real(XNUM(I))/real(NUMCL(I))
              YC=real(YNUM(I))/real(NUMCL(I))
              ZC=real(ZNUM(I))/real(NUMCL(I))
              DO J=1,NUMCL(I)
                AX(J)=real(X(I,J))-XC
                AY(J)=real(Y(I,J))-YC
                AZ(J)=real(Z(I,J))-ZC
              ENDDO
!
! HERE TO CALCULATE THE SQUARE RADIUS OF EVERY CLUSTER
              RS_SQ(I)=0.0
              DO J=1,NUMCL(I)
                RR=AX(J)*AX(J)+AY(J)*AY(J)+AZ(J)*AZ(J)
                RS_SQ(I)=RS_SQ(I)+RR/real(NUMCL(I))
              ENDDO
!
              DO J=1,NUMCL(I)
                TMT(1,1)=TMT(1,1)+AX(J)*AX(J)/real(NUMCL(I))
                TMT(1,2)=TMT(1,2)+AX(J)*AY(J)/real(NUMCL(I))
                TMT(1,3)=TMT(1,3)+AX(J)*AZ(J)/real(NUMCL(I))
                TMT(2,2)=TMT(2,2)+AY(J)*AY(J)/real(NUMCL(I))
                TMT(2,3)=TMT(2,3)+AY(J)*AZ(J)/real(NUMCL(I))
                TMT(3,3)=TMT(3,3)+AZ(J)*AZ(J)/real(NUMCL(I))
              ENDDO
!
! For small clusters with 0 value(s) in the dialogue of the orientation matrix
! This value may be not exactly 0.0 because of digital error
              IF (MIN(TMT(1,1),TMT(2,2),TMT(3,3)).LE.1.0E-8) THEN
                TMT(1,1)=TMT(1,1)+1.0e-4
                TMT(2,2)=TMT(2,2)+1.0e-4
                TMT(3,3)=TMT(3,3)+1.0e-4
              ENDIF
!
              TMT(2,1)=TMT(1,2)
              TMT(3,1)=TMT(1,3)
              TMT(3,2)=TMT(2,3)
              TMT0=TMT
              call ssyev('V','U',3,TMT,3,ev,work,8,info)
! The finite dimensions are Feret diameters in 3 principal axes. 
! When fdim(i) is very small (<=0.1), it means this direction there is only one voxel
! and the Feret diameter should be 1.
              fdim=4.0*sqrt(ev)
              do j=1,3
                if (fdim(j).lt.0.1) fdim(j)=1.0
              enddo
!
              WRITE(3,'(3I10,3x,3f10.2,2F12.3)') I,NUMCL(I),OUTSURF1(I),(fdim(J),J=1,3),&
                                                &fdim(1)/fdim(3),1.0-fdim(2)/fdim(3)
              WRITE(4,*) I,XC,YC,ZC
              WRITE(4,'(3e12.3)') (TMT0(1,J),J=1,3)
              WRITE(4,'(3e12.3)') (TMT0(2,J),J=1,3)
              WRITE(4,'(3e12.3)') (TMT0(3,J),J=1,3)
              WRITE(5,*) I
              WRITE(5,'(3f10.3,12X,3F8.1)') (TMT(1,J),J=1,3),(ACOSD(TMT(1,J)),J=1,3)
              WRITE(5,'(3f10.3,12X,3F8.1)') (TMT(2,J),J=1,3),(ACOSD(TMT(2,J)),J=1,3)
              WRITE(5,'(3f10.3,12X,3F8.1)') (TMT(3,J),J=1,3),(ACOSD(TMT(3,J)),J=1,3)
            ENDDO
!           write(*,*) '  Calculating orientation matrix of every cluster, finished. '
          DEALLOCATE (X,Y,Z,AX,AY,AZ)
          DEALLOCATE (NUMCL, OUTSURF1)
          DEALLOCATE (XNUM,YNUM,ZNUM)
!
          NTURN=NTURN+NGP
        ENDDO 
    RETURN
    END SUBROUTINE ORIAN_TMT   

