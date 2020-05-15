C
      SUBROUTINE RDCHG(R,MG,RV,FILTP)
C     SUBROUTINE RDCHG(CHMTYPE)
c
C        input file is modified, shortened version of charmm topology file
C
C        line definition:
C        First word: RESI, GROU, ATOM, PRES
C        Second entry: atom type
C        Third entry:  Internal atom type
C        fourth entry: partial charges, neutral groups
C        fifth entry: partial charges, charged groups
C        next entry: Rekker coefficients
C        remaining entries are Born radii.
C
C       read charges and Rekker hydrophobic fragmental constants.
C        For nonprotonatable residues charges are read into
C       R(4,I); for titratable residues, charges for charged state are
C       read into R(5,I) & charges for neutral state are read into R(4,I)
C        Hydrophobic constants are read into R(6,I)
c
c        e.l. mehler, January 30, 1997, dept. of pysiology/biophys
c        mt. sinai school of medicine, NY
c
      IMPLICIT REAL*8 (A-H,O-Z)
        character filnm*80,line*80,line2*80,atom*4,chatm(50)*4,res*3
        character lino*80,ans*1,patch*4,CHMTYPE*4,filtp*3
        logical parres
        integer*2 MG(*)
        real*8 R(6,*),RV(*)
        dimension ngp(50),chg(50),chgt(50),rkr(50),rsasa(50)
        data iout,inp/6,5/,lino/' '/,ncres/0/
        data parres/.false./
C       INCLUDE 'esp.inc'
C       do I=1,LATMX
C         MG(I)=NG(I)
C       enddo
        iatm=0
    1   read (1,'(a)',end=900,iostat=ierr) line
        if (line(1:1).eq.'*') go to 1
c       read(1,'(i5)')iatnum
    4   read (1,'(a)',end=900,iostat=ierr) line
C
C       iatm = atom count
C
        if (filtp.eq."CHM") then
          iatm=iatm+1
          read (line,'(5x,i5,1x,a3,2x,a4)') nres,res,atom
        else if (filtp.eq."PDB") then
c         if (line(:6).ne.'ATOM  '.and.line(:6).ne.'HETATM') go to 4
          if (line(:6).ne.'ATOM  ') go to 4
          iatm=iatm+1
          read (line,'(12x,a4,1x,a3,2x,i4)') atom,res,nres
c         write (*,'(12x,a4,1x,a3,2x,i4)') atom,res,nres
        else if (filtp.eq."EXT") then
          read(line,'(10x,i10,2x,a4,6x,a4)') nres,res,atom
        endif
C
c       if (atom.eq.'HT1 '.or.atom.eq.' H1 '.or. atom.eq.'HN1 '
c    *  .or. atom.eq.'HN2 ') then
        if (atom.eq.'HT1 ')  then
          patch='NTER'
          write(*,*)"patch=NTER"
          go to 13
        else if (atom.eq.'CAY ') then
          patch='ACE '
          go to 13
        else if (atom.eq.'NT  ') then
          patch='CT3 '
          go to 13
c       else if (atom.eq.'CAP ') then
c         patch='NTE2'
c         go to 13
c       else if (atom.eq.'CQ  ') then
c         patch='CTE2'
c         go to 13
        endif
        if (atom.eq.'OT2 '.or.atom.eq.' OXT') go to 12
        if (nres.eq.ncres) go to 201
        rewind 2
    2   read(2,'(a)',end=904,iostat=ierr) line2
         write(*,'(a)') line2
        if (line2(1:4).eq.'RESI'.and.line2(6:8).eq.res) then
          parres=.true.
          ic=0
          ncres=nres
    3     read (2,'(a)',end=1901,iostat=ierr) line2
          if (line2(:4).eq.'RESI') then
               ngp(ic)=1
               go to 200
          endif
C         READ IN CHARGES FOR RESIDUE
         if (line2(:4).eq.'GROU') then
             ng=1
             if(ic.ne.0) ngp(ic)=1
             go to 3
         elseif (line2(:4).eq.'ATOM') then
            ic=ic+1
            ng=0
            if (ic.gt.50) then
                write (*,*) '>50 atoms in residue'
                stop 132
            endif
            if(FILTP.eq.'CHM') 
     *      read (line2,'(5x,a4,8x,f6.3,1x,f6.3,2x,f6.3,1x,f6.4)') 
     *      chatm(ic),chg(ic),chgt(ic),rkr(ic),rsasa(ic)
            if (FILTP.eq.'PDB')
     *      read (line2,'(10x,a4,3x,f6.3,1x,f6.3,2x,f6.3,1x,f6.4)') 
     *      chatm(ic),chg(ic),chgt(ic),rkr(ic),rsasa(ic)
            ngp(ic)=ng
            go to 3
         endif
        else
c           write(*,*)"entering go to 2"
            go to 2
        endif
  200   icm=ic ! modified 11.07.14 DB
C       ALL CHARGES READ IN, SO PUT THEM IN ATOM ARRAYS
C       the charges and rekker coefficients are first read into local arrays
C       when done they are transferred to the global arrays
C
  201   if (icm.ne.0) then
         do ii=1,icm
          if (parres.and.chatm(ii).eq.atom) then
            R(4,iatm)=chg(ii)
            R(5,iatm)=chgt(ii)
            R(6,iatm)=rkr(ii)
            write(*,*)"R(6,atm) within loop in rdchg", R(6,iatm)
            RV(iatm)=rsasa(ii)
            MG(iatm)=ngp(ii)
            go to 4
          endif
c         write(*,*)"rsasa(ii)", rsasa(ii)
         enddo
          go to 4
        else 
          go to 4
        endif
   13   ic=0
        rewind 2
   11   read (2,'(a)',end=2901,iostat=ierr) line2
        if (line2(:4).eq.'PRES'.and.line2(6:9).eq.patch) then
   10     read (2,'(a)',end=3901,iostat=ierr) line2
c         write(*,*)"line2 in patch", line2,patch
190       read (2,'(a)',end=3901,iostat=ierr) line2
          write(*,*)"line2 in patch", line2,patch
          if (line2(:4).eq."ATOM") then
            ic=ic+1
            read (line2,'(17x,f6.3,1x,f6.3,2x,f6.3,1x,f6.4)') chg(ic),
     *      chgt(ic),rkr(ic),rsasa(ic)
            write (*,'(17x,f6.3,1x,f6.3,2x,f6.3,1x,f6.4,
     *      "patch within NTER")') chg(ic),chgt(ic),rkr(ic),rsasa(ic) 
c           read (line2,'(17x,f6.3,1x,f6.3,2x,f6.3,1x,f6.4)') 
c    *      R(4,iatm-1+ic),
c    *      R(5,iatm-1+ic),R(6,iatm-1+ic),RV(iatm-1+ic)
            MG(iatm+ic)=0
            go to 190
          else if (line2(:4).eq.'PRES') then
            MG(iatm)=1
            go to 190
          else
           do ii=1,ic
            R(4,iatm-2+ii)=chg(ii)
            R(5,iatm-2+ii)=chgt(ii)
            R(6,iatm-2+ii)=rkr(ii)
            RV(iatm-2+ii)=rsasa(ii)
            MG(iatm-2+ii)=0
c           write(*,*)"rsasa(ii)",RV(iatm-2+ii)
           enddo
c          MG(iatm-1)=1
            iadv=ic-1
            go to 14
          endif
        else
          go to 11
        endif
   14   do ii=1,iadv
   15   read (1,'(a)',end=900,iostat=ierr) line
        if (FILTP.EQ.'PDB'.AND.(line(:6).ne.'ATOM  '.and.line(:6)
     *  .ne.'HETATM')) go to 15 
        iatm=iatm+1
        enddo
        go to 4
   12   ic=0
C       READ IN Cterminal --COO- PATCHES
        rewind 2
   21   read (2,'(a)',end=4901,iostat=ierr) line2
        if (line2(:4).eq.'PRES'.and.line2(6:9).eq.'CTER') then
          read (2,'(a)',end=5901,iostat=ierr) line2
          read (2,'(a)',end=5901,iostat=ierr) line2
          write(*,*)"line2 in CTER",line2
   20     read (2,'(a)',end=5901,iostat=ierr) line2
          write(*,*)"line2 in CTER",line2
          if (line2(:4).eq.'ATOM') then
c          if (line2(6:6).ne.'C') then
            ic=ic+1
            read (line2,'(17x,f6.3,1x,f6.3,2x,f6.3,1x,f6.4)') chg(ic),
     *      chgt(ic),rkr(ic),rsasa(ic)
            write (*,'(17x,f6.3,1x,f6.3,2x,f6.3,1x,f6.4),"in CTER"') 
     *      chg(ic),chgt(ic),rkr(ic),rsasa(ic)
            ngp(ic)=0
            go to 20
           else if (line2(:4).eq.'PRES') then
            ngp(ic)=1
            go to 20
           else
            do ii=1,ic
             R(4,iatm-2+ii)=chg(ii)
             R(5,iatm-2+ii)=chgt(ii)
             R(6,iatm-2+ii)=rkr(ii)
             RV(iatm-2+ii)=rsasa(ii)
             MG(iatm-2+ii)=0
c            write(*,*)"rsasa(ii)",RV(iatm-2+ii)
            enddo
            MG(iatm)=1
16          read (1,'(a)',end=900,iostat=ierr) line
            if (FILTP.EQ.'PDB'.AND.(line(:6).ne.'ATOM  '.and.line(:6)
     *     .ne.'HETATM')) go to 16 
            iatm=iatm+1
             go to 4
c          endif
          endif
        else
          go to 21
        endif
  900   continue
       write (*,1000) (jj,R(4,jj),R(5,jj),MG(jj),R(6,jj)
     * ,RV(jj),jj=1,iatm)
       write(*,*)"R(6,JJ) in rdchg", (R(6,jj),jj=1,iatm)
 1000   format (i5,2f9.3,i6,2f9.3)
        rewind 1
        rewind 2
        return
c       call exit
1901   write (iout,*) ' end of file on unit 2, in 3'
2901   write (iout,*) ' end of file on unit 2, in 10'
3901   write (iout,*) ' end of file on unit 2, in 11'
4901   write (iout,*) ' end of file on unit 2, in 21'
5901   write (iout,*) ' end of file on unit 2, in 20'
        call exit
  904   write (iout,*) ' residue',res,nres,' not in par19'
        ncres=nres
        parres=.false.
        icm=0
        go to 200
        end

