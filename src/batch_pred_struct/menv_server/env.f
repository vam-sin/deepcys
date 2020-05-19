         subroutine env(inpt,iout,renv,debug,ngrp,ntatm,idb,nres,isgrd)
c-------------------------------------------------------------------------
c        declaration of statements
c-------------------------------------------------------------------------
         IMPLICIT REAL*8 (A-H,O-Z)
         dimension amax(1000000)
         CHARACTER ANS*1,CTRM*4,NTRM*4
         logical debug,grp
         INTEGER IATC(100),JATC(100),NGRDSTART
         REAL*8 KAPPA,RCUT(2),RCUT0(2),DR(3),VDW,HYDR
         PARAMETER (zro=0.,DST=78.4,RG=.0019872)
         include 'esp.inc'
         include 'chd.inc'
         include 'screen.inc'
         DATA KAPPA/.32915/,SR1,SR2/10.,30./
         DATA CTRM,NTRM/'CTRM','NTRM'/
c-------------------------------------------------------------------------
c        initialize few parameters
c-------------------------------------------------------------------------
         ist=1
         alpha=10.
C---------------------------------------------------------------------------
c        calculation per group
c------------------------------------------------------------------------------
         do 100 igrp=1,ngrp
            do iii=1,3
               rgrp(iii,igrp)=0.
c              write(*,*)"ngrps in env.f",ngrps(iii,igrp)
            enddo
            itiatm=NGRPS(2,igrp)
            istp=ist+NGRPS(2,igrp)-1
            jst=1
c           write(*,*)"igrp,ist,istp,ngrps(2,igrp)",igrp,ist,istp,
c    *      ngrps(2,igrp)
c--------------------------------------------------------------------------
c          defining mean position of a group, igrp
c---------------------------------------------------------------------------
            do iatm=ist,istp
              do iii=1,3
                rgrp(iii,igrp)=rgrp(iii,igrp)+r(iii,iatm)
c             if (igrp.eq.26) write(*,*)"x,y,z for igrp",rgrp(iii,igrp)
              enddo
            enddo
            do iii=1,3
c             if (igrp.eq.26) write(*,*)"x,y,z for igrp, in group",
c    *rgrp(iii,igrp)
                rgrp(iii,igrp)=rgrp(iii,igrp)/(istp-ist+1)
c             if (igrp.eq.26) write(*,*)"x,y,z for igrp in group/avg",
c    *rgrp(iii,igrp)
            enddo
c--------------------------------------------------------------------------
c           initializing variables
c--------------------------------------------------------------------------
            kres=ares(ist)
            ihb=0.
            tot=0.
            hydr=0. 
c           write(*,*)"hydr before 102 loop is initialized every time",hydr
c           write(*,*)"ntatm in env.f", ntatm
c---------------------------------------------------------------------------
c          find the neighbors around group igrp
c---------------------------------------------------------------------------
            do 102 jgrp=1,ngrp
c             jst=jst+ngrps(2,jgrp)
              do iii=1,3
                rgrp(iii,jgrp)=0.
              enddo
              jstp=jst+ngrps(2,jgrp)-1
              jtiatm=ngrps(2,jgrp)
c             if (igrp.eq.jgrp) go to 102
c----------------------------------------------
c            find mean position of group jgrp
c----------------------------------------------
              do jatm=jst,jstp
                do iii=1,3
                  rgrp(iii,jgrp)=rgrp(iii,jgrp)+r(iii,jatm)
                enddo
              enddo
              do iii=1,3
                rgrp(iii,jgrp)=rgrp(iii,jgrp)/(jstp-jst+1)
              enddo
c---------------------------------------------------------------------------
c            calculate the distance between mean position of grs igrp & jgrp
c---------------------------------------------------------------------------
              RD=0.0
              do iii=1,3
               rd=rd+((rgrp(iii,jgrp)-rgrp(iii,igrp))*
     *                (rgrp(iii,jgrp)-rgrp(iii,igrp)))
              enddo
c             write(*,*)"igrp, jgrrp,rd",igrp,jgrp,rd
c--------------------------------------------------------------
c         if distance greater than a threshold go to next jgrp
c--------------------------------------------------------------
              if ((rd.gt.49.0) .or. (ares(jst).eq.ares(ist))) then
c             if ((rd.gt.49.0) ) then
                 jst=jst+ngrps(2,jgrp)
                 go to 102
              endif
c            write(*,*)"jgrp,jst,jstp",jgrp,jst,jstp
C--------------------------------------------------------------------------
c----------------- loop over igrp------------------------------------------
C--------------------------------------------------------------------------
            do jatm=jst,jstp
                amax(jatm)=0.
            enddo
c---------------------------------------------------------------------------
c         if distance less than the threshold find the Rekkers constant from 
c         jgrp group; initialize variable for rekker's constant
c---------------------------------------------------------------------------
c--------------------------------------------------------------------------
c          loop over jgrp
c--------------------------------------------------------------------------
         do 101 iatm=ist,istp
          do 103 jatm=jst,jstp
              ires=ares(iatm)
              jres=ares(jatm)
c--------------------------------------------------------------------------
c          calculate distance between each iatm and all jatm, only once
c--------------------------------------------------------------------------
           if (ares(iatm).eq.ares(jatm)) then 
            goto 103
           else
             rd=0.0
             do jjj=1,3
              rd=rd+((r(jjj,jatm)-r(jjj,iatm))*
     *             (r(jjj,jatm)-r(jjj,iatm)))
             enddo
             rd=sqrt(rd)
c---------------------------------------------------------------------------
c          defining van der Waals radii)
c----------------------------------------------------------------
            if(atm(iatm)(:1).eq.'N'.or.atm(iatm)(:1).eq.'O') then
             vdw=rv(iatm)+2.275+0.2!angstrom; vdw ca=2.275+polar H radius, 0.2
            else
             vdw=rv(iatm)+2.275
            endif
c----------------------------------------------------------------------------
c          if distance rd is less than van der waals radii of iatm  - then
c----------------------------------------------------------------------------
            if (rd.le.vdw) then
             tot=1.
             if (tot.gt.amax(jatm)) amax(jatm)=tot        
c----------------------------------------------------------------------------
c           if distance between vdw and vdw+2 use guassian scale factor
c---------------------------------------------------------------------------
            elseif (rd.le.vdw+0.2) then
             rd=rd-vdw
             tot=exp(-alpha*rd*rd)
             if (tot.gt.amax(jatm)) amax(jatm)=tot
            endif
           endif
c---------------------------------------------------------------------------
c          end of the menv calculation
c---------------------------------------------------------------------------
103      continue 
c        jst=jst+ngrps(2,jgrp)
101      continue
         do jatm=jst,jstp
c          if(amax(jatm).gt.0.0.and. r(6,jatm).gt.0.0) 
           write(*,*)"amax(jatm)",jatm,amax(jatm)," r(6,jatm)",r(6,jatm)
     *      ," atm(jatm) ",atm(jatm),"hydr",hydr
           hydr=hydr+r(6,jatm)*amax(jatm)
           write(*,*)"hydr in side 102 loop of env.f",hydr
         enddo
         JST=JST+NGRPS(2,JGRP)
102      continue
         rkrh(igrp)=hydr
         write(*,*)"igrp,rkrh(igrp) in env.f",igrp, rkrh(igrp)
         ist=ist+ngrps(2,igrp)
100      continue
         return
         end
