              Character filename*20,res*3,atom*4,aothers*20,a*80
              character b(30)*80,atot*80,temp*80,a1(30)*80,arg(10)*32
              dimension x(1000000),y(1000000),z(1000000),num(1000000),
     *        iatm(1000000),ires(1000000),res(1000000),atom(1000000),
     *        aothers(10000000)

              do i = 1,iargc()
                call getarg(i, arg(i))
                  write(*,*)i,arg(i)
              enddo
              write(*,*)"input file:",arg(1)
              write(*,*)"outputfile:", arg(2)
              open(unit=10,file=arg(1))
              open(unit=4,file=arg(2))
              do i=1,5
                read(10,'(a80)')a
                write(*,'(a80)')a
                write(4,'(a80)')a
              enddo
              irescopy=1
              read(10,'(6x,i5)')iresnum
              if (iresnum.eq.1) then
               ii=1
              else 
               ii=0
              endif
              rewind(unit=10)
              do i=1,5
                read(10,'(a80)')a
              enddo
              iii=1
              do i=1,100000
               read(10,1,END=99)iatm(i),ires(i),res(i),atom(i),
     *         x(i),y(i),z(i),aothers(i)
c              if (ires(1).eq.1)  then 
c                    ii=1
c              else
c                    ii=0
c              endif
               if (ires(i).ne.irescopy)  then
                ii=ii+1
               endif 
               if (ii.eq.1) then
                 if (atom(i).eq."HT1 ") then 
                 write(*,*)"check HT1"
                    iii=1
                    iatm(i)=iii
                     write(3,2)iatm(i),ii,res(i),atom(i),
     *               x(i),y(i),z(i),aothers(i) 
                 elseif (atom(i).eq."N   ") then 
                 write(*,*)"check N"
                    iii=2
                    iatm(i)=iii
                     write(3,2)iatm(i),ii,res(i),atom(i),
     *               x(i),y(i),z(i),aothers(i) 
                 elseif (atom(i).eq."HT2 ") then 
                 write(*,*)"check HT2"
                    iii=3
                    iatm(i)=iii
                     write(3,2)iatm(i),ii,res(i),atom(i),
     *               x(i),y(i),z(i),aothers(i) 
                 elseif (atom(i).eq."HT3 ") then
                 write(*,*)"check HT3"
                    iii=4
                    iatm(i)=iii
                     write(3,2)iatm(i),ii,res(i),atom(i),
     *               x(i),y(i),z(i),aothers(i) 
                 elseif (atom(i).eq."CA  ") then 
                 write(*,*)"check CA"
                    iii=5
                    iatm(i)=iii
                     write(3,2)iatm(i),ii,res(i),atom(i),
     *               x(i),y(i),z(i),aothers(i) 
                 else
                     iii=iii+1
                    iatm(i)=iii
                     write(3,2)iatm(i),ii,res(i),atom(i),
     *               x(i),y(i),z(i),aothers(i) 
                    ilist=iii
                 endif
                 write(*,*)"iii",iii
               else
                     write(3,2)i,ii,res(i),atom(i),x(i),y(i),z(i),
     *               aothers(i)
               endif
               irescopy=ires(i)
              enddo
99            continue
              iserial=i
                 write(*,*)"ilist",ilist
              close(unit=3)
              open (unit=3,file='fort.3')
              do i5=1,ilist
               read(3,'(a80)')a1(i5)
               write(*,'(a80)')a1(i5)
                read(a1(i5)(:5),'(i5)')num(i5)
              enddo
               DO 15 I=1,ilist
                AMIN=num(i)
                temp=a1(i)
                DO 25 K=i+1,ilist
                 IF (AMIN.GE.num(K)) then
                  Zx=AMIN
                  AMIN=num(K)
                  num(K)=Zx
                  atot=temp
                  temp=a1(k)
                  a1(k)=atot
                 endif
25              CONTINUE
                num(I5)=AMIN
                a1(i5)=temp
                write(4,'(a80)')a1(i5)
15             continue 
               write(*,*)"ilist, iserial-1",ilist,iserial-1
               read (3,2) (iatm(i5),ires(i5),res(i5),atom(i5),
     *         x(i5),y(i5),z(i5),aothers(i5),I5=ilist+1,iserial-1)
               WRITE (4,2) (iatm(i5),ires(i5),res(i5),atom(i5),
     *         x(i5),y(i5),z(i5),aothers(i5),I5=ilist+1,iserial-1)
1             format(i5,3x,i3,1x,a3,2x,a4,f8.3,3x,f8.3,2x,f8.3,a20)
2             format(i5,1x,i4,1x,a3,2x,a4,3f10.5,a20)
              stop
              end
