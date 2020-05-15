              Character filename*20,res*3,atom*4,aothers*20,a*80
              character b(30)*80,atot*80,temp*80,a1(30)*80,arg(10)*320,
     *        chain(1000000)*1,atemp1*1,atemp2*1,atemp3*1
              dimension x(1000000),y(1000000),z(1000000),num(1000000),
     *        iatm(1000000),ires(1000000),res(1000000),atom(1000000),
     *        aothers(10000000)

              do i = 1,iargc()
                call getarg(i, arg(i))
                  write(*,*)i,arg(i)
              enddo
              write(*,*)"input file:",arg(1)
              write(*,*)"outputfile:", arg(2)
              write(*,*)"number of atoms",arg(3)
              read(arg(3),'(i5)')iatomnum
              open(unit=10,file=arg(1))
              open(unit=4,file=arg(2))
                write(4,'("*",2x,a32)') arg(1)
              do i=1,3
                write(4,'("*")')
              enddo
                write(4,'(i5)')iatomnum
              irescopy=-100
              ii=0
              is=0
              do i=1,100000
               read(10,'(a80)',END=99) a
               if(a(:4).eq."ATOM") then
                is=is+1
                read(a,1,END=99)iatm(is),atom(is),res(is),chain(is),
     *          ires(is),x(is),y(is),z(is)
                 if (atom(is)(1:1).eq." ") then
c                   write(*,*)"beginning of rearrangment",atom(is)
                    atemp1 = atom(is)(2:2) 
                    atom(is)(1:1) = atemp1
                    atemp2 = atom(is)(3:3) 
                    atom(is)(2:2) = atemp2
                    atemp3 = atom(is)(4:4) 
                    atom(is)(3:3) = atemp3
                    atom(is)(4:4)=" "
c                   write(*,*)"end of rearrangment",atom(is)
                 endif
                if (ires(is).ne.irescopy)  then
                 ii=ii+1
                endif
                write(4,2)is,ii,res(is),atom(is),x(is),y(is),z(is),
     *          arg(1)(:4),ii
                irescopy=ires(is)
               endif
              enddo
99            continue
1             format(6x,i5,1x,a4,1x,a3,1x,a1,i4,4x,3f8.3)
2             format(i5,1x,i4,1x,a3,2x,a4,3f10.5,2x,a4,1x,i4,2x,
     *        "0.00000")
              stop
              end
