      dimension x(100000),y(100000),z(100000),iss(100000),in(1000)
     *,im(1000),rss(1000),line(100000),is(100000),iab(100000)
      character line*80
      i1=0
      open(unit=3,file="ssbond.log",access='append')
      open (unit=4,file="num")
      do i=1,100000
       read(*,'(a80)',END=99)line(i)
       if (line(i)(:1).ne."*" .and. line(i)(17:18).eq."SG") then
c      if (line(i)(:4).eq."ATOM" .and. line(i)(14:15).eq."SG") then
        i1=i1+1
        iss(i1)=i
        read(line(i),1) is(i),x(i1),y(i1),z(i1)
c       write(*,1) is(i),x(i1),y(i1),z(i1)
        iab(i1)=is(i)
        else
c       write(*,'(a80)')line(i)
       endif
      enddo
99    continue
      n=0
      do k=1,i1
       do j=k+1,i1
         if (iab(k).eq.iab(j)) goto 101
         a=x(k)-x(j)
         b=y(k)-y(j)
         c=z(k)-z(j)
         r=a*a+b*b+c*c
         r=sqrt(r)
         if (r.lt.2.5.and. r .gt. 1.5) then
          n=n+1
          in(n)=iss(k)
          im(n)=iss(j)
          rss(n)=r
         endif
101    enddo
      enddo
      do i=1,n
       write(3,'(2i5,f8.2)')is(in(i)),is(im(i)),rss(i)
       write(*,'("PATCH DISU MAIN", i5," MAIN",i5)')is(in(i)),is(im(i))
c      write(*,'("sed s/i"",i5," CYS  SG"/"",i5," CYS  SM/")')is(in(i))
c      write(*,'("sed s/i"",i5," CYS  SG"/"",i5," CYS  SM/")')is(im(i))
c      write(*,'(a11,"CY2",a66)')line(in(i))(:11),line(in(i))(15:)
c      write(*,'(a11,"CY2",a66)')line(im(i))(:11),line(im(i))(15:)
      enddo
      write(4,*) (is(in(i)),is(im(i)),i=1,n)
1     format(6x,i4,10x,3f10.5)
c1     format(6x,i5,19x,3f8.3)
2     format(a6,i4,1x,"CY2",a66)
c1     format(22x,i4,4x,3f8.3)
      stop
      end
 
     

