c penalty2 calcula el castigo en la chi² por picos en las estratificaciones
c         calcula tambien dcastigo que es la derivada de castigo respecto
c         a la varicion de cada parametro del modelo
c         d2castigo es la derivada segunda.
c es llamada por marqcoef2 (linea 83)

	subroutine penalty2(mfit,mnodos,atry,castigo,dcastigo)

	implicit real*4 (a-h,o-z)
	include 'PARAMETER'   !por kld y mmx==kn que se tomaba 200

	dimension mnodos(*),atry(*),dcastigo(*)

	external signo

	castigo=0
	kred=0
	do j=1,18            !do en grupos de varibles (1=t,2=p,...etc)
	   if (mnodos(j) .gt. 2) then
	      c=1
	      c1=abs(atry(kred+1))
	      c2=c1+abs(atry(kred+2))
	      if(c1.ne. 0)then
	         c=c1                          !lo supondre constante
	      else if (c2 .ne. 0)then
                 c=c2/2.
	      endif
              suma=0.
	      do k=1,mnodos(j)-2
	         y1=atry(kred+k)
		 y2=atry(kred+k+1)
		 y3=atry(kred+k+2)
	         f=(y2-y1)*(y3-y2)
	         suma=suma+abs(f-abs(f))
	      end do
	      castigo=castigo+suma/c
c	      print*,'castigo sobre la variable',j,suma/c
	   endif
	   kred=kred+mnodos(j)
	end do

c	print*,'castigo total',castigo


	kred=0
	do j=1,18            !do en grupos de varibles (1=t,2=p,...etc)
	   if (mnodos(j) .gt. 2) then
	      c=1
	      c1=abs(atry(kred+1))
	      c2=c1+abs(atry(kred+2))
	      if(c1.ne. 0)then
	         c=c1                          !lo supondre constante
	      else if (c2 .ne. 0)then
                 c=c2/2.
	      endif

	      do k=1,mnodos(j)
		 y_1=0.
		 if(k.ge.2)y_1=atry(kred+k-1)
		 y_2=0.
	         if(k.ge.3)y_2=atry(kred+k-2)
		 y1=atry(kred+k)
	         y2=atry(kred+k+1)
	         y3=atry(kred+k+2)
	         fk=(y2-y1)*(y3-y2)
                 fk_1=(y1-y_1)*(y2-y1)
                 fk_2=(y_1-y_2)*(y1-y_1)
		 dc=0.
                 if(k.le.mnodos(j)-2)dc=signo(fk-abs(fk))*(signo(fk)-1.)*(y3-y2)
                 if(k.ge.2 .and.k.le.mnodos(j)-1)dc=dc+signo(fk_1-abs(fk_1))*(1.-signo(fk_1))*(y2-2*y1+y_1)
                 if(k.ge.3)dc=dc+signo(fk_2-abs(fk_2))*(1.-signo(fk_2))*(y_1-y_2)
		 dcastigo(kred+k)=dc/c

c		 if(k.ge.3)d2castigo(kred+k,kred+k-2)=-signo(fk_2-abs(fk_2))*(1.-signo(fk_2))/c
c		 if(k.ge.2)d2castigo(kred+k,kred+k-1)=(signo(fk_2-abs(fk_2))
c     &                              *(1.-signo(fk_2))+signo(fk_1-abs(fk_1))*(1.-signo(fk_1)))/c
c		 if(k.ge.3)d2castigo(kred+k,kred+k)=-2*signo(fk_1-abs(fk_1))*(1.-signo(fk_1))/c
c     		 d2castigo(kred+k,kred+k+1)=(signo(fk  -abs(fk  ))*(1.-signo(fk  )))/c
c                 if(k.ge.2 .and.k.le.mnodos(j)-1)d2castigo(kred+k,kred+k+1)=d2castigo(kred+k,kred+k+1)+
c     &                                      (signo(fk_1-abs(fk_1))*(1.-signo(fk_1)))/c
c     		 if(k.le.mnodos(j)-2)d2castigo(kred+k,kred+k+2)=-signo(fk  -abs(fk  ))*(1.-signo(fk  ))/c
	      enddo
	   endif
	   kred=kred+mnodos(j)
	end do

	return
	end

c-----------------------------------------------------------------------------
c signo: signo(x)=1 si x>=0 y signo(x)=-1 si x<0

	function signo(x)
	real*4 x

	signo=1
	if(x .lt. 0)signo=-1

	return
	end
c-----------------------------------------------------------------------------







