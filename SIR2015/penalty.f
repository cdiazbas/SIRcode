c penalty calcula el castigo en la chi² por picos en las estratificaciones
c         calcula tambien dcastigo que es la derivada de castigo respecto
c         a la varicion de cada parametro del modelo
c         ddcastigo es la derivada segunda. Como la unica componente
c         no nula es j=1 para todo k la guardo en un vector
c es llamada por marqcoef2 (linea 83)

	subroutine penalty(mfit,mnodos,atry,castigo,dcastigo,ddcastigo)

	implicit real*4 (a-h,o-z)
	include 'PARAMETER'   !por kld y mmx==kn que se tomaba 200

	dimension mnodos(*),atry(*),dcastigo(*),ddcastigo(*)

	external signo

	castigo=0
	kred=0
	do j=1,18            !do en grupos de varibles (1=t,2=p,...etc)
	   if (mnodos(j) .gt. 2) then
	      c=1
	      if(atry(kred+1) .ne. 0)c=atry(kred+1)
              suma=abs(atry(kred+2)-atry(kred+1))
	      do k=2,mnodos(j)-1
	         suma=suma+abs(atry(kred+k+1)-atry(kred+k))
	      end do
	      castigo=castigo+(suma-abs(atry(kred+mnodos(j))-atry(kred+1)))/c
c	      print*,'castigo sobre la variable',j,(suma-abs(atry(kred+mnodos(j))-atry(kred+1)))/c
	   endif
	   kred=kred+mnodos(j)
	end do

c	print*,'castigo total',castigo


	kred=0
	do j=1,18            !do en grupos de varibles (1=t,2=p,...etc)
	   if (mnodos(j) .gt. 2) then
	      c=1
	      fac=0  !para anular las derivadas si c=1
	      if(atry(kred+1) .ne. 0)then
	         c=atry(kred+1)
		 fac=1
	      end if
              dcastigo(kred+1)=-castigo/c*fac+( signo(atry(kred+mnodos(j))-atry(kred+1))
     &                                     -signo(atry(kred+2)-atry(kred+1))        )/abs(c)
              ddcastigo(kred+1)=(-dcastigo(kred+1)/c-castigo/c/c)*fac
c	      print*,'derivada castigo respecto a la variable',j,' en',kred+1,dcastigo(kred+1)
	      do k=2,mnodos(j)-1
	         dcastigo(kred+k)=( signo(atry(kred+k)-atry(kred+k-1))
     &                             -signo(atry(kred+k+1)-atry(kred+k)))/abs(c)
                 ddcastigo(kred+k)=-dcastigo(kred+k)/c*fac
c                 print*,'derivada castigo respecto a la variable',j,' en',kred+k,dcastigo(kred+k)
	      end do
	      dcastigo(kred+mnodos(j))=( signo(atry(kred+mnodos(j))-atry(kred+mnodos(j)-1))
     &                                  -signo(atry(kred+mnodos(j))-atry(kred+1))          )/abs(c)
              ddcastigo(kred+mnodos(j))=-dcastigo(kred+mnodos(j))/c*fac

c              print*,'derivada castigo respecto a la variable',j,' en',kred+mnodos(j),dcastigo(kred+mnodos(j))

	   endif
	   kred=kred+mnodos(j)
	end do

	return
	end








