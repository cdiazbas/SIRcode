c invierte una matriz a por medio del metodo svd
c y calcula la solucion x al problema a*x=b
c a la salida la solucion x se guardara en b
c la matriz a queda destruida a la salida
        subroutine svdmatriz2(a,n,mnodos,b,tol,v,w)

	implicit real*4 (a-h,o-z)
	include 'PARAMETER'  !para kt y mfitmax

        dimension a(mfitmax,mfitmax),b(*),w(mfitmax),v(mfitmax,mfitmax)
     &            ,mnodos(*),mm(18),tolv(mfitmax)
	real wt(kt),ww(kt),wi(2,kt)


c      print*,'mfitmax,n=',mfitmax,n
c      do i=1,mfitmax
c          do j=1,mfitmax
c             print*,'a(i,j)=',i,j,a(i,j)
c          enddo
c      enddo

c      do i=1,mfitmax
c         print*,'i,b=',i,b(i)
c      enddo

       call svdcmp2(a,n,n,n,n,w,v)

c contamos el numero de variables (no param. libres sino variables t,pe...)
	nvar=0
	do i=1,18
	   if(mnodos(i).gt.0)then
	       nvar=nvar+1
	       mm(nvar)=mnodos(i)
	       tolv(nvar)=tol
	       if(i .eq.2 .or. i.eq.10)tolv(nvar)=tol*1.e-8
c               print*,'en svdmatriz2 mm(',nvar,')=',mm(nvar),' n=',n
	   end if
	end do	

c el numero minimo de param. libres es el numero de variables.
c Si el num. de variables es > que el de param. libres estudiamos si
c se debe reducir el numero de variables

	if(nvar.ne.n)then
	  
c evidentemente si el num. de variables es 1 no hacemos nada

	   if(nvar.eq.1)then

c calculamos el maximo
	      wmax=0.d0
	      do j=1,n
	         if(w(j).gt.wmax)wmax=w(j)
	      end do

c anulamos los terminos menores que el max. por tol.
c	      wmin=wmax*tol
	      do j=1,n
	         wmin=wmax*tolv(j)
	         if(w(j).lt.wmin)w(j)=0
	      end do	
           else

c repito dos veces el proceso para optimizarlo
	      do ido=1,2

	         do i=1,n
	            ww(i)=0.
	         end do
	         i2=0
	         do i0=1,nvar 
                    i1=i2+1	         !indice inicial para el parametro i
		    i2=i2+mm(i0)         !indice final
                    do j=1,n
                       wt(j)=0           !inicializo las w para cada parametro
                    end do

c construyo la submatriz
	            do j=1,n	
                       do i=i1,i2        !w para cada parametro
	                  wt(j)=wt(j)+v(i,j)*v(i,j)*w(j)
	               end do
	            end do

c calculamos el maximo
	            wtx=0.
                    do j=1,n
	               if(wt(j).gt.wtx)wtx=wt(j)
	            end do

c anulamos los terminos menores que el max. por tol.
c	            wtn=wtx*tol
	            do j=1,n
	               wtn=wtx*tolv(j)
	               if(wt(j).lt.wtn)then
		          wt(j)=0.
	               else
	                  ww(j)=ww(j)+wt(j)
	               end if
	            end do
	         end do
	         do j=1,n
	            w(j)=ww(j)
	            wi(ido,j)=w(j)
	         end do
	      end do

c adopto la primera correccion modificada por la segunda
	      do j=1,n
	         if(abs(wi(2,j)).gt.1.e-10)then
	            w(j)=wi(1,j)*wi(1,j)/wi(2,j)
	         else
	            w(j)=0.e0
	         end if
	      end do
	   end if

        end if

c        print*,'en svdmatriz2 a',a(1,1),' b=',b(1)

	call svbksb(a,w,v,n,n,n,n,b,b)

c        print*,'en svdmatriz2 a*x=b x=',b(1)


	return
	end
