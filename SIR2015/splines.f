	subroutine splines(x,y,no,ntau,tau,yy,f)

	include 'PARAMETER'   !por kt
	parameter (kt2=kt*kn)
	implicit real*4 (a-h,o-z)
	real*4 a(kt,kt),x(*),y(*),tau(*),yy(*),f(kt,kt),ff(kt2)
	real*4 taure(kt)
	COMMON/FACTORSPLIN/FF
 	common/cambio/ncambno,ncambpre  !si 0,0 cambio nodos, cambio prec
c ncambno vale 0 la primera vez que se corre un ciclo
c      "       "   1 las restantes
	data iset/0/
	
 	iset=iset*ncambno

	n=no+2
	IF(ISET.EQ.0)THEN
           paso=x(2)-x(1)
	   iset=1
 	   ncambno=1
	   if(n.gt.2)then
	      CALL SPLINB(n,kt,a)

	      ks=0
	      do i=1,ntau
                j1=1+int((tau(i)-x(1))/paso)
	        j2=j1+1
	        f2=(tau(i)-x(j1))/paso
	        f1=1.d0-f2
	        f3=f1*(f1*f1-1.d0)
	        f4=f2*(f2*f2-1.d0)
	        do k=1,n
	           f(i,k)=f3*a(j1,k)+f4*a(j2,k)
	           if(k.eq.j1)f(i,k)=f(i,k)+f1
	           if(k.eq.j2)f(i,k)=f(i,k)+f2
	        end do
	        do k=1,n
	           ks=ks+1
	           ff(ks)=f(i,k)
	        end do
	     end do
	   else
	      do i=1,ntau
	         taure(i)=(tau(i)-x(1))/paso
	      end do
	      if(no.eq.0)then
	          ks=0
	          do i=1,ntau
	             f(i,1)=1.d0-taure(i)
	             f(i,2)=taure(i)
	             do k=1,n
		        ks=ks+1
                        ff(ks)=f(i,k)
	             end do
	           end do
	       else
	       ks=0
               do i=1,ntau
 	          if(taure(i).lt.1.d0)then
			f(i,1)=1.d0-taure(i)
			f(i,2)=taure(i)
			f(i,3)=0.d0
	           else
			f(i,1)=0.d0
			f(i,2)=2.d0-taure(i)
			f(i,3)=taure(i)-1.d0
		   end if
		   do k=1,n
		      ks=ks+1
		      ff(ks)=f(i,k)
		   end do
		end do
	      end if
	    end if
	END IF

c a partir de aqui lo hace siempre

	ks=0
	do i=1,ntau
	   yy(i)=0.d0
	   do k=1,n
	      ks=ks+1
	      f(i,k)=ff(ks)
	      yy(i)=yy(i)+ff(ks)*y(k)
	   end do
	end do

	return
	end

