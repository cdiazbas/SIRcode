	subroutine splines22(x,y,no,ntau,tau,yy,f)
	include 'PARAMETER'   !por kt

	parameter (kt2=kt*kn)  !lo pongo en ff, que antes tenia 5250
	implicit real*4 (a-h,o-z)
	real*4 a(kt,kt),x(*),y(*),tau(*),yy(*),f(kt,kt),ff(kt2)
	real*4 taure(kt),peso(5),sum(kt)

	COMMON/FACTORSPLIN/FF
        common/nspl_lin/nspl_lin !para seleccionar interp.lineal=(1,3 o 5) o splines=(0,2 o 4)


	n=no+2

        if(nspl_lin.eq.0 .or. nspl_lin.eq.2 .or. nspl_lin.eq.4) then !if 0
c Interpolamos por splines ---------------------------------
c           print*,'splines22 interpola splines'


           paso=x(2)-x(1)
           p2=2*paso*paso
	   if(n.gt.3)then    !if 1
	      CALL SPLINB(n,kt,a)

	      ks=0
	      do i=1,ntau
                j1=1+int((tau(i)-x(1))/paso)
	        j2=j1+1
	        f2=(tau(i)-x(j1))/paso
	        f1=1.0-f2
	        f3=f1*(f1*f1-1.0)
	        f4=f2*(f2*f2-1.0)
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
	   else    ! else (if 1)
                 do i=1,ntau
	            taure(i)=(tau(i)-x(1))/paso
	         end do
              if(n.eq.2)then !if 2
      		      ks=0
c                      print*,'trabajamos con 2 nodos bb'
		      do i=1,ntau
			 f(i,1)=1.0-taure(i)
			 f(i,2)=taure(i)
			    do k=1,n
			       ks=ks+1
			       ff(ks)=f(i,k)
			    end do
		      end do
		end if !fin (if 2)
              if(n.eq.3)then !if 2b
		  ks=0
		  do i=1,ntau
		     xx=(tau(i)-x(1))/(x(2)-x(1))
		     xx2=xx*xx
		     f(i,1)=1.0-1.5*xx+0.5*xx2
	             f(i,2)=    2.0*xx- xx2
	             f(i,3)=   -0.5*xx+0.5*xx2
		     do k=1,n
		        ks=ks+1
		        ff(ks)=f(i,k)
		     end do
	          end do
              end if !fin (if 2b)
           end if !fin (if 1)
       else
c Interpolamos por rectas ---------------------------------
c        print*,'splines22 interpola lineal'
	if(n.gt.2)then

              i_inter=(ntau-1)/(n-1) !numero de puntos entre 2 nodos
                                     !incluye uno de los nodos
              do j=1,n
                 do i=1,ntau
                    f(i,j)=0.
                 enddo
              enddo

              do j=1,n-1
                 j1=j
                 j2=j+1
                 i1=(j-1)*i_inter+1
                 pasoj=x(j2)-x(j1)
                 do i=i1,i1+i_inter-1
                    taure(i)=(tau(i)-x(j1))/pasoj
                    f(i,j1)=1.0-taure(i)
                    f(i,j2)=taure(i)
c           print*,j1,j2,x(j1),x(j2),'--',i,taure(i),'---',f(i,j1),f(i,j2)
                 enddo
              enddo
              f(ntau,n)=1.

c filtramos las esquinas
              nw=2                   !anchura del filtro
              peso(1)=.11111
              peso(2)=.22222
              peso(3)=.33334
              peso(4)=.22222
              peso(5)=.11111
              if(i_inter.eq.1)then
                 nw=1
                 peso(1)=0.
                 peso(2)=.25
                 peso(3)=.50
                 peso(4)=.25
                 peso(5)=0.
              end if

              do j=1,n
                 do i=nw+1,ntau-nw
                    sum(i)=0.
                    do ii=-nw,nw
                       sum(i)=sum(i)+peso(ii+3)*f(i+ii,j)
                    enddo
                    f(i,j)=sum(i)
                 enddo
              enddo

	      ks=0
	      do i=1,ntau
	         do k=1,n
		    ks=ks+1
                    ff(ks)=f(i,k)
	         end do
	      end do


	   else

              paso=x(2)-x(1)

	      do i=1,ntau
	         taure(i)=(tau(i)-x(1))/paso
	      end do
	      if(no.eq.0)then
	          ks=0
	          do i=1,ntau
	             f(i,1)=1.0-taure(i)
	             f(i,2)=taure(i)
	             do k=1,n
		        ks=ks+1
                        ff(ks)=f(i,k)
	             end do
	           end do
	       else
	       ks=0
               do i=1,ntau
 	          if(taure(i).lt.1.0)then
		        f(i,1)=1.0-taure(i)
			f(i,2)=taure(i)
			f(i,3)=0.0
	           else
			f(i,1)=0.0
			f(i,2)=2.0-taure(i)
			f(i,3)=taure(i)-1.0
		   end if
		   do k=1,n
		      ks=ks+1
		      ff(ks)=f(i,k)
		   end do
		end do
	      end if
	    end if

	end if


c a partir de aqui lo hace siempre

	ks=0
	do i=1,ntau
	   yy(i)=0.0
	   do k=1,n
	      ks=ks+1
	      f(i,k)=ff(ks)
	      yy(i)=yy(i)+f(i,k)*y(k)
	   end do	
	end do

	return
	end


	
