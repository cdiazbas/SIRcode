
C igual que automatico pero lleva el factor multiplicatico cte
c de la presion en el ultimo nodo (como nodosp que cuelga de blends2)

	subroutine automaticop(mi,mp,ntotal4,difer,npos,rt,t)

	implicit real*4 (a-h,o-z) 
	include 'PARAMETER' !por kt,kn,kl y kld
	parameter (kld4=4*kld)

        integer npos(*),ntotal4,nodosposibles(kt),icalerr
        real*4 difer(*),rt(*),derivada(kt),t(*)
        real*4 x(kt),z(kt),f(kt,kt),tau(kt)


        common/nodosposibles/nummax,nodosposibles
        common/ndata/ndata !para fperfil2 y automatico
        common/tau/tau
        common/calerr/icalerr !si calculo errores=1 else =0

	data iprimeraveza/0/

c calculamos los posibles nodos 
       ntau=mi
        if(iprimeraveza.eq.0)then 
           iprimeraveza=1
           k=2
           nodosposibles(1)=1
           nodosposibles(2)=2
           do j=3,ntau
              resto=(ntau-1)-(j-1)*((ntau-1)/(j-1))
              if(resto.eq.0)then
                 k=k+1
                 nodosposibles(k)=j
              end if
           end do 
           nummax=k
        end if

c calculamos la derivada de la chi^2
c        print*,'automatico nmax ',mi
        do i=1,mi
           s=0.
           do j=1,ndata
              k=(i-1)*ntotal4+npos(j)
              s=s+difer(j)*rt(k)
           end do
           derivada(i)=s
c           print*,i,derivada(i)
        end do

c buscamos el numero de nodos optimo
        if(icalerr.eq.0)then
           if(mp.gt.1)then
             call criterio(ntau,derivada,mi,mp)
           else
             mi=1
           end if
        else
           print*,'nodos (aut) para errores',mi
        end if

c        print*,'cuantos pongo'
c        read*,mi

        if(mi.gt.1)then
	    m=(ntau-1)/(mi-1)   !n es el numero de taus
                            !m es el numero de pasos en cada subintervalo
	    do i=1,mi        !mi es el numero de nodos total i el indice del nodo
	       j=(i-1)*m+1   !j el indice en tau
	       x(i)=tau(j)   !x(i) es el valor de tau en el nodo
            end do

	    call splines22(x,x,mi-2,ntau,tau,z,f)
        end if
           
        call nodos_subp(rt,ntau,mi,f,ntotal4,t)
  
	return
        end
c _____________________________________________________
	subroutine nodos_subp(rt,n,no,f,ntotal4,t)

	implicit real*4 (a-h,o-z) 
	include 'PARAMETER' !por kt,kn,kl y kld

        integer no
        real*4 rt(*),f(kt,kt),gr(kt),t(*),y(kt)


	if(no.le.0)return

	if(no.eq.1)then   !si perturbacion constante
           ymedio=0
           do i=1,n
              ymedio=ymedio+t(i)
           end do
           ymedio=ymedio/float(n)

           do ilam=1,ntotal4
	      sumart=0.
              do i=1,n
                 kk=(i-1)*ntotal4+ilam
                 sumart=sumart+rt(kk)
              end do
              rt(ilam)=sumart*ymedio 
           end do
	else
c          do j=1,n
c            print*,'automatico ',j,t(j)
c          end do

	   m=(n-1)/(no-1)   !n es el numero de taus
                            !m es el numero de pasos en cada subintervalo
	   do i=1,no        !no es el numero de nodos total i el indice del nodo
	      j=(i-1)*m+1   !j el indice en tau
	      y(i)=t(j)    !y(i) es el valor del parametro en el nodo
c              print*,'automatico n no i j y',n,no,i,j,y(i)
           end do

           do ilam=1,ntotal4
	      do i=1,no
	         ggg=0.
	         do j=1,n
                    kk=(j-1)*ntotal4+ilam
	            ggg=ggg+rt(kk)*f(j,i)
                 end do
	         gr(i)=ggg
	      end do
	      do i=1,no
                 kk=(i-1)*ntotal4+ilam
 	         rt(kk)=gr(i)*y(no)
	      end do
           end do
c           print*,'automatico rt(',kk,')=',rt(kk)
	end if 
	return
	end
        
c _________________________________________________________________

