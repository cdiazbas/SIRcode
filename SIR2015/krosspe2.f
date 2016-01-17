C KROSSPE evalua la opacidad de Rosseland a una T y Pe dadas
c escribe tambien la presion gaseosa, la densidad y el peso molecular medio   
c escribe tambien la presion de las particulas neutras (pn)
c y sqrt[(1.+M)/mui]=sqrt(1./mui + 1./mun) donde M=mui/mun con
c mui=peso del ion medio
c mun=peso de la particula neutra media


	implicit real*4 (a-h,o-z)

        real*4 t,pe,pp(10),pmusum,kappa,lambda,mui,mun,pn
        real*8 suma,kross,pesomedio,asum,cgases
        real*8 dt0,dt,dpe0,dpe,dpg,dpp(10),dro,dkappa
	character*100 fichabun
        common/fichabun/fichabun
        common/constantes/g,avog	!gravedad,n. avogadro

	g=2.7414e+4	!gravedad cm/s^2 en fotosfera solar   
       	avog=6.022045e23
        sigma=5.67e-5  !cte de Stefan Boltzmann erg /(s cm^2 K^4)
        cgases=83144087.d0
        pi=3.1415916

        fichabun='THEVENIN'

	pmusum=0.0
        asum=0.0
	do i=1,92
  	   ii=i
	   call neldatb(ii,0.,wgt,abu,ei1,ei2)
	   pmusum=pmusum+wgt*abu
           asum=asum+abu
	end do
c  pesomedio=pmusum/(asum+ne/nH) !ne/nH=pp(8)

        dpe0=1.e-6       			     
        dt0=2500.d0           			    
        open(2,file='fich_2500_22500.datos')         
 

         do k=1,151    !1,151 indices en pe  desde 10^-6 hasta 10^9 paso .1 
         
           dpe=dpe0*10**((k-1)*.1)
           pe=dpe
           print*,k		     
          do j=1,201 ! 1,201 T desde 2500 a 22500 con paso 100 
c                           do j=1,4
              dt=dt0+(j-1)*100.	
c                            dt=5040./(0.5+j*0.3)
              t=dt
              call gase(dt,dpe,dpg,dpp,mui,mun,pn)
	      
              pesomedio=pmusum/(asum+dpp(8))
	   
              dro=(dpg*pesomedio)/(dt*cgases)
              do iii=1,10
                 pp(iii)=dpp(iii)
              end do
              
  	      pasol=.05
              suma=0.
              kappa=1.e20
             do il=1,1001  !desde 50 A hasta 50000A con paso 50A
                lambda=(.05+(il-1.)*pasol)*1.e-5
                bt=dt_planck(t,lambda)  !derivada de B con respecto a T
                call k_nu(t,pe,pp,pmusum,lambda,kappa)
                dkappa=dble(kappa)
                if(dkappa.lt.1.d-40)dkappa=1.d-40
                if(dkappa.gt.1.d40)dkappa=1.d40
                suma=suma+(bt*pasol*1.d-5)/dkappa
             end do
             kross=(4.*sigma*(dt**3.)/pi)/suma
	     coef=sqrt(1./mui + 1./mun)
             write(2,*) dt,dpe,dpg,dro,pesomedio,kross,pn,coef
	      
cc	      print*,dlog10(dpg),5040./dt,mui,mui/mun,pg,15+0.93*dlog10(dpe/dpg)
c              if(j.eq.1. or. j.eq.131)then
c                 print*,k,j,t,pe,dpg,dro,pesomedio,kross
c              end if
           end do
	end do
        close(1)

        end


c________________________________________________________________

c K_NU evalua la opacidad a una T, ro y lambda(cm)  dadas

	subroutine k_nu(t,pe,pp,pmusum,lambda,kappa)

	implicit real*4 (a-h,o-z)

	real*4 pp(10),t,pe,d1(10),lambda,kappa,kac

        common/constantes/g,avog	!gravedad,n. avogadro
        
        do i=1,10
           d1(i)=0
        end do

	call kappach(lambda,t,pe,pp,d1,d1,kac,d2,d2)
        kappa=kac*avog/pmusum  !cm^2/gr

        return
        end

c________________________________________________________________

c	dt_planck es una funcion que calcula la derivada de la funcion
c	de planck con respecto a la temperatura

	real*4 function dt_planck(t,lambda)
	real*4 t,lambda,c1,c2,d4,b,ex,dplnck

	c1=1.1910627e-5		!1./(2*h*c**2) con lambda en cm
	c2=1.43879
	d4=lambda*lambda*lambda*lambda
	b=dplnck(t,lambda)
	ex=c2/(lambda*t)
        if(ex .lt. 9) then
      	   dt_planck=(c2*b*b*d4*exp(ex))/(c1*t*t)
        else
       	   dt_planck=(c2*b)/(lambda*t*t)
        end if

c        print*,lambda,ex,b,dt_planck

	return
	end
