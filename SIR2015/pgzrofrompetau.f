c PGZROFROMPE calcula pg, z y ro desde pe
c_______________________________________________________________
c Basilio Ruiz 8/2/01 
c _______________________________________________________________
	subroutine pgzrofrompetau(ntau,tau,t,pe,pg,z,ro)

	implicit real*4 (a-h,o-z)

	parameter (nex=28,cgases=83145100.,nkt=1000)
	real*4 tau(*),t(*),pe(*),pg(*),kac,d2,x(nkt),kappa(nkt),taue(nkt)
        real*4 z(*),z1(nkt),ro(*),y(nkt)
    
	real*4 wgt,abu,ei1,ei2,pp(10),tsi,psi,d1(10)
        common/constantes/g,avog	!gravedad,n. avogadro/pmu
c        common/mu/cth                              !esto esta puesto a 1 en sir
        common/anguloheliocent/xmu 
        common/preciso/prec
               
	g=xmu*2.7414e+4		!gravedad cm/s^2 en fotosfera solar   
       	avog=6.023e23

        tau0=1.
        imin=1
        do i=2,ntau
           if(abs(tau(i)).lt.tau0)then
              imin=i
              tau0=tau(i)
           end if
        end do

        do i=1,10
           d1(i)=0
        end do
c calculamos el peso molecular medio pmu
	pmusum=0.0
        asum=0.0
	do i=1,92
  	   ii=i
	   call neldatb(ii,0.,wgt,abu,ei1,ei2)
	   pmusum=pmusum+wgt*abu
           asum=asum+abu
	end do

	do i=1,ntau
           taue(i)=10.**tau(i)
        end do

	do i=1,ntau-1
           x(i+1)=g*(tau(i)-tau(i+1))*2.3025851  
        end do

        do i=1,ntau
  	   tsi=t(i)
           psi=pe(i)
           call gasc(tsi,psi,pg(i),pp)
	   call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d2)
           pesomedio=pmusum/(asum+pp(8)) !peso molec. medio
           kappa(i)=kac*avog/pmusum
           ro(i)=pesomedio*pg(i)/tsi/cgases
           y(i)=taue(i)/kappa(i)/ro(i)
        end do

        z1(1)=0.
        do i=2,ntau
           z1(i)=z1(i-1)+x(i)/g*(y(i-1)+y(i))/2.
        end do 
        do i=1,ntau
           z(i)=(z1(i)-z1(imin))*1.e-5
        end do
           
       return
       end 

