c PETAUFROMPGZRO calcula pe y mu desde tau y t, pg
c_______________________________________________________________
c Basilio Ruiz 8/2/01 
c _______________________________________________________________
	subroutine pemufrompgtaut(ntau,tau,t,pe,pg,z,ro,pmu)
    
	implicit real*4 (a-h,o-z)

	parameter (nex=28,cgases=83145100.,nkt=1000)
	real*4 tau(*),t(*),pe(*),pg(*),kac,d2,kappa(nkt)
        real*4 z(*),ro(*),pmu(*)
    
	real*4 wgt,abu,ei1,ei2,pp(10),tsi,psi,d1(10)
        common/constantes/g,avog	!gravedad,n. avogadro/pmu
c        common/mu/cth                              !esto esta puesto a 1 en sir
        common/anguloheliocent/xmu 
        common/precisoitera/precitera 
        common/nmaxitera/nmaxitera
        
        precitera=1.e-5
        nmaxitera=250
               
	g=xmu*2.7414e+4		!gravedad cm/s^2 en fotosfera solar   
       	avog=6.023e23

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
  	   tsi=t(i)
           psi=pe(i)
           call pefrompg11(tsi,pg(i),psi)
           pe(i)=psi
           call gasc(tsi,psi,pg(i),pp)
	   call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d2)
           pesomedio=pmusum/(asum+pp(8)) !peso molec. medio
           kappa(i)=kac*avog/pmusum
           pmu(i)=pesomedio
c          print*,'pesomedio=',pmu(i)
           ro(i)=pesomedio*pg(ntau)/tsi/cgases
        end do

c        tau(ntau)=10.**(tau(ntau)) !cond. de contorno la tau superf. (entrada)
c        do i=ntau-1,1,-1
c           tau(i)=tau(i+1)+(kappa(i)*ro(i)+kappa(i+1)*ro(i+1))/2.*
c     &             (z(i+1)-z(i))*1.e5
c        end do 
c        do i=1,ntau
c           tau(i)=alog10(tau(i))
c        end do
           
       return
       end 

