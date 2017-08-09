c   ztaufrompetro calcula z y tau desde pg, z y ro
c_______________________________________________________________
c Basilio Ruiz 8/2/01 
c _______________________________________________________________
	subroutine ztaufrompetro(ntau,tau,t,pe,pg,z,ro)

	implicit real*4 (a-h,o-z)

	parameter (nex=28,cgases=83145100.,nkt=1000)
	real*4 tau(*),t(*),pe(*),pg(*),kac,d2,kappa(nkt)
        real*4 z(*),ro(*)
    
	real*4 wgt,abu,ei1,ei2,pp(10),tsi,psi,d1(10)
        common/constantes/g,avog	!gravedad,n. avogadro/pmu
c        common/mu/cth                              !esto esta puesto a 1 en sir
        common/anguloheliocent/xmu 
        common/precisoitera/precitera 
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
c           ro(i)=pesomedio*pg(ntau)/tsi/cgases
        end do

        tau(ntau)=10.**(tau(ntau)) !cond. de contorno la tau superf. (entrada)
        do i=ntau-1,1,-1
           tau(i)=tau(i+1)+(kappa(i)*ro(i)+kappa(i+1)*ro(i+1))/2.*
     &             (z(i+1)-z(i))*1.e5
        end do 
        do i=1,ntau
           tau(i)=alog10(tau(i))
        end do
           
       return
       end 

