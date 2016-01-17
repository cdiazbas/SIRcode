c EQUISUBMU_cont evalua la presion en equilibrio hidrostatico
c_______________________________________________________________
c equisubmu rutina que evalua la presion en equilibrio hidrostatico
c teniendo en cuenta el numero de electrones para el calculo de mu
c Luis Bellot y Basilio Ruiz 27/7/95 
c Cambio del metodo de integracion Basilio Ruiz 3/9/96
c _______________________________________________________________
	subroutine equisubmu_cont(ntau,tau1,t,pe,pg0,pg,z,ro)

	implicit real*4 (a-h,o-z)

	include 'PARAMETER'   !solo por kt, que era igual a 1000
	parameter (nex=28,cgases=83145100.)
	real*4 tau(kt),t(kt),pe(kt),pg(kt),kac,d2,x(kt),kappa(kt),taue(kt)
        real*4 z(kt),z1(kt),ro(*),y(kt)
        real*4 tau1(kt)
  	real*4 mu
	integer*4 nmaxitera

 
	real*4 wgt,abu,ei1,ei2,pp(10),tsi,psi,d1(10)
        common/constantes/g,avog	!gravedad,n. avogadro/pmu
	common/mu/cth                   !esto esta deshabilitado desde sir (entra 1)
        common/precisoitera/precitera      
        common/anguloheliocent/mu
        common/nmaxitera/nmaxitera
        precitera=1.e-5
        nmaxitera=250
               
	g=mu*2.7414e+4		!gravedad cm/s^2 en fotosfera solar
       	avog=6.023e23

        do i=1,ntau
           tau(i)=tau1(i)+alog10(abs(cth))  !this is the optical depth in
        end do                        !the VERTICAL direction for
	                              !hydrostatic equilib. computations

        tau0=1.
        imin=1
        do i=2,ntau
           if(abs(tau1(i)).lt.tau0)then
              imin=i        !zero of z is at log tau (LOS) =0
              tau0=tau1(i)
           end if
        end do
c        print*,'imin tau(imin)=',imin,tau0,tau1(imin)

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

c inicializamos 
	tsi=t(ntau)
        psi=pe(ntau)

        call pefrompg11(tsi,pg0,psi)
        pe(ntau)=psi
        call gasc(tsi,psi,pg(ntau),pp)
        pg(ntau)=pg0
	call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d2)
        pesomedio=pmusum/(asum+pp(8)) !peso molec. medio
        kappa(ntau)=kac*avog/pmusum
        ro(ntau)=pesomedio*pg(ntau)/tsi/cgases
        y(ntau)=taue(ntau)/kappa(ntau)/ro(ntau)

c integramos
        do i=ntau-1,1,-1
           pgpr=pg(i+1)+x(i+1)*(taue(i+1)/kappa(i+1))
	   tsi=t(i)

           call pefrompg11(tsi,pgpr,psi)
           pe(i)=psi

           call gasc(tsi,psi,pgpr,pp)
	   call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d2)
           kappa(i)=kac*avog/pmusum

           pg(i)=pg(i+1)+x(i+1)*( taue(i+1)/kappa(i+1) + taue(i)/kappa(i) )/2.

           nit=1
           dif=2.*abs((pg(i)-pgpr))/(pg(i)+pgpr)

           do while (nit.lt.30.and.dif.gt.1.e-5)
              nit=nit+1
              pgpr=pg(i)
              call pefrompg11(tsi,pgpr,psi)
              pe(i)=psi
              call gasc(tsi,psi,pgpr,pp)
	      call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d2)
              kappa(i)=kac*avog/pmusum
              pg(i)=pg(i+1)+x(i+1)*(taue(i+1)/kappa(i+1)+taue(i)/kappa(i) )/2.
              dif=2.*abs((pg(i)-pgpr))/(pg(i)+pgpr)
            end do

            if(dif.gt.0.1) print*,'WARNING:
     & Hydrostatic equilibrium is resulting in inaccurate electronic pressures '
            call pefrompg11(tsi,pg(i),psi)

            pe(i)=psi
            ro(i)=pesomedio*pg(i)/tsi/cgases

            y(i)=taue(i)/kappa(i)/ro(i)
c	    print*,'equisubmu_cont 3',i,'pg(i) ro kappa= ',pg(i),ro(i),pesomedio,tsi,cgases,kappa(i),y(i)
	end do

        z1(1)=0.
        do i=2,ntau   !this is the z scale along the vertical direction,: NO MORE, NOW its along line of sight
           z1(i)=z1(i-1)+x(i)/g*(y(i-1)+y(i))/2.
        end do
	z00=z1(imin)
	
        do i=1,ntau
           z(i)=(z1(i)-z00)*1.e-5
c           print*,i,z(i)
        end do
c        print*,imin,abs(tau(imin)),z(imin)
c        do i=1,ntau   !we output geometrical heights ALONG the LOS
c          z(i)=z(i)/abs(cth)
c        enddo
         
       return
       end 

