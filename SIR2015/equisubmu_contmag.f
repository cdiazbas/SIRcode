c EQUISUBMU_contmag evalua la presion en equilibrio hidrostatico
c teniendo en cuenta la presion magnetica
c_______________________________________________________________
c equisubmu rutina que evalua la presion en equilibrio hidrostatico
c teniendo en cuenta el numero de electrones para el calculo de mu
c Luis Bellot y Basilio Ruiz 27/7/95 
c Cambio del metodo de integracion Basilio Ruiz 3/9/96
c _______________________________________________________________
	subroutine equisubmu_contmag(ntau,tau1,t,pe,pg0,pg,z,ro,b,gamma)

	implicit real*4 (a-h,o-z)

	include 'PARAMETER'  
	parameter (nex=28,cgases=83145100.)
	real*4 tau(kt),t(kt),pe(kt),pg(*),kac,d2,x(kt),kappa(kt),taue(kt)
        real*4 z(kt),z1(kt),ro(*),y(kt),b(kt),gamma(*),ptot(kt)
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
        epsilon=.95              !cota  epsilon=max(Pmag/Ptot)
        pi=3.1415916

        do i=1,ntau
           tau(i)=tau1(i)+alog10(abs(cth))  !this is the optical depth in
        end do                        !the VERTICAL direction for 
	                              !hydrostatic equilib. computations

        tau0=1.
        imin=1!esto esta deshabilitado desde sir (entra 1)
        do i=2,ntau
           if(abs(tau(i)).lt.tau0)then
              imin=i         !zero of z is at log tau (vertical) =0
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

c inicializamos (atencion gamma debe entrar en radianes)
	tsi=t(ntau)
        psi=pe(ntau)
        pmagsi=b(ntau)*b(ntau)/8./pi*sin(gamma(ntau))*sin(gamma(ntau))
        if(pmagsi.gt.epsilon*pg0)pmagsi=epsilon*pg0

        pg00=pg0-pmagsi !nueva presion gaseosa
c        print*,'pg0 pg00 pmagsi epsilon',pg0,pg00,pmagsi,epsilon,sin(gamma(ntau))

        call pefrompg11(tsi,pg00,psi)
        pe(ntau)=psi
        call gasc(tsi,psi,pg(ntau),pp)
        pg(ntau)=pg00
        ptot(ntau)=pg0
	call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d2)
        pesomedio=pmusum/(asum+pp(8)) !peso molec. medio
        kappa(ntau)=kac*avog/pmusum
        ro(ntau)=pesomedio*pg(ntau)/tsi/cgases
        y(ntau)=taue(ntau)/kappa(ntau)/ro(ntau)

c integramos
        do i=ntau-1,1,-1
           ptotpr=ptot(i+1)+x(i+1)*(taue(i+1)/kappa(i+1))
	   tsi=t(i)

           pmagsi=b(i)*b(i)/8./pi*sin(gamma(i))*sin(gamma(i))
           if(pmagsi.gt.epsilon*ptotpr)pmagsi=epsilon*ptotpr
           pgpr=ptotpr-pmagsi
           call pefrompg11(tsi,pgpr,psi)
           pe(i)=psi

           call gasc(tsi,psi,pgpr,pp)
	   call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d2)
           kappa(i)=kac*avog/pmusum

           ptot(i)=ptot(i+1)+x(i+1)*(taue(i+1)/kappa(i+1)+taue(i)/kappa(i))/2.

           nit=1
           dif=2.*abs((ptot(i)-ptotpr))/(ptot(i)+ptotpr)

           do while (nit.lt.30.and.dif.gt.7.e-6)
              ptotpr=ptot(i)
              nit=nit+1
              pmagsi=b(i)*b(i)/8./pi*sin(gamma(i))*sin(gamma(i))
              if(pmagsi.gt.epsilon*ptotpr)pmagsi=epsilon*ptotpr
              pgpr=ptotpr-pmagsi

              call pefrompg11(tsi,pgpr,psi)
              pe(i)=psi
              call gasc(tsi,psi,pgpr,pp)
	      call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d2)
              kappa(i)=kac*avog/pmusum
            ptot(i)=ptot(i+1)+x(i+1)*(taue(i+1)/kappa(i+1)+taue(i)/kappa(i) )/2.
              dif=2.*abs((ptot(i)-ptotpr))/(ptot(i)+ptotpr)
            end do

            if(dif.gt.0.1) print*,'WARNING:
     & Hydrostatic equilibrium results in inaccurate electron pressures '
            pg(i)=pgpr
            call pefrompg11(tsi,pg(i),psi)
            pe(i)=psi
            ro(i)=pesomedio*pg(i)/tsi/cgases
            y(i)=taue(i)/kappa(i)/ro(i)
	end do
        z1(1)=0.
        do i=2,ntau
           z1(i)=z1(i-1)+x(i)/g*(y(i-1)+y(i))/2.
        end do
	z00=z1(imin)

        do i=1,ntau
           z(i)=(z1(i)-z00)*1.e-5
        end do
           
c        do i=1,ntau   !we output geometrical heights ALONG the LOS
c          z(i)=z(i)/abs(cth)
c        enddo


       return
       end 

