c EQUISUBMU evalua la presion en equilibrio hidrostatico
c_______________________________________________________________
c equisubmu rutina que evalua la presion en equilibrio hidrostatico
c teniendo en cuenta el numero de electrones para el calculo de mu
c Luis Bellot y Basilio Ruiz 27/7/95 
c Cambio del metodo de integracion Basilio Ruiz 3/9/96
c _______________________________________________________________
	subroutine equisubmu_noiter(ntau,tau1,t,pe,pg,z,ro)

	implicit real*4 (a-h,o-z)

	include 'PARAMETER'   !solo por kt, que era igual a 1000
	parameter (nex=28,cgases=83145100.)
        real*4 tau1(kt)
	real*4 tau(kt),t(*),pe(*),pg(kt),kac,d2,x(kt),kappa(kt),taue(kt)
	real*4 p99(99),dp99(99),ddp99(99)  !NEW no iter
        real*4 z1(kt),z(kt),ro(kt),y(kt)
	real*4 mu
    
	real*4 wgt,abu,ei1,ei2,pp(10),tsi,psi,d1(10)
        common/constantes/g,avog	!gravedad,n. avogadro/pmu
	common/mu/cth                   !esto esta deshabilitado desde sir (entra 1)
        common/preciso/prec      
        common/anguloheliocent/mu
       
	g=mu*2.7414e+4		!gravedad cm/s^2 en fotosfera solar   
       	avog=6.023e23

        do i=1,ntau
           tau(i)=tau1(i)+alog10(cth)  !this is the optical depth in 
        end do                        !the VERTICAL direction for 
	                              !hydrostatic equilib. computations
        tau0=1.
        imin=1
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
            x(i+1)=g*(tau(i)-tau(i+1))*2.30259  
        end do

c inicializamos 
	tsi=t(ntau)
        psi=pe(ntau)
        call gasc(tsi,psi,pg(ntau),pp)
	call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d2)
c        theta=5040./tsi
c        call gasb(theta,pe,p99,dp99,ddp99)
c	call kappach_noiter(5.e-5,tsi,psi,p99,dp99,ddp99,kac,d2,d2)
	
        pesomedio=pmusum/(asum+pp(8)) !peso molec. medio
        kappa(ntau)=kac*avog/pmusum
        ro(ntau)=pesomedio*pg(ntau)/tsi/cgases
        y(ntau)=taue(ntau)/kappa(ntau)/ro(ntau)

c integramos
        do i=ntau-1,1,-1
           pgpr=pg(i+1)+x(i+1)/kappa(i+1)
	   tsi=t(i)

           call pefrompg10(tsi,pgpr,psi)
           pe(i)=psi

c          call gasc(tsi,psi,pgpr,pp)                                      !NEW NO iter
c	   call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d2)                  !NEW NO iter
c          pesomedio=pmusum/(asum+pp(8)) !peso molec. medio                !NEW NO iter
	   theta=5040./tsi                                                 !NEW NO iter
           call gasb(theta,pe,p99,dp99,ddp99)                              !NEW NO iter
	   call kappach_noiter(5.e-5,tsi,psi,p99,dp99,ddp99,kac,d2,d2)     !NEW NO iter
	   pesomedio=pmusum/(asum+p99(90)) !peso molec. medio              !NEW NO iter
           
           kappa(i)=kac*avog/pmusum

           pg(i)=pg(i+1)+x(i+1)*(taue(i+1)/kappa(i+1)+taue(i)/kappa(i))/2.

           nit=1
           deltaPg=pg(i)-pgpr
           dif=2.*abs(deltaPg)/(pg(i)+pgpr)
           pe_old=pe(i)

           do while (nit.lt.20.and.dif.gt.1.e-4)
              nit=nit+1
              pgpr=pg(i)
              
              psi=pe_old+deltaPg/ddp99(84)
              print*,nit,pe_old,psi,dif
              
c             call pefrompg10(tsi,pgpr,psi)
              pe(i)=psi
              pe_old=psi
              
c             call gasc(tsi,psi,pgpr,pp)                                      !NEW NO iter
c	      call kappach(5.e-5,tsi,psi,pp,d1,d1,kac,d2,d2)                  !NEW NO iter  
c             pesomedio=pmusum/(asum+pp(8)) !peso molec. medio                !NEW NO iter  
              call gasb(theta,pe,p99,dp99,ddp99)                              !NEW NO iter
	      call kappach_noiter(5.e-5,tsi,psi,p99,dp99,ddp99,kac,d2,d2)     !NEW NO iter
              pesomedio=pmusum/(asum+p99(90)) !peso molec. medio              !NEW NO iter

              kappa(i)=kac*avog/pmusum
              pg(i)=pg(i+1)+x(i+1)*(taue(i+1)/kappa(i+1)+taue(i)/kappa(i))/2.
              deltaPg=pg(i)-pgpr
              dif=2.*abs(deltaPg)/(pg(i)+pgpr)            
c              dif=2.*abs((pg(i)-pgpr))/(pg(i)+pgpr)
            end do

            if(dif.gt.0.1)print*,'WARNING: 
     & Hydrostatic equilibrium results in inaccurate electron pressures '
            call pefrompg10(tsi,pg(i),psi)
            pe(i)=psi
            ro(i)=pesomedio*pg(i)/tsi/cgases
            y(i)=taue(i)/kappa(i)/ro(i)
	end do
        z1(1)=0.
        do i=2,ntau  !this is the z scale along the vertical direction
           z1(i)=z1(i-1)+x(i)/g*(y(i-1)+y(i))/2.
        end do 
        do i=1,ntau
           z(i)=(z1(i)-z1(imin))*1.e-5
        end do

c        do i=1,ntau   !we output geometrical heights ALONG the LOS
c          z(i)=z(i)/mu
c        enddo


       return
       end 

