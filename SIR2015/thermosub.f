c thermo calcula los coeficientes termodinamicos:
c alpha,delta, cv
c tambien calcula pg, rho, pesomolecularmedio
c Basilio  11 abril 2013

       subroutine thermosub(t,pe,thermo)
       
       parameter (cgases=83145100.,boltz=1.38054e-16)
       
       real*4 theta,t,pe
       real p(99),dp(99),ddp(99)
       real*4 u_ion,du_ion,ddu_ion,thermo(14)
       real*4 pg,dpg,ddpg,grado1,grado,alpha,delta
       real*4 pmusum,asum,pesomedio,rho
       real*4 uint_neutro,uint,dluint,ddluint
       real*4 cv,ga31,ga1,nabad,cs,cp
       
       theta=5040./t
       call gasb_thermo(theta,pe,p,dp,ddp,u_ion,du_ion,ddu_ion)

       pg=p(84)      ! gas pressure
       dpg=dp(84)    ! dlog(pg)/dT
       ddpg=ddp(84)  ! ddlog(pg)/dpe
       
       if(pg-pe .le. 0)then
         print*,'WARNING: total pressure lower than pe'
         print*,'have been found in subroutine thermo.f'
         stop
       end if
       
       grado1=pg/(pg-pe)    !1+grado ionizacion
       grado=grado1-1.      !grado de ionizacion 
 
c calculamos la compresibilidad y el coef. de dilatacion
       alpha=grado1*(1.-1./(pg*ddpg))        !dlog(rho)/dlog(Pg) a T constante
       
       delta=1.-t*dpg*grado1/pg/ddpg          !-dlog(rho)/dlog(T) a Pg constante
       
c calculamos el peso molecular medio pmu
	pmusum=0.0
        asum=0.0
	do i=1,92
  	   ii=i
	   call neldatb(ii,0.,wgt,abu,ei1,ei2)
	   pmusum=pmusum+wgt*abu
           asum=asum+abu
	end do
	
		
c        pesomedio=pmusum/(asum+p(90)) !peso molec. medio
        pesomedio0=pmusum/asum        !peso molec. medio por atomo
        pesomedio=pesomedio0/grado1
        
c        print*,pmusum/(asum+p(90)),pesomedio,pmusum/(asum+p(90))/pesomedio

c calculamos la densidad y las derivadas de su logaritmo a T y Pe
        rho=pesomedio*pg/t/cgases
        drho=grado1*dpg-1./t
        ddrho=grado1*(ddpg-1./pg)
                       
c calculamos la energía interna por unidad de masa  y sus derivadas
c        uint_neutro=1.5*pg/rho
        
c calculamos la energia interna de ionizacion por unidad de masa
        u_ion_m=u_ion/rho

        du_ion_m=du_ion-ddu_ion*drho/ddrho
        
c calculamos el calor especifico a volumen constante        
        cv=1.5*pg*delta/(rho*t*alpha)+u_ion_m*du_ion_m
       
c calculamos ga31=gamma3-1, ga1=gamma1 , nabad=nabla adibatico , cs=velocidad sonido 
       ga31=(pg*delta)/(rho*T*alpha*cv)  ! gamma3-1= dlnT/dlnrho a S cte       
       chit=delta/alpha                  !dlnP/dlnT a rho cte
       chir=1./alpha                     !dlnP/dlnrho a T cte
       ga1=chir+chit*ga31                !dlnP/dlnrho a S cte

       nabad=ga31/ga1
       if(ga1 .lt.0 )print*,' gamma1=',ga1,' delta=',delta,' gamma3-1=',ga31,' alpha=',alpha
       cs=sqrt(ga1*pg/rho)
       
c Calculamos el calor especifico a presion cte 
       cp=pg/(rho*t)*chit*chit/chir+cv
c       gam=cp/cv
       
       thermo(1)=t
       thermo(2)=pg
       thermo(3)=rho
       thermo(4)=pesomedio
       thermo(5)=p(91)                        ! ne
       thermo(6)=p(85)*p(1)/(boltz*t)   ! nH neutro
       thermo(7)=p(85)*p(86)/(boltz*t)  ! nH+
       thermo(8)=grado                        !grado de ionizacion
       thermo(9)=alpha
       thermo(10)=delta
       thermo(11)=cv
       thermo(12)=cp
       thermo(13)=cs                          !velocidad sonido adibatica
       thermo(14)=nabad                       !nabla adibatica
       
       
	return
	end
       
       
