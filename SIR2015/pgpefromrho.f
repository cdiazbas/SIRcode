c PGPEFROMRHO evalua la presion electronica y gaseosa en un punto 
c a partir de la temperatura y la densidad
c_______________________________________________________________
c Basilio Ruiz 28/03/2012 
c _______________________________________________________________
	subroutine pgpefromrho(t,ro,pe,pg)

	implicit real*4 (a-h,o-z)
        real*4 t,ro,pe,pg,pmusum,asum,pp(10),change
        parameter (nex=28,cgases=83145100.,avog=6.023e23)
        common/precisoitera/precitera      
        common/nmaxitera/nmaxitera
        precitera=1.e-5
        nmaxitera=250

c calculamos el peso molecular medio pmu
	pmusum=0.0
        asum=0.0
	do i=1,92
  	   ii=i
	   call neldatb(ii,0.,wgt,abu,ei1,ei2)
	   pmusum=pmusum+wgt*abu
           asum=asum+abu
	end do

c inicializamos pp(8)=pe/ph'=0.1
c        pp(8)=asum-pmusum/1.302
c        if(pp(8) .lt. 1.e-3)pp(8)=1.e-3
c        pesomedio=pmusum/(asum+pp(8))
        pesomedio=1.302
        pg=ro*cgases*t/pesomedio
c        pgold=pg
        call pefrompg11(t,pg,pe)
        call gasc(t,pe,pg,pp)
        if(pp(8) .lt. 1.e-20)pp(8)=1.e-20 
        if(pg .lt. 1.e-8)pg=1.e-8
        pe=pg*1.e-6
        change=100.
        iv=0
        do while (iv .lt. 30 .and. change .gt. 1.e-5) 
           iv=iv+1 
           pgold=pg
c          print*,iv,pesomedio,ro,pg,pe,change,pe/pg
c            print*,pg,pe,change,pe/pg

           pesomedio=pmusum/(asum+pp(8))
           pg=ro*cgases*t/pesomedio
           call pefrompg11(t,pg,pe)
           call gasc(t,pe,pg,pp)
           if(pp(8) .lt. 1.e-20)pp(8)=1.e-20 
           if(pg .lt. 1.e-8)pg=1.e-8
           change=abs(1.-pg/pgold)*100.
        end do
        if(pg. lt. 1.1e-8 .or. iv. eq. 30)then 
           print*,'Bad boundary condition in density'
           print*,'We chose pg=',500.
           pg=500.
c           stop
        endif  
c        print*,iv,pesomedio,ro,pg,pe,change,pe/pg
c        print*,pg,pe,change,pe/pg
c
c        print*,'estoy en  linea 35 pgpefromrho'
       
        return
        end