c PefrompgTmult calcula pe desde pg y T

	implicit real*4 (a-h,o-z)
        parameter (kk=100000)

	real*4 t(kk),pe(kk),pg(kk),ti,pei,psi,pp(10),pgasi  
        character*100 fichabun
        common/fichabun/fichabun
        common/constantes/g,avog	!gravedad,n. avogadro/pmu
	common/mu/cth

        common/precisoitera/precitera     
        common/nmaxitera/nmaxitera
        
        precitera=1.e-5
        nmaxitera=250
        
        fichabun='THEVENIN'

        open(1,file='entrada')
           read(1,*),n,(pg(i),i=1,n),(T(i),i=1,n)
        close(1)

c        do i=1,n
c           print*,i,t(i),pg(i)
c        enddo
               
        cth=1.
	g=cth*2.7414e+4		!gravedad cm/s^2 en fotosfera solar   
       	avog=6.023e23

        do i=1,n
           pgi=pg(i) 
           ti=t(i)
           pei=1.e-6*pgi
           iter=0
           dif=1.e8
	   do while(dif .gt. 1.e-3 .and. iter.lt.50)
             psi=pei
             call pefrompg11(ti,pgi,pei)  !en equisubmu.f
             call gasc(ti,pei,pgasi,pp)
             dif=abs(pei-psi)/(pei+psi)
             iter=iter+1
c            print*,iter,psi,pei,dif
             pe(i)=pei
          enddo
c          print*,pe(i)
       enddo
       open(1,file='salida')
           write(1,*),(pe(i),i=1,n)
       close(1)


       end 

