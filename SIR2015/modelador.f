c************************************************
		program MODELador
c************************************************
c
c Se trata de un programa que es capaz de fabricar un fichero
c perturbado en T,Pe,H,V a partir de un modelo tau,t,pe,magfield
c definiciones
	parameter (np=50000)	!numero maximo de puntos en tau
	character nombre*20,nomin*20
	real ttau(np),tt(np),ppe(np),hh(np),MMIC(NP),VVZ(NP),
     &       GG(NP),FFI(NP)
	real tau(np),t(np),pe(np),h(np),vz(np),MIC(NP),G(NP),FI(NP)
	real xa(11),ya(11)
c	________________________________________________________________
c se abre el fichero a modificar
11	print*,'Input model?:'
	read(*,'(a)')nomin
c	print*,'se espera una linea de cabecera : Vmac,fill '
c	print*,'los datos estan dados en tau (1),o en logtau (2): '
c	read*,niflog
        niflog=2
 
	open(2,file=nomin)
c leemos la cabecera
	read(2,*)vmac,fill,stray

c se lee contando las lineas
	num=0
	do while (num.lt.499)
	   num=num+1
	   read(2,*,end=10,err=11)ttau(num),tt(num),ppe(num),MMIC(NUM),
     $                            hh(num),VVZ(NUM),gg(num),ffi(num)
	end do
10	num=num-1
	close(2)

	if(niflog.eq.1)then
	     do i=1,num
		ttau(i)=alog10(ttau(i))
	     end do
	end if
	
c Que vamos a variar?
	print*,'Do you want to modify the step in log tau (yes=1,no=0): '
	read*,nvtau
	if(nvtau.eq.0)then
	   do i=1,num
	      tau(i)=ttau(i)
	      t(i)=tt(i)
	      pe(i)=ppe(i)
	      h(i)=hh(i)
              vz(i)=vvz(i)
	      mic(i)=mmic(i)
              g(i)=gg(i)
              fi(i)=ffi(i)
	   end do
	   n=num
	else
c interpolaremos las presiones en logaritmos neperianos
	    do i=1,num
	       ppe(i)=alog(ppe(i))
	    end do

c     	print*,'paso igual a cero no equiespaciado'
    	print*,'Give the initial and the final log tau, and the step (1.2,-4,.1): '
	    read*,tau1,taun,paso
	    paso=-paso
	    if(paso.eq.0)then
	       print*,'Number of points in log tau'
	       read*,n
	       do i=1,n
		  print*,'tau numero ',i,' ='
		  read*,tau(i)
	       end do
	    else
	       n=nint((taun-tau1)/paso)+1
c definimos la red en tau
	       do i=1,n
	          tau(i)=tau1+(i-1)*paso
	       end do
	    end if
c interpolamos
	    print*,'Degree of the polynomial to be used for interpolating :'
	    read*,ngrado
c            ngrado=2
	    n2=int(ngrado/2)
	
	    do i=1,n
	       CALL LOCATE(TTAU,NUM,TAU(I),J)
	       n3=j-n2-1
               if(n3.lt.0)n3=0
               if(n3+ngrado+1.gt.num)n3=num-ngrado-1
	       do k=1,ngrado+1
		     xa(k)=ttau(n3+k)
	       end do
	       do k=1,ngrado+1
		     ya(k)=tt(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(I),T(I),ERROR)
	       do k=1,ngrado+1
		     ya(k)=ppe(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(I),pe(I),ERROR)
	       do k=1,ngrado+1
		     ya(k)=hh(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(I),h(I),ERROR)
	       do k=1,ngrado+1
	             ya(k)=vvz(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(I),vz(I),ERROR)
	       do k=1,ngrado+1
		     ya(k)=mmic(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(I),mic(I),ERROR)
	       do k=1,ngrado+1
		     ya(k)=gg(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(I),g(I),ERROR)
	       do k=1,ngrado+1
		     ya(k)=ffi(n3+k)
	       end do
	       CALL POLINT(XA,YA,NGRADO+1,TAU(I),fi(I),ERROR)
	    end do
	    do i=1,n
		pe(i)=exp(pe(i))
	    end do
	end if

        npar=-1
        do while(npar.ne.0)
	print*,'Select the physical quantity to be modified:'
	print*,'None=0,t=1,p=2,h=3,vz=4,mic=5,gamma=6,phi=7 :'
	read*,npar
        if(npar.eq.0)goto 101
c	print*,'entre que taus? (entre 1 y ',n,'=todos)'
c	read*,i1,i2
        i1=1
        i2=n
	print*,'xx=a+b*x+c*log10(tau) dame a , b ,c '
	read*,a,b,c
      
	do i=i1,i2
	   if(npar.eq.1)t(i)=a+b*t(i)+c*tau(i)
	   if(npar.eq.2)pe(i)=a+b*pe(i)+c*tau(i)
	   if(npar.eq.3)h(i)=a+b*h(i)+c*tau(i)
	   if(npar.eq.4)vz(i)=a+b*vz(i)+c*tau(i)
	   if(npar.eq.5)mic(i)=a+b*mic(i)+c*tau(i)
	   if(npar.eq.6)g(i)=a+b*g(i)+c*tau(i)
	   if(npar.eq.7)fi(i)=a+b*fi(i)+c*tau(i)
	end do
101     end do
c se abre el fichero en donde se escribira el modelo
	print*,'Output model? '
	read(*,'(a)')nombre
	open(1,file=nombre)
	
	write(1,*)vmac,fill,stray
	do i=1,n
      	    write(1,100)tau(i),t(i),pe(i),mic(i),h(i),vz(i),g(i),fi(i)
	end do
100     format(1x,f7.4,1x,f8.1,1x,1pe12.5,1x,e10.3,1x,4(e11.4,1x))

	close(1)
	end
