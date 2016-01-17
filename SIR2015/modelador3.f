c************************************************
		program MODELador3
c************************************************
c
c Se trata de un programa que es capaz de fabricar un fichero
c perturbado en T,Pe,H,V a partir de un modelo tau,t,pe,magfield
c definiciones
	parameter (np=1000)	!numero maximo de puntos en tau
	character nombre*20,nomin*20
        character nombreth*20,nombrethermo*24,chathermo*4
	real*4 ttau(np),tt(np),ppe(np),hh(np),MMIC(NP),VVZ(NP),
     &       GG(NP),FFI(NP),pmu(np)
	real*4 tau(np),t(np),pe(np),h(np),vz(np),MIC(NP),G(NP),FI(NP)
	real*4 z(np),pg(np),ro(np),zz(np),ppg(np),rro(np),gamma(np),nabla(np)
	real*4 an(10),bn(10),cn(10)
	real*4 xa(11),ya(11),xmu,cth
	real*4 thermo(14)
        real*4 mreadr2,mreadr3
        integer mreadi2,mreadi3,meves,i1,i2,inada
        character*100 mreadc2,fcontrol
	character*100 fichabun
	
	common/mu/cth
        common/anguloheliocent/xmu    
        common/fichabun/fichabun
        common/preciso/prec 
        

c       EXTERNAS
	external mreadi2,mreadr2,mreadc2,meves,mreadi3,mreadr3

c	________________________________________________________________

	pi=3.1415916
        ici=1
	ican=23	        !canal de lectura del fichero de control
        prec=1.e-6
        
        print*,'control file'
        read*,fcontrol  
	open(ican,file=fcontrol)
   
        nomin=mreadc2(ican,ici)   !Input model?:
        nombre=mreadc2(ican,ici)  !Ouput model?:


c leemos la cabecera; se lee contando las lineas
 	open(2,file=nomin)
	read(2,*)vmac,fill,stray
 	num2=0
	do while (num2.le.np-1)
	   num2=num2+1
	   read(2,*,end=110,err=111)ttau(num2),tt(num2),ppe(num2),mmic(num2),
     $              hh(num2),vvz(num2),gg(num2),ffi(num2),zz(num2),ppg(num2),rro(num2)
	end do
110	num2=num2-1
	close(2)

c comprobamos si el fichero contiene z,pg,ro
 	open(2,file=nomin)
	read(2,*)vmac,fill,stray
 	num=0
	do while (num.le.np-1)
	   num=num+1
	   read(2,*,end=10,err=11)ttau(num),tt(num),ppe(num),mmic(num),
     $                            hh(num),vvz(num),gg(num),ffi(num)
	end do
10	num=num-1
	close(2)

        longversion=0
        if(num2.eq.num)longversion=1  !si contiene pg,z,ro

        tau1i=ttau(1)
        tauni=ttau(num)
        pasoi=ttau(2)-ttau(1)

        call mreadr33(ican,tau1i,tauni,pasoi,tau1,taun,paso)

        nvtau=0
        if(tau1.ne.tau1i.or.taun.ne.tauni.or.paso.ne.pasoi)nvtau=1
        paso=-paso

	ngrado=mreadi3(ican,ici,2) !Degree of the polynomial to be used for 
                                   ! interpolating       
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

	   if(longversion.eq.1)then
	      do i=1,num
	         z(i)=zz(i)
	         pg(i)=ppg(i)
	         ro(i)=rro(i)
	      end do
           end if

	   n=num
	else
c interpolaremos las presiones (y/o densidades) en logaritmos neperianos
	   do i=1,num
	      if(ppe(i).gt.0)ppe(i)=alog(ppe(i))
	   end do
	   if(longversion.eq.1)then
	      do i=1,num
	         if(ppg(i).gt.0)ppg(i)=alog(ppg(i))
	         if(rro(i).gt.0)rro(i)=alog(rro(i))
	      end do
           end if

	   n=nint((taun-tau1)/paso)+1
           if(n.gt.np)then
             print*,' '
             print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
             print*,'number of grid points larger than ',np
             n=np
             num=n
             print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
             print*,' '
          end if
          npaso=nint(1000.*(taun-tau1)/(n-1))
          paso=npaso/1000.
          taun=tau1+(n-1)*paso

          print*,' '
          print*,'        the new grid has         ',n,' grid points'
          print*,'        the new values are tau1= ',tau1,' taun=',taun,' |step|=',-paso
          print*,' '


c definimos la red en tau
	   do i=1,n
	      tau(i)=tau1+(i-1)*paso
c	      print*,'nueva red en tau ',tau(i)
	   end do

	    n2=int(ngrado/2)
	
	    do i=1,n
	       call locate(ttau,num,tau(i),j)
	       n3=j-n2-1
               if(n3.lt.0)n3=0
               if(n3+ngrado+1.gt.num)n3=num-ngrado-1
	       do k=1,ngrado+1
		  xa(k)=ttau(n3+k)
		  ya(k)=tt(n3+k)
	       end do
	       call polint(xa,ya,ngrado+1,tau(i),t(i),error)
	       do k=1,ngrado+1
		  ya(k)=ppe(n3+k)
	       end do
	       call polint(xa,ya,ngrado+1,tau(i),pe(i),error)
	       do k=1,ngrado+1
		  ya(k)=hh(n3+k)
	       end do
	       call polint(xa,ya,ngrado+1,tau(i),h(i),error)
	       do k=1,ngrado+1
	          ya(k)=vvz(n3+k)
	       end do
	       call polint(xa,ya,ngrado+1,tau(i),vz(i),error)
	       do k=1,ngrado+1
		  ya(k)=mmic(n3+k)
	       end do
	       call polint(xa,ya,ngrado+1,tau(i),mic(i),error)
	       do k=1,ngrado+1
		  ya(k)=gg(n3+k)
	       end do
	       call polint(xa,ya,ngrado+1,tau(i),g(i),error)
	       do k=1,ngrado+1
		  ya(k)=ffi(n3+k)
	       end do
	       call polint(xa,ya,ngrado+1,tau(i),fi(i),error)
	       if(longversion.eq.1)then
	          do k=1,ngrado+1
		     ya(k)=zz(n3+k)
	          end do
	          call polint(xa,ya,ngrado+1,tau(i),z(i),error)
	          do k=1,ngrado+1
		     ya(k)=ppg(n3+k)
	          end do
	          call polint(xa,ya,ngrado+1,tau(i),pg(i),error)
	          do k=1,ngrado+1
		     ya(k)=rro(n3+k)
	          end do
	          call polint(xa,ya,ngrado+1,tau(i),ro(i),error)
               end if
	    end do
	    do i=1,n
	       pe(i)=exp(pe(i))
	    end do
	    if(longversion.eq.1)then
	       do i=1,n
	          pg(i)=exp(pg(i))
	          ro(i)=exp(ro(i))
	       end do
            end if
	end if

c i=1...10 corresponde a t,pe,mic,h,vz,g,f,z,pg,ro
c la transformacion es t(i)=an(1)+bn(1)*t(i)+cn(1)*tau(i)
c la transformacion es pe(i)=an(2)+bn(2)*pe(i)+cn(2)*tau(i) ,etc
c la transformacion se hace entre i1 e i2
        call mreadi33(ican,1,n,0,i1,i2,inada) 

        do j=1,10
           call mreadr33(ican,0.,1.,0.,an(j),bn(j),cn(j)) 
        end do
	do i=i1,i2
	   t(i)  =an(1) + bn(1)*t(i)  + cn(1)*tau(i)
	   pe(i) =an(2) + bn(2)*pe(i) + cn(2)*tau(i)
	   mic(i)=an(3) + bn(3)*mic(i)+ cn(3)*tau(i)
	   h(i)  =an(4) + bn(4)*h(i)  + cn(4)*tau(i)
	   vz(i) =an(5) + bn(5)*vz(i) + cn(5)*tau(i)
	   g(i)  =an(6) + bn(6)*g(i)  + cn(6)*tau(i)
	   fi(i) =an(7) + bn(7)*fi(i) + cn(7)*tau(i)
	end do
	if(longversion.eq.1)then
	   do i=i1,i2
	      z(i) =an(8) + bn(8)*z(i)  + cn(8)*tau(i)
	      pg(i)=an(9) + bn(9)*pg(i) + cn(9)*tau(i)
	      ro(i)=an(10)+ bn(10)*ro(i)+ cn(10)*tau(i)
           end do
        end if

c calculamos el equilibrio hidrostatico?
        iee=mreadi3(ican,ici,0)
        ipgmag  =mreadi3(ican,ici,0)      !terminos en P magnetica?
        if(ipgmag.eq.1)then
           do i=1,n
              gamma(i)=g(i)*pi/180.
           end do
        end if
        itautoz_or_ztotau=mreadi3(ican,ici,1)
        print*,'itautoz_or_ztotau=',itautoz_or_ztotau

        fichabun=mreadc2(ican,ici)   !Abundance file?

c        xmu=mreadr3(ican,ici,1.0)
        cth=mreadr3(ican,ici,1.0)
        xmu=cth


        if(iee.eq.1.or.longversion.eq.0)then
           pg0=mreadr3(ican,ici,pg(n))
           print*,'_________________________________________________________ '
           print*,'evaluation of pe pg z and rho from hydrostatic equilibrium'
	   if(pg0.eq.0)pg0=1.e2
           ro0=0.
           if(pg0.gt.0)then
              print*,'Boundary condition in gas pressure'
	      print*,'gas pressure boundary condition =',pg0
           end if
           if(pg0.lt.0.)then              
              ro0=abs(pg0)
              pg0=100.
              print*,'Boundary condition in density'
	      print*,'density boundary condition =',ro0
              ro(n)=ro0
              call pgpefromrho(t(n),ro(n),pe(n),pg(n))
              pg0=pg(n)
              print*,'Gas pressure at surface =',pg0
           end if
   
           if(ipgmag.ne.1)call equisubmu_cont(n,tau,t,pe,pg0,pg,z,ro)
	   
           if(ipgmag.eq.1)print*,'taking into account magnetic pressure'
           if(ipgmag.eq.1)call equisubmu_contmag(n,tau,t,pe,pg0,pg,z,ro,h,gamma)
           print*,'_________________________________________________________ '
           print*,' '
        else
          print*,'itautoz_or_ztotau=',itautoz_or_ztotau
          if(itautoz_or_ztotau.eq.2)then
              !calcular pe,tau desde pg,z
              print*,'_________________________________________ '
              print*,'evaluation of pe and tau from pg z and rho'
	      call petaufrompgzro(n,tau,t,pe,pg,z,ro)
              print*,'_________________________________________ '
              print*,' '
c         else if(itautoz_or_ztotau.eq.4)then
c              !calcular tau,z desde pe,rho,T
c              print*,'_________________________________________ '
c              print*,'evaluation of z and tau from T pe and rho'
c	      call ztaufrompetro(n,tau,t,pe,pg,z,ro)
c              print*,'_________________________________________ '
c              print*,' '
           else if (itautoz_or_ztotau.eq.1)then
              !calcular pg,z,ro desde pe,tau
              print*,'__________________________________________ '
              print*,'evaluation of pg z and rho from pe and tau'
	      call pgzrofrompetau(n,tau,t,pe,pg,z,ro)
              print*,'__________________________________________ '
              print*,' '
           else if (itautoz_or_ztotau.eq.3)then
           print*,'_________________________________________ '
              print*,'evaluation of mu pe from tau t pg'
	      call pemufrompgtaut(n,tau,t,pe,pg,z,ro,pmu)
              print*,'_________________________________________ '
              print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
              print*,' WARNING:in the output model MU is written instead micro'
              print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! '
       
              
              do i=1,n
c                print*,'pesomedio=',pmu(i)
                mic(i)=pmu(i)
              enddo   
c              do i=1,n
c      	         print*,tau(i),t(i),pe(i),pmu(i),h(i),vz(i),
c     &                      g(i),fi(i),z(i),pg(i),ro(i)
c             
c	      end do
c	      stop
           end if
        end if

        open(1,file=nombre)	
	write(1,*)vmac,fill,stray
        do i=1,n
      	       write(1,1000)tau(i),t(i),pe(i),mic(i),h(i),vz(i),
     &                      g(i),fi(i),z(i),pg(i),ro(i)
	end do	
        close(1)

	
c	Se define el nombre del fichero con los parametros termodinamicos :
	chathermo='.th'
        nombreth=nombre
	call quitaex(nombreth)
	call concatena(nombreth,chathermo,nombrethermo)
	
	print*,'fichero output con parametros thermodinamicos=',nombrethermo
	print*,' '
	print*,'logtau pesomedio ne nH nH+ grado_ioniz alpha delta c_v c_p c_s nabla_adiab'
   	print*,' '
   	
   	do i=2,n-1
	  nabla(i)=(alog(t(i+1))-alog(t(i-1)))/(alog(pg(i+1))-alog(pg(i-1)))
	enddo
	nabla(1)=nabla(2)
	nabla(n)=nabla(n-1)
   	
	open(2,file=nombrethermo)
	write(2,1002),'logtau','pesomed','  ne   ','  nH   ','   nH+  ',
     &	'   grd_ion','   alpha ','   delta ','   c_v   ','  c_p   ',
     &  '   c_s  ','  nabl_ad',' nabla'
	do i=1,n
	   call thermosub(t(i),pe(i),thermo)
	   write(2,1001)tau(i),thermo(4:14),nabla(i)
	end do
        close(2)
	
1000     format(1x,f7.4,1x,f8.1,1x,1pe12.5,1x,e10.3,1x,7(e11.4,1x))
1001     format(1x,f6.3,1x,12(1pe10.3,1x))
1002     format(1x,a5,1x,12(a10,1x))			
	
        go to 12
11      print*,'error opening input file'
        go to 12
111     print*,'error opening input file (longversion)'

12      close(ican)
	end
