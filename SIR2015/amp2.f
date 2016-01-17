c "amp2" construye la atmosfera completa a partir de la antigua y
c de las perturbaciones multiplicativas nuevas (atmosfera reducida) 
c 'atmos' es la atmosfera antigua completa a la salida es la atm. nueva
c 'pert'  es la atmosfera antigua reducida
c 'atmosr' es la atmosfera perturbada reducida
c a la salida atmos contendra la nueva atmosfera perturbada en todos los
c puntos, como una interpolacion por splines cubicos de la atmosfera
c perturbada en los nodos (reducida)
c
c    i variable   i  variable
c    -   tau1     -    tau2   profundidad optica a 5000 /AA
c    1   t1       9    t2     temperatura en ambos modelos (k)
c    2   p1      10    p2     presion electronica (dinas/cm**2)
c    3   mic1    11    mic2   microturbulencia (cm/s)
c    4   h1      12    h2     campo magnetico (G)
c    5   v1      13    v2     velocidad eje z (cm/s)
c    6   g1      14    g2     gamma (radianes)
c    7   f1      15    f2     fi (radianes)
c    8   mac1    16    mac2   macroturbulencia (en Km/s)
c    -   ff1     17    ff2    factor de llenado (ff1=1-ff2)
c                18    %      peso de la luz difusa
c trabajo siempre con el ff del segundo modelo
c m(i) es el numero de nodos de la varible i. 
c Asi si m(i)=0 la variable i no se modifica
c     si m(i)=1 la variable i se modifica mediante un factor mult. cte.
c     si m(i)=2 la variable i se modifica mediante un factor mult. lineal
c     .........
c     si m(i)=-1 la variable i se modifica igual que la misma variable 
c               de la otra atmosfera (es decir como i+/-8)
c En mdata(1-18) guardo los indices anteriores a la variable i (atm. ampliada)
c
c Basilio 22-3-93 
c Basilio 23-3-93 (modificacion para contemplar el caso m(i)=-1)
c Basilio y Jose Carlos 9-1-95 (modificacion perturbaciones aditivas en fi)
c Basilio y Jose Carlos 6-2-95 (modificacion perturbaciones aditivas en gamma)
c
c _______________________________________________________________________

	subroutine amp2(ntau,m,atmos,atmosr)

	implicit real*4 (a-h,o-z)
	include 'PARAMETER'

	integer m(*),mdata(18)
	real*4 atmos(*),atmosr(*)
	real*4 x(kt),y(kt),yy(kt),pert(14*kt+4),f(kt,kt)
	real*4 tau(kt),t1(kt),p1(kt),tnew1(kt),pnew1(kt)
	real*4 t2(kt),p2(kt),tnew2(kt),pnew2(kt)
        real*4 pg1(kt),z1(kt),ro1(kt),b1(kt),gam1(kt)
        real*4 pg2(kt),z2(kt),ro2(kt),b2(kt),gam2(kt)
        character*29 var(18) 
        character*100 control 
        integer icalerr
        common/preciso/prec      
        common/calerr/icalerr !si calculo errores=1 else =0
	common/canal/icanal
	common/nombrecontrol/control
        common/contornopg/ncontpg,pg01,pg02
        common/contornoro/ro01,ro02
        common/zetas/pg1,z1,ro1,pg2,z2,ro2
        common/pgmag/ipgmag

c	epsilon=1.d-2
	epsilon=1.e-8
	toffset=2500.
	pi=3.14159265

c precision equilibrio hidrostatico en tanto por uno (necesaria para equisubmu)
        prec=1.e-3  ! 
     
	call comprime2(ntau,m,atmos,pert) !atmosfera antigua reducida
	
        do i=1,16
	   mdata(i)=i*ntau+2*int(i/9)  ! indi. anteri. a la var. i (ampliada)
	end do
        mdata(17)=mdata(16)+1
        mdata(18)=mdata(17)+1
   
	do i=1,ntau
	   tau(i)=atmos(i)
           t1(i)=atmos(ntau+i)
	   p1(i)=atmos(2*ntau+i)	!inicializamos la presion
           t2(i)=atmos(9*ntau+2+i)
	   p2(i)=atmos(10*ntau+2+i)	
	end do

	kred=0		!indice reducido
        kamp=ntau	!indice ampliado (los ntau puntos de tau1)
        cota=.5
	cotapres=.1
        cotafi=.395      !pi/8  

	do i=1,18	!do en grupos de varibles (1=t,2=p,...etc)
           ntau2=ntau
           if(i.eq.8.or.i.eq.16.or.i.eq.17.or.i.eq.18)ntau2=1!mac1,mac2,ff2,%

           if(m(i).eq.1)then	            !si pert. constante sumo
              kred=kred+1
              if(i.eq.6.or.i.eq.14)then
                 y1=atmosr(kred)-pert(kred)
		 if(y1.lt.0. .and. atmosr(kred). lt.0)then
		    y1=-pert(kred)/2. 
		 endif
		 if(atmosr(kred). gt. pi)then
		   y1=(pi-pert(kred))/2. 
		 endif
		 if(y1.lt.-cotafi)y1=-cotafi   !acoto inferiormente
	         if(y1.gt.cotafi)y1=cotafi     !acoto superiormente
	      else if(i.eq.7.or.i.eq.15.)then
                 y1=atmosr(kred)-pert(kred)
  	         if(y1.lt.-cotafi)y1=-cotafi   !acoto inferiormente
	         if(y1.gt.cotafi)y1=cotafi     !acoto superiorment
	      else if (i.eq.2.or.i.eq.10) then
	         if(pert(kred).eq.0)goto 999
	         y1=(atmosr(kred)/pert(kred))-1.    !perturbacion multiplicativa 
c  	         if(y1.lt.-cota)y1=-cota   !acoto inferiormente
c	         if(y1.gt.cota)y1=cota     !acoto superiormente
c	         if(i .eq.2.)print*,'amp2 125 multiplicativa',y1
	         y1=y1*pert(kred) 	            !perturbaciona aditiva	         
c	         if(i .eq.2.)print*,'amp2 127 indoce - reducida- aditiva',kred,pert(kred),y1
      	      else
                 if(pert(kred).eq.0)goto 999
	         y1=(atmosr(kred)/pert(kred))-1.    !perturbacion multiplicativa 
  	         if(y1.lt.-cota)y1=-cota   !acoto inferiormente
	         if(y1.gt.cota)y1=cota     !acoto superiormente
	         y1=y1*pert(kred) 	            !perturbaciona aditiva
              end if

              if(i.eq.17)then        !si es el f.f
                 varfill=(1./atmos(kamp+1))-1.
                 varfill=varfill+y1  
	         atmos(kamp+1)=1./(1.+ varfill)
              else if(i.eq.18)then        !si es el %
                 varpercen=atmos(kamp+1)/(100.-atmos(kamp+1))
                 varpercen=varpercen+y1
	         atmos(kamp+1)=100.*varpercen/(1.+ varpercen)
	      else if(i.eq.2 .or. i.eq.10)then  
	         pmin=1.e15  
	         do j=1,ntau2
	            if(atmos(kamp+j) .lt. pmin)pmin=atmos(kamp+j)
	         end do
	         if(y1+pmin .lt. 0.05*pmin)y1=-0.95*pmin
	         do j=1,ntau2
	            presionelec=atmos(kamp+j)+y1
                    if(presionelec .lt. 1.e-6)presionelec=1.e-6
                    atmos(kamp+j)=presionelec
                 end do	         
	      else   !if(i.ne.5.and.i.ne.13)then    !para una cte en velocidad
                 do j=1,ntau2
                    atmos(kamp+j)=atmos(kamp+j)+y1
                 end do
              end if 

           else if(m(i).gt.1)then

              mm=(ntau-1)/(m(i)-1)   !espaciado entre nodos 
              kred1=kred+m(i)           !indice del ultimo nodo de cada variable
              do j=1,m(i)
                 kred=kred+1                
                 jj=(j-1)*mm+1	               !indice de tau en los nodos
              	 x(j)=atmos(jj)	               !tau en los nodos 
                 if(i.eq.7.or.i.eq.15)then
                    y(j)=atmosr(kred)-pert(kred)
  	            if(y(j).lt.-cotafi)y(j)=-cotafi   !acoto inferiormente
	            if(y(j).gt.cotafi)y(j)=cotafi     !acoto superiormente
		 else if(i.eq.6.or.i.eq.14)then
                    y1=atmosr(kred)-pert(kred)
		    if(y1.lt.0. .and. atmosr(kred). lt.0)then
		       y1=-pert(kred)/2. 
		    endif
		    if(atmosr(kred). gt. pi)then
		       y1=(pi-pert(kred))/2.
		    endif	    
  	            if(y1.lt.-cotafi)y1=-cotafi   !acoto inferiormente
	            if(y1.gt.cotafi)y1=cotafi     !acoto superiormente
		    y(j)=y1	
		 else if (i.eq.2.or.i.eq.10) then
		    y1=(atmosr(kred)/pert(kred1))-1.  !pert. multiplicativa
c                    y1=atmosr(kred)-pert(kred)        !perturbaciona aditiva en la presion e
 	            if(y1.lt.-cotapres)y1=-cotapres   !acoto inferiormente
	            if(y1.gt.cotapres)y1=cotapres     !acoto superiormente
                    y(j)=y1*pert(kred1)               !escala con la presion en el ultimo no
c                    print*,'amp2 171 =',j,y(j) 
c                    y(j)=y1*pert(kred)     !perturb. aditiva en los nodos
c                    print*,'amp2 178 ',pert(kred),y1,y(j)
      	         else
	            y1=(atmosr(kred)/pert(kred))-1.  !pert. multiplicativa
	            if(y1.lt.-5.d-1)y1=-5.d-1 !acoto inferiormente
	            if(y1.gt.1.5d+0)y1=1.5d0  !acoto superiormente
	            y(j)=y1*pert(kred)     !perturb. aditiva en los nodos
                 end if
	      end do
    
	      if(ntau2.ne.ntau)then
		  print*,'STOP: Are there more than one parameter for the macroturbulence or filling factor?'
	          stop
              end if

	      call splines22(x,y,m(i)-2,ntau,tau,yy,f)	  
 
              if(i.eq.2.or.i.eq.10)then
                 do j=1,ntau2
                    presionelec=atmos(kamp+j)+yy(j)
c                    print*,'amp2 203 atmos(kamp+j) yy presionelec',j,atmos(kamp+j),yy(j),presionelec
c                   if(presionelec .lt. atmos(kamp+j)*.25)print*,'amp2 204 atmos(kamp+j) yy presionelec',atmos(kamp+j),yy(j),presionelec
                    if(presionelec .lt. atmos(kamp+j)*.1)presionelec=.1*atmos(kamp+j)
                    if(presionelec .lt. 1.e-6)presionelec=1.e-6
c                    print*,'amp2 207 ',i,kamp+j,3*ntau,presionelec
                    atmos(kamp+j)=presionelec
                 end do                  
	      else
                 do j=1,ntau2
                    atmos(kamp+j)=atmos(kamp+j)+yy(j)
                 end do
              end if               
	   end if	     
	   kamp=kamp+ntau2
	   if(i.eq.8)kamp=kamp+ntau+1	!los ntau puntos de tau2 y el de ff1	
	end do
c	print*,'pe amp2 ',3*ntau,atmos(3*ntau)

c en caso de que no se corrija la presion
c ponemos las presiones en equilibrio hidrostatico con las temperaturas
	vart1=0.0
	vart2=0.0
        do i=1,ntau
           tnew1(i)=atmos(ntau+i)
           tnew2(i)=atmos(9*ntau+2+i)
	   avar1=abs((tnew1(i)-t1(i))/t1(i))
	   avar2=abs((tnew2(i)-t2(i))/t2(i))
	   vart1=amax1(vart1,avar1)           !maxima variacion de t1
	   vart2=amax1(vart2,avar2)           !maxima variacion de t2
	end do

	varp1=0.0
	varp2=0.0
        do i=1,ntau
           pnew1(i)=atmos(2*ntau+i)
           pnew2(i)=atmos(10*ntau+2+i)
           b1(i)=atmos(4*ntau+i)
           b2(i)=atmos(12*ntau+2+i)
           gam1(i)=atmos(6*ntau+i)
           gam2(i)=atmos(14*ntau+2+i)
	   pavar1=abs((pnew1(i)-p1(i))/p1(i))
	   pavar2=abs((pnew2(i)-p2(i))/p2(i))
	   varp1=amax1(varp1,pavar1)           !maxima variacion de p1
	   varp2=amax1(varp2,pavar2)           !maxima variacion de p2
	end do

c suavizo la variacion de la temperatura si es muy fuerte (mayor del 40%)
c	if(vart1.gt.0.50)then                !si la max. var. de t1
c           do i=1,ntau
c              tnew1(i)=(tnew1(i)+2.*t1(i))/3.
c	   end do
c        end if	
c	if(vart2.gt.0.50)then                !si la max. var. de t2
c           do i=1,ntau
c              tnew2(i)=(tnew2(i)+2.*t2(i))/3.
c	   end do
c        end if

	do i=1,ntau  !no dejamos que la temperatura baje de toffset
	      if(tnew1(i).lt.toffset)tnew1(i)=toffset
	      if(tnew2(i).lt.toffset)tnew2(i)=toffset
	enddo
	
c no dejamos que b sea negativo ni gamma
        do i=1,ntau
	    if( m(4).gt.0 .and. atmos(i+4*ntau).le.0.)atmos(i+4*ntau)=1.
	    if(m(12).gt.0 .and. atmos(i+12*ntau+2).le.0.)atmos(i+12*ntau+2)=1.
	end do
 
c suavizo la variacion de la presion si es muy fuerte (mayor del 40%)
c	if(varp1.gt.0.40 .and. m(2).eq.0 )then                !si la max. var. de p1
c           do i=1,ntau
c              dif=pnew1(i)-p1(i)
c              esc=1.2
c              if(dif .lt. 0.)esc=.8              
c              atmos(2*ntau+i)=esc*p1(i)
c	   end do
c        end if	
c	if(varp2.gt.0.40 .and. m(10).eq.0)then                !si la max. var. de p2
c           do i=1,ntau
c              dif=pnew2(i)-p2(i)
c              esc=1.2
c              if(dif .lt. 0.)esc=.8              
c              atmos(10*ntau+2+i)=esc*p2(i)
c	   end do
c        end if
 
c coloco la p1 en equilibrio hidrostatico si t1 cambio mas de .002% (epsilon) 
c y no estoy evaluando errores
	   if(m(2).eq.0) then
	        if(ncontpg.eq.0)call equisubmu(ntau,tau,tnew1,p1,pg1,z1,ro1)	                                      
                if(ncontpg.eq.-1)then
                   ro1(ntau)=ro01
                   call pgpefromrho(tnew1(ntau),ro1(ntau),p1(ntau),pg1(ntau))
                   pg01=pg1(ntau)
                endif   
	        if(ncontpg.ne.0)then
                   if(ipgmag.ne.1)call 
     & equisubmu_cont(ntau,tau,tnew1,p1,pg01,pg1,z1,ro1)
                   if(ipgmag.eq.1)call 
     & equisubmu_contmag(ntau,tau,tnew1,p1,pg01,pg1,z1,ro1,b1,gam1)
                end if
	        do i=1,ntau
                   atmos(ntau+i)=tnew1(i)
                   atmos(2*ntau+i)=p1(i)
	        end do
	     else
	        do i=1,ntau
                   atmos(ntau+i)=tnew1(i) 
                   p1(i)= atmos(2*ntau+i)                  
                end do
                call pgzrofrompetau(ntau,tau,tnew1,p1,pg1,z1,ro1)
	     end if

        if(m(9).ne.0.or.m(10).ne.0.or.m(12).ne.0.or.m(14).ne.0)then
	     if(m(10).eq.0) then
	        if(ncontpg.eq.0)call equisubmu(ntau,tau,tnew2,p2,pg2,z2,ro2)
                if(ncontpg.eq.-1)then
                   ro2(ntau)=ro02
                   call pgpefromrho(tnew2(ntau),ro2(ntau),p2(ntau),pg2(ntau))
                   pg02=pg2(ntau)
                endif   
	        if(ncontpg.ne.0)then
                   if(ipgmag.ne.1)call 
     & equisubmu_cont(ntau,tau,tnew2,p2,pg02,pg2,z2,ro2)
                   if(ipgmag.eq.1)call 
     & equisubmu_contmag(ntau,tau,tnew2,p2,pg02,pg2,z2,ro2,b2,gam2)
                end if
                do i=1,ntau
                   atmos(9*ntau+2+i)=tnew2(i)
                   atmos(10*ntau+2+i)=p2(i)
	        end do
	     else
                do i=1,ntau
                   atmos(9*ntau+2+i)=tnew2(i)
                   p2(i)= atmos(10*ntau+2+i)
	        end do
	        call pgzrofrompetau(ntau,tau,tnew2,p2,pg2,z2,ro2)
	     end if
	end if

c Si m(i)=-1 se toma la variable de la otra atmosfera
	do i=1,16   !el ff2 no puede tener -1 (en ese caso se toma como 0) 
           ntau2=ntau
           if(i.eq.8.or.i.eq.16)ntau2=1 !si mac1,mac2 o ff2
           if(m(i).eq.-1)then
              ii=i-8
              if(ii.le.0)ii=i+8
              do j=1,ntau2 
		 atmos(mdata(i)+j)=atmos(mdata(ii)+j)
              end do
           end if
        end do 

c en cualquier caso los ff tienen que ser complementarios
        atmos(8*ntau+2)=1.-atmos(16*ntau+4)
 
	return

999	var(1)='temperature (model 1)      '
        var(2)='electronic pressure (model 1)'
        var(3)='microturbulence (model 1)    '
        var(4)='magnetic field (model 1)     '
        var(5)='LOS velocity (model 1)       '
        var(6)='incl. magn. field  (model 1) '
        var(7)='azimuth magn. field (model 1)'
        var(8)='macroturbulence (model 1)    '
	var(9)='temperature (model 2)        '
        var(10)='electronic pressure (model 2)'
        var(11)='microturbulence (model 2)    '
        var(12)='magnetic field (model 2)     '
        var(13)='LOS velocity (model 2)       '
        var(14)='incl. magn. field  (model 2) '
        var(15)='azimuth magn. field (model 2)'
        var(16)='macroturbulence (model 2)    '
        var(17)='filling factor (model 2)     '
        var(18)='stray light factor           '
       
        

	print*,' '
	write(*,'(a10,a30,a18)') 'STOP: The ',var(i),'is zero somewhere.' 
        print*,'     Since RELATIVE perturbations are used in the inversion,'
        print*,'     non zero initial values for that parameter have to be provided.'
c	open(icanal,file=control,fileopt='eof')
c	write(icanal,*) ' '  
c	write(icanal,'(a10,a30,a18)') 'STOP: The ', var(i),'is zero somewhere.' 
c        write(icanal,*) '     Since RELATIVE perturbations are used in the inversion,'
c        write(icanal,*) '     non zero initial values for that parameter have to be provided.'
c	close(icanal)
        print*,' '
	print*,'__________________________________________________________________________________'
	stop

	end

