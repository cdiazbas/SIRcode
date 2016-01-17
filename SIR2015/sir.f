c 			SIR
c
c Stokes Inversion based on Response functions
c
c
c Modificacion inversion 2 componentes: Basilio 22-3-93 
c Modificacion calculo errores        : Basilio 23-5-94
c Modificacion malla variable         : Basilio 26-5-94
c Modificacion luz difusa             : Basilio 24-6-94
c Modificacion calculo errores        : Basilio y Jose Carlos 20-9-94
c Modificacion calculo presiones      : Basilio y Jose Carlos 6-2-95
c Modificacion pesos (escala automatica a los picos de I,Q,U,V)
c y modificacion del numero de iteraciones: Basilio  15-2-96
c Reforma optimizacion y user-friendly: Luis Bellot 1999
c Modificacion opcion splines/lineal  : Basilio 12-7-00 
c Busqueda autoamtica de nodos        : Basilio 28-11-00
c_______________________________________________________________________
c
c       DIMENSIONES:
	implicit real*4 (a-h,o-z)

	include 'PARAMETER'

	parameter (kt16=16*kt+5)         
	parameter (kl4=4*kl,kld4=4*kld)
	parameter (kldt=kld*kn,kldt4=4*kldt)
        parameter (mmax=16*kn+2,kldn=mmax*kld4)
c        parameter (mfitmax=200)  !incluido en PARAMETER
        parameter (kn16=16*kn+4)

c 	PARA LA MALLA
	integer nlin(kl),nlins(kl4),npas(kl),npass(kl4),npos(kld4),nposi(kld),indice(kl)
        integer nlincheck(kl),ihemi,imaya,ici,iauto
	integer ist(4),nble(kl),nlinsn(kl)
	integer il
	real*4 dlamda0(kl),dlamda(kld),pist(4)
        real*4 dlamdas(kld4),cth,vx
	integer m(18)
	integer icanal,ican,iciclo

c 	PARA LA ATMOSFERA
        integer ntau,ierror
	real*4 atmos(kt16),atmosr(kt16),atmosrlin(kt16)
	real*4 atmosout(kt16),atmoslin(kt16),tau(kt)
        real*4 pg1(kt),z1(kt),ro1(kt),gam1(kt)
        real*4 pg2(kt),z2(kt),ro2(kt),gam2(kt)
        real*4 pg1b(kt),z1b(kt),ro1b(kt)
        real*4 pg2b(kt),z2b(kt),ro2b(kt)

c 	PARA EL PERFIL
	real*4 stok(kld4),sig(kld4),stray(kld4),sigd(kld4) !sigd es copia de sig
	real*4 scal(kld4),smax(kl,4),vic(kl),perfil(kld4)
        integer npasobs(kl),numberfrec(4)
        real*4 dlamdaobs(kld),dlamdacheck(kld)

c 	PARA LA INVERSION
	integer mreadi2,mreadi3,meves,iratio(4)
	real*4 covar(mfitmax,mfitmax),alpha(mfitmax,mfitmax),beta(mfitmax)
        real*4 mreadr2,mreadr3

c       PARA LOS COEFICIENTES DE ALEJAMIENTO
        real*4 beta1fich(kl,300),beta2fich(kl,300),taufich(300),linfich(kl)
        real*4 beta1(kl,kt),beta2(kl,kt)
	real*4 xa(11),ya(11)
	integer*4 indexpar(kl)

c 	PARA LOS ERRORES
        real*4 dvcal(kldn),x(mfitmax)   !,taur(kt)
        real*4 errores(mfitmax)
        integer ifiable(kn16),indi(kn16)
	real*4 maximo

c	PARA EL FICHERO LOG
	integer	nguvec(500),posicionciclo(20),mvec(500,18),mvecmax(18)
	real*4  addvec(500),snchivec(500),chprintvec(500),amac1vec(500),amac2vec(500)
	real*4  amic1vec(500),amic2vec(500),fill2vec(500),contrastevec(500)
	real*4  porcienvector(500)
        real*4  chisnpr(20),snpr(20)

c 	CADENAS
        character*100 uveobs,uveout,modelout1,modelout2,modelin1,modelin2,mod1inic,mod2inic
        character*100 malla,control,controlb,fcontrol,modelerr1,modelerr2,difusa
        character*100 snychi
        character*100 filtro
	character cha1*2,cha2*2,chag*1,chamod*4,chaper*4,chaerr*4,chalog*4,ca(18)*3
	character chasnychi*4
        character*100 mreadc2
	character*100 nomlineas,fichabun,departfile
	character*80 men1,men2,men3
	character*120 linlarga
	character*80 cartel

c 	COMUNES
	common/responde2/ist,ntau,ntl,nlin,npas,nble
	common/ldeo/dlamda,dlamda0
        common/smalla/ntls,nlins,npass
        common/smalla1/dlamdas
	common/atmosfera/atmos
	common/atmosferaout/atmosout 
	common/tol/tol
	common/filfac/fill
	common/alamda0/alamda0
	common/uvesalida/scal
	common/observaciones/stok,sigd
c	common/mu/cth				     !para equisub en amp22			    
	common/sigrealchi/sig,chireal,sumsq
	common/canal/icanal
	common/nombrecontrol/control
	common/cambio/ncambiono,ncambioprec          !para splines
	common/ifiltro/ifiltro
	common/filtro/filtro
	common/nlte/nlte
	common/contraste/contr
        common/repeticion/factorrep !factor de diag para evitar repeticiones
        common/nciclos/nciclos   !para marquardt2, leeuve3, leemodi22
        common/derivchi/dvcal  ! deriv. chi^2 para los errores
        common/posiciones/npos 
        common/errores/errores
        common/calerr/icalerr     !si calculo errores=1 else =0 (para amp2)
        common/difusa/stray
        common/fiable/ifiable,indi
c       El common siguiente es para pasarle a marquardt2 los indices iniciales y 
c       finales de las perturbaciones a gamma y fi aditivas
c        common/ifies/iga1,ifi11,iga2,ifi22 
        common/ifies/ipa1,ipa11,iga1,ifi11,ipa2,ipa22,iga2,ifi22
        common/ieliminofrec/ielimino !para el print de la S/N en marquardt2
        common/ivez/ivez
        common/ficlineas/nomlineas   !para leelineasii
        common/fichabun/fichabun     !para atmdatb
        common/numberfrec/numberfrec,iratio     !para marqcoef2
	common/nommod/modelin1,modelin2       !para escribeFR
        common/atmosr/atmosr                  !para escribeFR
        common/tau/tau
        common/iautomatico/iauto ! se pasa a fperfil2
        common/iprimeravez/iprimeravez !se cambia en fperfil2
        common/nspl_lin/nspl_lin !para seleccionar interp.lineal=1 o splines=0
	                         !o penalty (sol. regularizada =3 o 2, 5 o 4)
        common/contornopg/ncontpg,pg01,pg02
        common/contornoro/ro01,ro02
        common/primerchi/chisn,snn
	common/departcoef/beta1,beta2
        common/zetas/pg1,z1,ro1,pg2,z2,ro2    !para el calculo de z,pg y ro en amp22 y amp22err
        common/pgmag/ipgmag
	common/mu/cthabs	
        common/anguloheliocent/xmu    !for equisubmu and related routines
                                     !to evaluate equ. hidr. in the z direction,
	                             !not along the LOS

c       EXTERNAS
	external mreadi2,mreadr2,mreadc2,meves,mreadi3,mreadr3

c ntau : numero de puntos en tau
c tau  : log10(tau)
c ntl  : numero total de lineas
c nlin : indice de cada linea
c npas : numero de puntos en cada linea
c dlamda:cada delta de l.d.o. en ma

c_______________________________________________________________________
        cotaminima=-0.99  !valor minimo de cualquier p. Stokes
                          !cualquier Stokes menor de cotaminima tendra una sigma=1.e15

        vx=0.d0         !  vx=1.95d5 velocidad rotacion fotosfera solar ecuador

        icanal=24	!canal del fichero log
	ican=23	        !canal de lectura del fichero de control

        print*,' '
	print*,'__________________________________________________________________________________'
	print*,' '
	CARTEL=" SIR VERSION (10/November/2015) "
	print*,CARTEL
	print*,'__________________________________________________________________________________'
	print*,' '
	print*,' '
c el cambio del 20 April 05 (respecto a 16 marzo 2005) consisitio en el cambio del calculo
c de la polaridad (ahora a partir del areaVroja y areaVazul de todas las lineas en lugar del
c maximo de V de la primera)lineas 808-809 842-843 y 878

	print*,'Control file: '
	read(*,'(a)')fcontrol

c	Se define el fichero log:
	chalog='.log'
        chasnychi='.chi'
	controlb=fcontrol

	call quitaex(controlb)
	call concatena(controlb,chalog,control)
        call concatena(controlb,chasnychi,snychi)

c	open(icanal,file=control,fileopt='eof')    
	linlarga='--------------------------------------------------------------------------------'

c 	Se leen las condiciones de ejecucion

	do i=1,18
           mvecmax(i)=0  !para los errores
        end do

        nciclos=1
	ici=0
	ll=0
	do while (ici.lt.nciclos)
        ivez=0
        iprimeravez=0
        ici=ici+1 
	ncambiono=0	!cambio el numero de nodos
	ncambioprec=0	!cambio la precision
	open(ican,file=fcontrol)
	nciclos =mreadi2(ican,ici)        !numero de ciclos (repeticiones)
	uveobs  =mreadc2(ican,ici)        !fichero de entrada perfiles
	difusa  =mreadc2(ican,ici)        !fichero entrada  luz difusa
	filtro  =mreadc2(ican,ici)        !fichero con transformada de PSF
        malla   =mreadc2(ican,ici)        !malla,(npas=n. puntos por linea)
	nomlineas=mreadc2(ican,ici)       !fichero con parametros atomicos
	fichabun =mreadc2(ican,ici)       !fichero con abundancias
        modelin1 =mreadc2(ican,ici)       !nombre del modelo 1 anterior 
        modelin2 =mreadc2(ican,ici)       !nombre del modelo 2 anterior
        pist(1)  =mreadr3(ican,ici,1.)        !i (peso de i en ajuste)
        pist(2)  =mreadr3(ican,ici,0.)        !q
        pist(3)  =mreadr3(ican,ici,0.)        !u
        pist(4)  =mreadr3(ican,ici,0.)        !v
        iauto   =mreadi3(ican,ici,0)        !seleccion automatica de nodos
	m(1)    =mreadi3(ican,ici,0)        !nodos en t1
	m(2)    =mreadi3(ican,ici,0)        !nodos en pe1
	m(3)    =mreadi3(ican,ici,0)        !nodos en micro1
	m(4)    =mreadi3(ican,ici,0)        !   "   "    h1
	m(5)    =mreadi3(ican,ici,0)        !   "   "    vz1
	m(6)    =mreadi3(ican,ici,0)        !   "   "    gamma1
	m(7)    =mreadi3(ican,ici,0)        !   "   "    fi1
	m(8)    =mreadi3(ican,ici,0)        !   "   "    macro1
	m(9)    =mreadi3(ican,ici,0)        !   "   "    t2
	m(10)    =mreadi3(ican,ici,0)       !   "   "    pe2
	m(11)    =mreadi3(ican,ici,0)       !   "   "    micro2
	m(12)    =mreadi3(ican,ici,0)       !   "   "    h2
	m(13)    =mreadi3(ican,ici,0)       !   "   "    vz2
	m(14)    =mreadi3(ican,ici,0)       !   "   "    gamma2
	m(15)    =mreadi3(ican,ici,0)       !   "   "    fi2
	m(16)    =mreadi3(ican,ici,0)       !   "   "    macro2
	m(17)    =mreadi3(ican,ici,0)       !   "   "    f.f.2
        m(18)    =mreadi3(ican,ici,0)       !luz difusa
	cth     =mreadr3(ican,ici,1.)        !coseno angulo helic. (- hem oeste)
	sn      =mreadr3(ican,ici,1000.)        !segnal ruido estimada para i
	contr    =mreadr3(ican,ici,-1.)       !contraste Ic1/Ic2
c        nlte    =mreadi2(ican,ici)	  !1=efectos de nlte , 0=lte
	tol     =mreadr3(ican,ici,1.e-4)        !tolerancia inversion 
	alamda0 =mreadr3(ican,ici,1.e-3)        !factor diagonal inicial 
        nspl_lin =mreadi3(ican,ici,0)	  !0=splines, 1=lineal (regularizada splines(2 o 4) o lineal(3 o 5)
        pg01    =mreadr3(ican,ici,0.)     !condi. cont. pg atm. 1/negative value= density
        pg02    =mreadr3(ican,ici,0.)     !condi. cont. pg atm. 2/negative value= density
        ipgmag  =mreadi3(ican,ici,0)      !terminos en P magnetica?
        departfile=mreadc2(ican,ici)	  !fichero con departure coef.

c ____________________________________________________________________
      	print*,'  '
        print*,linlarga
	
	xmu=abs(cth)
	if(cth .ne. 1.)then
	   print*,'Resulting models are stratified versus line of sight'
	   print*,'We will use heliocentric angle value only for hydrostatic equilibrium equation'
           cth=1.  !por defecto, no pasamos de la linea de vision a z
	endif
        cthabs=abs(cth)

        ncontpg=0
        ro01=0.
        ro02=0.
        if(pg01.gt.0..or.pg02.gt.0. .and. m(2)+m(10) .eq.0 )then
           print*,'Boundary condition in gas pressure'
	   if(pg02.eq.0.)pg02=pg01
	   if(pg01.eq.0.)pg01=pg02
           ncontpg=1
        end if
        if(pg01.lt.0..or.pg02.lt.0 .and. m(2)+m(10) .eq.0 )then
           print*,'Boundary condition in density'
	   if(pg02.eq.0.)pg02=pg01
	   if(pg01.eq.0.)pg01=pg02
           pg01=abs(pg01)
           pg02=abs(pg02)
           ro01=pg01
           ro02=pg02
           ncontpg=-1
        end if
	if(ncontpg.eq.0 .and. m(2)+m(10) .eq.0 )then
           print*,'Boundary condition in electronic pressure'
	   print*,'Better results are expected setting boundary condition in gas pressure/density'
	end if

      	nlte=meves(departfile,100)

c cadenas para el output
        ca(1)=' t1'      
        ca(2)=' p1'      
        ca(3)=' m1'      
        ca(4)=' h1'      
        ca(5)=' v1'      
        ca(6)=' g1'      
        ca(7)=' f1'      
        ca(8)=' M1'      
        ca(9)=' t2'      
        ca(10)=' p2'      
        ca(11)=' m2'      
        ca(12)=' h2'      
        ca(13)=' v2'      
        ca(14)=' g2'      
        ca(15)=' f2'      
        ca(16)=' M2'      
        ca(17)=' ff'      
        ca(18)=' st'     


	if(nciclos.eq.0)then   !para sintesis uso los 4 parametros (LRB)
          do i=1,4
            iratio(i)=0     !iratio(Q)=1 -> se invierte Q/I, etc
            ist(i)=1
            if(pist(i).lt.0)iratio(i)=1
            pist(i)=1
          enddo
          iratio(1)=0  !I no se puede sintetizar como I/I
	else
          do i=1,4
             ist(i)=0
             iratio(i)=0     !iratio(Q)=1 -> se invierte Q/I, etc
             if(pist(i).lt.0)iratio(i)=1
             pist(i)=abs(pist(i))
             if(pist(i).ne.0.)ist(i)=1
          end do
	  if(ist(1).eq.0)print*,'Stokes I is not considered in this inversion'
	  if(ist(2).eq.0)print*,'Stokes Q is not considered in this inversion'
	  if(ist(3).eq.0)print*,'Stokes U is not considered in this inversion'
	  if(ist(4).eq.0)print*,'Stokes V is not considered in this inversion'
	endif
        if(iratio(2).eq.1.or.iratio(3).eq.1.or.iratio(4).eq.1)iratio(1)=0

c definimos el fichero con los parametros atomicos
        ilineas=0
	ilineas=meves(nomlineas,100)
        if(ilineas.eq.0)then 
		men1='The default LINES file is being used for supplying atomic parameters.'
		print*,men1
                nomlineas='~/sir/default/LINES'
        endif

c definimos el fichero con las abundancias
        iabun=0
	iabun=meves(fichabun,100)
        if(iabun.eq.0)then 
		men1='The default ABUNDANCES file is being used for supplying abundances.'
		print*,men1
                fichabun='~/sir/default/ABUNDANCES'
        endif

c definimos los nombres de los ficheros de salida
 	
	chamod='.mod'
	chaerr='.err'
	chaper='.per'
        ncifno=1
        if(ici.ge.10)ncifno=2 
        chag='_'
	cha1(1:2)='  '
	cha2(1:2)='  '
	modelout1=modelin1
	modelout2=modelin2
  

	call quitaex2(modelout1,ixt1)
        call quitaex2(modelout2,ixt2)

        if(ici.eq.1.or.ici.eq.0)then
             mod1inic=modelout1
             mod2inic=modelout2
	     call chachi(1+ixt1,cha1,ncifno)
	     call chachi(1+ixt2,cha2,ncifno)
        else
             if(modelout1.eq.mod1inic)then 
	        call chachi(ici-1+ixt1,cha1,ncifno)
                call concatena(modelout1,chag,modelin1)
                call concatena(modelin1,cha1,modelin1)
                call concatena(modelin1,chamod,modelin1)
	        call chachi(ici+ixt1,cha1,ncifno)
             else
	        call chachi(1+ixt1,cha1,ncifno)
             endif

             if(modelout2.eq.mod2inic)then
	        call chachi(ici-1+ixt2,cha2,ncifno)
                call concatena(modelout2,chag,modelin2)
                call concatena(modelin2,cha2,modelin2)
                call concatena(modelin2,chamod,modelin2)
	        call chachi(ici+ixt2,cha2,ncifno)
             else
	        call chachi(1+ixt2,cha2,ncifno)
             endif
        endif

        call concatena(modelout1,chag,modelout1)
        call concatena(modelout1,cha1,modelout1)
        call concatena(modelout1,chaper,uveout)
        call concatena(modelout1,chaerr,modelerr1)
        call concatena(modelout1,chamod,modelout1) 
        call concatena(modelout2,chag,modelout2)
        call concatena(modelout2,cha2,modelout2)
        call concatena(modelout2,chaerr,modelerr2)
        call concatena(modelout2,chamod,modelout2)

	if(contr.lt.0)then
           contr=-1.e6
        else
           contr=(contr-1.)/(contr+1.)
        end if

	ifiltro=0
	ifiltro=meves(filtro,100)
        if(ifiltro.eq.0)then 
                men1='PSF NOT taken into account'
		print*,men1
        endif

        if(ifiltro.eq.1)then
		men1='The PSF will be read from: '//filtro(1:40)
		print*,men1
	endif

	istray=meves(difusa,100)
        if(istray.eq.0)then
		men1='NO stray light file has been specified.'
		print*,men1
        endif
 

        if(nspl_lin.eq. 0)print*,'Splines interpolation between nodes '
        if(nspl_lin.eq. 1)print*,'Linear interpolation between nodes '
        if(nspl_lin.eq. 2)print*,'Light regularized solution (Splines interpolation between nodes)'
	if(nspl_lin.eq. 3)print*,'Light regularized solution (Linear interpolation between nodes)'
	if(nspl_lin.eq. 4)print*,'Regularized solution (Splines interpolation between nodes)'
	if(nspl_lin.eq. 5)print*,'Regularized solution (Linear interpolation between nodes)'


	ierror=0

c       ------------------------------------------------------
c       leemos la malla, el perfil ,el modelo y la luz difusa:
c       ------------------------------------------------------

  	call leemodi22(0,modelin1,modelin2,atmos,ntau,ierror)
c  	call leemodi222(0,modelin1,modelin2,atmos,ntau,ierror,
c     &                  z1,pg1,ro1,z2,pg2,ro2)

        gam1med=0.
	gam2med=0.
        do i=1,ntau
           tau(i)=atmos(i)	
	   gam1(i)=atmos(6*ntau+i)
           gam2(i)=atmos(14*ntau+2+i)
	   gam1med=gam1med+gam1(i)
	   gam2med=gam2med+gam2(i)
        end do
	gam1med=gam1med/ntau  !gamma1 media
	gam2med=gam2med/ntau  !gamma1 media

	imodel2=meves(modelin2,100)
        if(imodel2.eq.1.and.ierror.eq.2)then
          men1='     WARNING: Model 2 does NOT exist. Only model 1 is considered'
          print*,men1
	  imodel2=0  !imodel2 =1 significa que existe modelo 2
	endif


	if(ierror.eq.2)then    !no podemos invertir el modelo 2!!
            do jj=9,17
               m(jj)=0
            end do
        endif

        if(atmos(16*ntau+5).gt.1e-5.or.m(18).ne.0)then
          if(istray.ne.0)then
            call leeuveobsindic(difusa,ist,ntl,nliobs,nlin,npasobs,dlamdaobs,nble,stray)  !(LRB)
	      
	    ntlcheck=ntl  !para comprobar que difusa tiene las mismas lambdas que perf. obs.
	    do l=1,nliobs
	      dlamdacheck(l)=dlamdaobs(l)
c	      print*,'SIR ',l,dlamdacheck(l)
	    enddo
	    do l=1,ntl
	      nlincheck(l)=nlin(l)
	    enddo

c	    -----------------------------------------------------------------
c           Por si no tenemos el fichero con la malla en sintesis:
c	    (si lo tenemos, estas variables se machacan luego y no pasa nada)

            nli=nliobs
            do i=1,ntl
               npas(i)=npasobs(i)
            end do
            do i=1,nli
               nposi(i)=i
               dlamda(i)=dlamdaobs(i)
            end do
c	    -----------------------------------------------------------------

           else
	     men1='     WARNING: The STRAY LIGHT FACTOR is being changed to ZERO'
             if(atmos(16*ntau+5).gt.0.)print*,men1
c	     if(ici.eq.1.and.atmos(16*ntau+5).gt.1e-5)write(icanal,*) men1
             atmos(16*ntau+5)=0.
             m(18)=0  !LB (si no casca, porque trato de invertir una variable = cero)
           end if
        end if

c -----------------------------------------------------------------
c Para el calculo de las FR
c -----------------------------------------------------------------

	  if(nciclos.lt.0)then !para calculo de las FR mnodos=ntau
            men1=' STOP: The number of depth points in the model' 
            men2=' is larger than the current kn value' 
            men3=' Decrease the number of depth points or change the PARAMETER file.' 

            if(kn.lt.ntau)call mensaje(4,men1,men2,men3)

            do i=1,15
               if(m(i).ne.0)m(i)=ntau
            enddo
            if(m(8).ne.0)m(8)=1 !la macro solo puede tener 1 nodo
            iauto=0
          end if
c -----------------------------------------------------------------

        imaya=meves(malla,100)       !(LRB)

     	if(nciclos.lt.1)then
           uveout=uveobs
	   if(imaya.eq.1)then   !si tenemos malla
	         call leemallab(malla,ntl,nli,nlin,npas,dlamda,nble)
                 nliobs=nli
                 do i=1,ntl
                    npasobs(i)=npas(i)
                 end do
                 do i=1,nli
                   nposi(i)=i
                   dlamdaobs(i)=dlamda(i)
                 end do
           else
              if(istray.eq.1)then
		 men1='The wavelength grid is being generated automatically from the'
		 men2='stray light profile. NO blends are considered.'
                 print*,men1
                 print*,men2
              else
	         men1='STOP: In synthesis mode, you MUST specify the wavelength grid'
		 call mensaje(1,men1,men2,men3)
	      endif
	   endif
	     
           print*,'Output profiles: ',uveout
c           print*,linlarga
           if(ici.eq.1)then
	     open(icanal,file=control,access='append')
	        write(icanal,*)linlarga 
	     close(icanal)
	   endif
       
	else 	
           if(imaya.eq.0)then
	     men1='The wavelength grid is being generated automatically from the'
	     men2='observed profiles. NO blends are considered.'
             print*,men1
             print*,men2
	           
             call leeuveobsindic(uveobs,ist,ntl,nliobs,nlin,npasobs,dlamdaobs,nble,stok)

             nli=nliobs
             do i=1,ntl
                npas(i)=npasobs(i)
             end do
             do i=1,nli
                nposi(i)=i
                dlamda(i)=dlamdaobs(i)
             end do
           else
             call leeuveobsindic(uveobs,ist,ntl,nliobs,indice,npasobs,dlamdaobs,nble,stok)
	     call leemalla2(malla,ntl,nli,nliobs,nlin,npasobs,npas,
     &                             dlamdaobs,dlamda,nble,nposi,indice)
           endif

	endif

	if(atmos(16*ntau+5).gt.0.)then   !comprobaciones varias !esto es si tenemos difusa

	  if(ntlcheck.ne.ntl)then
	    men1='STOP: The number of lines in the files containing the observed and'
	    men2='      stray light profiles are not equal.'
	    call mensaje(2,men1,men2,men3)
          endif

	  do l=1,ntl
	    if (l.eq.1)then
		numm=1
	    else
		numm=numm+nble(l-1)  
	    endif
	    if(nlin(numm).ne.nlincheck(l))then
	      men1='STOP: The order of the lines in the stray light file is not'
	      men2='      equal to that in the file containing the observed profiles.'
	      call mensaje(2,men1,men2,men3)
	    endif
	  enddo
	    
          do l=1,nli
	    if(abs(dlamda(l)-dlamdacheck(l)).gt.1.)then
	      men1='STOP: The wavelengths in the file containing the stray light profile differ by'
	      men2='      more than 1. mA from those generated by the wavelength grid:'
c	      men3='      Wavelength #: '//char(l)   !,':',dlamda(l),dlamdacheck(l)
	      print*,'      Wavelength #: ',l,':',dlamda(l),dlamdacheck(l)
	      call mensaje(2,men1,men2,men3)
	    endif
	  enddo
	
	endif


	if(nciclos.ge.1)then
           if(ierror.eq.0)then
	      men1='Output models  : '//modelout1(1:20)//' and '//modelout2(1:20)
	      men2='Uncertainties  : '//modelerr1(1:20)//' and '//modelerr2(1:20)
	      men3='Output profiles: '//uveout(1:20)//' (2 components)'
	      print*,men1
	      print*,men2
	      print*,men3

c              if(ici.eq.1)then
c                   write(icanal,*) men1
c                   write(icanal,*) men2
c                   write(icanal,*) men3
c	      endif
		
           else if (ierror.eq.2)then
	      men1='Output model   : '//modelout1(1:50)
	      men2='Uncertainties  : '//modelerr1(1:50)
	      men3='Output profiles: '//uveout(1:43)//' (1 component)'

	      print*,men1
	      print*,men2
	      print*,men3

c             if(ici.eq.1)then
c                  write(icanal,*) men1
c                  write(icanal,*) men2
c                  write(icanal,*) men3
c	      endif
	   else 
              men1=' STOP: Models not specified or unexistent. '
	      call mensaje(1,men1,men2,men3)
	      stop
           end if
        end if

        print*,linlarga
	if(ici.eq.1)then
	   open(icanal,file=control,access='append')
	   write(icanal,*) linlarga 
	endif 

c       if(istray.ne.0)call convierte(stray,nposi,nliobs,nli)
        ixx=0
        do l=1,ntl
	   do ible=1,nble(l)
	   ixx=ixx+1
              do i=1,ntau
c                 beta1(nlin(ixx),i)=1.
c                 beta2(nlin(ixx),i)=1.
                 beta1(ixx,i)=1.  !entran ordenados sobre todas lineas*blends
                 beta2(ixx,i)=1.
	      end do	 
           end do
        end do


	if (nlte.eq.1)then
            men1=' Departure coefficients are being taken into account from'
            print*,' '
            print*,men1,departfile
            print*,' '
c el fichero debe tener 4 columnas (indice linea, logtau,depat_nivel_inferior, 
c depart_nivel_superior)
c y ntaufich*nlinineas filas (ordenado: primero los ntau valores
c de la prim linea, then de la seg ...etc) 
c incluyendo blends. La red en tau no tiene por que ser la misma 
c del modelo utilizado. No tienen porque estar todas las lineas presentes
c diemnsionado a kl lineas y 300 puntos en tau 
c maximo una coleecion de coeficientes de alejameinto en cada linea (solo un blend puede llevarlos)
            open(14,file=departfile,status='old')
            num_lineas=0              !numero de lineas espectrales del fichero departfile
            linea1old=0 
            do kdep=1,1000 
               read(14,*,end=1414)rlinea1
	       linea1=nint(rlinea1)
               if(linea1.ne.linea1old)num_lineas=num_lineas+1
               linea1old=linea1
            end do
1414	    kdep=kdep-1               !numero de filas del fichero departfile 
            ntaufich=kdep/num_lineas  !numero de puntos en tau del fichero departfile
            close(14)
            open(14,file=departfile,status='old')

            do l=1,num_lineas
               indexpar(l)=1
               do i=1,ntaufich
                  read(14,*)linfich(l),taufich(i),beta1fich(l,i),beta2fich(l,i) 
               end do
               do jj=1,ixx  !ixx lleva el numero total de lineas*blends
                 if(nlin(jj) .eq. nint(linfich(l)))indexpar(l)=jj
               end do  
            end do
            close(14)

	    ngrado=1  !grado del polinomio interpolador

	    n2=int(ngrado/2)
	    do i=1,ntau
               stau=atmos(i) 
	       call locate(taufich,ntaufich,stau,j)
	       n3=j-n2-1
               if(n3.lt.0)n3=0
               if(n3+ngrado+1.gt.ntaufich)n3=ntaufich-ngrado-1
	       do k=1,ngrado+1
	          xa(k)=taufich(n3+k)
	       end do
               do l=1,num_lineas
	          do k=1,ngrado+1
	             ya(k)=beta1fich(l,n3+k)
	          end do
	          call polint(xa,ya,ngrado+1,stau,sa,error)
c                  beta1(nint(linfich(l)),i)=sa
                  lx=indexpar(l) 
                  beta1(lx,i)=sa
                  
	          do k=1,ngrado+1
	             ya(k)=beta2fich(l,n3+k)
	          end do
	          call polint(xa,ya,ngrado+1,stau,sa,error)
c                  beta2(nint(linfich(l)),i)=sa
                  beta2(lx,i)=sa
c                  print*,stau,beta1(nint(linfich(l)),i)/sa

	       end do
            end do
	end if

c	---------------------
c       Sentencias de control
c	---------------------

c       Expreso el contraste como diferencia:
 
        if(ierror.ne.0.and.contr.gt.-1.e-6)then
           contr=-1.e-6
	   men1=' No contrast calculation can be done without two models!'
           print*,men1
c	   if(ici.eq.1)then
c               write(icanal,*) men1
c	       close(icanal)
c	   endif
        end if

	nfrecmax=0
        ntlblends=0
	do i=1,ntl
	   nfrecmax=max0(nfrecmax,npas(i))  !numero max. de frec.
           nlinsn(i)=nlin(ntlblends+1)
           do j=1,nble(i)
              ntlblends=ntlblends+1
           end do
	end do


	if (ntlblends.gt.kl) then
             print*,' '
             write(*,'(a34,i2,a38,i2)') 'STOP: The number of lines (',ntlblends,') is larger than the current limit kl=',kl
	     print*,'      Decrease this number or change the PARAMETER file.'
	     print*,'hola',ntl,nble
             print*,' '
             print*,'_______________________________________________________________________________'
             stop 
        end if
	if (nli.gt.kld)then
           print*,' ' 
           write(*,'(a33,i4,a40,i4)') 'STOP: The number of wavelengths (',nli,') is larger than the current limit kld= ',kld
	   print*,'      Decrease the number of wavelengths or change the PARAMETER file.'
           print*,' '
           print*,'_______________________________________________________________________________'
           stop 
        end if
	if (nfrecmax.gt.kld) then
           men1='STOP: There is a line with more wavelengths than the current limit kld'
	   men2='      Decrease the number of wavelengths or change the PARAMETER file.'
	   call mensaje(2,men1,men2,men3)
        end if
        if (ntau.gt.kt) then
	   men1='STOP: The model has more depth points than the current limit kt.'
	   men2='      Decrease the number of grid points or change the PARAMETER file.'
	   call mensaje(2,men1,men2,men3)
        end if
	if (nli*ntau.gt.kldt) then
           print*,' ' 
           print*, 'STOP: The total number of wavelengths times the number of depth points ',nli*ntau
           print*, '      is larger than the current limit kldt=kld*kn',kldt
	   print*, '      Decrease this product or change the PARAMETER file.'
           print*,' ' 
           print*,'_______________________________________________________________________________'
           stop
        end if

c       Definimos nlins,ntls,npass,dlamda0s,dlamdas:

	ntls=0
        nb=0
	k4=0
	k4c=0
	do i=1,ntl
	   vic(i)=1.     !inicializo i del continuo
           do ii=1,4
	     smax(i,ii)=1./sn
           end do
	end do

	ii=0
	vmax=-2.
c	lanvmax=0
	vmin=2.
c	lanvmin=0
        areaVazul=0.
	areaVroja=0.
	ipolaridad=1
	do i=1,4
	   do j=1,ist(i)
              iinli=ii*nli
              ii=ii+1
	      k3=0
	      k3c=0
	      suma=0
	      jj=0
	      do k=1,ntl
                 do jble=1,nble(k)   
                    nb=nb+1
                    jj=jj+1
                    nlins(nb)=nlin(jj) 
	         end do
                 ntls=ntls+1
		 npass(ntls)=npas(k)
                 do l=1,npasobs(k)
	            k3=k3+1
	            k4=k4+1
                    npos(k4)=nposi(k3)+iinli
                    sk4=stok(k4)
c		    if(i.eq.4.and.nciclos.gt.0 .and. k.eq.1 )then  !calculamos la polaridad de V (primera linea)
c		       if(sk4.gt.vmax)then
c		          vmax=sk4
c			  lanvmax=l
c		       else if(sk4 .lt. vmin)then
c		          vmin=sk4
c			  lanvmin=l
c		       end if
		    if(i.eq.4.and.nciclos.gt.0)then
c		       print*,k,l,dlamdas(l)
		       if(dlamdas(l).le.0.and.(sk4.gt.cotaminima))areaVazul=areaVazul+sk4
		       if(dlamdas(l).gt.0.and.(sk4.gt.cotaminima))areaVroja=areaVroja+sk4
		    end if
  	            if(sk4.gt.cotaminima)then
    	               if(i.eq.1.and.nciclos.gt.0)then
                           xmax=-1.
                           do lll=1,npasobs(k)
                             if(xmax.le.stok(k4-l+lll))xmax=stok(k4-l+lll)
                           end do
	                   vic(k)=xmax	!i continuo
                           sss=abs( vic(k)-sk4 )
	               else
	                   sss=abs(sk4)
	               end if
	               if(sss.gt.smax(k,i))smax(k,i)=sss
                     end if

	         end do !fin del do en frecuencias obs
                 do l=1,npas(k)
	            k3c=k3c+1
	            k4c=k4c+1
	            dlamdas(k4c)=dlamda(k3c)
	         end do !fin del do en frecuencias malla

	      end do 	!fin del do en lineas
	   end do       !fin del do en j
	end do          !fin del do en i
c	print*,'(SIR) blue V Area=',areaVazul
c	print*,'(SIR) red  V Area=',areaVroja

c comprobamos la polaridad y cambiamos la inicializacion de gamma si polaridad
c negativa y gamma < 90 o si polridad positiva y gamma > 90
c        if(nciclos .gt. 0)then
c	   if(lanvmax .gt. lanvmin)then
c	      print*,'(SIR) Stokes V of first line has negative polarity'
c	      print*,vmin,lanvmin,vmax,lanvmax
c           if( areaVroja .gt. areaVazul)then
c	      ipolaridad=-1
c	    end if  
c	   if (imodel2 .eq. 0)then !o sea si no existe el modelo 2
c	      if(gam1med. lt. 90. .and. ipolaridad .eq. -1)then
c	         do i=1,ntau
c                    atmos(6*ntau+i)=90.
c                 end do 
c	         print*,' Negative polarity. Initial magnetic field inclination (model 1) set to 90 degrees'
c	      else if(gam1med. gt. 90. .and. ipolaridad .eq. 1)then
c     	         do i=1,ntau
c                    atmos(6*ntau+i)=90.
c                 end do 
c	         print*,' Positive polarity. Initial magnetic field inclination (model 1) set to 90 degrees '        
c	      end if
c	   else !o sea si existe el modelo 2
c	      if(gam1med. lt. 90. .and. gam2med .lt. 90. .and. ipolaridad .eq. -1)then
c	         do i=1,ntau
c                    atmos(6*ntau+i)=90.
c                    atmos(14*ntau+2+i)=180.-gam2(i)
c                end do 
c	         print*,' Negative polarity. Changing initial magnetic field inclination by 90 (model1)
c     & and by suplementary angle (model2) '
c	      else if(gam1med. gt. 90. .and. gam2med .gt. 90. .and. ipolaridad .eq. 1)then
c	         do i=1,ntau
c                    atmos(6*ntau+i)=90.
c                    atmos(14*ntau+2+i)=180.-gam2(i)
c                end do 
c	         print*,'Positive polarity. Changing initial magnetic field inclination by 90 (model1) 
c    & and by suplementary angle (model2) '
c	      end if
c	   end if
c	end if

	print*,'  '
	
	maximo=smax(1,1)
	do k=1,ntl
            if(smax(k,1).gt.maximo)maximo=smax(k,1)
        enddo

	do k=1,ntl
           if(smax(k,1).lt..01.and.nciclos.gt.0)then
              print*,'Weight is being changed for the line',k
              smax(k,1)=maximo
           endif
        enddo     

	nfrecs=k4

c Para controlar el maximo numero de nodos en modo automatico ("*" en m(i))
        if(m(8).gt.1)m(8)=1
        if(m(16).gt.1)m(16)=1
        if(m(17).gt.1)m(17)=1

        if(iauto.eq.1)then
           n_1000=0
           n_det_min=0
           do i=1,18
              if(m(i).eq.1000)then
                 n_1000=n_1000+1
              else 
                if(m(i).gt.0)n_det_min=n_det_min+m(i)
              end if
           end do
           ndd1=min(mfitmax,nfrecs)
           nmax1=nint((ndd1-n_det_min)/float((n_1000+2)))
           
           if(m(1).eq.1000)m(1)=2*nmax1 !temperatura (doble num de nodos)
           do i=2,7
              if(m(i).eq.1000)m(i)=nmax1
           end do
           if(m(9).eq.1000)m(9)=2*nmax1 !temperatura (doble num de nodos)
           do i=10,15
              if(m(i).eq.1000)m(i)=nmax1
           end do
           kv=0
           if(n_1000.gt.0.or.iauto.eq.1)then
             do i=1,18
                if(m(i).ne.0)then
                   print*,'maximun number of nodes for variable ',i,'= ',m(i)
                end if
             end do
           end if
	else
          do i=1,18
             if(m(i).eq.1000)then
               print*,'* is not allowed in non automatic mode'
               stop
             end if
          end do
   	end if
c ____________________________________________________________________

c	------------------------------------------------------------
c OLD VERSION
c	Seleccion automatica del numero de nodos para cada variable
c	o comprobacion/modificacion del numero de nodos:
c	------------------------------------------------------------
    
c	if(nciclos.gt.0)call nodos2aut(ntau,iauto,ierror,ici,m)  !en la version LBR estaba comentada
 	if(nciclos.gt.0)call nodos2aut(ntau,0,ierror,ici,m)      !en la version LBR estaba sin comentar

c	----------------------------------------------
c       Calculo de la atmosfera en la linea de vision:
c	----------------------------------------------
	ihemi=1
	if(cth.lt.0.)then
	   ihemi=-1
	   cth=-cth
	end if

        il=0
        call taulinea2(il,cth,ihemi,vx,atmos,ntau)


c 	if(nciclos.gt.0)call nodos2(ntau,m) !comprueba que el num. de nodos es valido
                                            ! y lo modifica caso de no serlo
c	???ES NECESARIO PASAR DE NUEVO POR NODOS????????? (YO creo que no)
	

	mfit=0
	do i=1,18
           if(m(i).gt.0)mfit=mfit+m(i)
	end do

	if(mfit.gt.mfitmax)then 
           print*,' ' 
           print*, 'STOP: The number of free parameters ',mfit
           print*, '      is larger than the current limit mfitmax',mfitmax
	   print*,'       Decrease the number of free parameters.'
           print*,' ' 
           stop
	endif

 
c       Escribe en atmosr los valores de atmos en los nodos:
c	if(nciclos.ne.0)call comprime2(ntau,m,atmos,atmosr)

	nfree=nfrecs-mfit

	print*,'Number-of-observ-freqs/fitable-params/free-params: ',nfrecs,mfit,nfree	
        print*,'________________________________________________________________________________'
	print*,' '
	
        if(mfit.gt.nfrecs.and.nciclos.ne.0)then
           print*,' ' 
           write(*,'(a37,i2,a27)') 'STOP: The number of free parameters (',mfit,') is larger than the number'
           write(*,'(a22,i3,a2)') '      of observables (',nfrecs,').'
           print*,' '
           print*,'_____________________________________________________________________________'
           stop
        end if

	k40=0
	ntotal4=0
        ielimino=0
	sigmamax=1.e15

	do i=1,4
	   numberfrec(i)=0

	   do j=1,ist(i)
              indicei=0
              sigc=1./(sn*sqrt(pist(i)))
	      do k=1,ntl
                 numberfrec(i)=numberfrec(i)+npasobs(k) !numero de ldo para
                                             !cada stokes invertido 
                                             !es cero si no se invierte 

 	         sigc2=sigc*sqrt(smax(k,i)/vic(k))
	         do l=1,npasobs(k)
	            k40=k40+1
                    indicei=indicei+1
		    sig(k40)=sigc2

                    if(stok(k40).lt.cotaminima)then
                       sigx=sigmamax
                       sig(k40)=sigx
                       ielimino=ielimino+1
                    end if
   
                    sigd(k40)=sig(k40)

	         end do
                 do l=1,npas(k)
                    ntotal4=ntotal4+1
                 end do
              end do
              if(contr.gt.-1.e5)then
                 if(i.eq.1)then
	             sig(k40)=4.*sig(k40-1)/sqrt(float(k3)*pist(i))
                 else
                     sig(k40)=sigmamax/sqrt(pist(i))
                 end if
                 sigd(k40)=sig(k40)
              end if
	   end do
	end do

        if(ielimino.gt.0)then 
            write(*,'(i3,a32,f5.1,a34)') ielimino,' data points, being smaller than',cotaminima,', are NOT considered for inversion'
	    if(ici.eq.1)then
c	      open(icanal,file=control,access='append')
c	      write(icanal,'(i3,a31,f5.1,a34)') ielimino,' data points, being smaller than',cotaminima,', are NOT considered for inversion' 
c	      close(icanal)
	    endif
	endif

	if(nciclos.gt.0)then
	     if(iauto.eq.0)write( *,1000)"it","DE","s/n","chi**2","Mac1","Mac2","mic1","mic2"
     &                ,"fill2","Ic1/Ic2","stray %"
             if(iauto.eq.1)write(*,2000)"ite","DE","s/neqv","chi**2",
     &                                   (ca(jj),jj=1,18)

        end if




c _____________________________________________________________________________
c                        empieza el ciclo iterativo
c _____________________________________________________________________________

        factorrep=1
	diag=-1.d0
	chi0=1.d20
	varchi=1.d0
	isigo=1
	it=0
	iamplio=0
	ngu=0
        icalerr=0   !=0 pues de momento no vamos a calcular errores

	do while(isigo.eq.1)
	    it=it+1
	    diag0=diag
c            print*,'SIR 1131',m(1),m(4),m(6)
	    call marquardt2(stok,sig,nfrecs,atmosr,m,mfit,
     &		         covar,alpha,chisq,diag,iamplio,beta)
c            print*,'SIR 1134',m(1),m(4),m(6)

	    if(nciclos.le.0)then
		call leeuve2(1,uveout,ist,ntl,nlinsn,npasobs,dlamdaobs,scal)
                stop
	    end if           

	    nfree=nfrecs-mfit
  
	    if(it.eq.1)diag0=alamda0
	
	    varchi=(chi0-chisq)/(chi0+chisq)

	    call diagonal(diag0,diag,isigo,iamplio,varchi,it)

	    if(iamplio.eq.1)then  !converge
 
		ngu=ngu+1
	        chi0=chisq

c copiamos la atmosfera de salida (via common) atmosout en atmos para no 
c machacar el common
c y la duplicamos en atmoslin para pasarla al centro del disco con taulinea2
c la inversion se realiza en la linea de vision pero los modelos de entrada 
c y salida estan en el centro del disco

c		print*,'sir ',atmosout(4*ntau+1)
	        do i=1,16*ntau+5
	           atmos(i)=atmosout(i)
	        end do

	        do i=1,ntau
	           pg1b(i)=pg1(i)
	           pg2b(i)=pg2(i)
	           z1b(i)=z1(i)
	           z2b(i)=z2(i)
	           ro1b(i)=ro1(i)
	           ro2b(i)=ro2(i)
	        end do

                do i=1,k40   !numero de longitudes de onda
                   perfil(i)=scal(i)
                end do

	        fill2=atmos(16*ntau+4)
	        amac1=atmos(8*ntau+1)
	        amac2=atmos(16*ntau+3)
	        amic1=atmos(3*ntau+1)*1.e-5
	        amic2=atmos(11*ntau+3)*1.e-5
                porcien=atmos(16*ntau+5)

                contraste=contr
	        if(nlin(ntlblends).eq.0)contraste=scal(nliobs)
                contraste=(1+contraste)/(1.-contraste)

c                print*,'SIR 1189',m(1),m(4),m(6)
               	call comprime2(ntau,m,atmos,atmosr)
c		print*,'SIR 1191',m(1),m(4),m(6)

c	        if(ngu.lt.31)then
c	           cotavar=10**(-15.+float(ngu))-1.e-4
c	           if(cotavar.gt.0.005)cotavar=0.005
c                else 
c                   cotavar=(ngu-30)*.01
c                   if(ngu.eq.31)then
c                      print*,'More iterations will not be permitted unless significant variations of chi**2 occur.'
c                      print*,'If you still want to use the same nodes, restart the inversion using the output models.'
c                   end if            
c                end if

	        if(ngu.lt.25)then
	           cotavar=exp((-42+float(ngu))/3.)-1.e-4
	           if(cotavar.gt.0.0033)cotavar=0.0033
                else 
                   if(ngu.lt.50)cotavar=(ngu-24)*.003
                   if(ngu.ge.50)cotavar=(ngu-49)*.08
                   if(ngu.eq.50)then
                      print*,'More iterations will not be permitted unless significant variations of chi**2 occur.'
                      print*,'If you still want to use the same nodes, restart the inversion using the output models.'
                   end if            
                end if
                
	        if(varchi.le.cotavar)isigo=0
                if(ngu.eq.100)then 
                   isigo=0
                   print*,'100 iterations. If you want to continue, restart the inversion and'
                   print*,'use the output models of this cycle'                 
	        end if      

	        chiw=sumsq/float(nfrecs-ielimino)
	        snchi=1./sqrt(chiw)	

	        add=-20.
	        if(diag0.gt.0.d0)add=alog10(diag0)
                chprint=chireal/float(nfree-ielimino)

		ll=ll+1
		nguvec(ll)=ngu
                addvec(ll)=add
		snchivec(ll)=snchi
		chprintvec(ll)=chprint
		amac1vec(ll)=amac1
		amac2vec(ll)=amac2
		amic1vec(ll)=amic1
		amic2vec(ll)=amic2
		fill2vec(ll)=fill2
		contrastevec(ll)=contraste
		porcienvector(ll)=porcien
		if(isigo.eq.1)posicionciclo(ici)=ll
                kv=0
                do i=1,18
                   mvec(ll,i)=m(i)
                   if(mvecmax(i).lt.m(i))mvecmax(i)=m(i)
                end do
c                  print*,'SIR 1247',m(1),m(4),m(6)

	        if(iauto.eq.0)write(*,1002)ngu,add,snchi,chprint,amac1,amac2,amic1,amic2,
     &                        fill2,contraste,porcien

	        if(iauto.eq.1)write(*,2002)ngu,add,snchi,chprint,(m(jj),jj=1,18)
	    end if



	end do 
c _____________________________________________________________________________
c                        acaba el ciclo iterativo
c _____________________________________________________________________________


	print*,' '
	write(*,423)'============================ End of cycle',ici,' ================================='
	print*,' '

	if(ngu.eq.0)then
 	     call leeuve2(0,uveobs,ist,ntl,nlinsn,npasobs,dlamdaobs,scal)
             ierror=0 
	     call leemodi222(0,modelin1,modelin2,atmoslin,ntau,ierror,
     &                       z1,pg1,ro1,z2,pg2,ro2)
	     call leemodi222(1,modelout1,modelout2,atmoslin,ntau,ierror,
     &                       z1,pg1,ro1,z2,pg2,ro2)

 	 end if
 
	 call leeuve2(1,uveout,ist,ntl,nlinsn,npasobs,dlamdaobs,perfil)

	 do i=1,16*ntau+5
	    atmoslin(i)=atmos(i)
	 end do
         icalerr=1   !=1 pues vamos a calcular errores
         iauto=0
      
         mfit=0
c	 do i=1,5
c           if(mvecmax(i).gt.0)mfit=mfit+mvecmax(i)
c	 end do
         if(mvecmax(1).gt.0)mfit=mfit+mvecmax(1)
         ipa1=mfit       !ipa1 es el indice anterior a la gamma comp. 1
         if(mvecmax(2).gt.0)mfit=mfit+mvecmax(2)
         ipa11=mfit
	 do i=3,5
           if(mvecmax(i).gt.0)mfit=mfit+mvecmax(i)
	 end do
         iga1=mfit    !iga1 es el indice anterior a la gamma comp. 1
         if(mvecmax(6).gt.0)mfit=mfit+mvecmax(6)
         if(mvecmax(7).gt.0)mfit=mfit+mvecmax(7)
         ifi11=mfit   !ifi11 es el indice ultimo de la fi comp. 1
c	 do i=8,13
c           if(mvecmax(i).gt.0)mfit=mfit+mvecmax(i)
c	 end do
	 do i=8,9
           if(mvecmax(i).gt.0)mfit=mfit+mvecmax(i)
	 end do
	 ipa2=mfit         !ipa2 es el indice anterior a la gamma comp. 2
	 if(mvecmax(10).gt.0)mfit=mfit+mvecmax(10)
	 ipa22=mfit
	 do i=11,13
           if(mvecmax(i).gt.0)mfit=mfit+mvecmax(i)
	 end do	 
         iga2=mfit    !iga2 es el indice anterior a la gamma comp. 2
         if(mvecmax(14).gt.0)mfit=mfit+mvecmax(14)
         if(mvecmax(15).gt.0)mfit=mfit+mvecmax(15)
         ifi22=mfit   !ifi22 es el indice ultimo de la fi comp. 1
	 do i=16,18
           if(mvecmax(i).gt.0)mfit=mfit+mvecmax(i)
	 end do


         call comprime2(ntau,mvecmax,atmos,atmosr)

	 do j=1,mfit
            atmosrlin(j)=atmosr(j)
	 end do

	 call marquarderr(stok,nfrecs,atmosr,mvecmax,mfit)

	 call taulinea2(1,cth,ihemi,vx,atmos,ntau)

c	 call leemodi22(1,modelout1,modelout2,atmos,ntau,ierror)
	 call leemodi222(1,modelout1,modelout2,atmos,ntau,ierror,
     &                       z1,pg1,ro1,z2,pg2,ro2)
                                                        
         do j=1,iga1
            x(j)=atmosrlin(j)*(1.+sqrt(abs(errores(j))))
	 end do
         do j=iga1+1,ifi11
            x(j)=atmosrlin(j)+sqrt(abs(errores(j)))
	 end do
         do j=ifi11+1,iga2
            x(j)=atmosrlin(j)*(1.+sqrt(abs(errores(j))))
	 end do
         do j=iga2+1,ifi22
            x(j)=atmosrlin(j)+sqrt(abs(errores(j)))
	 end do
         do j=ifi22+1,mfit
            x(j)=atmosrlin(j)*(1.+sqrt(abs(errores(j))))
	 end do

         do i=1,ntau                      !nuevo para el calculo errores
	    pg1b(i)=pg1(i)
	    pg2b(i)=pg2(i)
	    z1b(i)=z1(i)
	    z2b(i)=z2(i)
	    ro1b(i)=ro1(i)
	    ro2b(i)=ro2(i)
         end do

	 call amp2err(ntau,mvecmax,atmoslin,x)

	 call taulinea2(1,cth,ihemi,vx,atmoslin,ntau)

	 do i=ntau+1,16*ntau+5
	    atmos(i)=abs(atmos(i)-atmoslin(i))     !errores
	 end do
	 do i=1,ntau                      !nuevo para el calculo errores
	    z1(i)=abs(z1(i)-z1b(i))
	    z2(i)=abs(z2(i)-z2b(i))
            pg1(i)=abs(pg1(i)-pg1b(i))
	    pg2(i)=abs(pg2(i)-pg2b(i))
            ro1(i)=abs(ro1(i)-ro1b(i))
	    ro2(i)=abs(ro2(i)-ro2b(i))
         end do

c escribimos los errores
c	 call leemodi22(1,modelerr1,modelerr2,atmos,ntau,ierror)
	 call leemodi222(1,modelerr1,modelerr2,atmos,ntau,ierror,
     &                       z1,pg1,ro1,z2,pg2,ro2)

	     
         close(13)

	 close(ican)

         chisnpr(ici)=chisn
         snpr(ici)=snn

c   	 open(25,file='chi',access='append')
c          write(25,*) snchi,chprint
c	 close(25)
	
	 open(78,file=snychi,access='append')
           write(78,*) snchi,chprint
         close(78)
	 end do	!fin del do en el numero de ciclos (ici)
	 close(icanal)



	open(icanal,file=control,access='append')
        if(iauto.eq.0)write(icanal,1000)"ite","DE","s/neqv","chi**2","Mac1","Mac2","mic1","mic2"
     &                ,"fill2","Ic1/Ic2","stray %"



        if(iauto.eq.1)write(icanal,2000)"ite","DE","s/neqv","chi**2",(ca(kk),kk=1,18)


             write(icanal,786)0,chisnpr(1),snpr(1)
!esta es la correcta


        do i=1,ll
 	  if(iauto.eq.0)write(icanal,1002)nguvec(i),addvec(i),snchivec(i),chprintvec(i),amac1vec(i),
     &                      amac2vec(i),amic1vec(i),amic2vec(i),
     &                      fill2vec(i),contrastevec(i),porcienvector(i)

 	  if(iauto.eq.1) write(icanal,2002)nguvec(i),addvec(i),snchivec(i),chprintvec(i),
     & (mvec(i,jj),jj=1,18)

	  do j=1,nciclos-1   
             if(i.eq.posicionciclo(j))then         
	       if(iauto.eq.0)write(icanal,1000)"ite","DE","s/neqv","chi**2","Mac1","Mac2",
     &                "mic1","mic2","fill2","Ic1/Ic2","stray %"
               if(iauto.eq.1)write(icanal,2000)"ite","DE","s/neqv","chi**2",
     &                          (ca(kk),kk=1,18)

               write(icanal,786)0,chisnpr(j+1),snpr(j+1)

	     endif
	  enddo
	enddo
	close(icanal)

       
	print*,'__________________________________________________________________________________'
	print*,' '
	print*,CARTEL
	print*,'__________________________________________________________________________________'
	print*,' '


c formatos

786	format(1x,i3,6x,1pe9.2,1x,e10.3)

1000     format(2x,a2,2x, a2,5x,a3,6x,a6,2x,  a6,1x,a6,3x,  a4,3x,a4,3x,
     &           a5,1x,a7,1x,a6)
1002	 format(1x,i3,1x,f4.0,1x,1pe9.2,1x, e10.3,1x, 0pf6.3,5(1x, f6.3),2x,f6.3)

2000     format(2x,a2,2x, a2,5x,a3,6x,a6,3x,18(a3))

2002	 format(1x,i3,1x,f4.0,1x,1pe9.2,1x, e10.3, 18(i3))


999	 format(3x,a4,1x,1pe9.2,1x)
423	format(a42,i2,a36)
	end
c ...................................................................

	subroutine diagonal(diag0,diag,isigo,iamplio,varchi,it)

	implicit real*4(a-h,o-z)
        real*4 diagon(100) 
	character diver*10,espacio*36

	common/iteradiagonal/itt
        common/repeticion/factorrep !factor de diag para evitar repeticiones
	data ir/1/

	diver=' increases '
	espacio=' _________________________________ '
  
        if(it.eq.1)then
           ir=1  
           itt=0       
	endif
      
	if(ir.eq.1)diagon(ir)=diag0
        ir=ir+1
	if(ir.eq.100)ir=6
        diagon(ir)=diag
                                                                
	if(diag.eq.0)then
	    isigo=0
	    iamplio=0
	    return
	else if(diag.gt.diag0)then
            isigo=1
	    iamplio=0
	    itt=itt+1
	else if(diag.le.diag0)then
            isigo=1
	    iamplio=1
	    itt=0
	end if

	if(ir.gt.5.)then
           if(diagon(ir).eq.diagon(ir-2).or.
     &                  diagon(ir-1).eq.diagon(ir-3))then
	     factorrep=0.
           else
	     factorrep=1.
           end if 
        end if

	if(ir.gt.7)then
          if(diagon(ir).eq.diagon(ir-1).and.
     &                  diagon(ir).eq.diagon(ir-2).and.
     &                  diagon(ir).eq.diagon(ir-3).and.
     &                  diagon(ir).eq.diagon(ir-4).and.
     &                  diagon(ir).eq.diagon(ir-5))factorrep=1.
        end if
       		
	if(diag0.gt.0)then
	     if(diag.gt.diag0)then
	          varchi=1.0
	          add=-20.
	          if(diag0.gt.0.0)add=alog10(diag0)
	          write(*,1100)add,diver,espacio
	          if(itt.ge.7)then
	               isigo=0
	               iamplio=0
	          end if
	     end if
	end if
	return
1100	format(5x,f4.0,1x,a10,a53)
	end

c _______________________________________________________________
c rutina convierte
c reescribe el perfil de luz difusa en las longitudes de onda de la 
c malla
c entradas: s     - perfil de luz difusa en las ldo observadas
c           nposi - array de enteros que contiene las posiciones
c	            de la malla correspondientes a las de las ldo
c                   observadas
c           nliobs- numero total de longitudes de onda observadas
c           nli   - numero total de longitudes de onda de la malla
c salidas:  stray - (via common) perfil de luz difusa en las ldo de la malla
c Basilio 24 Junio 1994
c ................................................................

        subroutine convierte(s,nposi,nliobs,nli)


	include 'PARAMETER' !para kld
	parameter (kld4=4*kld)           
	real*4 stray(kld4),s(*)
        integer nposi(*)
        common/difusa/stray

        do i=1,nliobs
           stray(nposi(i))=s(i)
  	end do

	return
        end 
c _______________________________________________________________
c "amp2err" construye la atmosfera completa a partir de la antigua y 
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
c Basilio y Jose Carlos 9-1-96 (modificacion eliminacion de cotas para errores)
c
c _______________________________________________________________________

	subroutine amp2err(ntau,m,atmos,atmosr)

	implicit real*4 (a-h,o-z)

	include 'PARAMETER'  !para kt
	integer m(*),mdata(18)
	real*4 atmos(*),atmosr(*)
	real*4 x(kt),y(kt),yy(kt),pert(14*kt+4),f(kt,kt)
	real*4 tau(kt),t1(kt),p1(kt),tnew1(kt),pnew1(kt)
	real*4 t2(kt),p2(kt),tnew2(kt),pnew2(kt)
	real*4 pg1(kt),z1(kt),ro1(kt),b1(kt),gam1(kt)
        real*4 pg2(kt),z2(kt),ro2(kt),b2(kt),gam2(kt)
        character*26 var(18) 
        integer icalerr
        common/preciso/prec      
        common/calerr/icalerr !si calculo errores=1 else =0
	common/zetas/pg1,z1,ro1,pg2,z2,ro2
        common/contornopg/ncontpg,pg01,pg02
        common/contornoro/ro01,ro02
        common/pgmag/ipgmag

	epsilon=1.d-2

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

	cota=1.e10
	cotapres=2.
        cotafi=1.e10 
        cota1=-1.e10
        cota2=1.e10

	do i=1,18	!do en grupos de varibles (1=t,2=p,...etc)

           ntau2=ntau
           if(i.eq.8.or.i.eq.16.or.i.eq.17.or.i.eq.18)ntau2=1!mac1,mac2,ff2,%

           if(m(i).eq.1)then	            !si pert. constante sumo
              kred=kred+1
              if(i.eq.7.or.i.eq.15.or.i.eq.6.or.i.eq.14)then
                 y1=atmosr(kred)-pert(kred)
  	         if(y1.lt.-cotafi)y1=-cotafi   !acoto inferiormente
	         if(y1.gt.cotafi)y1=cotafi     !acoto superiormente
c	      else if (i.eq.2.or.i.eq.10) then
c	         y1=atmosr(kred)-pert(kred)        !perturbaciona aditiva en el logaritmo
c  	         if(y1.lt.-cotapres)y1=-cotapres   !acoto inferiormente
c	         if(y1.gt.cotapres)y1=cotapres     !acoto superiormente
      	      else
                 if(pert(kred).eq.0)goto 999   !por aqui nunca va a pasar!!!!
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
c	      else if(i.eq.2 .or. i.eq.10)then   
c	         do j=1,ntau2
c                    atmos(kamp+j)=atmos(kamp+j)*exp(y1)
c                 end do	     
	      else
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
                 if(i.eq.7.or.i.eq.15.or.i.eq.6.or.i.eq.14)then
                    y(j)=atmosr(kred)-pert(kred)
  	            if(y(j).lt.-cotafi)y(j)=-cotafi   !acoto inferiormente
	            if(y(j).gt.cotafi)y(j)=cotafi     !acoto superiormente
	         else if (i.eq.2.or.i.eq.10) then
                    y1=atmosr(kred)-pert(kred)        !perturbaciona aditiva en la presion e
  	            if(y1.lt.-cotapres)y1=-cotapres   !acoto inferiormente
	            if(y1.gt.cotapres)y1=cotapres     !acoto superiormente
                    y(j)=y1*pert(kred1)               !escala con la presion en el ultimo nodo
      	         else
	            y(j)=(atmosr(kred)/pert(kred))-1.  !pert. multiplicativa
	            if(y(j).lt.cota1)y(j)=cota1 !acoto inferiormente
	            if(y(j).gt.cota2)y(j)=cota2  !acoto superiormente
	            y(j)=y(j)*pert(kred)     !perturb. aditiva en los nodos
                 end if
	      end do
    
	      if(ntau2.ne.ntau)then
		  print*,'a la macro o al ff se les asigna mas de 1 variable? '
	          stop
              end if

	      call splines22(x,y,m(i)-2,ntau,tau,yy,f)
              
c              if(i.eq.2.or.i.eq.10)then
c                  do j=1,ntau2
c                    yyy=yy(j)
c                    if(yyy .gt. cotapres)yy=cotapres
c                    if(yyy .lt. -cotapres)yy=-cotapres
c                    atmos(kamp+j)=atmos(kamp+j)*exp(yyy)
c                  end do
c	      else
                 do j=1,ntau2
                    atmos(kamp+j)=atmos(kamp+j)+yy(j)
                 end do    
c              end if                 
	   end if	     
	   kamp=kamp+ntau2
	   if(i.eq.8)kamp=kamp+ntau+1	!los ntau puntos de tau2 y el de ff1	
	end do

c en caso de que no se corrija la presion
c ponemos las presiones en equilibrio hidrostatico con las temperaturas
c	vart1=0.0
c	vart2=0.0
        do i=1,ntau
           tnew1(i)=atmos(ntau+i)
           tnew2(i)=atmos(9*ntau+2+i)
c	   avar1=abs((tnew1(i)-t1(i))/t1(i))
c	   avar2=abs((tnew2(i)-t2(i))/t2(i))
c	   vart1=amax1(vart1,avar1)           !maxima variacion de t1
c	   vart2=amax1(vart2,avar2)           !maxima variacion de t2
	end do

c	varp1=0.0
c	varp2=0.0
        do i=1,ntau
           pnew1(i)=atmos(2*ntau+i)
           pnew2(i)=atmos(10*ntau+2+i)
	   b1(i)=atmos(4*ntau+i)
           b2(i)=atmos(12*ntau+2+i)
	   gam1(i)=atmos(6*ntau+i)
           gam2(i)=atmos(14*ntau+2+i)
c	   pavar1=abs((pnew1(i)-p1(i))/p1(i))
c	   pavar2=abs((pnew2(i)-p2(i))/p2(i))
c	   varp1=amax1(varp1,pavar1)           !maxima variacion de p1
c	   varp2=amax1(varp2,pavar2)           !maxima variacion de p2
	end do


c ----------------------------------------------------------------------------
c coloco la p1 y p2 en equilibrio hidrostatico

	   if(m(2).eq.0) then	   
	        if(ncontpg.eq.0)call equisubmu(ntau,tau,tnew1,p1,pg1,z1,ro1)
                if(ncontpg.eq.-1)then
                   pg01=ro01*83145100.*tnew1(ntau)/1.302
                   pg1(ntau)=pg01
                   ro1(ntau)=ro01
                   call pgpefromrho(tnew1(ntau),ro1(ntau),p1(ntau),pg1(ntau))
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

	     if(m(10).eq.0) then
	        if(ncontpg.eq.0)call equisubmu(ntau,tau,tnew2,p2,pg2,z2,ro2)

                if(ncontpg.eq.-1)then
                   pg02=ro02*83145100.*tnew2(ntau)/1.302
                   pg2(ntau)=pg02
                   ro2(ntau)=ro02
                   call pgpefromrho(tnew2(ntau),ro2(ntau),p2(ntau),pg2(ntau))
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
c ----------------------------------------------------------------------------


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

999	var(1)='Temperatura (modelo 1)    '
        var(2)='Presion elec.  (modelo 1) '
        var(3)='V. microturb.   (modelo 1)'
        var(4)='Campo magnetico (modelo 1)'
        var(5)='Velocidad  (modelo 1)     '
        var(6)='Inclinacion mag.(modelo 1)'
        var(7)='Azimuth campo (modelo 1)  '
        var(8)='V. macroturb.   (modelo 1)'
        var(9)='Temperatura (modelo 2)    '
        var(10)='Presion elec.  (modelo 2) '
        var(11)='V. microturb.   (modelo 2)'
        var(12)='Campo magnetico (modelo 2)'
        var(13)='Velocidad  (modelo 2)     '
        var(14)='Inclinacion mag.(modelo 2)'
        var(15)='Azimuth campo (modelo 2)  '
        var(16)='V. macroturb.   (modelo 2)'
        var(17)='Factor llenado (modelo 2) '
        var(18)='Peso % de la luz difusa '
       
        print*,' '
        print*,'*********************************************************************'
	print*,'se esta tratando de invertir la variable '
        print*,var(i)
        print*,'y su valor inicial es 0 en algun punto'
        print*,'o inicializas de otra manera o no inviertas dicha variable'
        print*,'*********************************************************************'
        print*,' '
        stop 'en amp2err'

	end
c _______________________________________________________________________

