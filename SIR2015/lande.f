c_______________________________________________________________
c	Lande
c	Calcula los factores de Lande efectivos
        integer mc
	parameter (mc=20)	!numero maximo de componentes zeeman
c	parameter (zeff=5.)	!semiempirico para el damping gamma6
	integer mult(2),ji(2),jf(2),np,nl,nr
	real tam(2),G(2),loggf,wlengt,zeff,alfa,sigma
	real dlp(mc),dll(mc),dlr(mc),sp(mc),sl(mc),sr(mc)
	character design(2)*1,atom*2,linea*100	!,multno*6

	character*100 nomlineas
	character*100 fichabun
        common/ficlineas/nomlineas
        common/fichabun/fichabun

c	nomlineas='LINEAS'
	print*,'atomic parameter filename:'
	read*,nomlineas
	fichabun='THEVENIN'
	PRINT*,'NUMERO DE LA LINEA'
	READ*,ILN

	CALL LEELINEASii(ILN,ATOM,ISTAGE,WLENGT,zeff,ENERGY,
     &                 LOGGF,MULT,DESIGN,TAM,alfa,sigma)




c	____________________________________________________________
c	PARAMETROS ATOMICOS
c	llamo a ATMDAT que devuelve weight (peso molecular),abu (abunda
c	cia), chi1,chi2 (pot.de ionizacion del atomo neutro e ion),
c	u1,u2,u3 (funciones de particion atomo neutro,ion,ion2).
C
c	PRINT*,'TEMPERATURA =5000 '
C	READ*,TT
	TT=5000.
c	PRINT*,'MULTIPLETE ',MULTNO
	call ATMDATB(atom,TT,nel,weight,abu,chi10,chi20,u1,u2,u3,du1,du2,du3)
	PRINT*,' ELEM=',NEL,' WEIGHT=',WEIGHT,' ABU=',ABU
	PRINT*,'CHI1=',CHI10,' CHI2=',CHI20
	PRINT*,'U1=',U1, 'U2= ',U2,' U3=',U3
	PRINT*,'DU1=',DU1, 'DU2= ',DU2,' DU3=',U3
c	________________________________________________________________
c	@6)SUBNIVELES ZEEMAN
c	calculo la parte entera "ji" y la parte fraccionaria "jf" del
c	momento angular total (tam), necesarios para Zeeman
	do i=1,2
	   ji(i)=int(tam(i))
	   jf(i)=int(10*(tam(i)-ji(i)))
c	   print*,'ji(',i,')=',ji(i),'  jf(',i,')=',jf(i)
	end do

c	llamo a SUBLANDE para que calcule el numero de transiciones PI,
c	L,R (NP,NL,NR),desplazamientos de cada componente (dlp,dll,dlr)
c	e intensidad (sp,sl,sr).dlo se calculo en el apartado @5

	call SUBLANDE2(mc,mult,design,tam,ji,jf,dlo,np,nl,nr,dlp,dll,dlr,
     &	sp,sl,sr,G)
	
	T1=TAM(1)
	T2=TAM(2)

	FACTOR=(T1*(T1+1)-T2*(T2+1))*(G(1)-G(2))/4.
	FACTOR=FACTOR+(G(2)+G(1))/2.

	PRINT*,' FACTOR DE LANDE EFECTIVO (LS) = ',FACTOR
	print*,'lower and upper g values:',g(1),g(2)
	END


