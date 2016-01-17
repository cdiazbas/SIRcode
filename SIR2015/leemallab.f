c rutina leemallab
c lee el fichero de control de lineas y longitudes de onda:'mallaobs'
c ntau : numero de puntos en tau
c tau  : log10(tau)
c ntl  : numero total de lineas
c nlin : indice de cada linea
c npas : numero de puntos en cada linea
c dlamda:cada delta de l.d.o. en ma
c nble :numero de blends de cada linea

	subroutine leemallab(mallaobs,ntl,nli,nlin,npas,dlamda,nble)

	include 'PARAMETER'

	integer nlin(*),npas(*),nble(*),nd(kl),blanco,ntl,nli,ierror
	real*4 dlamda(*),dini,dipa,difi
        character mallaobs*(*)
	character*100 control
	character*80 men1,men2,men3
	character*31 mensajito
	common/canal/icanal
	common/nombrecontrol/control

	men1=' '
	men2=' '
	men3=' '

        
	ican=57
	mensajito=' containing the wavelength grid'
	call cabecera(ican,mallaobs,mensajito,ifail)
	if(ifail.eq.1)goto 999

        jj=0
        k=0 
        nli=0
	numlin=0
	do i=1,1000   !Vamos a tener mas de 1000 lineas???? Pues entonces....  
	  
	  call mreadmalla(ican,ntonto,nd,dini,dipa,difi,ierror,blanco)
	  if(ierror.eq.1)goto 8
	  if(blanco.ne.1)then

	    numlin=numlin+1   
	    if(numlin.gt.kl.or.jj+ntonto.gt.kl)then !comprobamos numero de lineas
	     men1='STOP: The number of lines in the wavelength grid is larger than the '
	     men2='      current limit. Decrease this number or change the PARAMETER file.'
	     call mensaje(2,men1,men2,men3)
            end if
	    nble(numlin)=ntonto

	    do j=1,nble(numlin)
	      jj=jj+1
              nlin(jj)=nd(j)
c	      print*,'el indice de la linea es', nlin(jj)
	    end do

	    npas(numlin)=(difi-dini)/dipa+1
            if(10*( (difi-dini)/dipa+1 -int( (difi-dini)/dipa+1 ) ).gt..5) npas(numlin)=npas(numlin)+1
            if(nli+npas(numlin).gt.kld)then  !comprobamos numero longitudes onda
	       men1='STOP: The wavelength grid has more points than the current limit kld.'
	       men2='      Decrease the number of wavelengths or change the PARAMETER file.'
	       call mensaje(2,men1,men2,men3)
            end if
	    do j=1,npas(numlin)
	       nli=nli+1
	       dlamda(nli)=dini+(j-1)*dipa
	    end do
         endif
	enddo
8	ntl=numlin

	print*,'Number of wavelengths in the wavelength grid: ',nli
	close(ican)

	return

c       -------------------------------------------------------------------------
c	Mensajes de error:

999	men1='STOP: The file containing the wavelength grid does NOT exist:'
	men2=mallaobs
	call mensaje(2,men1,men2,men3)

992	men1='STOP: Incorrect format in the file containing the wavelength grid:'
        men2=mallaobs
	call mensaje(2,men1,men2,men3)


	end
