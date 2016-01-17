c leeuveobs lee un fichero de datos de parametros de stokes :vobs
	
	subroutine leeuveobsindic(vobs,ist,ntl,nli,nlin,npas,dl,nble,stok)

	include 'PARAMETER'  !por kld
	real*4 stok(*),dl(*),si(kld),sq(kld),su(kld),sv(kld)
	integer npas(*),ist(*),nble(*),nlin(*),nli,icanal
	character vobs*(*)
	character*100 control
	character*80 men1,men2,men3
c	character mensajito*45
	common/canal/icanal
	common/nombrecontrol/control

        common/contraste/contr
	men1=' '
	men2=' '
	men3=' '

	ican=56

	open(ican,file=vobs,status='old',err=991)
	k=0
	n0=-1
        i=0
	do while(k.lt.kld)
	   k=k+1
	   read(ican,*,err=992,end=990)a,dl(k),si(k),sq(k),su(k),sv(k)
           n1=nint(a)
           if(n1.ne.n0)then
             j=1
             i=i+1

	     if(i.gt.kl) then  !Si el numero de lineas es mayor que kl
	         men1='STOP: The number of lines in the observed/stray light profiles is larger than'
	         men2='      the current limit. Decrease this number or change the PARAMETER file.'       
	         call mensaje(2,men1,men2,men3)
	     endif
             
	     nlin(i)=n1     !para eliminar la malla si es necesario
           endif
	   n0=n1
           npas(i)=j
           j=j+1
	end do
990     ntl=i
        nli=k-1     !este es el numero de longitudes de onda
	
	if(nli.gt.kld)then  
	   men1='STOP: The number of wavelengths in the observed/stray light profiles is larger than'
	   men2='      the current limit kld. Decrease this number or change the PARAMETER file.'
	   call mensaje(2,men1,men2,men3)
	endif

	do j=1,ntl
	    nble(j)=1      !si no tenemos malla, no tenemos blends
	    do jj=j+1,ntl
               if(nlin(j).eq.nlin(jj+1))goto 993
	    enddo
       	enddo


	print*,'Number of wavelengths in the observed profiles: ',nli 
c	open(icanal,file=control,fileopt='eof')
c	write(icanal,*)'Number of wavelengths in the observed profiles: ',nli 
c	close(icanal)

        if(contr.gt.-1.e5)then
           nli=nli+1
           ntl=ntl+1
           npas(ntl)=1
           si(nli)=contr
           sq(nli)=0.
           su(nli)=0.
           sv(nli)=0.
           dl(nli)=0.
        end if

	call sfromiquv(ist,nli,si,sq,su,sv,stok) !cuelga de leeuve2
	close(ican)
	return


c	Mensajes de error:

991	men1='STOP: The file containing the observed/stray light profiles does NOT exist:'
	men2=vobs
	call mensaje(2,men1,men2,men3)

992	men1='STOP: Incorrect format in the file containing the observed/stray light profiles:'
	men2=vobs
	call mensaje(2,men1,men2,men3)

993	men1='STOP: The wavelength samples of a given line have to be written contiguously'
	men2='      in the file containg the observed/stray light profiles:'
	men3=vobs
	call mensaje(3,men1,men2,men3)

	end
