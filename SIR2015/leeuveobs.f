c leeuveobs lee un fichero de datos de parametros de stokes :vobs
	
	subroutine leeuveobs(vobs,ist,ntl,nli,npas,dl,stok)

	include 'PARAMETER'  !por kld
	real*4 stok(*),dl(*),si(kld),sq(kld),su(kld),sv(kld)
	integer npas(*),ist(*)
	character vobs*(*)
	character*100 control
	common/canal/icanal
	common/nombrecontrol/control

        common/contraste/contr

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
           endif
	   n0=n1
           npas(i)=j
           j=j+1
	end do
990     ntl=i
        nli=k-1

	print*,'Number of wavelengths in the observed/stray light profiles: ',nli 
c	open(icanal,file=control,fileopt='eof')
c	write(icanal,*)'Number of wavelengths in the observed/stray light profiles: ',nli 
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

	call sfromiquv(ist,nli,si,sq,su,sv,stok)  !cuelga de leeuve2
	close(ican)
	return

991	print*,' '
        print*,'STOP: The file containing the observed or stray light profiles does NOT exist:'
        print*,vobs
        print*,'_______________________________________________________________________________'
c	open(icanal,file=control,fileopt='eof')
c        write(icanal,*) 'STOP: The file containing the observed or the stray light profiles does NOT exist:'
c        write(icanal,*) vobs
c	close(icanal)
	stop

992	print*,' '
	print*,'STOP: Incorrect format in the file containing the observed or stray light profiles:'
        print*,vobs
        print*,'_______________________________________________________________________________'
c	open(icanal,file=control,fileopt='eof')
c	write(icanal,*) 'STOP: Incorrect format in the file containing the observed or stray light profiles'
c        write(icanal,*) vobs
c	close(icanal)

	stop


	end
