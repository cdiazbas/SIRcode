c leeuve2 lee un fichero de datos de parametros de stokes :vobs
	
	subroutine leeuve2(idis,vobs,ist,ntl,nlin,npas,dl,stok)

	include 'PARAMETER' !por kld
	real*4 stok(*),dl(*),si(kld),sq(kld),su(kld),sv(kld)
	integer npas(*),nlin(*),ist(*),idis
	character*100 vobs
	character*100 control
	character*80 men1,men2,men3
	common/canal/icanal
	common/nombrecontrol/control

	men3=' '
        men1=' '
	men2=' '

	ican=56

	if(idis.eq.0)then
c	   ican=23
c	   mensajito=' containing the observed/stray light profiles'
c	   call cabecera(ican,vobs,mensajito,ifail)
c	   if(ifail.eq.1)goto 991
		
	        open(ican,file=vobs,status='old',err=991)
		k=0
		do i=1,ntl
		   do j=1,npas(i)
			k=k+1
			read(ican,*,err=992)a,dl(k),si(k),sq(k),su(k),sv(k)
	                nlin(i)=nint(a)
		   end do
		end do
	        call sfromiquv(ist,k,si,sq,su,sv,stok)
	else
c           print*,'estoy en leeuve2 ntl=',ntl
		open(ican,file=vobs)
		k=0
		do i=1,ntl
c                       print*,'estoy en leeuve2 i npas =',i,npas(i)
		   do j=1,npas(i)
			k=k+1
	           end do
	        end do
c                print*,'estoy en leeuve2  voy a llamar a iquvfroms'
	        call iquvfroms(ist,k,si,sq,su,sv,stok)
c                print*,'estoy en leeuve2  voy a escribir',vobs


	        k=0
		do i=1,ntl
		   do j=1,npas(i)
	              k=k+1
c                      print*,        nlin(i),dl(k),si(k),sq(k),su(k),sv(k)
                      write(ican,993)nlin(i),dl(k),si(k),sq(k),su(k),sv(k)
c                      print*,nlin(i),dl(k),si(k),sq(k),su(k),sv(k)
		   end do
		end do
c                print*,'estoy en leeuve2  voy a salir',vobs
c                stop

	end if
	close(ican)
	return


c	Mensajes de error:
991	men1='STOP: The file containing the observed/stray light profiles does NOT exist:'
	men2=vobs
	call mensaje(2,men1,men2,men3)

992	men1='STOP: Incorrect format in the file containing the observed/stray light profiles:'
	men2=vobs
	call mensaje(2,men1,men2,men3)

c formato de escritura
993     format(1x,i5,1x,f11.4,1x,4(e14.7,1x))

	end



c______________________________________________________________________
c sfromiquv guarda en stok los vectores si,sq,su,sv segun ist(4)
c si por ejemplo ist=1,0,1,0 stok contiene si y a continuacion su

	subroutine sfromiquv(ist,n,si,sq,su,sv,stok)

	integer ist(*),n
	real*4 si(*),sq(*),su(*),sv(*),stok(*)

	ks=0
	if(ist(1).eq.1)then
	   do i=1,n
	      ks=ks+1
	      stok(ks)=si(i)
	   end do
	end if
	if(ist(2).eq.1)then
	   do i=1,n
	      ks=ks+1
	      stok(ks)=sq(i)
	   end do
	end if
	if(ist(3).eq.1)then
	   do i=1,n
	      ks=ks+1
	      stok(ks)=su(i)
	   end do
	end if
	if(ist(4).eq.1)then
	   do i=1,n
	      ks=ks+1
	      stok(ks)=sv(i)
	   end do
	end if

	return
	end
c______________________________________________________________________
c iquvfroms divide stok en los vectores si,sq,su,sv segun ist(4)
c si por ejemplo ist=1,0,1,0 sq=0 sv=0

	subroutine iquvfroms(ist,n,si,sq,su,sv,stok)

	integer ist(*),n
	real*4 si(*),sq(*),su(*),sv(*),stok(*)

	        ks=0
	        if(ist(1).eq.1)then
	           do i=1,n
	              ks=ks+1
	              si(i)=stok(ks)
	           end do
	        else
	           do i=1,n
	              si(i)=0.
	           end do
	        end if
	        if(ist(2).eq.1)then
	           do i=1,n
	              ks=ks+1
	              sq(i)=stok(ks)
	           end do
	        else
	           do i=1,n
	              sq(i)=0.
	           end do
	        end if
	        if(ist(3).eq.1)then
	           do i=1,n
	              ks=ks+1
	              su(i)=stok(ks)
	           end do
	        else
	           do i=1,n
	              su(i)=0.
	           end do
	        end if
	        if(ist(4).eq.1)then
	           do i=1,n
	              ks=ks+1
	              sv(i)=stok(ks)
	           end do
	        else
	           do i=1,n
	              sv(i)=0.
	           end do
	        end if
	
	return
	end
