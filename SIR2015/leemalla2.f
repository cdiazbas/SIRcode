c rutina leemallab
c lee el fichero de control de lineas y longitudes de onda:'mallaobs'
c ntau : numero de puntos en tau
c tau  : log10(tau)
c ntl  : numero total de lineas
c nlin : indice de cada linea y de cada blend
c npas : numero de puntos en cada linea
c dlamda:cada delta de l.d.o. en ma
c nble :numero de blends de cada linea

	subroutine leemalla2(mallaobs,ntl,nli,nliobs,nlin,npasobs,
     &                         npas,dlamdaobs,dlamda,nble,nposi,indice)

	include 'PARAMETER'  !por kl
	integer npasobs(*),npas(kl),nble(kl),nd(kl),nposi(*),nlin(kl),ndd(kl,kl),nblee(kl)
	real*4 dlamdaobs(*),dlamda(kld),dini,dipa,difi,dinii(kl),dipaa(kl),difii(kt)
	integer indice(*),nddd(50),blanco
        character mallaobs*(*),mensajito*31
	character*80 men1,men2,men3
	character*100 control
	common/canal/icanal
	common/nombrecontrol/control
        common/nciclos/nciclos   
        common/contraste/contr

	men1=' '
	men2=' '
	men3=' '

        if(contr.gt.-1.e5.and.nciclos.gt.0)then
           ntl=ntl-1
           nliobs=nliobs-1
        end if 

	ican=57
	mensajito=' containing the wavelength grid'
	call cabecera(ican,mallaobs,mensajito,ifail)
	if(ifail.eq.1)goto 999
        
	jj=0
        k=0 
        nli=0
	numlin=0
	do i=1,1000   
	  call mreadmalla(ican,ntonto,nddd,ddinii,ddipaa,ddifii,ierror,blanco)
	  if(ierror.eq.1)goto 8   !si he llegado al final del fichero
	  if(blanco.ne.1)then 
	    numlin=numlin+1
	    if(numlin.gt.kl.or.ntonto.gt.kl)then 
              men1='STOP: The number of lines in the wavelength grid is larger than the '
	      men2='      current limit. Decrease this number or change the PARAMETER file.'
	      call mensaje(2,men1,men2,men3)
            end if
	     
	    nblee(numlin)=ntonto
            dinii(numlin)=ddinii
            dipaa(numlin)=ddipaa
            difii(numlin)=ddifii

	    do j=1,nblee(numlin)
              ndd(numlin,j)=nddd(j)
            enddo
          endif   !el de blanco.ne.1
	enddo
8	close(ican)

	ntl=numlin

c	Ahora las ordenamos como en los perfiles:

	do i=1,ntl     !indice en los perfiles 
	   icheck=0
           do l=1,numlin   !indice en la malla
c	       print*,'indice(i) vale',i,indice(i)
               if(indice(i).eq.ndd(l,1).and.icheck.lt.1)then   !he encontrado una
		   icheck=icheck+1
                   nble(i)=nblee(l)
                   do ll=1,nblee(l)
                       nd(ll)=ndd(l,ll)
                   enddo
                   dini=dinii(l)
                   dipa=dipaa(l)
                   difi=difii(l)
	 
	           do j=1,nble(i)
	              jj=jj+1
                      nlin(jj)=nd(j)
c		      print*,'nlin,jj',jj,nlin(jj)
	           end do
                   ntlblends=jj
	           pasmin=1.e30
                   primero=dlamdaobs(k+1) 
                   ultimo=dlamdaobs(k+npasobs(i)) 
                   do j=1,npasobs(i)-1
                        k=k+1
                        pasoobs=dlamdaobs(k+1)-dlamdaobs(k) 
                        pasmin=min(pasmin,pasoobs) !cambio amin0 por min
	           end do

                   k=k+1
c                  if(dipa.gt.pasmin)dipa=pasmin      !provisional!!!!
                  if(dini.gt.primero)dini=primero
                  if(difi.lt.ultimo)difi=ultimo
                  ntimes=int(pasmin/dipa)	 !numero de veces que la red fina divide a la obs 
c                  dipa=pasmin/ntimes
 	          
                  n1=nint((primero-dini)/dipa)  !numero de puntos anteriores
                  n2=nint((difi-ultimo)/dipa)   !nuero de puntos posteriores
                  dini=primero-n1*dipa
                  difi=ultimo+n2*dipa

c		  print*,'primero,ultimo,n1,n2,dini,difi',primero,ultimo,n1,n2,dini,difi

	          npas(i)=(difi-dini)/dipa+1
		  if(10*( (difi-dini)/dipa+1 -int( (difi-dini)/dipa+1 ) ).gt..5) npas(i)=npas(i)+1
        	  do j=1,npas(i)
	             nli=nli+1
	             dlamda(nli)=dini+(j-1)*dipa
c	print*,'en la malla,la serie es',nli,dlamda(nli),dipa
	          end do
              endif
	   end do  !fin del do en l
	   if(icheck.eq.0)then
	      men1='STOP: There are lines in the observed/stray light profiles which'
	      men2='      do not appear in the wavelength grid.'
	      men3='      Check also the wavelength grid for missing : symbols.'
	      call mensaje(3,men1,men2,men3)
	   endif 
9	enddo      !fin del do en i
c		print*,'pasos 1 y 2',npas(1),npas(2)

        k=1
	epsilon=dipa/2.
        do i=1,nliobs
           dd=dlamdaobs(i)
           do while(k.lt.nli.and.dd.gt.dlamda(k)+epsilon)
              k=k+1
           end do          
           do while(k.lt.nli.and.dd.lt.dlamda(k)-epsilon)
              k=k+1
           end do
           do while(k.lt.nli.and.dd.gt.dlamda(k)+epsilon)
              k=k+1
           end do  
           nposi(i)=k
        end do  

        if(contr.gt.-1.e5.and.nciclos.gt.0)then
           ntl=ntl+1
           ntlblends=ntlblends+1	
           nlin(ntlblends)=0
           npas(ntl)=1
           nble(ntl)=1
           nli=nli+1
	   nliobs=nliobs+1 
           dlamda(nli)=0.
           nposi(nliobs)=nli
        end if 

212	print*,'Number of wavelengths in the wavelength grid  : ',nli
	close(ican)

	return

c	------------------------------------------------------------------------
c	Mensajes de error:

999	men1='STOP: The file containing the wavelength grid does NOT exist:'
	men2=mallaobs
	call mensaje(2,men1,men2,men3)

992	men1='STOP: Incorrect format in the file containing the wavelength grid:'
        men2=mallaobs
	call mensaje(2,men1,men2,men3)



	end
