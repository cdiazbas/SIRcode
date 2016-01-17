c fperfil2
c dada una atmosfera (2 componentes) "atmos" que entra por el common
c y dada una perturbacion "atmosr" esta rutina construye la nueva 
c atmosfera (mediante amp2), calcula los nuevos perfiles "scal" y las nuevas 
c fr (blends0), convoluciona los perfiles y las fr con la macro y 
c la psf (deconv y deconv2) y calcula las derivadas de la chi^2 que 
c salen por "dscal"
c
c Basilio 23-3-93
c Modificacion luz difusa Basilio 20-6-94
c ________________________________________________________________________
	subroutine fperfil2err(m,atmosr,scal,dscal)

	implicit real*4 (a-h,o-z)
 
	include 'PARAMETER'   !por kt,kn,kl,kld
	parameter (kt8=8*kt+2,kt16=16*kt+5)         !
	parameter (kl4=4*kl)    !numero maximo de lineas
	parameter (kld4=4*kld)
	parameter (kldt=kld*kn,kldt4=4*kldt)

c para la malla
	integer nlin(kl),npas(kl),ist(4),m(*),mdata(18),mnod1(8),mnod2(8)
	integer nlins(kl4),npass(kl4),nble(kl)
	real*4 dlamda0(kl),dlamda0s(kl4),dlamda(kld),dlamdas(kld4)
c para la atmosfera
	real*4 atmosr(*),atmos(kt16),atmostry(kt16)
	real*4 atmos1(kt8),atmos2(kt8)
c para los perfiles y f. respuesta
	real*4 scal(*),scal1(kld4),scal2(kld4),dscal(*),sin0(kld4),stray(kld4)
	real*4 rt(kldt4),rh(kldt4),rv(kldt4),rg(kldt4),rf(kldt4)
	real*4 rp(kldt4),rm(kldt4)
c comunes
	common/responde2/ist,ntau,ntl,nlin,npas,nble
        common/ldeo/dlamda,dlamda0
        common/smalla/ntls,nlins,npass
        common/smalla1/dlamdas
        common/atmosfera/atmos
        common/atmosferaout/atmostry
	common/nohaycampo/nohaycampo1,nohaycampo2
        common/primera2/ntotal,ntotal4,ists
        common/nciclos/nciclos   !del principal
        common/difusa/stray
	common/numtotblends/ntlblends

	data iprimeravez/0/

	epsilon=2.e-5 !precision para comparar reales


	nohaycampo=nohaycampo1+nohaycampo2


	if(iprimeravez.eq.0)then
           iprimeravez=1

	   if(nohaycampo.eq.0)then
              ntls=ntl
              do i=1,ntl
                 npass(i)=npas(i)
              end do
	      do i=2,4
                 ist(i)=0
              end do
              ist(1)=1
           end if
	   ntotal=0
           iii=0
   	   do i=1,ntl
	      ntotal=ntotal+npas(i)
              do ii=1,nble(i)
                 iii=iii+1 
              end do
              ntlblends=iii
	   end do
	   ists=ist(1)+ist(2)+ist(3)+ist(4)
	   ntotal4=ntotal*ists
	end if

        nli=0
	do i=1,ntl
           nli=nli+npas(i)
	end do
        

c copiamos la atmosfera antigua en atmostry; "amp2" escribira en 
c atmostry la nueva atmosfera
	do i=1,16*ntau+5
	   atmostry(i)=atmos(i)
	end do


	if(nciclos.ne.0)call amp2(ntau,m,atmostry,atmosr)

        do i=1,16
	   mdata(i)=i*ntau+2*int(i/9)  ! indi. anteri. a la var. i (ampliada)
	end do
        mdata(17)=mdata(16)+1
        mdata(18)=mdata(17)+1

c la macro y los f.f son
	vmac1=atmostry(mdata(8)+1)
	vmac2=atmostry(mdata(16)+1)
	fill2=atmostry(mdata(17)+1)
	fill1=1.0-fill2
        peso2=atmostry(mdata(18)+1)/100.
        peso1=1.-peso2
c	peso1=1.     !luz difusa com velo

	do i=1,ntotal4
           scal1(i)=0.
           scal2(i)=0.
           sin0(i)=0.
	end do

	do i=1,8
           mnod1(i)=m(i)
           mnod2(i)=m(i+8)
	end do

c dividimos la atmosfera en las 2 componentes
        do i=1,8*ntau+2
           atmos1(i)=atmostry(i)
           atmos2(i)=atmostry(i+8*ntau+2)
	end do

c calculamos las funciones respuesta y el perfil observado para cada atmosfera
c solo en el caso de que fill1.ne.0 se hace lo siguiente
c ************************* atmosfera 1 ***************************************

	k=0	!inicializo el indice de las derivadas
	if(abs(fill1).gt.epsilon)then
	   if(nohaycampo1.eq.0.and.ist(1).ne.0)then
	      call blendscon2(atmos1,scal1,rt,rp,rv,rm,mnod1)
	      if(nlin(ntlblends).eq.0)then
                 do i=nli,m(1)*ntotal,ntotal               
	            rt(i)=rt(i)/2./fill1
	         end do
                 do i=nli,m(2)*ntotal,ntotal               
	            rp(i)=rp(i)/2./fill1
                 end do
              end if
	      if(vmac1.gt.0)then
	         do j=1,ntotal
	            sin0(j)=scal1(j)
	         end do
	         call deconv(scal1,1,ntl,npas,dlamda0,dlamda,vmac1)
                 call deconv2(sin0,1,ntl,npas,dlamda0,dlamda,vmac1)
              end if 

              fill=fill1*peso1
	      call sub1err(m(1),vmac1,fill,rt,k,dscal)      !temperatura 
	      call sub1err(m(2),vmac1,fill,rp,k,dscal)      !presion  
	      call sub1err(m(3),vmac1,fill,rm,k,dscal)      !micro 
              call cero1err(m(4),ntotal4,k,dscal)           !campo 
	      call sub1err(m(5),vmac1,fill,rv,k,dscal)      !velocidad 
              call cero1err(m(6),ntotal4,k,dscal)           !inclinacion 
              call cero1err(m(7),ntotal4,k,dscal)           !azimuth 
              call mult1err(m(8),ntotal4,k,fill,vmac1,sin0,dscal) !macro
   	   else  !o sea si hay campo1 

       	      call blends2(atmos1,scal1,rt,rp,rh,rv,rg,rf,rm,mnod1)

c Dado que dlamda0 entra en el common via blends2 la sentencia siguiente 
c no puede colocarse antes de la llamada a blends2  
	      k1=0
	      do i=1,4
	         do j=1,ist(i)
	            do klin=1,ntl
	               k1=k1+1
	               dlamda0s(k1)=dlamda0(klin)
	            end do
	         end do
              end do

	      if(nlin(ntlblends).eq.0)then
                 do i=nli,m(1)*ntotal4,ntotal4               
	            rt(i)=rt(i)/2./fill1
	         end do
                 do i=nli,m(2)*ntotal4,ntotal4               
	            rp(i)=rp(i)/2./fill1
                 end do
              end if
c calculamos la convolucion con el perfil gaussiano y la psf
c calculamos la convolucion con la derivada del perfil gaussiano y la psf
	      if(vmac1.gt.0)then
	         do j=1,ntotal4
	            sin0(j)=scal1(j)
	         end do
	         call deconv(scal1,1,ntls,npass,dlamda0s,dlamdas,vmac1)
                 call deconv2(sin0,1,ntls,npass,dlamda0s,dlamdas,vmac1)
              end if 
        
c calculamos la convolucion de las funciones respuesta
              fill=fill1*peso1
	      call sub2err(m(1),vmac1,fill,dlamda0s,ist,rt,k,dscal) !tempe.
	      call sub2err(m(2),vmac1,fill,dlamda0s,ist,rp,k,dscal) !presi.
	      call sub2err(m(3),vmac1,fill,dlamda0s,ist,rm,k,dscal) !micro.
	      call sub2err(m(4),vmac1,fill,dlamda0s,ist,rh,k,dscal) !campo
	      call sub2err(m(5),vmac1,fill,dlamda0s,ist,rv,k,dscal) !veloc.
	      call sub2err(m(6),vmac1,fill,dlamda0s,ist,rg,k,dscal) !incli.
	      call sub2err(m(7),vmac1,fill,dlamda0s,ist,rf,k,dscal) !azimuth.
              call mult1err(m(8),ntotal4,k,fill,vmac1,sin0,dscal)   !macro
           end if
	end if

c solo en el caso de que fill2.ne.0 se hace lo siguiente
c ************************* atmosfera 2 ***************************************

	if(abs(fill2).gt.epsilon)then
	   if(nohaycampo2.eq.0.and.ist(1).ne.0)then
              call blendscon2(atmos2,scal2,rt,rp,rv,rm,mnod2)

	      if(nlin(ntlblends).eq.0)then
                 do i=nli,m(9)*ntotal,ntotal               
	            rt(i)=-rt(i)/2./fill2
	         end do
                 do i=nli,m(10)*ntotal,ntotal               
	            rp(i)=-rp(i)/2./fill2
                 end do
              end if
	      if(vmac2.gt.0)then
	         do j=1,ntotal
	            sin0(j)=scal2(j)
	         end do
	         call deconv(scal2,1,ntl,npas,dlamda0,dlamda,vmac2)
	         call deconv2(sin0,1,ntl,npas,dlamda0,dlamda,vmac2)
              end if
              fill=fill2*peso1
	      call sub1err(m(9),vmac2,fill,rt,k,dscal)       !temperatura  
	      call sub1err(m(10),vmac2,fill,rp,k,dscal)      !presion   
	      call sub1err(m(11),vmac2,fill,rm,k,dscal)      !micro
              call cero1err(m(12),ntotal4,k,dscal)            !campo  
	      call sub1err(m(13),vmac2,fill,rv,k,dscal)      !velocidad
              call cero1err(m(14),ntotal4,k,dscal)            !inclinacion  
              call cero1err(m(15),ntotal4,k,dscal)            !azimuth 
              call mult1err(m(16),ntotal4,k,fill,vmac2,sin0,dscal) !macro

	   else  !o sea si hay campo2 

       	      call blends2(atmos2,scal2,rt,rp,rh,rv,rg,rf,rm,mnod2)

c Dado que dlamda0 entra en el common via blends2 la sentencia siguiente 
c no puede colocarse antes de la llamada a blends2  
	      k1=0
	      do i=1,4
	         do j=1,ist(i)
	            do klin=1,ntl
	               k1=k1+1
	               dlamda0s(k1)=dlamda0(klin)
	            end do
	         end do
              end do

	      if(nlin(ntlblends).eq.0)then
                 do i=nli,m(9)*ntotal4,ntotal4               
	            rt(i)=-rt(i)/2./fill2
	         end do
                 do i=nli,m(10)*ntotal4,ntotal4               
	            rp(i)=-rp(i)/2./fill2
                 end do
              end if

c calculamos la convolucion con la derivada del perfil gaussiano y la psf
	      if(vmac2.gt.0)then
	         do i=1,ntotal4
	            sin0(i)=scal2(i)
	         end do

                 call deconv(scal2,1,ntls,npass,dlamda0s,dlamdas,vmac2)
	         call deconv2(sin0,1,ntls,npass,dlamda0s,dlamdas,vmac2)
              end if

c calculamos la convolucion de las funciones respuesta
              fill=fill2*peso1
	      call sub2err(m( 9),vmac2,fill,dlamda0s,ist,rt,k,dscal) !tempe.
	      call sub2err(m(10),vmac2,fill,dlamda0s,ist,rp,k,dscal) !presi.
	      call sub2err(m(11),vmac2,fill,dlamda0s,ist,rm,k,dscal) !micro
	      call sub2err(m(12),vmac2,fill,dlamda0s,ist,rh,k,dscal) !campo
	      call sub2err(m(13),vmac2,fill,dlamda0s,ist,rv,k,dscal) !veloc.
	      call sub2err(m(14),vmac2,fill,dlamda0s,ist,rg,k,dscal) !inclin
	      call sub2err(m(15),vmac2,fill,dlamda0s,ist,rf,k,dscal) !azimut
              call mult1err(m(16),ntotal4,k,fill,vmac2,sin0,dscal)   !macro

	   end if

c la variable es (1-fill2)/fill2 . Su derivada respecto a fill2 es -1/fill2^2
c como tenemos modf. multiplicativas el factor es fill2*(fill2-1)
   
           if(m(17).eq.1)then
              factor=fill2*(fill2-1.)*peso1
              do j=1,ntotal4
                 k=k+1 
	         dscal(k)=factor*(scal2(j)-scal1(j)) 
	      end do
	   end if

	end if

c *********************mezclamos las dos atmosferas****************************
c y los perfiles de salida seran (teniendo en cuenta la luz difusa)

        do j=1,ntotal4
	   scal(j)=peso1*(scal1(j)*fill1+scal2(j)*fill2)+peso2*stray(j)
        end do

c calculamos la funcion respuesta al peso de la luz difusa
c la variable es peso2/(1.-peso2). La derivada de peso2 respecto 
c a la variable es (1.-peso2)^2 Y teniendo en cuenta pert. mult.
c el factor es peso2*(1.-peso2)
        if(m(18).eq.1)then
           factor=peso1*peso2
c           factor=(1.-peso2)*peso2   !luz difusa com velo
           do j=1,ntotal4
              k=k+1 
	      dscal(k)=factor*(stray(j)-scal1(j)*fill1-scal2(j)*fill2)
c	      dscal(k)=factor*stray(j)  !luz difusa com velo
	   end do
	end if


	if(nlin(ntlblends).eq.0)then
	   scal(nli)=(scal1(nli)-scal2(nli))/2.
        end if
    
	return
	end
c _____________________________________________________________________________
	subroutine sub1err(mi,vmac,fill,rt,k,dscal)

	implicit real*4 (a-h,o-z) 
	include 'PARAMETER' !por kt,kn,kl y kld

	integer ist(4),nlin(kl),npas(kl),nble(kl)
	real*4 dlamda(kld),dlamda0(kl)
	real*4 rt(*),dscal(*)
	common/responde2/ist,ntau,ntl,nlin,npas,nble
        common/ldeo/dlamda,dlamda0
        common/primera2/ntotal,ntotal4,ists

	do i=1,mi
	   npun=ntotal*(i-1)+1
           if(vmac.gt.0)call deconv(rt(npun),1,ntl,npas,dlamda0,dlamda,vmac)
           do j=npun,npun+ntotal-1
              k=k+1 
	      dscal(k)=rt(j)*fill      !t1(l1,l2...lntot4),t2(l1,l2...
	   end do
           do j=ntotal+1,ntotal4
              k=k+1 
	      dscal(k)=0.
	   end do
        end do		   

	return
	end
c _____________________________________________________________________________
	subroutine sub2err(mi,vmac,fill,dlamda0s,ist,rt,k,dscal)

	implicit real*4 (a-h,o-z) 
	include 'PARAMETER' !por kt,kn,kl y kld

	parameter (kl4=4*kl)    !numero maximo de lineas
	parameter (kld4=4*kld)

	integer ist(4),nlins(kl4),npass(kl4)
	real*4 dlamda0s(*),dlamdas(kld4)
	real*4 rt(*),dscal(*)
        common/smalla/ntls,nlins,npass
        common/smalla1/dlamdas
        common/primera2/ntotal,ntotal4,ists

	do i=1,mi
	   npun=ntotal4*(i-1)+1
           if(vmac.gt.0)call deconv(rt(npun),1,ntls,npass,dlamda0s,dlamdas,vmac)
           do j=npun,npun+ntotal4-1
              k=k+1 
	      dscal(k)=rt(j)*fill      !t1(l1,l2...lntot4),t2(l1,l2...
	   end do
        end do		
	return
	end
c _____________________________________________________________________________
	subroutine cero1err(mi,ntotal4,k,dscal)

	real*4 dscal(*)

	do i=1,mi
           do j=1,ntotal4
              k=k+1 
	      dscal(k)=0.
	   end do
        end do		   

	return
	end
c _____________________________________________________________________________
	subroutine mult1err(mi,ntotal4,k,fill,vmac,sin,dscal)
        
	real*4 fill,vmac,sin(*),dscal(*) 

        if(mi.eq.1)then
           prod=fill*vmac
           do j=1,ntotal4
              k=k+1 
	      dscal(k)=sin(j)*prod           
	   end do
	end if

        return
	end 
c _____________________________________________________________________________
