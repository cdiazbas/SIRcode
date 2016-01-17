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
	subroutine fperfil2(m,atmosr,scal,dscal)

	implicit real*4 (a-h,o-z)
 
	include 'PARAMETER'   !por kt,kn,kl,kld
	parameter (kt8=8*kt+2,kt16=16*kt+5)         !
	parameter (kl4=4*kl)    !numero maximo de lineas
	parameter (kld4=4*kld)
	parameter (kldt=kld*kn,kldt4=4*kldt)

c para la malla
	integer nlin(kl),npas(kl),ist(4),m(*),mdata(18) !,mnod1(8),mnod2(8)
	integer mnod1tot(8),mnod2tot(8),npos(kld4)
	integer mprimera(18)
	integer nlins(kl4),npass(kl4),nble(kl),ifiltro
	real*4 dlamda0(kl),dlamda0s(kl4),dlamda(kld),dlamdas(kld4)
c para la atmosfera
	real*4 atmosr(*),atmos(kt16),atmostry(kt16)
	real*4 atmos1(kt8),atmos2(kt8),vof(kt),gam1(kt),fi1(kt)
c para los perfiles y f. respuesta
	real*4 scal(kld4),scal1(kld4),scal2(kld4),dscal(*),sin01(kld4),sin02(kld4),stray(kld4)
        real*4 stok(kld4),difer(kld4),sig(kld4)
	real*4 rt1(kldt4),rh1(kldt4),rv1(kldt4),rg1(kldt4),rf1(kldt4)
	real*4 rp1(kldt4),rm1(kldt4)
	real*4 rt2(kldt4),rh2(kldt4),rv2(kldt4),rg2(kldt4),rf2(kldt4)
	real*4 rp2(kldt4),rm2(kldt4)
        integer icalerr
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
        common/iautomatico/iautomatico
        common/observaciones/stok,sig
        common/posiciones/npos 
        common/mprimera/mprimera ! para guardar m inicial para iteaciones posteriores 
        common/ndata/ndata !para fperfil2 y automatico
        common/offset/voffset  !para respuestas
        common/iprimeravez/iprimeravez
        common/calerr/icalerr !si calculo errores=1 else =0
	common/ifiltro/ifiltro

	epsilon=2.e-5 !precision para comparar reales

	nohaycampo=nohaycampo1+nohaycampo2

c	print*,'fperfil2 nohaycampo=',nohaycampo1,nohaycampo2,nohaycampo
	if(iprimeravez.eq.0)then

           do i=1,18
	      mprimera(i)=m(i) ! guardamos los nodos de entrada
           end do 
           
c           print*,'fperfil 72',m(4),m(6)
           
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

c duplicado del calculo en comprime2
           vmin=7.e5
           if(m(13).ne.0.or.m(5).ne.0)vmin=1.e20

           do i=1,ntau              
	      if(atmos(2+13*ntau+i).lt.vmin.and.m(13).ne.0)
     &                              vmin=atmos(2+13*ntau+i)
              if(atmos(5*ntau+i).lt.vmin.and.m(5).ne.0)vmin=atmos(5*ntau+i)
              voffset=vmin-7.e5    !cm/s
           end do

	end if
        iprimeravez=iprimeravez+1


        nli=0
	do i=1,ntl
           nli=nli+npas(i)
	end do
        
c copiamos la atmosfera antigua en atmostry; "amp2" escribira en 
c atmostry la nueva atmosfera
	do i=1,16*ntau+5
	   atmostry(i)=atmos(i)
	end do

c	print*,'fperfil2 122',m(1),m(4),m(6)
	if(nciclos.ne.0.and.iprimeravez.ne.1)call amp2(ntau,m,atmostry,atmosr)
c 	if(nciclos.ne.0)call amp2(ntau,m,atmostry,atmosr)
c	print*,'fperfil2 125',m(1),m(4),m(6)
	
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
           sin01(i)=0.
           sin02(i)=0.
	end do
        
        if(icalerr.eq.1)then
	   do i=1,8
              mnod1tot(i)=m(i)
c              print*,'fperfil2 m err',m(i)
              mnod2tot(i)=m(i+8)
           end do
        else
  	   do i=1,8
c             mnod1(i)=m(i)
c             mnod2(i)=m(i+8)
              if(mprimera(i).gt.1.and.iautomatico.eq.1)then
                 mnod1tot(i)=ntau
              else
                 mnod1tot(i)=mprimera(i)
              end if
              if(mprimera(i+8).gt.1.and.iautomatico.eq.1)then
                 mnod2tot(i)=ntau
              else
                 mnod2tot(i)=mprimera(i+8)
              end if

	   end do
        end if

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
	      call blendscon2(atmos1,scal1,rt1,rp1,rv1,rm1,mnod1tot)

c	      print*,'mnod1tot=',mnod1tot
              
c	      do llk=0,kldt4
c                 print*,'i,rt1=',llk,rt1(llk)
c		  if(rt1(kldt4).ne.0)print*,'encontre uno no negativo=',llk,rt1(llk)
c	      enddo

	      if(vmac1.gt.0. .or. ifiltro.eq.1)then
	         do j=1,ntotal
	            sin01(j)=scal1(j)
	         end do
		 do j=ntotal+1,ndata
                   scal1(j)=0.
                   sin01(j)=0.
	         end do
	         if(vmac1 .lt. 1.e-7)vmac1=1.e-7
	         call deconv(scal1,1,ntl,npas,dlamda0,dlamda,vmac1)
                 call deconv2(sin01,1,ntl,npas,dlamda0,dlamda,vmac1)
              end if 
   	   else  !o sea si hay campo1 
       	      call blends2(atmos1,scal1,rt1,rp1,rh1,rv1,rg1,rf1,rm1,mnod1tot)


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

c calculamos la convolucion con el perfil gaussiano y la psf
c calculamos la convolucion con la derivada del perfil gaussiano y la psf
	      if(vmac1.gt.0 .or. ifiltro.eq.1)then
	         do j=1,ntotal4
	            sin01(j)=scal1(j)
	         end do
	         if(vmac1 .lt. 1.e-7)vmac1=1.e-7
	         call deconv(scal1,1,ntls,npass,dlamda0s,dlamdas,vmac1)
                 call deconv2(sin01,1,ntls,npass,dlamda0s,dlamdas,vmac1)
              end if 
           end if
        end if
c ************************* atmosfera 2 ***************************************
	if(abs(fill2).gt.epsilon)then
	   if(nohaycampo2.eq.0.and.ist(1).ne.0)then
              call blendscon2(atmos2,scal2,rt2,rp2,rv2,rm2,mnod2tot)
	      if(vmac2.gt.0 .or. ifiltro.eq.1)then
	         do j=1,ntotal
	            sin02(j)=scal2(j)
	         end do
		 do j=ntotal+1,ndata
                   scal2(j)=0.
                   sin02(j)=0.
	         end do
	         if(vmac2 .lt. 1.e-7)vmac2=1.e-7
	         call deconv(scal2,1,ntl,npas,dlamda0,dlamda,vmac2)
	         call deconv2(sin02,1,ntl,npas,dlamda0,dlamda,vmac2)
              end if
   	   else  !o sea si hay campo2
       	      call blends2(atmos2,scal2,rt2,rp2,rh2,rv2,rg2,rf2,rm2,mnod2tot)

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


c calculamos la convolucion con la derivada del perfil gaussiano y la psf
	      if(vmac2.gt.0 .or. ifiltro.eq.1)then
	         do i=1,ntotal4
	            sin02(i)=scal2(i)
	         end do
                 if(vmac2 .lt. 1.e-7)vmac2=1.e-7
                 call deconv(scal2,1,ntls,npass,dlamda0s,dlamdas,vmac2)
	         call deconv2(sin02,1,ntls,npass,dlamda0s,dlamdas,vmac2)
              end if
           end if
	end if 
c *********************mezclamos las dos atmosferas****************************
c y los perfiles de salida seran (teniendo en cuenta la luz difusa)

        do j=1,ntotal4
	   scal(j)=peso1*(scal1(j)*fill1+scal2(j)*fill2)+peso2*stray(j)
c	   print*,'fperfil2(276) ',j,scal(j)
        end do
	if(nlin(ntlblends).eq.0)then
	   scal(nli)=(scal1(nli)-scal2(nli))/2.
        end if
c        print*,'fperfil2(282) ndata ntotal4 nli',ndata,ntotal4,nli,npos(nli)+1,' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
c        do j=npos(nli)+1,ndata
c           scal(npos(j))=0.
c	end do
c        do j=ntotal4,ndata
c           print*,'fperfil2(286) ',scal1(j),scal2(j),scal(j)
c        end do
        do j=1,ndata
c           print*,'fperfil2(284) ',j,npos(j),stok(j),scal(npos(j)),sig(j)
           difer(j)=(stok(j)-scal(npos(j)))/sig(j)/sig(j)
        end do

c *******************reevaluamos el numero de nodos y escribimos el vector derivada

	if(abs(fill1).gt.epsilon)then
	   if(nohaycampo1.eq.0.and.ist(1).ne.0)then

	      if(nlin(ntlblends).eq.0)then
                 do i=nli,mnod1tot(1)*ntotal,ntotal               
	            rt1(i)=rt1(i)/2./fill1
	         end do
                 do i=nli,mnod1tot(2)*ntotal,ntotal               
	            rp1(i)=rp1(i)/2./fill1
                 end do
              end if

            
              fill=fill1*peso1
	      call sub1(mnod1tot(1),difer,npos,vmac1,fill,rt1,k,
     &                  dscal,atmos1(ntau+1),mprimera(1)) !tempe.
	      call sub1p(mnod1tot(2),difer,npos,vmac1,fill,rp1,k,             !ojo llamo a sub2p
     &                  dscal,atmos1(2*ntau+1),mprimera(2)) !presi. 
	      call sub1(mnod1tot(3),difer,npos,vmac1,fill,rm1,k,
     &                  dscal,atmos1(3*ntau+1),mprimera(3)) !micro.
              call cero1(mnod1tot(4),ntotal4,k,dscal)           !campo 

              if(mnod1tot(5).gt.0)then
                 do i=1,ntau
                    vof(i)=atmos1(5*ntau+i)-voffset
                 end do
              end if

	      call sub1(mnod1tot(5),difer,npos,vmac1,fill,rv1,k,
     &                  dscal,vof,mprimera(5)) !veloc.
              call cero1(mnod1tot(6),ntotal4,k,dscal)           !inclinacion 
              call cero1(mnod1tot(7),ntotal4,k,dscal)           !azimuth 
              call mult1(mnod1tot(8),ntotal4,k,fill,vmac1,sin01,dscal) !macro
                  do i=1,8
                     m(i)=mnod1tot(i)   !nuevo numero de nodos
                  end do

   	   else  !o sea si hay campo1 

	      if(nlin(ntlblends).eq.0)then
                 do i=nli,mnod1tot(1)*ntotal4,ntotal4               
	            rt1(i)=rt1(i)/2./fill1
	         end do
                 do i=nli,mnod1tot(2)*ntotal4,ntotal4               
	            rp1(i)=rp1(i)/2./fill1
                 end do
              end if
        
c calculamos la convolucion de las funciones respuesta
              fill=fill1*peso1


	      call sub2(mnod1tot(1),difer,npos,vmac1,fill,dlamda0s,ist,rt1,k,
     &                  dscal,atmos1(ntau+1),mprimera(1)) !tempe.
	      call sub2p(mnod1tot(2),difer,npos,vmac1,fill,dlamda0s,ist,rp1,k,  !ojo llamo a sub2p
     &                  dscal,atmos1(2*ntau+1),mprimera(2)) !presi.
	      call sub2(mnod1tot(3),difer,npos,vmac1,fill,dlamda0s,ist,rm1,k,
     &                  dscal,atmos1(3*ntau+1),mprimera(3)) !micro.
	      call sub2(mnod1tot(4),difer,npos,vmac1,fill,dlamda0s,ist,rh1,k,
     &                  dscal,atmos1(4*ntau+1),mprimera(4)) !campo



              if(mnod1tot(5).gt.0)then
                 do i=1,ntau
                    vof(i)=atmos1(5*ntau+i)-voffset
                 end do
              end if
              
	      call sub2(mnod1tot(5),difer,npos,vmac1,fill,dlamda0s,ist,rv1,k,
     &                  dscal,vof,mprimera(5)) !veloc.


              if(mnod1tot(6).gt.0)then
                 do i=1,ntau
                    gam1(i)=1.
                 end do
              end if
c                   print*,'fperfil 379',mnod1tot(6),mprimera(6)
	      call sub2(mnod1tot(6),difer,npos,vmac1,fill,dlamda0s,ist,rg1,k,
     &                  dscal,gam1,mprimera(6)) !incli.
     
c	      print*,'fperfil 379',mnod1tot(6),mprimera(6)
              if(mnod1tot(7).gt.0)then
                 do i=1,ntau
                    fi1(i)=1.
                 end do
              end if

	      call sub2(mnod1tot(7),difer,npos,vmac1,fill,dlamda0s,ist,rf1,k,
     &                  dscal,fi1,mprimera(7)) !azimuth.
              call mult1(mnod1tot(8),ntotal4,k,fill,vmac1,sin01,dscal)   !macro

 
c              open(15,file='escfperfil2')
c              do i=1,k
c                 write(15,*),dscal(i)
c              end do
c              close(15)
c              print*,'fperfil2 acabe de escribir fichero'

                  do i=1,8
                     m(i)=mnod1tot(i)   !nuevo numero de nodos
                  end do
           end if
	end if

c solo en el caso de que fill2.ne.0 se hace lo siguiente
c ************************* atmosfera 2 ***************************************

	if(abs(fill2).gt.epsilon)then
	   if(nohaycampo2.eq.0.and.ist(1).ne.0)then

	      if(nlin(ntlblends).eq.0)then
                 do i=nli,mnod2tot(1)*ntotal,ntotal               
	            rt2(i)=-rt2(i)/2./fill2
	         end do
                 do i=nli,mnod2tot(2)*ntotal,ntotal               
	            rp2(i)=-rp2(i)/2./fill2
                 end do
              end if
              fill=fill2*peso1

	      call sub1(mnod2tot(1),difer,npos,vmac2,fill,rt2,k,
     &                  dscal,atmos2(ntau+1),mprimera(9)) !tempe.
	      call sub1p(mnod2tot(2),difer,npos,vmac2,fill,rp2,k,               !ojo llamo a sub2p
     &                  dscal,atmos2(2*ntau+1),mprimera(10)) !presi.
	      call sub1(mnod2tot(3),difer,npos,vmac2,fill,rm2,k,
     &                  dscal,atmos2(3*ntau+1),mprimera(11)) !micro.
              call cero1(mnod2tot(4),ntotal4,k,dscal)           !campo 

              if(mnod2tot(5).gt.0)then
                 do i=1,ntau
                    vof(i)=atmos2(5*ntau+i)-voffset
                 end do
              end if

	      call sub1(mnod2tot(5),difer,npos,vmac2,fill,rv2,k,
     &                  dscal,vof,mprimera(13)) !veloc.


              call cero1(mnod2tot(6),ntotal4,k,dscal)            !inclinacion  
              call cero1(mnod2tot(7),ntotal4,k,dscal)            !azimuth 
              call mult1(mnod2tot(8),ntotal4,k,fill,vmac2,sin02,dscal) !macro
                  do i=1,8
                     m(i+8)=mnod2tot(i)   !nuevo numero de nodos
                  end do

	   else  !o sea si hay campo2 


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
                 do i=nli,mnod2tot(1)*ntotal4,ntotal4               
	            rt2(i)=-rt2(i)/2./fill2
	         end do
                 do i=nli,mnod2tot(2)*ntotal4,ntotal4               
	            rp2(i)=-rp2(i)/2./fill2
                 end do
              end if


c calculamos la convolucion de las funciones respuesta
              fill=fill2*peso1

	      call sub2(mnod2tot(1),difer,npos,vmac2,fill,dlamda0s,ist,rt2,k,
     &                  dscal,atmos2(ntau+1),mprimera(9)) !tempe.
	      call sub2p(mnod2tot(2),difer,npos,vmac2,fill,dlamda0s,ist,rp2,k, !ojo llamo a sub2p
     &                  dscal,atmos2(2*ntau+1),mprimera(10)) !presi.
	      call sub2(mnod2tot(3),difer,npos,vmac2,fill,dlamda0s,ist,rm2,k,
     &                  dscal,atmos2(3*ntau+1),mprimera(11)) !micro.
	      call sub2(mnod2tot(4),difer,npos,vmac2,fill,dlamda0s,ist,rh2,k,
     &                  dscal,atmos2(4*ntau+1),mprimera(12)) !campo

              if(mnod2tot(5).gt.0)then
                 do i=1,ntau
                    vof(i)=atmos2(5*ntau+i)-voffset
                 end do
              end if

	      call sub2(mnod2tot(5),difer,npos,vmac2,fill,dlamda0s,ist,rv2,k,
     &                  dscal,vof,mprimera(13)) !veloc.

              if(mnod2tot(6).gt.0)then
                 do i=1,ntau
                    gam1(i)=1.
                 end do
              end if
	      call sub2(mnod2tot(6),difer,npos,vmac2,fill,dlamda0s,ist,rg2,k,
     &                  dscal,gam1,mprimera(14)) !incli.

              if(mnod2tot(7).gt.0)then
                 do i=1,ntau
                    fi1(i)=1.
                 end do
              end if
	      call sub2(mnod2tot(7),difer,npos,vmac2,fill,dlamda0s,ist,rf2,k,
     &                  dscal,fi1,mprimera(15)) !azimuth.
              call mult1(mnod2tot(8),ntotal4,k,fill,vmac2,sin02,dscal)   !macro


                  do i=1,8
                     m(i+8)=mnod2tot(i)   !nuevo numero de nodos
                  end do

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

       	if(nciclos.ne.0.)
     &         call comprime2(ntau,m,atmostry,atmosr)

    
	return
	end
c _____________________________________________________________________________
	subroutine sub1(mi,difer,npos,vmac,fill,rt,k,dscal,t,mp)

	implicit real*4 (a-h,o-z) 
	include 'PARAMETER' !por kt,kn,kl y kld
	parameter (kl4=4*kl)    !numero maximo de lineas
	parameter (kld4=4*kld)

	integer ist(4),nlin(kl),npas(kl),nble(kl),npos(*),ifiltro
	real*4 dlamda(kld),dlamda0(kl)
	real*4 rt(*),dscal(*),difer(*),t(*)
	common/responde2/ist,ntau,ntl,nlin,npas,nble
        common/ldeo/dlamda,dlamda0
        common/primera2/ntotal,ntotal4,ists
        common/iautomatico/iautomatico
	common/ifiltro/ifiltro
c !!!!!!!!!esto seria lo correcto pero es mas lento y no se nota la mejora
c	do i=1,mi
c	   npun=ntotal*(i-1)+1
c           if(vmac.gt.0)call deconv(rt(npun),1,ntl,npas,dlamda0,dlamda,vmac)
c        end do

        if(iautomatico.eq.1.and.mp.gt.1.and.mi.gt.1)call automatico(mi,mp,ntotal,difer,npos,rt,t) 

c TRABAJO con las FR sin convolucionar en todo tau !!!!!!!!!!
	do i=1,mi
	   npun=ntotal*(i-1)+1
           if(vmac.gt.0 .or. ifiltro .eq. 1)then
              if(vmac.lt.1.e-7)vmac=1.e-7 
              call deconv(rt(npun),1,ntl,npas,dlamda0,dlamda,vmac)
           end if
        end do 
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	do i=1,mi
	   npun=ntotal*(i-1)+1
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
C igual que sub1 pero llamando a automaticop (escala con el ultimo nodo de la presion)
	subroutine sub1p(mi,difer,npos,vmac,fill,rt,k,dscal,t,mp)

	implicit real*4 (a-h,o-z) 
	include 'PARAMETER' !por kt,kn,kl y kld
	parameter (kl4=4*kl)    !numero maximo de lineas
	parameter (kld4=4*kld)

	integer ist(4),nlin(kl),npas(kl),nble(kl),npos(*),ifiltro
	real*4 dlamda(kld),dlamda0(kl)
	real*4 rt(*),dscal(*),difer(*),t(*)
	common/responde2/ist,ntau,ntl,nlin,npas,nble
        common/ldeo/dlamda,dlamda0
        common/primera2/ntotal,ntotal4,ists
        common/iautomatico/iautomatico
	common/ifiltro/ifiltro
c !!!!!!!!!esto seria lo correcto pero es mas lento y no se nota la mejora
c	do i=1,mi
c	   npun=ntotal*(i-1)+1
c           if(vmac.gt.0)call deconv(rt(npun),1,ntl,npas,dlamda0,dlamda,vmac)
c        end do

        if(iautomatico.eq.1.and.mp.gt.1.and.mi.gt.1)call automaticop(mi,mp,ntotal,difer,npos,rt,t) 

c TRABAJO con las FR sin convolucionar en todo tau !!!!!!!!!!
	do i=1,mi
	   npun=ntotal*(i-1)+1
           if(vmac.gt.0 .or. ifiltro .eq. 1)then
              if(vmac.lt.1.e-7)vmac=1.e-7 
              call deconv(rt(npun),1,ntl,npas,dlamda0,dlamda,vmac)
           end if
        end do 
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	do i=1,mi
	   npun=ntotal*(i-1)+1
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
	subroutine sub2(mi,difer,npos,vmac,fill,dlamda0s,ist,rt,k,dscal,t,mp)

	implicit real*4 (a-h,o-z) 
	include 'PARAMETER' !por kt,kn,kl y kld

	parameter (kl4=4*kl)    !numero maximo de lineas
	parameter (kld4=4*kld)

	integer ist(4),nlins(kl4),npass(kl4),npos(*),ifiltro
	real*4 dlamda0s(*),dlamdas(kld4)
	real*4 rt(*),dscal(*),difer(*),t(*)
        common/smalla/ntls,nlins,npass
        common/smalla1/dlamdas
        common/primera2/ntotal,ntotal4,ists
        common/iautomatico/iautomatico
	common/ifiltro/ifiltro
c !!!!!!!!!esto seria lo correcto pero es mas lento y no se nota la mejora
c	do i=1,mi
c	   npun=ntotal4*(i-1)+1
c          if(vmac.gt.0)call deconv(rt(npun),1,ntls,npass,dlamda0s,dlamdas,vmac)
c        end do 

         if(iautomatico.eq.1.and.mp.gt.1.and.mi.gt.1)call automatico(mi,mp,ntotal4,difer,npos,rt,t) 

c TRABAJO con las FR sin convolucionar en todo tau !!!!!!!!!!
	do i=1,mi
	   npun=ntotal4*(i-1)+1
	   if(vmac.gt.0 .or. ifiltro .eq. 1)then
              if(vmac.lt.1.e-7)vmac=1.e-7 
              call deconv(rt(npun),1,ntls,npass,dlamda0s,dlamdas,vmac)
           end if   
        end do 
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        do i=1,mi
	   npun=ntotal4*(i-1)+1
           do j=npun,npun+ntotal4-1
              k=k+1 
	      dscal(k)=rt(j)*fill      !t1(l1,l2...lntot4),t2(l1,l2...
	   end do
        end do	
	return
	end
c___________________________________________________________________________________________

	subroutine sub2p(mi,difer,npos,vmac,fill,dlamda0s,ist,rt,k,dscal,t,mp)

	implicit real*4 (a-h,o-z) 
	include 'PARAMETER' !por kt,kn,kl y kld

	parameter (kl4=4*kl)    !numero maximo de lineas
	parameter (kld4=4*kld)

	integer ist(4),nlins(kl4),npass(kl4),npos(*),ifiltro
	real*4 dlamda0s(*),dlamdas(kld4)
	real*4 rt(*),dscal(*),difer(*),t(*)
        common/smalla/ntls,nlins,npass
        common/smalla1/dlamdas
        common/primera2/ntotal,ntotal4,ists
        common/iautomatico/iautomatico
	common/ifiltro/ifiltro
c !!!!!!!!!esto seria lo correcto pero es mas lento y no se nota la mejora
c	do i=1,mi
c	   npun=ntotal4*(i-1)+1
c          if(vmac.gt.0)call deconv(rt(npun),1,ntls,npass,dlamda0s,dlamdas,vmac)
c        end do 

         if(iautomatico.eq.1.and.mp.gt.1.and.mi.gt.1)call automaticop(mi,mp,ntotal4,difer,npos,rt,t) 

c TRABAJO con las FR sin convolucionar en todo tau !!!!!!!!!!
	do i=1,mi
	   npun=ntotal4*(i-1)+1
	   if(vmac.gt.0 .or. ifiltro .eq. 1)then
              if(vmac.lt.1.e-7)vmac=1.e-7 
              call deconv(rt(npun),1,ntls,npass,dlamda0s,dlamdas,vmac)
           end if   
        end do 
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        do i=1,mi
	   npun=ntotal4*(i-1)+1
           do j=npun,npun+ntotal4-1
              k=k+1 
	      dscal(k)=rt(j)*fill      !t1(l1,l2...lntot4),t2(l1,l2...
	   end do
        end do	
	return
	end
	
	
c _____________________________________________________________________________
	subroutine cero1(mi,ntotal4,k,dscal)

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
	subroutine mult1(mi,ntotal4,k,fill,vmac,sinp,dscal)
        
	real*4 fill,vmac,sinp(*),dscal(*) 

        if(mi.eq.1)then
           prod=fill*vmac
           do j=1,ntotal4
              k=k+1 
	      dscal(k)=sinp(j)*prod           
	   end do
	end if

        return
	end 
c _____________________________________________________________________________
