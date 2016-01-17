        subroutine marqcoef2(y,sig,ndata,a,m,mfit,alpha,beta,chisq)

	implicit real*4 (a-h,o-z)

	include 'PARAMETER'  !por kn,kl,kld,mfitmax

	parameter (kn16=16*kn+4)         !
        parameter (mmax=16*kn+2,kld4=4*kld,kldn=mmax*kld4)        
        dimension y(ndata),sig(ndata),alpha(mfitmax,mfitmax),beta(*)
        dimension dyda(mfitmax),a(*),ymod(kld4),ymodobs(kld4),dvcal(kldn)
	integer m(*),npos(kld4),numberfrec(4),iratio(4)
	real*4 sigreal(kld4),stokesratio(kld4)
	real*4 dcastigo(mfitmax),ddcastigo(mfitmax)
c	real*4 d2castigo(mfitmax,mfitmax)
c	common/uvesalida/ymodobs
	common/uvesalida/stokesratio
	common/sigrealchi/sigreal,chireal,sumsq
        common/nciclos/nciclos   !del principal
        common/derivchi/dvcal  ! deriv. chi^2 para los errores
        common/primera2/ntotal,ntotal4,ists
        common/posiciones/npos 
        common/numberfrec/numberfrec,iratio     !para marqcoef2
        common/ndata/ndatosobs !para fperfil2 y automatico
        common/ifies/ipa1,ipa11,iga1,ifi11,ipa2,ipa22,iga2,ifi22
	common/nspl_lin/nspl_lin !para seleccionar penalty2 (4,5); penalty1(2 ,3) o no(0,1)
     
        ndatosobs=ndata
        icastigo=0
	if(nspl_lin .eq. 2 .or. nspl_lin .eq. 3)icastigo=1
	if(nspl_lin .eq. 4 .or. nspl_lin .eq. 5)icastigo=2

        chisq=0.
        chireal=0.
        sumsq=0.

        call fperfil2(m,a,ymod,dvcal)

         mfit=0
c	 do i=1,5
c           if(m(i).gt.0)mfit=mfit+m(i)
c	 end do
         if(m(1).gt.0)mfit=mfit+m(1)
         ipa1=mfit       !ipa1 es el indice anterior a la gamma comp. 2 
 	 if(m(2).gt.0)mfit=mfit+m(2)
         ipa11=mfit
	 do i=3,5
           if(m(i).gt.0)mfit=mfit+m(i)
	 end do
         iga1=mfit    !iga1 es el indice anterior a la gamma comp. 1
         if(m(6).gt.0)mfit=mfit+m(6)
         if(m(7).gt.0)mfit=mfit+m(7)
         ifi11=mfit   !ifi11 es el indice ultimo de la fi comp. 1
c	 do i=8,13
c           if(m(i).gt.0)mfit=mfit+m(i)
c	 end do
	 do i=8,9
           if(m(i).gt.0)mfit=mfit+m(i)
	 end do
	 ipa2=mfit         !ipa2 es el indice anterior a la gamma comp. 2
	 if(m(10).gt.0)mfit=mfit+m(10)
	 ipa22=mfit
	 do i=11,13
           if(m(i).gt.0)mfit=mfit+m(i)
	 end do	
         iga2=mfit    !iga2 es el indice anterior a la gamma comp. 2
         if(m(14).gt.0)mfit=mfit+m(14)
         if(m(15).gt.0)mfit=mfit+m(15)
         ifi22=mfit   !ifi22 es el indice ultimo de la fi comp. 1
	 do i=16,18
           if(m(i).gt.0)mfit=mfit+m(i)
	 end do
	 	

      do 12 j=1,mfit
         do 11 k=1,j
            alpha(j,k)=0.
11       continue
         beta(j)=0.
12    continue

c escribimos la FR
c =======================================================================
      if(nciclos .eq. -1)then
	i=0
        iratio(1)=0. !no permitimos escribir la FR de I/I 
        do ii=1,4
           do iii=1,numberfrec(ii)
              i=i+1
	      do j=1,mfit
	         jj=(j-1)*ntotal4+i
	         jji=(j-1)*ntotal4+iii
	         if(iratio(ii).eq.1)dvcal(jj)=
     &             ( dvcal(jj) - dvcal(jji)*ymod(i)/ymod(iii) ) /ymod(iii)
	       end do
	    end do
	 end do
         call escribeFR(m,ndata,ntotal4,npos,dvcal)
      end if
c =======================================================================
c calculamos el castigo (y sus derivadas) por picos en las estratificaciones
      if(icastigo.eq.1) call penalty(mfit,m,a,castigo,dcastigo,ddcastigo)
      if(icastigo.eq.2) call penalty2(mfit,m,a,castigo,dcastigo)
      wcast=1.  !pesocastigo

      do i=1,ndata
         ymodobs(i)=ymod(npos(i))
c	 print*,'marqcoef2(94) ',i,npos(i),ymod(npos(i)),ymodobs(i)
      end do

	i=0
      do ii=1,4
c        print*,'marqcoef2 (99) numberfrec(',ii,')=',numberfrec(ii)
        do iii=1,numberfrec(ii)
        i=i+1
c      do 15 i=1,ndata

         sss=sig(i)
	 do j=1,mfit
	    jj=(j-1)*ntotal4+npos(i)
	    jji=(j-1)*ntotal4+npos(iii)
	    if(iratio(ii).eq.1)then
               dyda(j)=( dvcal(jj)-dvcal(jji)*ymodobs(i)/ymodobs(iii) )/ymodobs(iii)
            else
               dyda(j)=dvcal(jj)
            end if
	 end do

         sig2i=1./sss/sss
         sig2reali=1./(sigreal(i)*sigreal(i))  !sin comentarios
         if(iratio(ii).ne.1)stokesratio(i)=ymodobs(i)
c         print*,'marqcoef2 ',iratio(ii),ii,iii,i,ymodobs(i),ymodobs(iii),y(i),stokesratio(i)
         if(iratio(ii).eq.1)stokesratio(i)=ymodobs(i)/ymodobs(iii)
         dy=y(i)-stokesratio(i)

         do 14 j=1,mfit
            wt=dyda(j)*sig2i
            do 13 k=1,j
               alpha(j,k)=alpha(j,k)+wt*dyda(k)
13          continue

c            print*,'marqcoef2 i j dy dyda',i,j,dyda(j),wt,dy
	    beta(j)=beta(j)+dy*wt

14       continue

         chisq=chisq+dy*dy*sig2i
         chireal=chireal+dy*dy*sig2reali       !sin comentarios
c	  print*,'marqcoef2 (134)',dy,sig2reali,dy*dy*sig2reali,ii,iii

         if(sss.lt.1.e14)sumsq=sumsq+dy*dy

	end do
	end do
c	print*,'marqcoef2(139) chireal=',chireal,sumsq

c 15    continue

c------------------------------castigo 1------------------------------------------------
        if(icastigo .eq.1)then
           castigo=1.+castigo*wcast   !multiplicativamente
c	   print*,'castigo',castigo

           do j=1,mfit
	      alpha(j,1)=alpha(j,1)*castigo+beta(j)*dcastigo(1)*wcast+
     &                 beta(1)*dcastigo(j)*wcast+chisq*ddcastigo(j)*wcast
	      do k=2,j
	         alpha(j,k)=alpha(j,k)*castigo+beta(j)*dcastigo(k)*wcast+
     &                                      beta(k)*dcastigo(j)*wcast
              end do
           end do
	   do j=1,mfit
              beta(j)=beta(j)*castigo+chisq*dcastigo(j)*wcast
	   end do
	   chisq=chisq*castigo
	   chireal=chireal*castigo
        end if
c------------------------------castigo 2------------------------------------------------
        if(icastigo .eq.2)then
           castigo=1.+castigo*wcast   !multiplicativamente
c	   print*,'castigo',castigo

           do j=1,mfit
	      do k=1,j
c	         d2castigo(j,k)=dcastigo(j)*dcastigo(k)
	         alpha(j,k)=alpha(j,k)*castigo+beta(j)*dcastigo(k)*wcast+
     &                      beta(k)*dcastigo(j)*wcast+chisq*dcastigo(j)*dcastigo(k)*wcast
              end do
           end do
	   do j=1,mfit
              beta(j)=beta(j)*castigo+chisq*dcastigo(j)*wcast
	   end do
	   chisq=chisq*castigo
	   chireal=chireal*castigo
        end if
c----------------------------------------------------------------------------------------

        do 17 j=2,mfit
           do 16 k=1,j-1
              alpha(k,j)=alpha(j,k)
16         continue
17      continue

      return
      end
