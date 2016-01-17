      subroutine marquardt2(y,sig,ndata,a,mnodos,mfit,
     *    covar,alpha,chisq,alamda,iamplio,beta)

	implicit real*4 (a-h,o-z)

	include 'PARAMETER'   !por kld y mmx==kn que se tomaba 200
	parameter (kld4=4*kld,nodos=18)
c        dimension y(ndata),sig(ndata),a(*),mnodos(nodos),v(mmx,mmx),w(mmx),
c     *  covar(mfit,mfit),alpha(mfit,mfit),atry(mmx),beta(mfit),da(mmx)

        dimension y(ndata),sig(ndata),a(*),mnodos(*),v(mfitmax,mfitmax),
     *  w(mfitmax),covar(mfitmax,mfitmax),alpha(mfitmax,mfitmax),
     *  atry(mfitmax),beta(mfitmax),da(mfitmax)
	real*4 sigreal(kld4)

        integer mnodosold(nodos)

	character*20 control
	common/tol/tol
        common/alamda0/alamda0
	common/allc/all
	common/ochisq/ochisq
	common/sigrealchi/sigreal,chireal,sumsq
	common/canal/icanal
	common/nombrecontrol/control
        common/repeticion/factorrep !factor de lambda para evitar repeticiones
        common/nciclos/nciclos   !del principal
c el common siguiente viene del principal y es para pasarle a marquardt2 
c los indices iniciales y finales de las perturbaciones a gamma y fi aditivas
c        common/ifies/iga1,ifi11,iga2,ifi22 
        common/ifies/ipa1,ipa11,iga1,ifi11,ipa2,ipa22,iga2,ifi22
        common/ieliminofrec/ielimino !(del ppal) para el print de la S/N
        common/primerchi/chisn,snn

	all=alamda

        xpi=3.14159265

c        print*,'marquardt2  38',mnodos(1),mnodos(4),mnodos(4)
	
      if(alamda.lt.0.)then

        alamda=alamda0
c        print*,'marquardt2  43',mnodos(1),mnodos(4),mnodos(4)
        call marqcoef2(y,sig,ndata,a,mnodos,mfit,alpha,beta,chisq)
c        print*,'marquardt2  45',mnodos(1),mnodos(4),mnodos(4)
	if(nciclos.lt.1)return


        ochisq=chisq
        chip=sumsq/float(ndata-ielimino)
	chisn=1./sqrt(chip)
        snn=chireal/float(ndata-ielimino-mfit)
c	print*,'marquardt2(52) ',chireal,sumsq
c        print*,'marquardt2(53) ',chisn,snn,chip,ndata,ielimino,mfit,float(ndata-ielimino-mfit)
	write(*,786)0,chisn,snn  !esta es la correcta
c	open(icanal,file=control,fileopt='eof')
c	write(icanal,786)0,chisn,snn !esta es la correcta
c	close(icanal)
786	format(1x,i3,6x,1pe9.2,1x,e10.3)
      endif
c      print*,'marquardt2  61',mnodos(1),mnodos(4),mnodos(4) 
c ::::::::::::::::::::::::hasta aqui solo para la primera iteracion

c        print*,'marquardt2 mfit alamda=',mfit,alamda
	
        do j=1,mfit
           atry(j)=a(j)
c           print*,'marquardt2 j a beta=',j,atry(j),beta(j),alamda,mnodos(j),mfit,tol

           do k=1,mfit
              covar(j,k)=alpha(j,k)
c              print*,'marquardt2 (1) j k covar=',j,k,covar(j,k)
	   end do
           covar(j,j)=alpha(j,j)*(1.+alamda)
           da(j)=beta(j)
 	end do

	call svdmatriz2(covar,mfit,mnodos,da,tol,v,w)

        do j=1,ipa1                       !nodos en T atm 1
           atry(j)=a(j)*(1.+da(j))
	end do
	do j=ipa1+1,ipa11                 !nodos en Pe atm 1 
           atry(j)=a(ipa11)*(1.+da(j))   !escala con la Pe en el ultimo nodo
	end do
        do j=ipa11+1,iga1                 !nodos en mic,H,Vz atm 1
           atry(j)=a(j)*(1.+da(j))
	end do
        do j=iga1+1,ifi11                 !nodos en gamma y fi atm 1
           atry(j)=a(j)+da(j)
        end do
        do j=ifi11+1,ipa2                 !resto nodos atm1 y nodos T atm2
           atry(j)=a(j)*(1.+da(j))
        end do
        do j=ipa2+1,ipa22                 !nodos en Pe atm 2
           atry(j)=a(ipa22)*(1.+da(j))        !escala con la Pe en el ultimo nodo
        end do
       	do j=ipa22+1,iga2                 !nodos en mic,H,Vz atm 2
           atry(j)=a(j)*(1.+da(j))
	end do
        do j=iga2+1,ifi22                 !nodos en gamma y fi atm 2
           atry(j)=a(j)+da(j)            
	end do
        do j=ifi22+1,mfit                 !resto nodos atm2
           atry(j)=a(j)*(1.+da(j))
	end do

        mfitold=mfit
        do j=1,18
           mnodosold(j)=mnodos(j)
        end do
        ipa1old=ipa1
        ipa2old=ipa2
        ipa11old=ipa11
        ipa22old=ipa22
        iga1old=iga1
        iga2old=iga2
        if11old=ifi11
        if22old=ifi22


        call marqcoef2(y,sig,ndata,atry,mnodos,mfit,covar,da,chisq)

c        print*,'marquardt2 (20) mfit chisq ochisq', mfit, chisq ,ochisq
 
        if(chisq.lt.ochisq)then

c           if(factorrep.ne.0)alamda=0.1*alamda	

            if(factorrep.ne.0)then
                 if(alamda.gt.1.e-4)then
                    alamda=0.1*alamda
                 else	
                    alamda=0.5*alamda
                 end if	
            end if
           ochisq=chisq

           do j=1,mfit
              do k=1,mfit
                 alpha(j,k)=covar(j,k)
	      end do
c              beta1(j)=beta(j)
              beta(j)=da(j)
              a(j)=atry(j)
c              print*,'marquardt2 (3) j atrynew =',j,atry(j)
	      iamplio=1
	   end do

        else

c           alamda=10.*alamda
c                 if(alamda.lt.1.e3)then
c                    alamda=10.*alamda
c                 else	
c                    alamda=2.*alamda
c                 end if	
		
                 if(alamda.le.1.e-3)then
                    alamda=100.*alamda
                 else 
                    if(alamda.lt.1.e3)then
                       alamda=10.*alamda
                    else	
                       alamda=2.*alamda
                    end if
                 end if	

           chisq=ochisq
	   iamplio=0

           mfit=mfitold
           do j=1,18
              mnodos(j)=mnodosold(j)
           end do
           ipa1=ipa1old
           ipa2=ipa2old
           ipa11=ipa11old
           ipa22=ipa22old
           iga1=iga1old
           iga2=iga2old
           if11=ifi11old
           if22=ifi22old

        endif

        return

        end
c _____________________________________________________________________
