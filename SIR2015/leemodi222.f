c leemodi222
c iesc=0 :lee dos atmosferas modelo contando puntos en tau de la 1
c iesc=1 :escribe dos atmosferas modelo
c ntau : numero de puntos en tau
c tau  : log10(tau)
c ierror=0 O.K. ierror=1 el modelo 1 no existe, ierror=2 el modelo 2 no existe
c ierror=3 no existen los modelos 1 ni 2

	subroutine leemodi222(iesc,model1,model2,atmos,ntau,ierror,
     &                        z1,pg1,ro1,z2,pg2,ro2)

	include 'PARAMETER'  !por kt
	real*4 atmos(*),a(8)
        real*4 pg1(*),z1(*),ro1(*)
        real*4 pg2(*),z2(*),ro2(*)
        character model1*(*),model2*(*)
	character*100 control
	character*80 men1
	common/canal/icanal
	common/nombrecontrol/control
	common/nohaycampo/nohaycampo1,nohaycampo2
        common/nciclos/nciclos   !para marquardt, leeuve3, leemodi2

	epsilon=2.e-5

	ican=52
	if(iesc.eq.0)then

	

c contamos las lineas del modelo 1
           if(ierror.ne.1)then  
	     open(ican,file=model1,status='old',err=991)
	     ntau=0
	     read(ican,*,err=999)a(1),fill1,peso1
	     do while(ntau.lt.kt)
	        ntau=ntau+1
	        read(ican,*,end=221,err=999)a
	     end do
221	     ntau=ntau-1
	     close(ican)
c            ahora leemos el modelo 1
	     open(ican,file=model1)
	     read(ican,*,err=999)atmos(8*ntau+1),atmos(8*ntau+2),pesostray
	     do i=1,ntau
	        read(ican,*,err=999)(atmos(i+j*ntau),j=0,7),z1(i),pg1(i),ro1(i)
	     end do
	     close(ican)
           end if


c contamos las lineas del modelo 2
           if(ierror.ne.2)then  
	     open(ican,file=model2,status='old',err=992)
	     ntau2=0
	     read(ican,*,err=999)a(1),fill2
	     do while(ntau2.lt.kt)
	       ntau2=ntau2+1
	       read(ican,*,end=222,err=999)a,z2(i),pg2(i),ro2(i)
	     end do
222	     ntau2=ntau2-1
	     close(ican)
	
	     go to 993    !el fichero del modelo 2 no existe

992          ierror=2
             if(abs(fill1-1.0).ne.0)then 
	        men1='     WARNING: The FILLING FACTOR is being changed to 1.0'
	        print*,men1
c	        open(icanal,file=control,fileopt='eof')
c                write(icanal,*) men1
c	        close(icanal)
	     endif
             fill1=1.
             fill2=0.
             ntau2=ntau
	     do i=1,ntau
                do j=0,2 
	           atmos(i+(8+j)*ntau+2)=atmos(i+j*ntau)
                end do
                do j=3,7 
	           atmos(i+(8+j)*ntau+2)=0.
                end do
	     end do
	     nohaycampo2=0
             nohaycampo1=0
             do i=1,ntau
                if(abs(atmos(4*ntau+i)).gt.1.)nohaycampo1=1
             end do 

	     if(pesostray.lt.0..or.pesostray.gt.100.) then
	      men1='STOP: In model 1, the stray light factor is outside the interval [0,100]'
	      call mensaje(1,men1,men1,men1)
             endif

             atmos(16*ntau+5)=pesostray  !% de luz difusa
	
c            Comprobamos que el modelo esta equiespaciado y con tau decreciente
	     paso=atmos(2)-atmos(1)
	     if(paso.ge.0)then
		men1='STOP: The model supplied is NOT ordered with DECREASING log(tau).'
   	        call mensaje(1,men1,men1,men1)
	     end if

             do i=2,ntau-1
	      paso1=atmos(i+1)-atmos(i)
               if(abs(paso1-paso).gt.epsilon)then
		 if(nciclos.gt.0)then
		    men1='STOP: For inversion, an EQUALLY SPACED grid is required to discretize the model.'
		    call mensaje(1,men1,men1,men1)
		  endif
              end if
   	     end do


             return   

993	     if(ntau2.ne.ntau)then
  	       men1='STOP: The initial models are discretized in DIFFERENT SPATIAL GRIDS.'             
	       call mensaje(1,men1,men1,men1)
             end if
           end if

c comprobamos que fill1 esta en el intervalo [0,1]
           if(fill1.lt.0.or.fill1.gt.1.)then
	      men1='STOP: For model 1, the filling factor is outside the interval [0,1] '
	      call mensaje(1,men1,men1,men1)
           end if
	   if(abs(1-(fill1+fill2)).gt.epsilon)then
              fill2=1.-fill1
              print*,'     WARNING: For model 2, the filling factor is taken to be ',fill2
	   end if

c ahora leemos los dos modelos
	   open(ican,file=model1)
	   read(ican,*,err=999)atmos(8*ntau+1),atmos(8*ntau+2),peso1
	   do i=1,ntau
	      read(ican,*,err=999)(atmos(i+j*ntau),j=0,7),z1(i),pg1(i),ro1(i)
	   end do
	   close(ican)
	   open(ican,file=model2)
	   read(ican,*,err=999)atmos(16*ntau+3),atmos(16*ntau+4),peso2
	   do i=1,ntau
	      read(ican,*,err=999)(atmos(i+j*ntau+2),j=8,15),z2(i),pg2(i),ro2(i)
	   end do
	   close(ican)

	   if(peso1.lt.0..or.peso1.gt.100.) then
	      men1='STOP: In model 1, the stray light factor is outside the interval [0,100] '
	      call mensaje(1,men1,men1,men1)
           endif

   	   if(peso1.ne.peso2)then
              print*,'     WARNING: The stray light factor is DIFFERENT in models 1 and 2.'
              print*,'              Its value in model 1 is being adopted:',peso1  
c	      open(icanal,file=control,fileopt='eof')
c              write(icanal,*)'     WARNING: The stray light factor is DIFFERENT in models 1 and 2.'
c              write(icanal,*)'              Its value in model 1 is being adopted:',peso1  
c	      close(icanal)
	   endif
           atmos(16*ntau+5)=peso1  !% de luz difusa
	   atmos(8*ntau+2)=fill1

c comprobamos que estan equiespaciados y con tau decreciente
	   paso=atmos(2)-atmos(1)
	   if(paso.ge.0)then
	      men1='STOP: The models supplied are NOT ordered with DECREASING log(tau).'    
	      call mensaje(1,men1,men1,men1)
	   end if

           do i=2,ntau-1
	      paso1=atmos(i+1)-atmos(i)
              paso2=atmos(8*ntau+i+3)-atmos(8*ntau+i+2)
              if(abs(paso1-paso).gt.epsilon)then
	       if(nciclos.gt.0)then
		 men1='STOP: Model 1 is NOT equally spaced. For inversion, an evenly spaced grid is required'
	         call mensaje(1,men1,men1,men1)
	       endif
              end if
              if(abs(paso2-paso1).gt.epsilon)then
	         men1='STOP: The same spatial grid has to be used to discretize models 1 and 2.'       
	         call mensaje(1,men1,men1,men1)
              end if
	   end do


           nohaycampo1=0
           nohaycampo2=0
           do i=1,ntau
              if(abs(atmos(4*ntau+i)).gt.1.)nohaycampo1=1
              if(abs(atmos(12*ntau+i+2)).gt.1.)nohaycampo2=1
           end do            
	else
c escribimos los modelos

	   peso=atmos(16*ntau+5)

	   open(ican,file=model1,err=990)
           goto 800 
990        print*,' '
	   print*,'WARNING: The file containing model 1 does NOT exist'
           ierror=1
	   goto 989

800	   write(ican,*)atmos(8*ntau+1),atmos(8*ntau+2),peso
	   do i=1,ntau
	      if(atmos(i+ntau) .lt. 0)atmos(i+ntau)=0.            !T
              if(atmos(i+ntau) .gt. 9.999e4)atmos(i+ntau)=9.999e4 !T
	      write(ican,100)(atmos(i+j*ntau),j=0,7),z1(i),pg1(i),ro1(i)
	   end do
	   close(ican)
          	
989	   if(ierror.eq.0)then
	      open(ican,file=model2,err=988)
	      write(ican,*)atmos(16*ntau+3),atmos(16*ntau+4),peso
	      do i=1,ntau
	      	 if(atmos(i+9*ntau+2) .lt. 0)atmos(i+9*ntau+2)=0.            !T
             if(atmos(i+9*ntau+2) .gt. 9.999e4)atmos(i+9*ntau+2)=9.999e4 !T
	         write(ican,100)(atmos(i+j*ntau+2),j=8,15),z2(i),pg2(i),ro2(i)
	      end do
	      close(ican)
	   end if
        end if 

100     format(1x,f7.4,1x,f8.1,1x,1pe12.5,1x,e10.3,1x,7(e11.4,1x))
	return

991     men1='STOP: The file containing model 1 does NOT exist'
	call mensaje(1,men1,men1,men1)


988 	print*,'WARNING: Is the file containing model 2, ',model2,', unexistent?? '

        if(ierror.eq.0)then
           ierror=2
        else
           ierror=3
        end if 
	return	


999     men1='STOP: Incorrect format in the file(s) containing the model(s).'
	call mensaje(1,men1,men1,men1)

	end
