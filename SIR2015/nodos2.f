c rutina nodos2

	subroutine nodos2(ntau,mnodos)

	include 'PARAMETER'   !por kt
	integer ndiv(kt),mnodos(*)
	character*4 vari(18)
      data vari/'t1  ','p1  ','mic1','h1  ','vz1 ','gam1','fi1 ','mac1',
     &          't2  ','p2  ','mic2','h2  ','vz2 ','gam2','fi2 ','mac2',
     &          'f.f.','  % '/

c calculamos todos los divisores de ntau
	ndiv(1)=-1	!el primer divisor (perturbacion = otra variable)
	ndiv(2)=0	!segundo divisor   (no hay perturbacion)
	ndiv(3)=1	!el tercero 	   (perturbacion constante)
	ndiv(4)=2	!los 2 extremos    (perturbacion lineal)
	k=4	        !ya llevamos 4 
	do i=3,ntau     !es i-1 divisor de ntau-1?
	   m=( (ntau-1)/(i-1) )*(i-1)
	   if(m.eq.(ntau-1))then
	      k=k+1
	      ndiv(k)=i
           end if
	end do
        num=k

        if(mnodos(8).gt.1)then
           mnodos(8)=1 
	   print*,'el num. de nodos de mac1. debe ser menor/igual a 1 (tomo 1)'
	end if
        if(mnodos(16).gt.1)then
           mnodos(16)=1 
	   print*,'el num. de nodos de mac2. debe ser menor/igual a 1 (tomo 1)'
	end if
        if(mnodos(17).gt.1.or.mnodos(17).lt.0)then
           mnodos(17)=1 
	   print*,'el num. de nodos del ff debe ser 0 o 1 (tomo 1)'
	end if
        if(mnodos(18).gt.1.or.mnodos(18).lt.0)then
           mnodos(18)=1 
	   print*,'el num. de nodos del % debe ser 0 o 1 (tomo 1)'
	end if

c ahora buscamos el divisor menor o igual mas cercano a cada mnodos
	do j=1,18
	   mold=mnodos(j)
	   if(mnodos(j).lt.-1)mnodos(j)=-1
           k=1
           do while( ndiv(k).le.mnodos(j) .and. k.le.num)
              k=k+1
           end do
           mnodos(j)=ndiv(k-1)
	   if(mnodos(j).ne.mold)print*,'num. de nodos de ',vari(j),'=',mnodos(j)
	end do
	return
	end

