c taulinea2
c il=0 :pasa de una atmosfera modelo en direccion z a una
c       atmosfera en la linea de vision
c il=1 :transforma una atmosfera en la linea de vision a una atmosfera
c       en direccion z
c ntau : numero de puntos en tau
c tau  : log10(tau)
	subroutine taulinea2(il,cth,ihe,vx,atmos,ntau)

	include 'PARAMETER'   !para kt
	parameter (kt8=8*kt+2)
	implicit real*4 (a-h,o-z)
	integer ihe,il,ntau
	real*4 atmos(*),atmos1(kt8),atmos2(kt8),cth,vx
c	real*4 gammal1(kt),gammal2(kt),phil1(kt),phil2(kt)
c	common/LOS/ gammal1,phil1,gammal2,phil2  !para ASP: subsir

	do i=1,8*ntau+2
           atmos1(i)=atmos(i)
           atmos2(i)=atmos(i+8*ntau+2)
	end do
        
	if(cth.ne.1.)then
	   print*,'The atmospheres are assumed to be line of sight!!'
           print*,'However, hydrostatic equilibrium is computed along'
           print*,'the vertical (using the heliocentric angle provided)',cth
	end if

c	call taulinea(il,cth,ihe,vx,atmos1,ntau,gammal1,phil1)  ! gammal1???
c	call taulinea(il,cth,ihe,vx,atmos2,ntau,gammal2,phil2)

        call taulinea(il,cth,ihe,vx,atmos1,ntau)
        call taulinea(il,cth,ihe,vx,atmos2,ntau)

	do i=1,8*ntau+2
           atmos(i)=atmos1(i)
           atmos(i+8*ntau+2)=atmos2(i)
	end do

	return
	end

