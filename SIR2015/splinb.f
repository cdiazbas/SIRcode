
c SPLINB
	SUBROUTINE SPLINB(N,nx,A)
	implicit real*4 (a-h,o-z)

	include 'PARAMETER'   !por nmax===kn
c	PARAMETER (nmax=1280)
c	dimension a(nx,nx),b(nmax,nmax),c(nmax,nmax),indx(nmax)
	dimension a(nx,nx),b(kn,kn),c(kn,kn),indx(kn)

	do i=1,n		!inicializo a cero
	   do j=1,n
	      a(i,j)=0.d0
	      b(i,j)=0.d0
           end do
	end do
	
	do i=2,n-1
	   a(i,i)=4.d0          !diagonales
	   b(i,i)=-2.d0
	   a(i,i-1)=1.d0        !bajo la diagonal
	   b(i,i-1)=1.d0
	   a(i,i+1)=1.d0	!sobre la diagonal
	   b(i,i+1)=1.d0
	end do

	a(1,1)=1.d0             !escribo las esquinas
	a(1,2)=-5.d0
	a(1,3)=1.d0
	a(n,n)=1.d0
	a(n,n-1)=-5.d0
	a(n,n-2)=1.d0
	b(1,1)=-5.d-1           !escribo las esquinas
	b(1,2)=1.d0
	b(1,3)=-5.d-1
	b(n,n)=-5.d-1
	b(n,n-1)=1.d0
	b(n,n-2)=-5.d-1
	
	do i=1,n                !matriz diagonal
	   do j=1,n
	      c(i,j)=0.d0
	   end do
	   c(i,i)=1.d0
	end do

	CALL LUDCMP(A,N,NX,INDX,D)
	do j=1,n
	   CALL LUBKSB(A,N,NX,INDX,c(1,j))
	end do

	do i=1,n
	   do j=1,n
	      a(i,j)=0.d0
              do k=1,n
                 a(i,j)=a(i,j)+c(i,k)*b(k,j)
	      end do
	   end do
	end do

	return
	end
