	subroutine alphasub(tau,n,alpha,depar)


	include 'PARAMETER'
        parameter (num=19)
	real*4 tau(*),alpha(*),depar(*),x(num),y(num),z(num)
        real xr(num)
	real xa(11),ya(11)

	data x/ 
     &    0.5235 ,   -0.0020 ,   -0.3570 ,   -0.6574 ,   -0.9519 ,
     &   -1.2627 ,   -1.5735 ,   -1.9115 ,   -2.2498 ,   -2.6212 ,
     &   -2.9927 ,   -3.4148 ,   -3.8368 ,   -4.1535 ,   -4.3514 ,
     &   -4.4831 ,   -4.5957 ,   -4.7022  ,  -8.0 /

	data y/ 
     &    1.00000e+00 ,    1.00000e+00 ,    9.89000e-01 ,    9.79185e-01 ,
     &    9.63803e-01 ,    9.34267e-01 ,    8.81333e-01 ,    8.11867e-01 ,
     &    7.05453e-01 ,    5.86343e-01 ,    4.67177e-01 ,    3.55368e-01 ,
     &    2.61316e-01 ,    1.87804e-01 ,    1.29386e-01 ,    8.72174e-02 ,
     &    5.19049e-02 ,    2.72269e-02 ,     0./

	data z/ 
     &    1.00000e+00 ,    1.01600e+00 ,    1.04200e+00 ,    1.04200e+00 ,
     &    1.02100e+00 ,    9.95000e-01 ,    9.60000e-01 ,    9.21000e-01 ,
     &    8.98000e-01 ,    8.98000e-01 ,    9.35000e-01 ,    1.10800e+00 ,
     &    1.60500e+00 ,    3.03500e+00 ,    6.90700e+00 ,    1.63790e+01 ,
     &    3.80500e+01 ,    9.21060e+01 ,    9.2106e1   /
   
	ngrado=1  !grado del polinomio interpolador

	n2=int(ngrado/2)

c        print*,'entro en alphasub'
        do i=1,num
	   xr(i)=x(i)
        end do
 	
	do i=1,n
           stau=tau(i) 
	   call locate(xr,num,stau,j)
	   n3=j-n2-1
           if(n3.lt.0)n3=0
           if(n3+ngrado+1.gt.num)n3=num-ngrado-1
	   do k=1,ngrado+1
	      xa(k)=x(n3+k)
	   end do
	   do k=1,ngrado+1
	      ya(k)=y(n3+k)
	   end do
	   call polint(xa,ya,ngrado+1,stau,sa,error)
           alpha(i)=dble(sa)

	   do k=1,ngrado+1
	      ya(k)=z(n3+k)
	   end do
	   call polint(xa,ya,ngrado+1,stau,da,error)
           depar(i)=dble(da) 
c         print*,tau(i),alpha(i),depar(i)
	end do

c         print*,'salgo de alphasub'


	return
	end

	 
