      subroutine marquarderr(y,ndata,a,mnodos,n)

	implicit real*4 (a-h,o-z)

	include 'PARAMETER'  !por mmx==kn, que aqui se tomaba igual a 200
c        parameter (mmx=200)
c        dimension y(ndata),a(*),mnodos(*),beta(200),errcv(mmx)  

        dimension y(ndata),a(*),mnodos(*),beta(kn),errcv(kn)  

        common/errores/errcv

	nn=n 
        call marqcoeferr(y,ndata,a,mnodos,nn,beta)

	 rn=1./float(nn) 
c	 rn=1.      
	 do i=1,nn
            b=beta(i)
            if(b.lt.1.e-10)b=1.e-10
	    errcv(i)=rn/b
	 end do 
       
        return
	end
