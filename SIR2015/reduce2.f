c ...................................................................
	subroutine reduce2(ntau,no,atmos,atmosr)

	implicit real*4 (a-h,o-z)
	real*4 atmos(*),atmosr(*)

	m=(ntau-1)/(no+1)
	ik=0
	do k=1,8
	   do i=1,no+2
	      ik=ik+1
	      j=(k-1)*ntau+(i-1)*m+1
	      jk=(i-1)*m+1
	      if(k.eq.7)then
	         atmosr(ik)=tan(atmos(j)/2.0)
	      else if(k.eq.8)then
	         atmosr(ik)=tan(atmos(j)/4.0)
	      else
	         atmosr(ik)=atmos(j)
	      end if
	   end do
	end do

	atmosr(8*no+17)=atmos(8*ntau+1)
	atmosr(8*no+18)=atmos(8*ntau+2)
	return
	end
c ......................................................................
