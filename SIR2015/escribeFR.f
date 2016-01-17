	subroutine escribeFR(m,ndata,ntotal4,npos,dscal)

	implicit real*4 (a-h,o-z)

	include 'PARAMETER'  !por kn,kl,kld

	parameter (kn16=16*kn+4,kt16=16*kt+5)         !
        parameter (mmax=16*kn+2,kld4=4*kld,kldn=mmax*kld4) 
	parameter (kldt=kld*kn,kldt4=4*kldt)
 
	integer m(*),ndata,ntotal4,npos(*)  
        real*4 dscal(*)
	real*4 FR(kldt4),atmosr(kt16),norma(kt16)

        character*100 modelin1,modelout1,modelFR
        character*100 modelin2,modelout2,modelFR2
	character*5 cha(18)

        common/nommod/modelin1,modelin2       !de SIR
        common/atmosr/atmosr                  !de SIR


	modelout1=modelin1
	modelout2=modelin2

	call quitaex0(modelout1) !quita la extension detras del punto
	call quitaex0(modelout2) 

        cha(1)='.rt'
        cha(2)='.rpe'
        cha(3)='.rmic'
        cha(4)='.rh'
        cha(5)='.rvz'
        cha(6)='.rinc'
        cha(7)='.raz'
        cha(8)='.rmac'

        cha(9)='.rt'
        cha(10)='.rpe'
        cha(11)='.rmic'
        cha(12)='.rh'
        cha(13)='.rvz'
        cha(14)='.rinc'
        cha(15)='.raz'
        cha(16)='.rmac'
        cha(17)='.rff'

        cha(18)='.rstl'
      
        ki=0
        do i=1,5
           do j=1,m(i)
              ki=ki+1
              norma(ki)=atmosr(ki)
c              if(i.eq.2)norma(ki)=exp(atmosr(ki)) !presion
           end do
        end do
        do i=6,7   !gamma y fi son pert. absolutas
           do j=1,m(i)
              ki=ki+1
              norma(ki)=1.
           end do
        end do
        do i=8,13
           do j=1,m(i)
              ki=ki+1
              norma(ki)=atmosr(ki)
c              if(i.eq.10)norma(ki)=exp(atmosr(ki)) !presion
           end do
        end do
        do i=14,15   !gamma y fi son pert. absolutas
           do j=1,m(i)
              ki=ki+1
              norma(ki)=1.
           end do
        end do
        do i=16,18
           do j=1,m(i)
              ki=ki+1
              norma(ki)=atmosr(ki)
           end do
        end do

	k=0
        k1=0
        do i=1,8
           if(m(i).ne.0)then
	      call extraeFR(m(i),FR,k,dscal) 
              call concatena(modelout1,cha(i),modelFR)
              print*,'Output RF file ',modelFR
              open(19,file=modelFR)
                 write(19,*) m(i),ntotal4
                 write(19,*) ((FR((j-1)*ntotal4+jj)/norma(k1+j),jj=1,ntotal4)
     &                      ,j=1,m(i))
              close(19)
           endif
           k1=k1+m(i)
        enddo
        do i=9,17
           if(m(i).ne.0)then
	      call extraeFR(m(i),FR,k,dscal) 
              call concatena(modelout2,cha(i),modelFR2)
              print*,'Output RF file  ',modelFR2
              open(19,file=modelFR2)
                 write(19,*) m(i),ntotal4
                 write(19,*)
     &           ((FR((j-1)*ntotal4+jj)/norma(k1+j),jj=1,ntotal4),j=1,m(i))
              close(19)
           endif
           k1=k1+m(i)
        enddo
	i=18
        if(m(i).ne.0)then
	   call extraeFR(m(i),FR,k,dscal) 
           call concatena(modelout1,cha(i),modelFR)
           print*,'Output RF file  ',modelFR
           open(19,file=modelFR)
                 write(19,*) m(i),ntotal4
                 write(19,*)
     &           ((FR((j-1)*ntotal4+jj)/norma(k1+j),jj=1,ntotal4),j=1,m(i))
           close(19)
        endif

       
        return
        end 

c _____________________________________________________________________________

	subroutine extraeFR(mi,rt,k,dscal)

	implicit real*4 (a-h,o-z) 

	real*4 rt(*),dscal(*)
        common/primera2/ntotal,ntotal4,ists

	do i=1,mi
	   npun=ntotal4*(i-1)+1
           do j=npun,npun+ntotal4-1
              k=k+1 
	      rt(j)=dscal(k)         !t1(l1,l2...lntot4),t2(l1,l2...
	   end do
        end do		
	return
	end
c _____________________________________________________________________________
    





