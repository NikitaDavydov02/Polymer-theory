subroutine enter 
	Implicit none
	Real*8 rN,rNA,rNB,sigma,R,aA,nu,coe,BA,pi,y_min,y_max,yacc,Lamb_Pol
	Real*8 chi(5,5),Nal(5),lambda(5), Fibulk(5),Lagrbulk(5),osmbulk,Osmmix
	integer i,j
	common/elastic/ R,BA,y_min,y_max,yacc,aA,rNA,rNB,sigma
	common/chilam/chi,Nal,lambda,osmbulk,Lamb_Pol

	! give all parameters and constants:
	pi = 3.1415926d0
	coe = (3.0/8.0)*pi**2
	aA = 6.8d-9 ! in decimeters Kuhn segment length
	rN =  80.0d0 !total number of polymer segments per chain
	rNA = 60.0d0 !number of segments in A-subchain
	rNB = rN - rNA
	nu = 2.0d0	                             !spherical micelle
	sigma=(7.0d-09)**2/0.12
	sigma=(7.0d-09)**2/0.6d0
	R = rNB*(nu+1.)*aA**3/sigma				 ! core radius	 ((nu+1.0)*rNB/SIGMA)**2	LAGRANGE_STR1 = BA*(y_try-1.0)**2*((nu+1.0)*rNB/SIGMA)**2*aA**3  CHECK THIS
	!write(*,*) R,sigma
	!stop
	BA = coe/(rNA*aA)**2 

	y_min = 1.0d0 + 0.1*aA/R
	y_max = 1.00001d0*(1.0d0 +1.0d0*aA*rNA/R)                   !y_max = 1.0d0*(1.0d0 +2.0*aA*rNA/R)
	yacc = 1.0d-8  

	! sigma min for all morphologies:
	!sigma_MIN=10.0*(nu+1.1)*aA**2
	!sigma_MAX=3.0d0*aA**2*(3.1415926*4.0*rNB**2/3.0d0)**(1.0/3.0)
	!write(*,*) BA,R,rNB*aA,sigma,10.0*(nu+1.1)*aA**2,y_max,BA*(R*(y_max-1.0))**2,
	!stop

	 Nal(1)=1.0d0    ! solvent
	 Nal(2)=3.0d0	 ! bioadditive
	 Nal(3)=rNA		 ! polymer

	 chi(1,2)= 0.3d0		! solv-bio
	 chi(1,3)=0.0d0			! solv-polym
	 chi(2,3)=0.05          ! bio-polym

	 chi(1,2)=1.0d0		! solv-bio
	 chi(1,3)=0.0d0			! solv-polym
	 chi(2,3)=-0.8d0          ! bio-polym
 
	 do i=1,5
		chi(i,i)=0.0d0
		do j=i+1,5
			chi(j,i)=chi(i,j)
		enddo
	 enddo

	! give bulk composition and calculate Lagr.multipliers and osm pressure in the bulk:
	Fibulk(1) =0.98d0
	Fibulk(2)=1.0d0-Fibulk(1)
	FiBulk(3)=0.0d0

	call  Lagrmix(2,Fibulk,Lagrbulk)
	lambda(2)=Lagrbulk(2)		  !  this is for biocomponent 
	osmbulk=Osmmix(2,Fibulk)  	  !  this is for solvent
	lambda(1)=Lagrbulk(1)         !  this should be identiacally zero

	! now calculate  mixing part of Lagrange multiplier at the edge of the brush, where  FiA=0:
	FiBulk(3)=0.0d0
	call  Lagrmix_PolA(3,Fibulk,Lagrbulk) ! NOTE  that  Lagr multipliers are spoiled in  Lagrbulk   but stored in common
	lambda(3)=Lagrbulk(3)
	return
end subroutine
!****************************************************************************************
subroutine Lagrmix(K,X,DfmixDfi)
	!input: K -number of components    X-volume fractions;   output:   DfmixDfi    dfmix/dFI(I) 
 	use GuggenheimModel 
	implicit none
	real*8 	chi(5,5),Nal(5),lambda(5),X(5),DfmixDfi(5),sum,osmbulk,Lamb_Pol, AlternativeDfmixDfi(5)
	integer K,i,j
	common/chilam/chi,Nal,lambda,osmbulk,Lamb_Pol

	do i=1,K
		sum=0.0d0
		do j=1,K
			sum=sum+X(j)*(chi(i,j)-chi(1,j))
		enddo
		DfmixDfi(i)=(dlog(X(i))+1.0d0)/Nal(i)-(dlog(X(1))+1.0d0)/Nal(1)+sum   ! dummy for solvent  identically 0
		!AlternativeDfmixDfi(i)=CalculateExchangeChemialPotentialOfComponent(3,X,i)
		!DfmixDfi(i)=AlternativeDfmixDfi(i)
	enddo
	return
end subroutine
!********************************************************************************************
subroutine Lagrmix_PolA(K,X,DfmixDfi)
	!input: K -number of components    X-volume fractions;   output:   DfmixDfi    dfmix/dFI in the brush 
 	use GuggenheimModel 
	implicit none
	real*8 	chi(5,5),Nal(5),lambda(5),X(5),DfmixDfi(5),sum,osmbulk,Lamb_Pol, AlternativeDfmixDfi(5)

	integer K,i,j,n
	common/chilam/chi,Nal,lambda,osmbulk,Lamb_Pol

	do i=1,K
		if (X(i).LT.0) then
			do n=1,K
				DfmixDfi(n)=50.0d0*((X(n)-0.5d0)/(X(n)-0.5d0))**2 +1.0d+03
			enddo
			return
		endif

		sum=0.0d0
		do j=1,K
			sum=sum+X(j)*(chi(i,j)-chi(1,j))
		enddo
		if(i.eq.K)then
			DfmixDfi(i)=-(dlog(X(1))+1.0d0)/Nal(1)+sum   ! for immobile polymer A
		else
			DfmixDfi(i)=(dlog(X(i))+1.0d0)/Nal(i)-(dlog(X(1))+1.0d0)/Nal(1)+sum   ! dummy for solvent  identically 0
		endif
		!AlternativeDfmixDfi(i)=CalculateExchangeChemialPotentialOfComponent(3,X,i)
		!DfmixDfi(i)=AlternativeDfmixDfi(i)
	enddo
	return
end subroutine
!******************************************************************************************
real*8 function Osmmix(K,X)
	 !input: K -number of components    X-volume fractions;   output:   mu solvent 
	implicit none
	real*8 	chi(5,5),Nal(5),lambda(5),X(5),sum,sum1,osmbulk,Lamb_Pol
	integer K,i,j
	common/chilam/chi,Nal,lambda,osmbulk,Lamb_Pol

	sum=0.0d0
	sum1=0.0d0
	do i=1,K
		sum1=sum1+X(i)*(1.0d0-1.0d0/Nal(i))	                 !  solvent does not contribute as long as Nal(1)=1
		do j=1,K
			sum=sum+X(i)*X(j)*(chi(1,j)-chi(i,j)/2.0d0)
		enddo
	enddo
	Osmmix=dlog(X(1))+sum1+sum
	return
end
!********************************************************************************************


