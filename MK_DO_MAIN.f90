
Program NumSOL
	!use GuggenheimModel

	Implicit None
	Integer L
	Parameter (L=2)
	Real*8 fipolimer(5),XBrush(L)
	Real*8 rNA,rNB,sigma,R,aA,BA,y_min,y_max,yacc,rtbis3,y_cur,y_try,y_edge,P_AB,fiav,Osmmix
	Real*8 chi(5,5),Nal(5),lambda(5),osmbulk,Lamb_Pol
	Real*8 Fmix, u_pol,u_sol,u_bio

	common/elastic/ R,BA,y_min,y_max,yacc,aA,rNA,rNB,sigma
	common/chilam/chi,Nal,lambda,osmbulk,Lamb_Pol

	External rtbis3,normal,FI_POLI,fiav

	open (unit=1,file = "profile.txt")
	write(1,*) '   y_cur      solvent   bio     polymer    osm pressure'

	call enter

	y_try=y_min
	y_edge = rtbis3(normal,y_min,y_max,yacc)

	y_cur=y_min

	!Calculate polymer concentration at the start position (A/B boundary)
	P_AB=fiav(y_cur,y_edge)      !   Polymer concentration at A/B boundary (highest)  must work before finding profile to set  Lagr multipliers in common block

	do while (y_cur.LT.y_edge)
		y_cur = y_cur + aA/R  
		!write(*,*) 'Phi poly subroutine is called'


		call FI_POLI(Xbrush,y_cur)	   ! calculates concentration profile in the brush after the solution is found
		fipolimer(2) =Xbrush(1)		   ! local compsition of the brush at point y_cur     (1)-solvent   (2)-biocomponent  (3)-polymer
		fipolimer(3) =Xbrush(2)
		fipolimer(1) =1.0d0- fipolimer(2)- fipolimer(3)

		!Fmix=0
		
		!write(*,*) 'phi_solvent=', fipolimer(1)
		!write(*,*) 'phi_conteminant=', fipolimer(2)
		!write(*,*) 'phi_polymer=', fipolimer(3)


		!Fmix= CalculateMixingFreeEnergy(3,fipolimer)
		!u_pol=CalculateExchangeChemialPotentialOfComponent(3,fipolimer,3)
		!u_sol=CalculateExchangeChemialPotentialOfComponent(3,fipolimer,1)
		!u_bio=CalculateExchangeChemialPotentialOfComponent(3,fipolimer,2)


		write(*,*) 'Fmix=', Fmix
		write(*,*) 'u_sol=', u_sol
		write(*,*) 'u_bio=', u_bio
		write(*,*) 'u_pol=', u_pol

		!Output resulting concentration profile
		write(1,'(F10.4,3(2X,F8.4),2X,E12.4)') y_cur, fipolimer(1),fipolimer(2),fipolimer(3), (Osmmix(3,fipolimer)-osmbulk)
	enddo
 
	write(*,*) 'beta',y_edge
	write(*,*) 'Calculation is done!'
	stop
END Program

!**********************************************************************************************
real*8 function normal(y_try)
	!Normalization equation for polymer profile
	implicit none
	real*8 nu,y_try
	real*8 aINT,bINT,s,R,BA,y_min,y_max,yacc,aA,rNA,rNB,sigma
	common/elastic/ R,BA,y_min,y_max,yacc,aA,rNA,rNB,sigma

	nu =2.0d0
	aINT = 1.0d0
	bINT = y_try
	call qtrap(aINT,bINT,s)
	normal = s - rNA/(rNB*(nu+1.0d0))
	return
end function
!************************************************************************************************
subroutine trapzd(aINT,bINT,s,n)
	!Integrate fiav (profile) function from aINT to bINT
	INTEGER n
	double precision aINT,bINT,s,fiav

	EXTERNAL fiav
	INTEGER it,counter
	double precision del,sum,tnm,x

	if (n.eq.1) then
		s=0.5d0*(bINT-aINT)*(fiav(aINT,bINT)+fiav(bINT,bINT))
		
	else
		it=2**(n-2)
		tnm=it
		del=(bINT-aINT)/tnm
		x=aINT+0.5d0*del
		sum=0.0d0
		do counter=1,it
			sum=sum+fiav(x,bINT)
			x=x+del
		enddo
		s=0.5d0*(s+(bINT-aINT)*sum/tnm)
	endif
	return
END
!*****************************************************************************************
subroutine qtrap(aINT,bINT,s)
	integer JMAX
	double precision aINT,bINT,fiav,s,EPS

	EXTERNAL fiav
	PARAMETER (EPS=1.e-1,JMAX=8)
	integer n
	double precision olds
	olds=-1.e-30
		do n=1,JMAX
			call trapzd(aINT,bINT,s,n)
			if (n.gt.5) then
				if (dabs(s-olds).lt.EPS*dabs(olds).or.(s.eq.0.and.olds.eq.0.)) return
			end if
			olds=s
		
		end do
		pause "too many steps in qtrap"
	return
END

!***************************************************************************************
Real*8 Function fiav (y_cur,bINT)
	Implicit None
	Integer L
	Parameter (L=2)
	Real*8 y_cur,y_try,fay,Xbrush(L),bINT
	Real*8 rNA,rNB,sigma,R,aA,nu,BA,y_min,y_max,yacc
	Real*8 chi(5,5),Nal(5),lambda(5), osmbulk,Lamb_Pol
	common/elastic/ R,BA,y_min,y_max,yacc,aA,rNA,rNB,sigma
	common/chilam/chi,Nal,lambda,osmbulk,Lamb_Pol
	 nu=2.0d0

	 y_try=bINT

	Lamb_Pol=lambda(3)+BA*(R*(y_try-1.0))**2			! form Lagr multiplier for polymer A at current thickness of the brush (y_try=BINT)

	call FI_POLI(Xbrush,y_cur)

	fay=Xbrush(2)
	fiav=fay*y_cur**2
	return 
end function
!******************************************************************************************
!real*8 function fiavZul(y_cur,bINT)
! b - is y_try 
!double precision y_try,nu,BA,aA,rNA,rNB,y_cur,bINT
!double precision z1,z2,zacc,fiBn,fay, rtbis2,R
!integer iflag
!common/elastic/ R,BA,y_min,y_max,yacc,aA,rNA,rNB
!common/flagg/ iflag
!external rtbis2
!
!y_try = bINT
!
!fiBn=0.0
!call eptit(fiBn,lambda,y_try)	 
!z1 = 5.0d-9
!
!z1 = 5.0d-13
!
!z2 = 1.0d0-z1
!zacc = 1.0d-9
!
!fay = rtbis2(z1,z2,zacc,lambda,y_cur)	
!if (iflag.LT.0) then
!write (*,*) z1,z2,y_cur,y_try,fay
!endif
!fiav = fay*y_cur**nu
!return
!end

 