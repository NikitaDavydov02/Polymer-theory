subroutine FI_POLI(XBrush,y_cur)
	!Solving (localy) system of non-linear equations 
	Implicit none
	Integer L,ITMAX
	Parameter (L=2)
	Real*8 ERRREL,FNORM,XBrush(L),XBrushGUESS(L),y_cur,y
	common/point/y
	external EU_LA

	 y=y_cur						  ! transmitt the current coordinate to EU_LA
	 ITMAX=600
	 ERRREL=1.0d-4
	 XBrushGUESS(1)=	1.0d-8		  ! this is the fraction of biocomponent in the brush
	 XBrushGUESS(2)=	0.97d0		  ! this is the fraction of polymer in the brush

	 CALL DNEQNF (EU_LA, ERRREL, 2, ITMAX, XBrushGUESS, XBrush, FNORM)

	 return
 end subroutine
!************************************************************************************************************************
!Calculating local Lagrangian multipliers		   
Subroutine EU_LA (XBrush,F,L)
	implicit none
	Integer L
	Real*8 fipolimer(5),DfmixDfi(5), XBrush(L),F(L)
	Real*8 rNA,rNB,sigma,R,aA,nu,BA,y_min,y_max,yacc,Lamb_Pol,y,y_cur
	Real*8 chi(5,5),Nal(5),lambda(5), osmbulk
	common/elastic/ R,BA,y_min,y_max,yacc,aA,rNA,rNB,sigma
	common/chilam/chi,Nal,lambda,osmbulk,Lamb_Pol
	common/point/y
	external Lagrmix_PolA
	nu=2.0d0
	y_cur=y
	fipolimer(2) =Xbrush(1)		   ! local compsition of the brush at point y_cur     (1)-solvent   (2)-biocomponent  (3)-polymer
	fipolimer(3) =Xbrush(2)
	fipolimer(1) =1.0d0- fipolimer(2)- fipolimer(3) !solvent


	!Calculate values that in ideal case must be equal to Lagrangian multipliers based on current concentrations
	call Lagrmix_PolA(3,fipolimer,DfmixDfi)


	F(1) = (DfmixDfi(2)-lambda(2))**2  !bio contaminant error
	F(2) = (DfmixDfi(3)+BA*(R*(y_cur-1.0d0))**2-Lamb_Pol)**2  !polymer error

	write(*,*) F,Xbrush

	return
end subroutine