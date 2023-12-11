module GuggenheimModel
contains
	real*8 function CalculateMixingFreeEnergy(K,X)
		!input: K -number of components    X-volume fractions;   output:   DfmixDfi    dfmix/dFI(I) 
		!1-solvent
		!2-contaminant
		!3-polymer
		!Calculates Fmix/kT in Guggenheim model with correlated strong interaction and non-correlated other interactions
		!but with taking into account that strong interaction occupies some amount of pairs from non-correlated pairs
		implicit none
		real*8 	chi(5,5),Nal(5),lambda(5),X(5),sum,osmbulk,Lamb_Pol,Fmix,beta
		integer K,z
		common/chilam/chi,Nal,lambda,osmbulk,Lamb_Pol
		common/additionalParameters/z


		!z=4	!Coordination number
		beta=sqrt(1+4*X(2)*X(3)*(exp(2*chi(2,3)/z)-1)) !Correlator for strong interaction

		Fmix=X(1)*dlog(X(1))+X(2)*dlog(X(2))/Nal(2)
		Fmix=Fmix+chi(2,3)*X(2)*X(3)*2/(beta+1)
		Fmix=Fmix+chi(1,2)*X(1)*(X(2)-(1-2/z)*X(3)*X(2)*2/(beta+1))/(X(1)+X(2)-(1-2/z)*X(3)*X(2)*2/(beta+1))
		Fmix=Fmix+chi(1,3)*X(1)*((1-2/z)*X(3)-(1-2/z)*X(3)*X(2)*2/(beta+1))/(X(1)+(1-2/z)*X(3)-(1-2/z)*X(3)*X(2)*2/(beta+1))
		CalculateMixingFreeEnergy=Fmix

		return
	end function
	!*********************************************************************************************************************
	real*8 function CalculateExchangeChemialPotentialOfComponent(numberOfComponents,volumeFractions,component)
		!1-solvent
		!2-contaminant
		!3-polymer
		!Calculates exchange chemical potential of component in the mixture based on Fmix energy
		implicit none
		real*8 	chi(5,5),Nal(5),volumeFractions(5),beta,lambda(5),osmbulk,Lamb_Pol
		real*8 x,x_dx,f,f_df
		integer numberOfComponents,component,z
		common/chilam/chi,Nal,lambda,osmbulk,Lamb_Pol
		common/additionalParameters/z

		CalculateExchangeChemialPotentialOfComponent=0

		return
	end function

end module





