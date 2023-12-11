Real*8 Function rtbis3(normal,y_min,y_max,yacc)
	!Solving non-linear equation with bisection method
	INTEGER JMAX
	Real*8 y_min,y_max,yacc,normal
	EXTERNAL normal
	PARAMETER (JMAX=100) !Maximum allowed number of bisections.
	!Using bisection, find the root of a function func known to lie between x1 and x2. The
	!root, returned as rtbis, will be refined until its accuracy is �xacc.
	INTEGER j
	Real*8 dx,f,fmid,xmid

	!calculating function value at amx
	fmid = normal(y_max)
	f = normal(y_min) !R2(beta2,lambda,itog2)

	if(f*fmid.ge.0.) pause "root must be bracketed in rtbis3"
	if(f.lt.0.)then !Orient the search so that f>0 lies at x+dx.
		rtbis3=y_min
		dx=y_max-y_min
	else
		rtbis3=y_max
		dx=y_min-y_max
	endif
	do j=1,JMAX !Bisection loop.
		dx=dx*.5d0
		xmid=rtbis3+dx

		fmid = normal(xmid) !call R3(beta3,lambda,itog3)

		if(fmid.le.0.)rtbis3=xmid
		if(abs(dx).lt.yacc .or. fmid.eq.0.) return
	enddo 
	pause "too many bisections in rtbis3"
	return
End Function




