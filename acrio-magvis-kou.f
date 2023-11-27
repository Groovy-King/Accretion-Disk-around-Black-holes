      implicit real*8(a-h,o-z)
c  This programe gives the solution of accretion flows from inner sonic
c  point to outside in case of pseudo-Kerr-potential for magnetized viscous
c  flow in 1.5D ----- v_z is chosen zero --- 
c  We are now including dvdx_c (of critical point) self-consistently.

      common/aaa/alfa,xf,gama,xmdot,xmu,xmui,xmue,xmp,xk,xo,an,xlc,beta
      common/ddd/b,g1,g3,xIfac,xInp1,xIn,xmas
	common/bbb/xmx,xmy,xmz,s1,s2,pi,xfm,l,ll
c good result for l=0: xo=5.8, xlc=3.2 alfa=.01 xf=.5 b=0. uptil sonic point
c with   x4=xo
c        v4=v-1.d-7
c        sou4=sou-1.d-7
c        xl4=xlc+1.e-7
c        rho4=rho-1.d-24
c        P4=P-1.e-26

	write(*,*)'give xo,xlc,alfa,xf,b,ll'
	read(*,*)xo,xlc,alfa,xf,b,ll
	l=1
	xfm=xf
	s3=5.d+1
	pi=acos(-1.)
	xmx=sqrt(3.)*s3
	xmy=sqrt(3.)*s3
	xmz=sqrt(3.)*s3
	s1=-0.5
	s2=0.5
	xmas=1.d+1
	xmdot=3.55e-22*1.d-1*xmas
	gama=1.335
	xmu=0.5
	xmui=1.23
	xmue=1.14
	xmp=1.67e-24
	xk=1.38e-16
	ome=xmu*xmp/gama/xk
	an=6.023e+23
	beta=(6.*gama-8.)/3./(gama-1.)
        g1=beta+(4.-3.*beta)**2*(gama-1.)/(beta+12.*(gama-1.)*(1.-beta))
        g3=1.+(g1-beta)/(4.-3.*beta)
c	print*,gama,g1,g3
        xIfac=2.*gama/(3.*gama-1.)
      x1=0.0000000001
      x2=x1
10    x2=x2+0.1
      if(x2.gt.1.0) then
c      write(*,*) 'no solution'
c      go to 11
	stop
      endif
      if((f(x1)*f(x2)).gt.0.0) go to 10
      call bisec(x1,x2,sou)
      write(*,*) sou

c	pause
c	if(sou.ge.0.)then
        amc=1./sqrt(2.)*sqrt((alfa*xf*gama**2*(-1 + g3)*l**2*xmx+ 
     &gama*(gama*(1.+g1)*l**2+2.*g1*xmx**2)*xmy+sqrt(gama**2*(4.*(-1.-g1
     &+alfa**2*xf*(-1.+g3))*xmx**2*xmy*(-alfa*xf*gama*(-1.+g3)*l**2*xmx+ 
     &2.*gama*g1*l**2*xmy+alfa**2*xf*(-1.+g3)*xmx**2*xmy)+(alfa*xf*gama*
     &(-1.+g3)*l**2*xmx+(gama*(1.+g1)*l**2+2.*g1*xmx**2)*xmy)**2)))/
     &(gama**2*(1.+g1-alfa**2*xf*(-1.+g3))*xmx**2*xmy))
	print*,amc
	v=amc*sou
	go to 12
c	endif
11    continue
c      x1=sou+1.e-10
c      x2=x1
c      go to 10
12	t9=ome*sou**2*9.e+20*beta*1.e-9
	fx=(b**2-2.*b*dsqrt(xo)+xo**2)**2/xo**2/(xo**2-2.*xo+b
     &*dsqrt(xo))**2
        xmac=v/sou
	xlkep=sqrt(xo**3*fx)

	rho=xmdot/v/sou/xo**1.5*sqrt(fx)/4./pi
	P=rho*sou**2/gama
       abx=sqrt(4.*pi*rho)*sou*xo/xmx*l
	bx=sqrt(4.*pi*rho)*sou/xmx*l
	by=sqrt(4.*pi*rho)*sou/xmy*l
	bz=sqrt(4.*pi*rho)*sou/xmz*l
c	abx=bx*xo

	xh=sou*sqrt(xo)/sqrt(fx)
	eta=rho*alfa*sou*xh
	ss=(P+rho)/t9
	etas=eta/ss*1.e+9*(6.67*2.d+25*xmas/3**3/1.d+30)


	write(2,9)xo,v/sou,v,sou,xlc,t9,rho,P,bx,by,bz,xlc/xlkep,etas
     &,dsqrt(bx**2+by**2+bz**2)
c	write(*,*)xo,v/sou,v,sou,xlc,t9,rho,P,bx,by,bz,xlc/xlkep
9       format(14e13.5)

c TO BE AMENDED FROM HERE................


	 x3=xo
	 v3=v
	 sou3=sou
         xl3=xlc
	P3=P
	rho3=rho
       bx3=bx
       by3=by
       bz3=bz

         dx=ll*1.e-4
cc         dx=1.e-3*x3
         v=v3
         x=x3
	 sou=sou3
         xl=xl3
	P=P3
	rho=rho3
	bx=bx3
        by=by3
        bz=bz3
c        print*,x,sou

	 call  f1(x,v,sou,xl,P,rho,bx,by,bz,dvdx,t9)
c         call  ff1(x,v,sou,xl,rho,bx,by,bz,dPdx)
         call  bb1(x,v,sou,xl,rho,bx,by,bz,dPdx)
c         call  fff1(x,v,sou,xl,rho,bx,by,bz,dxldx)
         call  bbb1(x,v,sou,xl,rho,bx,by,bz,dxldx)
c         call  fff1(x,v,sou,xl,dvdx,dadx,dxldx)
         call  c1(x,v,sou,xl,rho,bx,by,bz,drhodx)
         call  cc1(x,v,sou,xl,rho,bx,by,bz,dbydx)
         call  ccc1(x,v,sou,xl,rho,bx,by,bz,dbzdx)

         h1=dx*dvdx
	 h11=dx*dPdx
         h111=dx*dxldx
	 hc1=dx*drhodx
	hc11=dx*dbydx
	hc111=dx*dbzdx

         x=x3+dx/2.
         v=v3+h1/2.
	 P=P3+h11/2.
         xl=xl3+h111/2.
	rho=rho3+hc1/2.
	by=by3+hc11/2.
	bz=bz3+hc111/2.

	 call  f1(x,v,sou,xl,P,rho,bx,by,bz,dvdx,t9)
c         call  ff1(x,v,sou,xl,rho,bx,by,bz,dPdx)
         call  bb1(x,v,sou,xl,rho,bx,by,bz,dPdx)
c         call  fff1(x,v,sou,xl,rho,bx,by,bz,dxldx)
         call  bbb1(x,v,sou,xl,rho,bx,by,bz,dxldx)
c         call  fff1(x,v,sou,xl,dvdx,dadx,dxldx)
         call  c1(x,v,sou,xl,rho,bx,by,bz,drhodx)
         call  cc1(x,v,sou,xl,rho,bx,by,bz,dbydx)
         call  ccc1(x,v,sou,xl,rho,bx,by,bz,dbzdx)

         h2=dx*dvdx
         h12=dx*dPdx
         h122=dx*dxldx
	hc2=dx*drhodx
	hc12=dx*dbydx
	hc122=dx*dbzdx

         x=x3+dx/2.
         v=v3+h2/2.
         P=P3+h12/2.
         xl=xl3+h122/2.
         rho=rho3+hc2/2.
	 by=by3+hc12/2.
	 bz=bz3+hc122/2.

	 call  f1(x,v,sou,xl,P,rho,bx,by,bz,dvdx,t9)
c         call  ff1(x,v,sou,xl,rho,bx,by,bz,dPdx)
         call  bb1(x,v,sou,xl,rho,bx,by,bz,dPdx)
c         call  fff1(x,v,sou,xl,rho,bx,by,bz,dxldx)
         call  bbb1(x,v,sou,xl,rho,bx,by,bz,dxldx)
c         call  fff1(x,v,sou,xl,dvdx,dadx,dxldx)
         call  c1(x,v,sou,xl,rho,bx,by,bz,drhodx)
         call  cc1(x,v,sou,xl,rho,bx,by,bz,dbydx)
         call  ccc1(x,v,sou,xl,rho,bx,by,bz,dbzdx)

         h3=dx*dvdx
         h13=dx*dPdx
         h133=dx*dxldx
	hc3=dx*drhodx
	hc13=dx*dbydx
	hc133=dx*dbzdx

         x=x3+dx
         v=v3+h3
         P=P3+h13
         xl=xl3+h133
	rho=rho3+hc3
	by=by3+hc13
	bz=bz3+hc133

	 call  f1(x,v,sou,xl,P,rho,bx,by,bz,dvdx,t9)
c         call  ff1(x,v,sou,xl,rho,bx,by,bz,dPdx)
         call  bb1(x,v,sou,xl,rho,bx,by,bz,dPdx)
c         call  fff1(x,v,sou,xl,rho,bx,by,bz,dxldx)
         call  bbb1(x,v,sou,xl,rho,bx,by,bz,dxldx)
c         call  fff1(x,v,sou,xl,dvdx,dadx,dxldx)
         call  c1(x,v,sou,xl,rho,bx,by,bz,drhodx)
         call  cc1(x,v,sou,xl,rho,bx,by,bz,dbydx)
         call  ccc1(x,v,sou,xl,rho,bx,by,bz,dbzdx)

         h4=dx*dvdx
         h14=dx*dPdx
         h144=dx*dxldx
c	ht4=dx*dtedx
	hc4=dx*drhodx
	hc14=dx*dbydx
	hc144=dx*dbzdx

         v4=v3+h1/6.+h2/3.+h3/3.+h4/6.
         P4=P3+h11/6.+h12/3.+h13/3.+h14/6.
         xl4=xl3+h111/6.+h122/3.+h133/3.+h144/6.
	rho4=rho3+hc1/6.+hc2/3.+hc3/3.+hc4/6.
	by4=by3+hc11/6.+hc12/3.+hc13/3.+hc14/6.
	bz4=bz3+hc111/6.+hc122/3.+hc133/3.+hc144/6.

        sou4=sqrt(gama*P4/rho4)

         x4=x3+dx

        bx4=abx/x4


      fx4=x4**3*(-4.*b**2+x4**3)**2/(-2.*b**3+b*x4**3*(-4.+3.*x4)+x4**8
     &*sqrt((5.*b**2+x4**3)*(2.*b**2+x4**3*(-2.+x4))**2/x4**16))**2
        xmac=v4/sou4 
	xlkep=sqrt(x4**3*fx4)
      write(2,9)x4,v4/sou4,v4,sou4,xl4,t9,rho4,P4,bx4,by4,bz4,xl4/xlkep

c	x4=xo
c	v4=v-1.d-7
cc	sou4=sou-1.d-1
c	sou4=sou-1.d-7
c	xl4=xlc+1.e-7
c	rho4=rho-1.d-24
c	P4=P-1.e-26
c	bx4=bx-1.d-23*l
c	by4=by+1.d-23*l
c	bz4=bz+1.d-23*l
	

	dx=ll*1.e-3
         do 40 j=1,900000000
         if (dx.gt.0.) dvdxo=dvdx
         if (dx.lt.0.) dvdxo=0.
         x3=x4
	sou3=sou4
	v3=v4
        xl3=xl4
	rho3=rho4
	P3=P4
	bx3=bx4
	by3=by4
	bz3=bz4

c	dx=1.e-5*x3
c	dx=1.e-5
         v=v3
         x=x3
	 sou=sou3
         xl=xl3
	rho=rho3
	P=P3
	bx=bx3
        by=by3
        bz=bz3

         call  b1(x,v,sou,xl,P,rho,bx,by,bz,dvdx,t9)
c	print*,x,v,sou,xl,P,rho,bx,by,bz,dvdx,dPdx,t9
         call  bb1(x,v,sou,xl,rho,bx,by,bz,dPdx)
c	print*,x,v,sou,xl,P,rho,bx,by,bz,dvdx,dPdx,t9
c	pause
c	print*,dxldx,xl
        call  bbb1(x,v,sou,xl,rho,bx,by,bz,dxldx)
c	print*,dxldx,xl
c	pause
         call  c1(x,v,sou,xl,rho,bx,by,bz,drhodx)
         call  cc1(x,v,sou,xl,rho,bx,by,bz,dbydx)
         call  ccc1(x,v,sou,xl,rho,bx,by,bz,dbzdx)
c	print*,dxldx,xl,dvdx,v,dPdx,P,drhodx,rho,dbydx,by,dbzdx,bz

         h1=dx*dvdx
	 h11=dx*dPdx
	 h111=dx*dxldx
	 hc1=dx*drhodx
	hc11=dx*dbydx
	hc111=dx*dbzdx

	ht1=dx*dtedx
         x=x3+dx/2.
         v=v3+h1/2.
	 P=P3+h11/2.
	 xl=xl3+h111/2.
	rho=rho3+hc1/2.
	by=by3+hc11/2.
	bz=bz3+hc111/2.

         call  b1(x,v,sou,xl,P,rho,bx,by,bz,dvdx,t9)
         call  bb1(x,v,sou,xl,rho,bx,by,bz,dPdx)
         call  bbb1(x,v,sou,xl,rho,bx,by,bz,dxldx)
         call  c1(x,v,sou,xl,rho,bx,by,bz,drhodx) 
	 call  cc1(x,v,sou,xl,rho,bx,by,bz,dbydx)
	 call  ccc1(x,v,sou,xl,rho,bx,by,bz,dbzdx)

         h2=dx*dvdx
         h12=dx*dPdx
         h122=dx*dxldx
	hc2=dx*drhodx
	hc12=dx*dbydx
	hc122=dx*dbzdx
	
         x=x3+dx/2.
         v=v3+h2/2.
         P=P3+h12/2.
         xl=xl3+h122/2.
         rho=rho3+hc2/2.
	 by=by3+hc12/2.
	 bz=bz3+hc122/2.

         call  b1(x,v,sou,xl,P,rho,bx,by,bz,dvdx,t9)
         call  bb1(x,v,sou,xl,rho,bx,by,bz,dPdx)
         call  bbb1(x,v,sou,xl,rho,bx,by,bz,dxldx)
         call  c1(x,v,sou,xl,rho,bx,by,bz,drhodx) 
	 call  cc1(x,v,sou,xl,rho,bx,by,bz,dbydx)
	 call  ccc1(x,v,sou,xl,rho,bx,by,bz,dbzdx)

         h3=dx*dvdx
         h13=dx*dPdx
         h133=dx*dxldx
	hc3=dx*drhodx
	hc13=dx*dbydx
	hc133=dx*dbzdx

         x=x3+dx
         v=v3+h3
         P=P3+h13
         xl=xl3+h133
	rho=rho3+hc3
	by=by3+hc13
	bz=bz3+hc133

         call  b1(x,v,sou,xl,P,rho,bx,by,bz,dvdx,t9)
         call  bb1(x,v,sou,xl,rho,bx,by,bz,dPdx)
         call  bbb1(x,v,sou,xl,rho,bx,by,bz,dxldx)
         call  c1(x,v,sou,xl,rho,bx,by,bz,drhodx) 
	 call  cc1(x,v,sou,xl,rho,bx,by,bz,dbydx)
	 call  ccc1(x,v,sou,xl,rho,bx,by,bz,dbzdx)

         h4=dx*dvdx
         h14=dx*dPdx
         h144=dx*dxldx
	hc4=dx*drhodx
	hc14=dx*dbydx
	hc144=dx*dbzdx

         v4=v3+h1/6.+h2/3.+h3/3.+h4/6.
         P4=P3+h11/6.+h12/3.+h13/3.+h14/6.
         xl4=xl3+h111/6.+h122/3.+h133/3.+h144/6.
	rho4=rho3+hc1/6.+hc2/3.+hc3/3.+hc4/6.
	by4=by3+hc11/6.+hc12/3.+hc13/3.+hc14/6.
	bz4=bz3+hc111/6.+hc122/3.+hc133/3.+hc144/6.

	sou4=sqrt(gama*P4/rho4)

         x4=x3+dx

	bx4=abx/x4

c	print*,dbydx,dbzdx,dvdx,bx4,by4,bz4,x4
c	pause

	if(mod(j-1,100).eq.0) then
      fx4=x4**3*(-4.*b**2+x4**3)**2/(-2.*b**3+b*x4**3*(-4.+3.*x4)+x4**8
     &*sqrt((5.*b**2+x4**3)*(2.*b**2+x4**3*(-2.+x4))**2/x4**16))**2
        xmac=v4/sou4
c	t9=sou4**2*xmui*xmp*9.e+20/xk/gama/1.e+9
	t9=ome*sou4**2*9.e+20*beta*1.e-9
	xlkep=sqrt(x4**3*fx4)

c	print*,xlkep

	xh4=sou4*sqrt(x4)/sqrt(fx4)
	eta=rho4*alfa*sou4*xh4
	ss=(P4+rho4)/t9
	etas=eta/ss*1.d+9*(6.67*2.d+25*xmas/3**3/1.d+30)

        write(2,9)x4,v4/sou4,v4,sou4,xl4,t9,rho4,P4,bx4,by4,bz4,
     &xl4/xlkep,etas,dsqrt(bx4**2+by4**2+bz4**2)

	call ccc1(x4,v4,sou4,xl4,rho4,bx4,by4,bz4,dbzdx)
	write(3,9)x4,1./4./pi/rho4*(bx4*dbzdx+0.5*bz4**2/xh4),xh4/x4*fx4
	call cc1(x4,v4,sou4,xl4,rho4,bx4,by4,bz4,dbydx)
         call  bb1(x4,v4,sou4,xl4,rho4,bx4,by4,bz4,dP4dx)
	call  c1(x4,v4,sou4,xl4,rho4,bx4,by4,bz4,drho4dx)
	write(4,9)x4,(bx4*dbydx+bx4*by4/x4+s2*bz4*by4/xh4)*x4/rho4/v4,
     &dbzdx,xh4,bx4,3.*(bx4**2+by4**2+bz4**2)*v4/16./pi/x4,v4/(g3-1.)
     &*(dP4dx-g1*p4/rho4*drho4dx)

c	print*,eta/ss*1.d+9*(6.67*2.d+25*xmas/3**3/1.d+30)
c	pause

	endif

	if(abs(xl4).gt.xlkep)stop
	if(abs(v4).ge.1.)stop
	if(abs(sou4).ge.1.)stop
c	if(x.ge.1000.0.or.x.le.(1.001+sqrt(1-b**2)))stop
	if(x.le.(1.001+sqrt(1-b**2)))stop
        if(dvdxo.gt.0..and.dvdx.lt.0.)dx=-dx
40	continue
     	end
       
      subroutine bisec(x1,x2,sou)
      implicit real*8(a-h,o-z)
      common/aaa/alfa,xf,gama,xmdot,xmu,xmui,xmue,xmp,xk,xo,an,xlc,beta
      common/ddd/b,g1,g3,xIfac,xInp1,xIn,xmas
	common/bbb/xmx,xmy,xmz,s1,s2,pi,xfm,l,ll
      do 10 i=1,1000000
      xm=0.5*(x1+x2)
c      write(*,*) x1,x2,f(x1),f(xm)
      if((f(x1)*f(xm)).gt.0.0) then
      x1=xm
      else
      x2=xm
      endif
      check=abs((x1-x2)/x1)
      if(check.lt.1.e-10) goto20
10    continue
20    continue
      sou=x1
      return
      end

      function f(a)
      implicit real*8(a-h,o-z)
      common/aaa/alfa,xf,gama,xmdot,xmu,xmui,xmue,xmp,xk,xo,an,xlc,beta
      common/ddd/b,g1,g3,xIfac,xInp1,xIn,xmas
	common/bbb/xmx,xmy,xmz,s1,s2,pi,xfm,l,ll
        amc=1./sqrt(2.)*sqrt((alfa*xf*gama**2*(-1 + g3)*l**2*xmx+ 
     &gama*(gama*(1.+g1)*l**2+2.*g1*xmx**2)*xmy+sqrt(gama**2*(4.*(-1.-g1
     &+alfa**2*xf*(-1.+g3))*xmx**2*xmy*(-alfa*xf*gama*(-1.+g3)*l**2*xmx+ 
     &2.*gama*g1*l**2*xmy+alfa**2*xf*(-1.+g3)*xmx**2*xmy)+(alfa*xf*gama*
     &(-1.+g3)*l**2*xmx+(gama*(1.+g1)*l**2+2.*g1*xmx**2)*xmy)**2)))/
     &(gama**2*(1.+g1-alfa**2*xf*(-1.+g3))*xmx**2*xmy))
	v=amc*a
	ome=xmu*xmp/gama/xk
	t=ome*a**2
	t9=ome*a**2*9.e+20*1.e-9*beta
        fx=(b**2-2.*b*dsqrt(xo)+xo**2)**2/xo**2/(xo**2-2.*xo+b
     &*dsqrt(xo))**2
	xh=a*dsqrt(xo)/dsqrt(fx)
	rho=xmdot/v/a/xo**1.5*sqrt(fx)/4./pi
	P=rho*a**2/gama
	dfxdx=(-3.*b**5-6.*b**4*(xo-3.)*sqrt(xo)-2.*xo**(11./2.)+4.*
     &b**3*xo*(5.*xo-9.)+3.*b*xo**3*(5.*xo-4.)-2.*b**2*xo**(3./2.)*
     &(4.*xo**2+5.*xo-12.))/xo**4/(b+(xo-2.)*sqrt(xo))**3     
	x=xo
	xl=xlc
	W=P+rho*v**2
	bx=sqrt(4.*acos(-1.)*rho)*a/xmx*l
	by=sqrt(4.*acos(-1.)*rho)*a/xmy*l
	bz=sqrt(4.*acos(-1.)*rho)*a/xmz*l
      f=(-1.+g3)*P*(-3.*l*x**2*(bx**2+by**2+bz**2)*v*(l**2*bx**2-4.
     &*pi*rho*v**2)-32.*alfa*xf*l**2*pi*x*bx**2*xl*W+1./fx*16.*alfa*xf
     &*pi*x**2*v*(l*bx*by*fx+l*s2*x*by*bz*fx-2.*alfa*pi*(x*dfxdx-5.*fx)
     &*(P+rho*v**2))*W+(4.*(l*x**2*(s1*x*bx*bz-(bx**2+by**2)*xh)+4.*pi
     &*xh*(-x**3*fx+xl**2)*rho)*v*((1.+g1)*P*(l**2*bx**2-4.*pi*rho*v**2
     &)+4.*alfa**2*xf*(-1.+g3)*pi*(P-rho*v**2)*W+4.*alfa**2*xf*(-1.+g3)
     &*pi*W**2))/((-1.+g3)*xh*P)-1./((-1.+g3)*fx)*16.*pi*x**2*(x*dfxdx 
     &-3.*fx)*v*(g1*l**2*bx**2*P-2.*pi*(2.*g1*P*rho*v**2+alfa**2*xf*
     &(-1.+g3)*(2.*rho*v**2-W)*W)))
        return
        end

         subroutine f1(x,v,a,xl,P,rho,bx,by,bz,dvdx,t9)
      implicit real*8(a-h,o-z)
      common/aaa/alfa,xf,gama,xmdot,xmu,xmui,xmue,xmp,xk,xo,an,xlc,beta
      common/ddd/b,g1,g3,xIfac,xInp1,xIn,xmas
	common/bbb/xmx,xmy,xmz,s1,s2,pi,xfm,l,ll
	ome=xmu*xmp/gama/xk
        fx=(b**2-2.*b*dsqrt(x)+x**2)**2/x**2/(x**2-2.*x+b
     &*dsqrt(x))**2
	xh=a*sqrt(x)/sqrt(fx)
c	print*,rho
	rho=xmdot/v/a/x**1.5*sqrt(fx)/4./pi
c	print*,rho
c	pause
        dfxdx=(-3.*b**5-6.*b**4*(x-3.)*sqrt(x)-2.*x**(11./2.)+4.*
     &b**3*x*(5.*x-9.)+3.*b*x**3*(5.*x-4.)-2.*b**2*x**(3./2.)*
     &(4.*x**2+5.*x-12.))/x**4/(b+(x-2.)*sqrt(x))**3
	d2f=(24.*b**6+12.*x**7+3.*b**5*sqrt(x)*(-64.+27.*x)+4.*b**4
     &*x*(144.-116.*x+21.*x**2)-3.*b*x**(7./2.)*(40.-74.*x+47.
     &*x**2)-2.*b**3*x**(3./2.)*(384.-425.*x+126.*x**2)+8.*b**2
     &*x**2*(48.-54*x+7.*x**2+10.*x**3))/2./x**5/(b+(x-2.)*
     &sqrt(x))**4

c Koushik you wirte here...........

	l1=0 ! magnetic part of radial Euler equation's switch
        bsq= bx**2+by**2+bz**2
	  oa= v*x*(2.*alfa*bx-by)
	  ob= alfa*x*(5.-1./g1)/(2.*rho)
        og= 8*v**2*pi*alfa*rho-l*bx*by
	  xW= alfa*(P+rho*(v**2))	
	  qp= (xf*xW*xl/x**2)+3.*l*xfm*v*bsq/(16.*pi*x)
	  oc= (g3-1.)*qp/(v*g1*P)*(x*xW/(2.*rho)-alfa*x*v**2)+2.5*xW/
     &rho+l*x/(4.*pi*rho)*(s2*by*bz/xh+by*bx/x)-0.5*x*xW*dfxdx/(rho       
     &*fx)-2.*xl*v/x
        xd= 4.*pi*rho*v**2-l*bx**2
	  xe= s1*bx*bz/xh-(by**2+bx**2)/x
c	print*,P,rho,a,v,xl,bx,by,bz
c	print*,oa,ob,xW,qp,oc,xd,xe
c	pause
       xk1=(2*xl*og/(x**2*xd)-l1*s1*bx*bz
     &/(4.*pi*v*rho*xh)-2.*l1*oa*by/(x**2*xd))/v**2+(g3-1.)/(2.*v
     &*g1*P)*(xf*(alfa*2.*v*rho*xl/x**2+xW*og/(x*xd)
     &)+3.*l*xfm/(16.*pi*x)*(bsq-2.*bz**2+8.*pi*rho
     &*v*oa*by/(x*xd)))-0.5*(g3-1.)*qp/(g1*P*v**2)-rho*((8.*pi*xl*
     &ob*rho*v/(x**3*xd)-l1*xe/(4.*pi*g1*rho*P)-l1*s1*bx*bz*(1.-1./g1)/     
     &(8.*pi*xh*P*rho)-2.*l1*ob*bx*by/(x**2*xd))/v+(xl**2/x**3-fx+l1*xe
     &/(4.*pi*rho))*(1./g1-1.)/(v*P)+(g3-1.)/(2.*g1*P)*(xf*alfa*xl*(g1+
     &3.)/(x**2*(1.+g1))+4.*pi*ob*xf*xW*rho*v/(x**2*xd)+1.5*ob*bx
     &*by*rho*v*l*xfm/(x**2*xd))+(g3-1.)*qp/(2.*g1*P**2))
        xk2=-(8.*pi*rho*xl*v*oc/(x**3*xd)+xl**2/x**4-dfxdx+l1*xe*(g3-1.
     &)*qp/(4.*pi*rho*g1*v*P)-l1/(4.*pi*rho)*(1.5*s1*bx*bz/(x*xh)+s1*bx
     &*bz/(2.*xh)*((g3-1.)*qp/(g1*P*v)-dfxdx/fx)+xe/x+8.*pi*oc*rho*bx
     &*by/(x**2*xd)-2.*bx**2/x**2))/v**2+(xl**2/x**3-fx+l1*xe/(4.*pi*
     &rho))*(g3-1.)*qp/(g1*P*v**3)+(g3-1.)/(2.*g1*v*P*x**2)*(4.*pi*xf* 
     &oc*rho*v*xW/xd-alfa*xf*xl*v*rho*(g3-1.)*qp/(g1*P)+3.*l*xfm/(16.*
     &pi)*(2.*v*by**2-3.*v*bsq+8.*pi*rho*v*oc*bx*by/xd))+rho*(xl**2/x**  
     &3-fx+l1*xe/(4.*pi*rho))*((8.*pi*xl*ob*rho*v/(x**3*xd)-l1*xe/(4.*
     &pi*g1*rho*P)-l1*s1*bx*bz*(1.-1./g1)/(8.*pi*xh*P*rho)-2.*l1*ob*bx
     &*by/(x**2*xd))/v**2+(xl**2/x**3-fx+l1*xe/(4.*pi*rho))*(1./g1-
     &1.)/(v**2*P)+(g3-1.)/(2.*g1*P*v)*(xf*alfa*xl*(g1+3.)/(x**2*(1.
     &+g1))+4.*pi*xf*ob*xW*rho*v/(x**2*xd)+1.5*ob*bx*by*rho*v*l*xfm/(x*
     &*2*xd))-(g3-1.)*qp/(2.*g1*v*P**2))+d2f/fx-(dfxdx/fx)**2+
     &1.5/x**2
        xk3=-2./v**2-rho*(1.-1./g1)/P
        xk4=(xl**2/x**3-fx+l1*xe/(4.*pi*rho))*rho*(1./g1-1.)/(v*P)+
     &(g3-1.)*qp/(g1*P*v**2) 
       	dvdx=(xk1-xk4+sqrt((xk1-xk4)**2+4.*xk3*xk2)/(2.*xk3))
c	print*,(xk1-xk4)**2,4.*xk3*xk2,xk1,xk2
c	pause



         return
         end

         subroutine ff1(x,v,a,xl,rho,bx,by,bz,dPdx)
      implicit real*8(a-h,o-z)
      common/aaa/alfa,xf,gama,xmdot,xmu,xmui,xmue,xmp,xk,xo,an,xlc,beta
      common/ddd/b,g1,g3,xIfac,xInp1,xIn,xmas
	common/bbb/xmx,xmy,xmz,s1,s2,pi,xfm,l,ll
	ome=xmu*xmp/gama/xk
        fx=(b**2-2.*b*sqrt(x)+x**2)**2/x**2/(x**2-2.*x+b
     &*sqrt(x))**2
	xh=a*sqrt(x)/sqrt(fx)
	rho=xmdot/v/a/x**1.5*sqrt(fx)/4./pi
        dfxdx=(-3.*b**5-6.*b**4*(x-3.)*sqrt(x)-2.*x**(11./2.)+4.*
     &b**3*x*(5.*x-9.)+3.*b*x**3*(5.*x-4.)-2.*b**2*x**(3./2.)*
     &(4.*x**2+5.*x-12.))/x**4/(b+(x-2.)*sqrt(x))**3
         return
         end

         subroutine fff1(x,v,a,xl,rho,bx,by,bz,dxldx)
      implicit real*8(a-h,o-z)
      common/aaa/alfa,xf,gama,xmdot,xmu,xmui,xmue,xmp,xk,xo,an,xlc,beta
      common/ddd/b,g1,g3,xIfac,xInp1,xIn,xmas
	common/bbb/xmx,xmy,xmz,s1,s2,pi,xfm,l,ll
	ome=xmu*xmp/gama/xk
        fx=(b**2-2.*b*sqrt(x)+x**2)**2/x**2/(x**2-2.*x+b
     &*sqrt(x))**2
	xh=a*sqrt(x)/sqrt(fx)
	rho=xmdot/v/a/x**1.5*sqrt(fx)/4./pi
        return
        end

         subroutine b1(x,v,a,xl,P,rho,bx,by,bz,dvdx,t9)
      implicit real*8(a-h,o-z)
      common/aaa/alfa,xf,gama,xmdot,xmu,xmui,xmue,xmp,xk,xo,an,xlc,beta
      common/ddd/b,g1,g3,xIfac,xInp1,xIn,xmas
	common/bbb/xmx,xmy,xmz,s1,s2,pi,xfm,l,ll
	ome=xmu*xmp/gama/xk
	t=ome*a**2
	t9=ome*a**2*9.e+20*1.e-9*beta
        fx=(b**2-2.*b*sqrt(x)+x**2)**2/x**2/(x**2-2.*x+b
     &*sqrt(x))**2
	xh=a*sqrt(x)/sqrt(fx)
        dfxdx=(-3.*b**5-6.*b**4*(x-3.)*sqrt(x)-2.*x**(11./2.)+4.*
     &b**3*x*(5.*x-9.)+3.*b*x**3*(5.*x-4.)-2.*b**2*x**(3./2.)*
     &(4.*x**2+5.*x-12.))/x**4/(b+(x-2.)*sqrt(x))**3
	bt=sqrt(bx**2+by**2+bz**2)
	W=P+rho*v**2
	hh=3./2./x-dfxdx/2./fx
	aa=xl**2/x**3-fx-1./(4.*pi*rho)*((bx**2+by**2)/x-bz*bx*s1/xh)
	xj=hh+rho*aa/2./P
	xxk=v*rho/2./P-1./v
	if(l.eq.1)then
        f1=(-1.+g3)*P*(-3.*l*x**2*(bx**2+by**2+bz**2)*v*(l**2*bx**2-4.
     &*pi*rho*v**2)-32.*alfa*xf*l**2*pi*x*bx**2*xl*W+1./fx*16.*alfa*xf
     &*pi*x**2*v*(l*bx*by*fx+l*s2*x*by*bz*fx-2.*alfa*pi*(x*dfxdx-5.*fx)
     &*(P+rho*v**2))*W+(4.*(l*x**2*(s1*x*bx*bz-(bx**2+by**2)*xh)+4.*pi
     &*xh*(-x**3*fx+xl**2)*rho)*v*((1.+g1)*P*(l**2*bx**2-4.*pi*rho*v**2
     &)+4.*alfa**2*xf*(-1.+g3)*pi*(P-rho*v**2)*W+4.*alfa**2*xf*(-1.+g3)
     &*pi*W**2))/((-1.+g3)*xh*P)-1./((-1.+g3)*fx)*16.*pi*x**2*(x*dfxdx
     &-3.*fx)*v*(g1*l**2*bx**2*P-2.*pi*(2.*g1*P*rho*v**2+alfa**2*xf*
     &(-1.+g3)*(2.*rho*v**2-W)*W)))
        f2=16.*pi*x**3*(P*(2.*g1*P-(1.+g1)*rho*v**2)*(-l**2*bx**2+ 
     &4.*pi*rho*v**2)+alfa*xf*(-1.+g3)*(l**2*bx*by*P+4.*alfa*pi*rho*v**2
     &*(P-rho*v**2))*W-4.*alfa**2*xf*(-1.+g3)*pi*(P-rho*v**2)*W**2)
	dvdx=f1/f2
	else
	f1=v*aa*rho+2.*P*v*g1*xj-xf*alfa*W*xl*(g3-1.)/x**3-3.*l*v*bt**2
     &/16./pi/x*(g3-1.)*xfm
	f2=v**2*rho+2.*P*v*g1*xxk
	dvdx=f1/f2
c	write(*,*)'dv'
	endif
c	write(*,*)f1,v*bt**2/16./pi/x*(g3-1.)*xfm,v*aa*rho,P*v*g1*xj
c	pause
         return
         end

         subroutine bb1(x,v,a,xl,rho,bx,by,bz,dPdx)
      implicit real*8(a-h,o-z)
      common/aaa/alfa,xf,gama,xmdot,xmu,xmui,xmue,xmp,xk,xo,an,xlc,beta
      common/ddd/b,g1,g3,xIfac,xInp1,xIn,xmas
	common/bbb/xmx,xmy,xmz,s1,s2,pi,xfm,l,ll
	ome=xmu*xmp/gama/xk
        fx=(b**2-2.*b*sqrt(x)+x**2)**2/x**2/(x**2-2.*x+b
     &*sqrt(x))**2
	P=rho*a**2/gama
	xh=a*sqrt(x)/sqrt(fx)
        dfxdx=(-3.*b**5-6.*b**4*(x-3.)*sqrt(x)-2.*x**(11./2.)+4.*
     &b**3*x*(5.*x-9.)+3.*b*x**3*(5.*x-4.)-2.*b**2*x**(3./2.)*
     &(4.*x**2+5.*x-12.))/x**4/(b+(x-2.)*sqrt(x))**3
	aa=xl**2/x**3-fx-1./(4.*pi*rho)*((bx**2+by**2)/x-bz*bx*s1/xh)
	call b1(x,v,a,xl,P,rho,bx,by,bz,dvdx,t9)
        dPdx=(aa-v*dvdx)*rho
         return
         end

         subroutine bbb1(x,v,a,xl,rho,bx,by,bz,dxldx)
      implicit real*8(a-h,o-z)
      common/aaa/alfa,xf,gama,xmdot,xmu,xmui,xmue,xmp,xk,xo,an,xlc,beta
      common/ddd/b,g1,g3,xIfac,xInp1,xIn,xmas
	common/bbb/xmx,xmy,xmz,s1,s2,pi,xfm,l,ll
	ome=xmu*xmp/gama/xk
        fx=(b**2-2.*b*sqrt(x)+x**2)**2/x**2/(x**2-2.*x+b
     &*sqrt(x))**2
        dfxdx=(-3.*b**5-6.*b**4*(x-3.)*sqrt(x)-2.*x**(11./2.)+4.*
     &b**3*x*(5.*x-9.)+3.*b*x**3*(5.*x-4.)-2.*b**2*x**(3./2.)*
     &(4.*x**2+5.*x-12.))/x**4/(b+(x-2.)*sqrt(x))**3
	xh=a*sqrt(x)/sqrt(fx)
	P=rho*a**2/gama
	W=P+rho*v**2
	bt=sqrt(bx**2+by**2+bz**2)
	aa=xl**2/x**3-fx-1./(4.*pi*rho)*((bx**2+by**2)/x-bz*bx*s1/xh)
	hh=3./2./x-dfxdx/2./fx
	xj=hh+rho*aa/2./P
	xxk=v*rho/2./P-1./v
	xx=alfa*W/rho/x+alfa*W*hh/rho+l*(bx*by/4./pi/rho/x+bz/4./pi/rho
     &*by*s2)
	xxl=x/xf/alfa/W*(3.*l*v/2./x*bt**2/8./pi-v*aa*rho/(g3-1.)
     &-2.*P*v*g1*xj/(g3-1.))
	xxn=x/xf/W/alfa*(2.*P*v*g1*xxk/(g3-1.)+v**2*rho/(g3-1.))
	xe=l*(by+bx*xxn/x)
	zz=(2.*xl*bx/x**2+bx*xxl/x)*l
	T1=4.*pi*rho*(2.*alfa*(v**2-W/2./rho)*xj-alfa*aa*(1.+W/2./P)
     &-xx)+2.*xl*bx**2*l/x**2/v
	T2=4.*pi*rho*(alfa*v*(1.+W/2./P)-2.*alfa*(v**2-W/2./rho)*xxk
     &-2.*alfa*v)+bx*by*l/v
	call b1(x,v,a,xl,P,rho,bx,by,bz,dvdx,t9)
	if(l.eq.0)then
        dxldx=x/v*((alfa*(1.+W/2./P)*aa-2.*alfa*(v**2-W/2./rho)*xj-
     &l*bx*zz/4./pi/rho/v+xx)+dvdx*(2.*alfa*v-alfa*v*(1.+W/2./P)+2.*
     &xxk*alfa*(v**2-W/2./rho)-l*bx/4./pi/rho*xe/v))
	else
	dxldx=(T1+T2*dvdx)/(l*bx**2/v/x-v/x*4*pi*rho)
c	write(*,*)'dxl'
	endif
        return
        end

         subroutine c1(x,v,a,xl,rho,bx,by,bz,drhodx)
      implicit real*8(a-h,o-z)
      common/aaa/alfa,xf,gama,xmdot,xmu,xmui,xmue,xmp,xk,xo,an,xlc,beta
      common/ddd/b,g1,g3,xIfac,xInp1,xIn,xmas
	common/bbb/xmx,xmy,xmz,s1,s2,pi,xfm,l,ll
	ome=xmu*xmp/gama/xk
        fx=(b**2-2.*b*sqrt(x)+x**2)**2/x**2/(x**2-2.*x+b
     &*sqrt(x))**2
	xh=a*sqrt(x)/sqrt(fx)
        dfxdx=(-3.*b**5-6.*b**4*(x-3.)*sqrt(x)-2.*x**(11./2.)+4.*
     &b**3*x*(5.*x-9.)+3.*b*x**3*(5.*x-4.)-2.*b**2*x**(3./2.)*
     &(4.*x**2+5.*x-12.))/x**4/(b+(x-2.)*sqrt(x))**3
	P=rho*a**2/gama
	W=P+rho*v**2
	bt=sqrt(bx**2+by**2+bz**2)
	aa=xl**2/x**3-fx-1./(4.*pi*rho)*((bx**2+by**2)/x-bz*bx*s1/xh)
	hh=3./2./x-dfxdx/2./fx
	xj=hh+rho*aa/2./P
	xxk=v*rho/2./P-1./v
	call b1(x,v,a,xl,P,rho,bx,by,bz,dvdx,t9)
        drhodx=(-xj+xxk*dvdx)*2.*rho
         return
         end

         subroutine cc1(x,v,a,xl,rho,bx,by,bz,dbydx)
      implicit real*8(a-h,o-z)
      common/aaa/alfa,xf,gama,xmdot,xmu,xmui,xmue,xmp,xk,xo,an,xlc,beta
      common/ddd/b,g1,g3,xIfac,xInp1,xIn,xmas
	common/bbb/xmx,xmy,xmz,s1,s2,pi,xfm,l,ll
	ome=xmu*xmp/gama/xk
        fx=(b**2-2.*b*sqrt(x)+x**2)**2/x**2/(x**2-2.*x+b
     &*sqrt(x))**2
	xh=a*sqrt(x)/sqrt(fx)
        dfxdx=(-3.*b**5-6.*b**4*(x-3.)*sqrt(x)-2.*x**(11./2.)+4.*
     &b**3*x*(5.*x-9.)+3.*b*x**3*(5.*x-4.)-2.*b**2*x**(3./2.)*
     &(4.*x**2+5.*x-12.))/x**4/(b+(x-2.)*sqrt(x))**3
	P=rho*a**2/gama
	W=P+rho*v**2
	bt=sqrt(bx**2+by**2+bz**2)
	aa=xl**2/x**3-fx-1./(4.*pi*rho)*((bx**2+by**2)/x-bz*bx*s1/xh)
	hh=3./2./x-dfxdx/2./fx
	xxk=v*rho/2./P-1./v
	xj=hh+rho*aa/2./P
	xxl=x/xf/alfa/W*(xf*3.*l*v/2./x*bt**2/8./pi-v*aa*rho/(g3-1.)
     &-2.*P*v*g1*xj/(g3-1.))
	xxn=x/xf/W/alfa*(2.*P*v*g1*xxk/(g3-1.)+v**2*rho/(g3-1.))
	xe=l*(by+bx*xxn/x)
	zz=(2.*xl*bx/x**2+bx*xxl/x)*l
	xx=alfa*W/rho/x+alfa*W*hh/rho+l*(bx*by/4./pi/rho/x+bz/4./pi/rho
     &*by*s2)
	T1=4.*pi*rho*(2.*alfa*(v**2-W/2./rho)*xj-alfa*aa*(1.+W/2./P)
     &-xx)+2.*xl*bx**2*l/x**2/v
	T2=4.*pi*rho*(alfa*v*(1.+W/2./P)-2.*alfa*(v**2-W/2./rho)*xxk
     &-2.*alfa*v)+bx*by*l/v
	call b1(x,v,a,xl,P,rho,bx,by,bz,dvdx,t9)
	if(l.eq.0)then
        dbydx=-(zz/v+xe/v*dvdx)
	else
	dbydx=l*bx*T1/v/x/(l*bx**2/v/x-v*4.*pi*rho/x)-l*2.*xl*bx/v/x**2
     &+dvdx*(l*bx/v/x*T2/(l*bx**2/v/x-v/x*4.*pi*rho)-l*by/v)
c	dbydx=0.d-30
c	write(*,*)'dby'
	endif
         return
         end

         subroutine ccc1(x,v,a,xl,rho,bx,by,bz,dbzdx)
      implicit real*8(a-h,o-z)
      common/aaa/alfa,xf,gama,xmdot,xmu,xmui,xmue,xmp,xk,xo,an,xlc,beta
      common/ddd/b,g1,g3,xIfac,xInp1,xIn,xmas
	common/bbb/xmx,xmy,xmz,s1,s2,pi,xfm,l,ll
	P=rho*a**2/gama
	call b1(x,v,a,xl,P,rho,bx,by,bz,dvdx,t9)
        dbzdx=-l*(bz/x+bz/v*dvdx)
c        dbzdx=0.d-30
         return
         end



