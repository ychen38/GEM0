MODULE gem_equil
  IMPLICIT NONE
  integer :: itube,ibase,iperi,iperidf,ibunit,icandy=1,isprime=0,ildu=0,eldu=0
  real :: mimp=2,mcmp=12,chgi=1,chgc=6
  real :: elon0=1.0,tria0=0.0,rmaj0=500.0,r0,a=180.0,selon0=0.0,&
             stria0=0.0,rmaj0p=-0.0,q0p=0.006,q0=1.4, elonp0=0.,triap0=0.,erp=0.01,er0=0.,q0abs
  real :: beta,Rovera,shat0,teti,tcti,rhoia,Rovlni,Rovlti,Rovlne,Rovlte,Rovlnc,Rovltc,ncne,nuacs
  real :: gamma_E,mach
  real :: f0, f0p,bunit
  real :: rin,rout,dr,dth,delz,jacmax,eadj
  real :: cn0e,cn0i,cn0b,cn0c,n0emax,n0imax,n0bmax,n0cmax
  real :: r0a,lxa,lymult,delra,delri,delre,delrn,rina,routa,betai, &
               tir0,xnir0,xu,frequ,vu,eru

  integer :: nr=256,nr2=150,ntheta=100,idiag=0
  real,dimension(:,:),allocatable :: bfld,qhat,radius,gr,gth,grdgt,grcgt, &
                                        gxdgy,dydr,dbdr,dbdth,dqhdr,jacob, &
                                        yfn,hght,thflx
  real,dimension(:),allocatable :: rmaj,rmajp,elon,selon,tria,stria, psi,&
                                      f,psip,sf,jacoba,jfn,zfnth,thfnz,&
                                      t0i,t0e,t0b,t0c,t0ip,t0ep,t0bp,t0cp,&
                                      xn0i,xn0e,xn0c,xn0b,xn0ip,xn0ep,xn0bp,&
                                      xn0cp,vpari,vparc,vparb,&
                                      vparip,vparcp,vparbp, &
                                      capti,capte,captb,captc,capni,capne,&
                                      capnb,capnc,zeff,nue0,phinc,phincp,&
                                      er,upari,&
                                      dldth,sinu,cosu,dudl,dzdl,bps,&
                                      grr,grz,gtr,gtz, &
                                      grdgl,grdgrho,gtdgl,gtdgrho, &
                                      dldr,dldt,drhdr,drhdt,dbdl,dbdrho, &
                                      db2dl,db2drho,dbpsdl,dipdr, &
                                      rdtemp
!for Miller local flux-tube
  real :: candyf0p
  real,dimension(:),allocatable :: candyd0,candyd1,candyd2,candynus,candynu1,candydr

!for including bstar effects
  real,dimension(:),allocatable :: psip2
  real,dimension(:,:),allocatable :: curvbz,srbr,srbz,thbr,thbz,prsrbr,prsrbz,pthsrbr,pthsrbz,bdcrvb

  real,dimension(:,:),allocatable :: t0s,xn0s,capts,capns,vpars,vparsp
  real,dimension(:),allocatable :: cn0s,n0smax,tgis
  real :: tge
  character(len=32) :: trflnm ! profile-data-file name
!  real,external :: erf

contains
  subroutine new_equil()
      
      real :: pi,r,th,s,ss,c1,c2,lti,s0,s1,delsi,lte,delse,ln,deln
      parameter(c1=0.43236,c2=2.33528,lti=144.9,s1=0.5,delsi=0.6, &
                lte=144.9,delse=0.6,deln=0.6,ln=454.5)
      integer :: i,j,k
      real :: dum,x,denom,dum1,dum2,dum3,dum4,xp,elonp,triap
      real :: t0i0,t0e0,t0c0,ni0,nc0,ne0,t0i0p,t0e0p,t0c0p,ni0p,nc0p,ne0p
      real :: dr24,cutoff,cutoffti,cutoffte,cutoffi,cutoffe,cutoffn
      character(len=100) :: header,header2

      real :: e=1.6e-19,eps0=8.85e-12,me=0.9e-30,vte,pprime

      allocate(bfld(0:nr,0:ntheta),qhat(0:nr,0:ntheta),radius(0:nr,0:ntheta), &
               gr(0:nr,0:ntheta),gth(0:nr,0:ntheta),grdgt(0:nr,0:ntheta), &
               grcgt(0:nr,0:ntheta),gxdgy(0:nr,0:ntheta),dydr(0:nr,0:ntheta),&
               dbdr(0:nr,0:ntheta),dbdth(0:nr,0:ntheta),dqhdr(0:nr,0:ntheta),&
               jacob(0:nr,0:ntheta), jfn(0:ntheta), zfnth(0:ntheta),&
               thfnz(0:ntheta),yfn(0:nr,0:ntheta),hght(0:nr,0:ntheta),thflx(0:nr,0:ntheta))
               
      allocate(rmaj(0:nr), elon(0:nr), tria(0:nr), sf(0:nr), psi(0:nr), &
               rmajp(0:nr),selon(0:nr),stria(0:nr),psip(0:nr),&
               f(0:nr),jacoba(0:nr),t0i(0:nr),t0e(0:nr),t0b(0:nr),t0c(0:nr),&
               t0ip(0:nr),t0ep(0:nr),t0bp(0:nr),t0cp(0:nr),xn0i(0:nr),xn0e(0:nr),&
               xn0b(0:nr),xn0c(0:nr),xn0ip(0:nr),xn0ep(0:nr),xn0bp(0:nr),xn0cp(0:nr),&
               capti(0:nr),capte(0:nr),captb(0:nr),captc(0:nr),capni(0:nr),&
               capne(0:nr),capnb(0:nr),capnc(0:nr),zeff(0:nr),nue0(0:nr),&
               vpari(0:nr),vparc(0:nr),vparb(0:nr),phinc(0:nr), &
               vparip(0:nr),vparcp(0:nr),vparbp(0:nr),phincp(0:nr),er(0:nr), &
               upari(0:nr),dipdr(0:nr),rdtemp(0:nr))

      allocate(curvbz(0:nr,0:ntheta),srbr(0:nr,0:ntheta),srbz(0:nr,0:ntheta),&
               thbr(0:nr,0:ntheta),thbz(0:nr,0:ntheta), psip2(0:nr),bdcrvb(0:nr,0:ntheta),&
               prsrbr(0:nr,0:ntheta),prsrbz(0:nr,0:ntheta), &
               pthsrbr(0:nr,0:ntheta),pthsrbz(0:nr,0:ntheta))


      allocate(cn0s(1:5),n0smax(1:5),t0s(1:5,0:nr),xn0s(1:5,0:nr),&
               capts(1:5,0:nr),capns(1:5,0:nr),vpars(1:5,0:nr),&
               vparsp(1:5,0:nr),tgis(1:5))
      allocate(dldth(0:ntheta),sinu(0:ntheta),cosu(0:ntheta), &
               dudl(0:ntheta),dzdl(0:ntheta),bps(0:ntheta), &
               grr(0:ntheta),grz(0:ntheta),gtr(0:ntheta),gtz(0:ntheta))
      allocate(grdgl(0:ntheta),grdgrho(0:ntheta), &
               gtdgl(0:ntheta),gtdgrho(0:ntheta))
      allocate(dldr(0:ntheta),dldt(0:ntheta),drhdr(0:ntheta),drhdt(0:ntheta))
      allocate(dbdl(0:ntheta),dbdrho(0:ntheta))
      allocate(db2dl(0:ntheta),db2drho(0:ntheta),dbpsdl(0:ntheta))
      allocate(candyd0(0:ntheta),candyd1(0:ntheta),candyd2(0:ntheta),candynus(0:ntheta),candynu1(0:ntheta),candydr(0:ntheta))
     
     
      !write(*,*),'Rmaj0=',Rmaj0,'Rovera=',Rovera,'shat0=',shat0,'rmaj0p',rmaj0p
! if analytical equilibrium
      a = Rmaj0/Rovera
      r0=r0a*a
      q0p = shat0*q0/r0
      rin=rina*a
      rout=routa*a
      if(itube==1)then
         rin = r0-0.00001
         rout = r0+0.00001
      end if
      pi = atan(1.0)*4.
      dr = (rout-rin)/nr
      dth = pi*2/ntheta
      if(idiag==1)write(*,*)pi,dr,dth,nr,ntheta

! specify f(r),q(r)
      do i = 0,nr
         r = rin+i*dr
         f(i) = rmaj0
         s = r/a
         ss = s*s
         !sf(i) = q0+(r-r0)*q0p
         sf(i)=0.86-0.16*s+2.52*ss
      end do

! specify rmaj(r)
      do i = 0,nr
         r = rin+i*dr
         s=r/a
         rmaj(i) = rmaj0+(r-r0)*rmaj0p
         rmajp(i) = rmaj0p
      end do

! specify elon(r) and compute selon(r)
      elonp0 = selon0*elon0/r0
      do i = 0,nr
         r = rin+i*dr
         s=r/a
         elon(i) = elon0+(r-r0)*elonp0
         selon(i) = r*elonp0/elon(i)
      end do
      triap0 = stria0/r0 !*sqrt(1-tria0**2)/r0
      do i = 0,nr
         r = rin+i*dr
         s=r/a
         tria(i) = tria0+(r-r0)*triap0
         stria(i) = r*triap0/sqrt(1-tria(i)**2)
      end do

200   continue
! compute radius(r,theta)
      do i = 0,nr
         r = rin+i*dr
         do j = 0,ntheta
            th = -pi+dth*j
            radius(i,j) = rmaj(i)+r*cos(th+asin(tria(i))*sin(th))
            hght(i,j) = r*elon(i)*sin(th)
         end do
      end do

!define sinu, R_c on r0 surface ! for calculating of f^prime
      r = r0
      x = asin(tria0)
      xp = stria0/r0
      elonp = elonp0
      triap = triap0
      do j = 0,ntheta
         th = -pi+j*dth
         dum1 = elon0*r*cos(th) !dZ d theta
         dum2 = -r*sin(th+x*sin(th))*(1+x*cos(th)) !dR d theta
         dum3 = elonp*r*sin(th)+elon0*sin(th) !dZ d r
         dum4 = rmaj0p+cos(th+x*sin(th))-r*sin(th+x*sin(th))* &
                xp*sin(th)             !dR d r
         denom = dum4*dum1-dum3*dum2
         grr(j) = dum1/denom
         grz(j) = -dum2/denom
         gtr(j) = -dum3/denom
         gtz(j) = dum4/denom

         gr(nr2,j) = sqrt(dum1**2+dum2**2)/denom
         gth(nr2,j) = sqrt(dum3**2+dum4**2)/denom
         grdgt(nr2,j) = (-dum1*dum3-dum2*dum4)/denom**2
         grcgt(nr2,j) = 1/denom

         dldth(j) = sqrt(dum1**2+dum2**2) !Miller paper
         sinu(j) = dum1/dldth(j)
         cosu(j) = dum2/dldth(j)
         dzdl(j) = dum1/dldth(j)
         grdgl(j) = grr(j)*cosu(j)+grz(j)*sinu(j)
         grdgrho(j) = grr(j)*sinu(j)-grz(j)*cosu(j)
         gtdgl(j) = gtr(j)*cosu(j)+gtz(j)*sinu(j)
         gtdgrho(j) = gtr(j)*sinu(j)-gtz(j)*cosu(j)

         dldt(j) = 1./gtdgl(j) !"happens" to be equal to dldth
         dldr(j) = -dldt(j)*gtdgrho(j)/grdgrho(j) !verified to be the same as below
         drhdt(j) = 0. !verified with shaped parameters
         drhdr(j) = 1./grdgrho(j) !obvious if drhdt=0

         drhdt(j) = (grdgrho(j)*grdgt(nr2,j)-gtdgrho(j)*gr(nr2,j)**2)/(grdgt(nr2,j)**2-gr(nr2,j)**2*gth(nr2,j)**2)
         drhdr(j) = (grdgrho(j)*gth(nr2,j)**2-gtdgrho(j)*grdgt(nr2,j))/(gr(nr2,j)**2*gth(nr2,j)**2-grdgt(nr2,j)**2)

         dldt(j) = (grdgl(j)*grdgt(nr2,j)-gtdgl(j)*gr(nr2,j)**2)/(grdgt(nr2,j)**2-gr(nr2,j)**2*gth(nr2,j)**2)
         dldr(j) = (grdgl(j)*gth(nr2,j)**2-gtdgl(j)*grdgt(nr2,j))/(gr(nr2,j)**2*gth(nr2,j)**2-grdgt(nr2,j)**2)

      end do
      do j = 1,ntheta-1
         dudl(j) = 1/(cosu(j)+1.e-8)*(dzdl(j+1)-dzdl(j-1))/(2*dth)/dldth(j)
      end do
      dudl(0) = 1/(cosu(0)+1.e-8)*(dzdl(1)-dzdl(ntheta-1))/(2*dth)/dldth(0)
      dudl(ntheta) = dudl(0)

! compute grad r,grad theta, grdgt,grcgt
      do i = 0,nr
         r = rin+i*dr
         x = asin(tria(i))
         xp = stria(i)/r
         elonp = selon(i)*elon(i)/r
         triap = stria(i)*sqrt(1-tria(i)**2)/r
         do j = 0,ntheta
            th = -pi+j*dth
            dum1 = elon(i)*r*cos(th)                   !dZ d theta
            dum2 = -r*sin(th+x*sin(th))*(1+x*cos(th))  !dR d theta
            dum3 = elonp*r*sin(th)+elon(i)*sin(th)     !dZ d r 
            dum4 = rmajp(i)+cos(th+x*sin(th))-r*sin(th+x*sin(th))* & 
                    xp*sin(th)                         !dR d r
            denom = dum4*dum1-dum3*dum2
            srbr(i,j) = dum1/denom                     !component of (grad r) along R  
            srbz(i,j) = -dum2/denom                    !component of (grad r) along Z
            thbr(i,j) = -dum3/denom                    !component of (grad theta) along R
            thbz(i,j) = dum4/denom                     !component of (grad theta) along Z
            if(denom<0)write(*,*)elonp,triap,xp,denom
            gr(i,j) = sqrt(dum1**2+dum2**2)/denom
            gth(i,j) = sqrt(dum3**2+dum4**2)/denom
            grdgt(i,j) = (-dum1*dum3-dum2*dum4)/denom**2
            grcgt(i,j) = 1/denom
         end do
      end do

!!!!! ADD f^prime from Eq.(21) of Miller eq, for flux tube ONLY
      f0 = rmaj0
! compute psip(r0) on r0
      dum = 0.
      do j = 0,ntheta-1
         dum = dum+dth/radius(nr2,j)/grcgt(nr2,j)
      end do
      psip(nr2) = f0/2/pi/sf(nr2)*dum
      bunit = q0/r0*psip(nr2)

!convert input betai and rhostar into definition with GEM's B0
!This is the right position because pprime is needed for f0p
      if(ibunit==1)then
         betai = betai*bunit**2
         rhoia = rhoia*bunit
      end if
      !t0e0 = mimp*(rhoia*a*chgi/mimp)**2
      t0e0=1.0
      t0i0 = t0e0/teti
      t0c0 = t0i0*tcti
      ne0 = 1.0 !n_u
      nc0 = ncne*ne0
      ni0 = (ne0-chgc*nc0)/chgi
      t0i0p = -t0i0*Rovlti/rmaj0
      t0e0p = -t0e0*Rovlte/rmaj0
      t0c0p = -t0c0*Rovltc/rmaj0
      ni0p = -ni0*Rovlni/rmaj0
      ne0p = -ne0*Rovlne/rmaj0
      nc0p = -nc0*Rovlnc/rmaj0 !(ne0p-chgi*ni0p)/chgc !-nc0*Rovlnc/rmaj0
      betai = betai/(ne0*t0e0)/2
      nuacs = nuacs*sqrt(t0e0/mimp)/a

      do i = 0,nr
         r = rin+i*dr
         s=r/a
         xn0i(i) = ni0*exp(-Rovlni*0.36*0.3*tanh((s-0.5)/0.3))
         capni(i) = -ni0p/ni0*(1-tanh((s-0.5)/0.3)*tanh((s-0.5)/0.3))
         xn0e(i) = ne0*exp(-Rovlne*0.36*0.3*tanh((s-0.5)/0.3))
         capne(i) = -ne0p/ne0*(1-tanh((s-0.5)/0.3)*tanh((s-0.5)/0.3))
         xn0c(i) = nc0
         capnc(i) = -nc0p/nc0
         t0i(i) = t0i0*exp(-Rovlti*0.36*0.3*tanh((s-0.5)/0.3))
         capti(i) = -t0i0p/t0i0*(1-tanh((s-0.5)/0.3)*tanh((s-0.5)/0.3))
         t0e(i) = t0e0*exp(-Rovlte*0.36*0.3*tanh((s-0.5)/0.3))
         capte(i) = -t0e0p/t0e0*(1-tanh((s-0.5)/0.3)*tanh((s-0.5)/0.3))
         t0c(i) = t0c0
         captc(i) = -t0c0p/t0c0
         phincp(i) = 0.005*sin((r-rin)/(rout-rin)*2*pi)
         nue0(i) = 1.
         zeff(i) = 1.
      end do
      q0abs = q0
      
! compute B_p on r0
      do j = 0,ntheta
         bps(j) = psip(nr2)/radius(nr2,j)*gr(nr2,j)
      end do
!compute dbpsdl
      do j = 1,ntheta-1
         dbpsdl(j) = (bps(j+1)-bps(j-1))/(2*dth)/dldth(j)
      end do
      dbpsdl(0) = (bps(1)-bps(ntheta-1))/(2*dth)/dldth(j)
      dbpsdl(ntheta) = dbpsdl(0)

!compute term1 in (21) 
      dum1 = 0.
      do j = 0,ntheta-1
         dum1 = dum1+dth*dldth(j)/radius(nr2,j)**3/bps(j)**2*2*dudl(j)
      end do
      dum1 = dum1*f0/(2*pi)

!compute term2 in (21)
      dum2 = 0.
      do j = 0,ntheta-1
         dum2 = dum2+dth*dldth(j)/radius(nr2,j)**3/bps(j)**2*(-2)*sinu(j) &
                     /radius(nr2,j)
      end do
      dum2 = dum2*f0/(2*pi)

!compute term3 in (21) ... WWan
      pprime =  (t0i0p*ni0+t0i0*ni0p + t0e0p*ne0+t0e0*ne0p)/psip(nr2)*isprime
      dum3 = 0.
      do j = 0,ntheta-1
         dum3 = dum3+dth*dldth(j)/radius(nr2,j)**3/bps(j)**2*betai*radius(nr2,j)/bps(j)* &
              (pprime)
      enddo
      dum3 = dum3*f0/(2*pi)
!compute term4 in (21)
      dum4 = 0.
      do j = 0,ntheta-1
         dum4 = dum4+dth*dldth(j)/radius(nr2,j)**3/bps(j)**2*f0 &
                     /(radius(nr2,j)*bps(j))
      end do
      dum4 = dum4*f0/(2*pi)

!      f0p = psip(nr2)*(q0p/psip(nr2)-dum1-dum2)/(q0/f0+dum4)
      f0p = psip(nr2)*(q0p/psip(nr2)-dum1-dum2-dum3)/(q0/f0+dum4)!*0d0
!      write(*,*) 'f0p = ', f0p

!compute page 10 of Candy09
      do j = 0,ntheta
         bfld(nr2,j) = sqrt((f0/radius(nr2,j))**2+(psip(nr2)/radius(nr2,j)*gr(nr2,j))**2)
      end do
      candynus(ntheta/2) = 0
      candyd0(ntheta/2) = 0
      candyd1(ntheta/2) = 0
      candyd2(ntheta/2) = 0
      do j = ntheta/2+1,ntheta
         candyd0(j) = candyd0(j-1)+dth*dldth(j)/(radius(nr2,j)**2*bps(j))*(dudl(j)/(radius(nr2,j)*bps(j))-sinu(j)/(radius(nr2,j)**2*bps(j)))*f0 &
                                +dth*dldth(j-1)/(radius(nr2,j-1)**2*bps(j-1))*(dudl(j-1)/(radius(nr2,j-1)*bps(j-1))-sinu(j-1)/(radius(nr2,j-1)**2*bps(j-1)))*f0
         candyd1(j) = candyd1(j-1)+0.5*dth*dldth(j)/(radius(nr2,j)**2*bps(j))*bfld(nr2,j)**2/(bps(j)**2*f0) &
                                +0.5*dth*dldth(j-1)/(radius(nr2,j-1)**2*bps(j-1))*bfld(nr2,j-1)**2/(bps(j-1)**2*f0)
         candyd2(j) = candyd2(j-1)+0.5*dth*dldth(j)/(radius(nr2,j)**2*bps(j))*betai*f0/bps(j)**2 &
                                  +0.5*dth*dldth(j-1)/(radius(nr2,j-1)**2*bps(j-1))*betai*f0/bps(j-1)**2
         candyd0(ntheta-j) = -candyd0(j)
         candyd1(ntheta-j) = -candyd1(j)
         candyd2(ntheta-j) = -candyd2(j)
      end do
!compute f0p according Eq.(86)
      candyf0p = (pi*2*q0p/psip(nr2)- ((candyd0(ntheta)-candyd0(0))+(candyd2(ntheta)-candyd2(0))*pprime))/((candyd1(ntheta)-candyd1(0))*f0)*psip(nr2)
      f0p = candyf0p
      do j = 0,ntheta
         candynu1(j) = radius(nr2,j)*bps(j)*(candyd0(j)+candyd1(j)*f0*f0p/psip(nr2)+candyd2(j)*pprime)
         candydr(j) = r0/q0*(f0/(radius(nr2,j)**2*bps(j))*dldr(j)+candynu1(j)*drhdr(j))
      end do

!compute db2dl,db2drho
      do j = 0,ntheta
         db2dl(j) = -2*f0**2/radius(nr2,j)**3*cosu(j)+2*bps(j)*dbpsdl(j)
         db2drho(j) = f0**2/radius(nr2,j)**2*(2*f0p/psip(nr2)/f0*radius(nr2,j)*bps(j)-2*sinu(j)/radius(nr2,j)) &
                      -2*bps(j)**2*(dudl(j)+f0*f0p/psip(nr2)/(radius(nr2,j)*bps(j))+betai*radius(nr2,j)*(pprime)/bps(j))
      end do

!set f(i)
      do i = 0,nr
         r = rin+i*dr
         f(i) = f0+(r-r0)*f0p
         dipdr(i) = f0p
      end do

!!!!! end of ADD f^prime 

! compute psip(r), from f and q
      do i = 0,nr
         dum = 0.
         do j = 0,ntheta-1
            dum = dum+dth/radius(i,j)/grcgt(i,j)
         end do
         psip(i) = f(i)/2/pi/sf(i)*dum
      end do

!compute psi(r)
      dum = 0.
      psi(0) = 0.
      do i = 1,nr
         psi(i) = psi(i-1)+dr*(psip(i-1)+psip(i))/2
      end do

! compute B(r,theta),qhat(r,theta),dbdr,dbdth --modified due to f^prime
      do i = 0,nr
         r = rin+i*dr
         x = asin(tria(i))
         xp = stria(i)/r
         elonp = selon(i)*elon(i)/r
         triap = stria(i)*sqrt(1-tria(i)**2)/r
         do j = 0,ntheta
            th = -pi+j*dth
            bfld(i,j) = sqrt((f(i)/radius(i,j))**2+ &
                     (psip(i)/radius(i,j)*gr(i,j))**2)
            dum2 = -r*sin(th+x*sin(th))*(1+x*cos(th))  !dR d theta
            dum4 = rmajp(i)+cos(th+x*sin(th))-r*sin(th+x*sin(th))* & 
                    xp*sin(th)                         !dR d r
!            dbdr(i,j) = -f(i)/radius(i,j)**2*dum4      !dB_t/dr
!            dbdth(i,j) = -f(i)/radius(i,j)**2*dum2     !dB_t/dth
            qhat(i,j) = f(i)/radius(i,j)/(psip(i)*grcgt(i,j))
         end do
      end do

!compute dbdr,dbdth from dbdl,dbdrho --- added due to f^prime  
!(1) finite-difference 
      do i = 1,nr-1
         do j = 1,ntheta-1
            dbdr(i,j) = (bfld(i+1,j)-bfld(i-1,j))/(2*dr)
            dbdth(i,j) = (bfld(i,j+1)-bfld(i,j-1))/(2*dth)
         end do
         dbdr(i,0) = (bfld(i+1,0)-bfld(i-1,0))/(2*dr)
         dbdr(i,ntheta) = dbdr(i,0)
         dbdth(i,0) = (bfld(i,1)-bfld(i,ntheta-1))/(2*dth)
         dbdth(i,ntheta) = dbdth(i,0)
      end do
      do j = 0,ntheta
         dbdr(0,j) = dbdr(1,j)
         dbdr(nr,j) = dbdr(nr-1,j)
         dbdth(0,j) = dbdth(1,j)
         dbdth(nr,j) = dbdth(nr-1,j)
      end do
!(2) use analytic formulae
      if(icandy==1)then
         do j = 0,ntheta
            dbdl(j) = db2dl(j)/(2*bfld(nr2,j))
            dbdrho(j) = db2drho(j)/(2*bfld(nr2,j))
            dbdr(nr2,j) = dbdl(j)*dldr(j)+dbdrho(j)*drhdr(j)
            dbdth(nr2,j) = dbdl(j)*dldt(j)+dbdrho(j)*drhdt(j)
         end do
         do i = 0,nr
            do j = 0,ntheta
               dbdr(i,j) = dbdl(j)*dldr(j)+dbdrho(j)*drhdr(j)
               dbdth(i,j) = dbdl(j)*dldt(j)+dbdrho(j)*drhdt(j)
            end do
         end do
      end if


!compute dydr(r,theta)
      do i = 1,nr-1
         do j = 0,ntheta
            dqhdr(i,j) = (qhat(i+1,j)-qhat(i-1,j))/(2*dr)
         end do
      end do
      do j = 0,ntheta
         dqhdr(0,j) = (qhat(1,j)-qhat(0,j))/dr
         dqhdr(nr,j) = (qhat(nr,j)-qhat(nr-1,j))/dr
      end do

      do i = 0,nr
         yfn(i,ntheta/2) = 0.
         dydr(i,ntheta/2) = 0.
         dum = 0.
         dum1 = 0.
         do j = ntheta/2+1,ntheta
            dum = dum+(dqhdr(i,j-1)+dqhdr(i,j))*dth/2
            dydr(i,j) = r0/q0*dum
            dydr(i,ntheta-j) = -dydr(i,j)
            dum1 = dum1+r0/q0*(qhat(i,j-1)+qhat(i,j))*dth/2
            yfn(i,j) = dum1
            yfn(i,ntheta-j) = -yfn(i,j)
         end do
      end do

    !compute the flux coordinate theta                                                                                                                                                        \
                                                                                                                                                                                               
    do i = 0,nr
       thflx(i,ntheta/2) = 0.
       dum = 0.
       do j = ntheta/2+1,ntheta
          dum = dum+(qhat(i,j-1)+qhat(i,j))*dth/2
          thflx(i,j) = dum/sf(i)
          thflx(i,ntheta-j) = -thflx(i,j)
       end do
    end do

! use candydr for dydr()
      if(icandy==1)then
         do i = 0,nr
            do j = 0,ntheta
               dydr(i,j) = candydr(j)
            end do
         end do
      end if
            
! compute gxdgy
      jacmax = 0.
      do i = 0,nr
         dum = 0.
         do j = 0,ntheta
            gxdgy(i,j) = dydr(i,j)*gr(i,j)**2+r0/q0*qhat(i,j)*grdgt(i,j)
            jacob(i,j) = 1./(r0*rmaj0/radius(i,j)*grcgt(i,j))
            if(jacob(i,j)>jacmax)jacmax = jacob(i,j)
            if(j<ntheta)dum = dum+jacob(i,j)
         end do
         jacoba(i) = dum/ntheta
      end do

!adjust ion profiles
      eadj = 0.0
      do i = 0,nr-1
         xn0i(i) = xn0i(i)+xn0c(i)*eadj*6
         xn0c(i) = xn0c(i)-xn0c(i)*eadj
      end do

!compute cn0e,cn0i,cn0b,cn0e
      cn0e = 0.
      cn0i = 0.
      cn0b = 0.
      cn0c = 0.
      dum = 0.
      do i = 0,nr-1
         dum = dum+(jacoba(i)+jacoba(i+1))/2.0
      end do
      do i = 0,nr-1
         cn0e = cn0e+(xn0e(i)+xn0e(i+1))/2*(jacoba(i)+jacoba(i+1))/2.0
         cn0i = cn0i+(xn0i(i)+xn0i(i+1))/2*(jacoba(i)+jacoba(i+1))/2.0
         cn0b = cn0b+(xn0b(i)+xn0b(i+1))/2*(jacoba(i)+jacoba(i+1))/2.0
         cn0c = cn0c+(xn0c(i)+xn0c(i+1))/2*(jacoba(i)+jacoba(i+1))/2.0
      end do
      cn0e = cn0e/dum
      cn0i = cn0i/dum
      cn0b = cn0b/dum
      cn0c = cn0c/dum
      n0emax = 0.
      n0imax = 0.
      n0bmax = 0.
      n0cmax = 0.
      do i = 0,nr
         xn0e(i) = xn0e(i)/cn0e
         xn0i(i) = xn0i(i)/cn0i
         xn0b(i) = xn0b(i)/(cn0b+1e-10)
         xn0c(i) = xn0c(i)/(cn0c+1e-10)
         if(xn0e(i)>n0emax)n0emax=xn0e(i)
         if(xn0i(i)>n0imax)n0imax=xn0i(i)
         if(xn0b(i)>n0bmax)n0bmax=xn0b(i)
         if(xn0c(i)>n0cmax)n0cmax=xn0c(i)
      end do


!assign value to xn0s...
      xn0s(1,:) = xn0i(:)
      xn0s(2,:) = xn0c(:)
      xn0s(3,:) = xn0b(:)

      t0s(1,:) = t0i(:)
      t0s(2,:) = t0c(:)
      t0s(3,:) = t0b(:)

      tgis(1) = t0s(1,1)
      tgis(2) = t0s(2,1)
      tgis(3) = t0s(3,1)
      tge = t0e(1)

      capts(1,:) = capti(:)
      capts(2,:) = captc(:)
      capts(3,:) = captb(:)

      capns(1,:) = capni(:)
      capns(2,:) = capnc(:)
      capns(3,:) = capnb(:)

      cn0s(1) = cn0i
      cn0s(2) = cn0c
      cn0s(3) = cn0b

      vpars(1,:) = vpari(:)
      vpars(2,:) = vparc(:)
      vpars(3,:) = vparb(:)

      vparsp(1,:) = vparip(:)
      vparsp(2,:) = vparcp(:)
      vparsp(3,:) = vparbp(:)

!set kapn=0,n=1, etc.

!      t0i = 1.
!      capne = 0.
!      capte = 0.
!      xn0e = 1.
!      t0e = 1.
!      capns(1,:) = 0.
!      capts(1,:) = 0.
!      t0s(1,:) = 1.
!      xn0s(1,:) = 1.
!      cn0e = 1.
!      cn0s(1) = 1.

!compute jfn(theta)
      do j = 0,ntheta
         dum = 0.
         do i = 0,nr-1
            dum = dum+(jacob(i,j)+jacob(i+1,j))/2.0
         end do
         jfn(j) = dum/nr
      end do
      dum = 0.
      do j = 0,ntheta-1
         dum = dum+(jfn(j)+jfn(j+1))/2
      end do
      dum = dum/ntheta
      do j = 0,ntheta
         jfn(j) = dum/jfn(j)
      end do
      jfn = 1.
 300  continue

!bstar effects
! compute psip2(r), from psip(r)
      do i = 1,nr-1
          psip2(i) = (psip(i+1)-psip(i-1))/(2*dr)
      end do
      psip2(0) = psip2(1)
      psip2(nr) = psip2(nr-1)

! compute prsrbr, prsrbz, pthsrbr,pthsrbz
      do j = 0,ntheta
         do i = 1,nr-1
            prsrbr(i,j) = (srbr(i+1,j)-srbr(i-1,j))/(2.*dr)
            prsrbz(i,j) = (srbz(i+1,j)-srbz(i-1,j))/(2.*dr)
         end do
         prsrbr(0,j) = (srbr(1,j)-srbr(0,j))/dr
         prsrbr(nr,j) = (srbr(nr,j)-srbr(nr-1,j))/dr
         prsrbz(0,j) = (srbz(1,j)-srbz(0,j))/dr
         prsrbz(nr,j) = (srbz(nr,j)-srbz(nr-1,j))/dr         
      end do
      do i = 0,nr
         do j = 1,ntheta-1
            pthsrbr(i,j) = (srbr(i,j+1)-srbr(i,j-1))/(2.*dth)
            pthsrbz(i,j) = (srbz(i,j+1)-srbz(i,j-1))/(2.*dth)
         end do
         pthsrbr(i,0) = (srbr(i,1)-srbr(i,0))/dth
         pthsrbr(i,ntheta) = (srbr(i,ntheta)-srbr(i,ntheta-1))/dth
         pthsrbz(i,0) = (srbz(i,1)-srbz(i,0))/dth
         pthsrbz(i,ntheta) = (srbz(i,ntheta)-srbz(i,ntheta-1))/dth
      end do

! compute curvbz
      do i = 0,nr
         do j = 0,ntheta
            dum1 = prsrbz(i,j)*srbz(i,j)+pthsrbz(i,j)*thbz(i,j)
            dum2 = prsrbr(i,j)*srbr(i,j)+pthsrbr(i,j)*thbr(i,j)
            curvbz(i,j) = psip(i)*(dum1/radius(i,j)-srbr(i,j)/radius(i,j)**2 &
                                   +dum2/radius(i,j))
            bdcrvb(i,j) = f(i)/(bfld(i,j)**2*radius(i,j))*(psip2(i)*(gr(i,j))**2/radius(i,j)+curvbz(i,j)) &
                          -1./bfld(i,j)**2/radius(i,j)**2*psip(i)*dipdr(i)*gr(i,j)**2
         end do
      end do

!ExB flow for flux tube
      gamma_E = gamma_E*sqrt(t0e0/mimp)/a
      er0 = mach*r0*bunit/q0*sqrt(t0e0/mimp)/rmaj0
!      erp = -gamma_E*bunit-(q0p/q0-1./r0)*er0 !3/12/2013 note
      erp = -gamma_E*bunit+psip2(nr2)/psip(nr2)*er0 !8/8/2013 note
      if(itube==1)then
         lxa = 1./(q0p*lymult*a) 
         rina = r0a-lxa/2
         routa = r0a+lxa/2
         rin=rina*a
         rout=routa*a
         dr = (rout-rin)/nr
         do i = 0,nr
            r = rin+i*dr
            sf(i) = q0+q0p*(r-r0)
            er(i) = er0+erp*(r-r0)
            phincp(i) = -er(i) !/gr(i,ntheta/2)
         end do
      end if
!compute the parallel flow term to kap
      do i = 0,nr
         upari(i) = -phincp(i)*rmaj0/(psip(i)*bfld(i,ntheta/2))
      enddo
      do i = 1, nr-1
         vparip(i)=(upari(i+1)-upari(i-1))/(2.d0*dr)
      enddo
      vparip(0) = vparip(1)
      vparip(nr) = vparip(nr-1)
      vparsp(1,:) = vparip(:)
      vparsp(2,:) = vparip(:)

      dipdr=0
      bdcrvb=0
      curvbz=0
      psip2=0
      
      rdtemp=0

  end subroutine new_equil

END MODULE gem_equil
