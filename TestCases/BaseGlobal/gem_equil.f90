MODULE gem_equil
  IMPLICIT NONE
  integer :: itube,ibase,iperi,iperidf,ibunit,icandy=0,isprime=0,ildu=0,eldu=0
  real :: mimp=2,mcmp=12,chgi=1,chgc=6
  real :: elon0=1.0,tria0=0.0,rmaj0=1000.0,r0,a=360.0,selon0=0.0,&
             stria0=0.0,rmaj0p=-0.0,q0p=0.006,q0=1.4, elonp0=0.,triap0=0.,erp=0.01,er0=0.,q0abs
  real :: beta,Rovera,shat0,teti,tcti,rhoia,Rovlni,Rovlti,Rovlne,Rovlte,Rovlnc,Rovltc,ncne,nuacs
  real :: gamma_E,mach
  real :: f0, f0p,bunit
  real :: rin,rout,dr,dth,delz,jacmax,eadj
  real :: cn0e,cn0i,cn0b,cn0c,n0emax,n0imax,n0bmax,n0cmax
  real :: r0a,lxa,lymult,delra,delri,delre,delrn,rina,routa,betai, &
               tir0,xnir0,xu,frequ,vu,eru

  integer :: nr=300,nr2=150,ntheta=100,idiag=0,isgnf=1, isgnq=1
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
                                      er,upari,dipdr

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
      implicit none
      real(8) :: pi,pi2,r,th,s,ss,c1,c2,lti,s0,delsi,lte,delse,teti,lnh,delsh,sh
      real(8) :: lref,kappat,kappan,wt,wn
      parameter(kappat=6.96,kappan=2.23,wt=0.3,wn=0.3)
      parameter(c1=0.43236,c2=2.33528,lti=1200.0,s0=0.5,delsi=0.9, &
                lte=1200.0,delse=0.9,lnh=600,delsh=0.2,sh=0.5)
      integer :: i,j,k,m,myid,j1,j2,ileft(0:nr)
      real(8) :: dum,x,denom,dum1,dum2,dum3,dum4,xp,xpp,cutoff,delrai=0.01,delrao=0.01,t,psitmax,dt,beamconst,t0
      real(8) :: vpar,vperpc,vperp,v,plam,vpatcl,auglam,dplam,dlam2,dist,beamconst1,preconst1 
      real(8) :: e,proton,xu,vu,omegau,nu,Tu,Bu
      real(8) :: r1,r2,th1,th2,z1,z2,rdum,rdum1,thdum,thdum1,wx0,wx1,wy0,wy1,psitrin,psita,psitrout
      real(8) :: btor(0:nr,0:ntheta),rhogem(0:nr),rgemplot(0:ntheta-1),zgemplot(0:ntheta-1)
      real(8) :: elonp(0:nr),triap(0:nr)

      allocate(bfld(0:nr,0:ntheta),qhat(0:nr,0:ntheta),radius(0:nr,0:ntheta), &
               gr(0:nr,0:ntheta),gth(0:nr,0:ntheta),grdgt(0:nr,0:ntheta), &
               grcgt(0:nr,0:ntheta),gxdgy(0:nr,0:ntheta),dydr(0:nr,0:ntheta),&
               dbdr(0:nr,0:ntheta),dbdth(0:nr,0:ntheta),dqhdr(0:nr,0:ntheta),&
               jacob(0:nr,0:ntheta), jfn(0:ntheta), zfnth(0:ntheta),thfnz(0:ntheta),&
               yfn(0:nr,0:ntheta),hght(0:nr,0:ntheta),thflx(0:nr,0:ntheta))

      allocate(rmaj(0:nr), elon(0:nr), tria(0:nr), sf(0:nr), psi(0:nr),&
               rmajp(0:nr),selon(0:nr),stria(0:nr),psip(0:nr),&
               f(0:nr),jacoba(0:nr),t0i(0:nr),t0c(0:nr),t0e(0:nr),t0ip(0:nr),&
               t0ep(0:nr),capti(0:nr),captc(0:nr),capte(0:nr),&
               xn0e(0:nr),xn0i(0:nr),xn0c(0:nr),capne(0:nr),capni(0:nr),capnc(0:nr),zeff(0:nr),nue0(0:nr),&
               vpari(0:nr),vparc(0:nr),vparb(0:nr),phinc(0:nr), &
               vparip(0:nr),vparcp(0:nr),vparbp(0:nr),phincp(0:nr),er(0:nr), &
               upari(0:nr),dipdr(0:nr))

               
      allocate(psip2(0:nr), curvbz(0:nr,0:ntheta),srbr(0:nr,0:ntheta),srbz(0:nr,0:ntheta),&
               thbr(0:nr,0:ntheta),thbz(0:nr,0:ntheta),bdcrvb(0:nr,0:ntheta), &
               prsrbr(0:nr,0:ntheta),prsrbz(0:nr,0:ntheta), &
               pthsrbr(0:nr,0:ntheta),pthsrbz(0:nr,0:ntheta))

      allocate(cn0s(1:5),n0smax(1:5),t0s(1:5,0:nr),xn0s(1:5,0:nr),&
               capts(1:5,0:nr),capns(1:5,0:nr),vpars(1:5,0:nr),&
               vparsp(1:5,0:nr),tgis(1:5))

!Normalization
      e = 1.6e-19
      proton = 1.67e-27
      Bu = 2.0
      Tu = 2.14*1000*e  !1KeV
      omegau = e*Bu/proton
      vu = sqrt(Tu/proton)
      xu = proton*vu/(e*Bu)
      nu = 4.66e19*3.0

      rmaj0 = 1.67/xu
      a = rmaj0*0.36
      lref = rmaj0

      beta = 4*3.14159*1e-7*nu*Tu/Bu**2

      rin = rina*a
      rout = routa*a
      r0 = r0a*a
      elonp0 = elon0*selon0/r0
      triap0 = sqrt(1-tria0**2)/r0*stria0
      pi = atan(1.0)*4
      dr = (rout-rin)/nr
      dth = pi*2/ntheta
!      write(*,*)pi,dr,dth,nr,ntheta

      teti = 1.0
! specify f(r),q(r)                                                                                                                          
      do i = 0,nr
         r = rin+i*dr
         f(i) = isgnf*rmaj0
      end do

! specify rmaj(r)                                                                                                                            
      do i = 0,nr
         r = rin+i*dr
         rmaj(i) = rmaj0-rmaj0p/2*(1.-(r/a)**2)*a
         rmajp(i) = rmaj0p*r/a
      end do

! specify elon(r) and compute selon(r)                                                                                                       
      do i = 0,nr
         r = rin+i*dr
         elon(i) = (elon0-0.1)+0.1*(r/a)**4
         selon(i) = r*0.4*(r/a)**3/a/elon(i)
      end do
      if(elon0==1.0)then
         elon = 1.0
         selon = 0.
      end if

! specify tria(r) and compute stria(r)                                                                                                       
      do i = 0,nr
         r = rin+i*dr
         tria(i) = tria0*(r/a)**2
         stria(i) = r*tria0*2.*(r/a)/a/sqrt(1-tria(i)**2)
      end do

! compute radius(r,theta)                                                                                                                    
      do i = 0,nr
         r = rin+i*dr
         do j = 0,ntheta
            th = -pi+dth*j
            radius(i,j) = rmaj(i)+r*cos(th+asin(tria(i))*sin(th))
            hght(i,j) = r*elon(i)*sin(th)
         end do
      end do

! compute grad r,grad theta, grdgt,grcgt
      do i = 0,nr
         r = rin+i*dr
         x = asin(tria(i))
         xp = stria(i)/r
         elonp(i) = selon(i)*elon(i)/r
         triap(i) = stria(i)*sqrt(1-tria(i)**2)/r
         do j = 0,ntheta
            th = -pi+j*dth
            dum1 = elon(i)*r*cos(th)                   !dZ d theta
            dum2 = -r*sin(th+x*sin(th))*(1+x*cos(th))  !dR d theta
            dum3 = elonp(i)*r*sin(th)+elon(i)*sin(th)     !dZ d r 
            dum4 = rmajp(i)+cos(th+x*sin(th))-r*sin(th+x*sin(th))* & 
                    xp*sin(th)                         !dR d r
            denom = dum4*dum1-dum3*dum2

            srbr(i,j) = dum1/denom
            srbz(i,j) = -dum2/denom
            thbr(i,j) = -dum3/denom
            thbz(i,j) = dum4/denom
            if(denom<0)write(*,*)elonp(i),triap(i),xp,denom
            gr(i,j) = sqrt(dum1**2+dum2**2)/denom
            gth(i,j) = sqrt(dum3**2+dum4**2)/denom
            grdgt(i,j) = (-dum1*dum3-dum2*dum4)/denom**2
            grcgt(i,j) = 1/denom
         end do
      end do


!assign q, T, n profiles 
      do i = 0,nr
         r = rin+i*dr
         s = r/a
         sf(i) = 2.52*s**2-0.16*s+0.86
         sf(i) = sf(i)*isgnq
         t0i(i) = exp(-kappat*wt*a/lref*tanh((r-r0)/(wt*a)))
         t0e(i) = t0i(i)
         xn0e(i) = exp(-kappan*wn*a/lref*tanh((r-r0)/(wn*a)))
         xn0i(i) = xn0e(i)
         phincp(i) = 0.
         nue0(i) = 0.
         zeff(i) = 0.         
      end do
      q0 = sf(nr/2)
      q0abs = abs(q0)


      do i = 1,nr-1
         r = rin+i*dr
         s = r/a
         cutoff=1.-exp(-((routa-s)/delrao)**2)-exp(-((s-rina)/delrai)**2)
         ss = s*s
         t0ip(i) = (t0i(i+1)-t0i(i-1))/(2*dr)
         t0ep(i) = (t0e(i+1)-t0e(i-1))/(2*dr)
         capti(i) = -t0ip(i)/t0i(i)*cutoff
         capte(i) = -t0ep(i)/t0e(i)*cutoff
         dum = (xn0e(i+1)-xn0e(i-1))/(2*dr)
         capne(i) = -dum/xn0e(i)*cutoff
         dum = (xn0i(i+1)-xn0i(i-1))/(2*dr)
         capni(i) = -dum/xn0i(i)*cutoff
      end do

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

! compute B(r,theta),qhat(r,theta),dbdr,dbdth
      do i = 0,nr
         r = rin+i*dr
         x = asin(tria(i))
         xp = stria(i)/r
         elonp(i) = selon(i)*elon(i)/r
         triap(i) = stria(i)*sqrt(1-tria(i)**2)/r
         do j = 0,ntheta
            th = -pi+j*dth
            bfld(i,j) = sqrt((f(i)/radius(i,j))**2+ &
                     (psip(i)/radius(i,j)*gr(i,j))**2)
            dum2 = -r*sin(th+x*sin(th))*(1+x*cos(th))  !dR d theta
            dum4 = rmajp(i)+cos(th+x*sin(th))-r*sin(th+x*sin(th))* & 
                    xp*sin(th)                         !dR d r
            dbdr(i,j) = -abs(f(i))/radius(i,j)**2*dum4      !dB_t/dr
            dbdth(i,j) = -abs(f(i))/radius(i,j)**2*dum2     !dB_t/dth
            qhat(i,j) = f(i)/radius(i,j)/(psip(i)*grcgt(i,j))
!            qhat(i,j) = f(i)/rmaj0/(psip(i)*grcgt(i,j))
         end do
      end do

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

!compute the flux coordinate theta
      do i = 0,nr
         thflx(i,ntheta/2) = 0.
         dum = 0.
         do j = ntheta/2+1,ntheta
            dum = dum+(qhat(i,j-1)+qhat(i,j))*dth/2
            thflx(i,j) = dum/sf(i)
            thflx(i,ntheta-j) = -thflx(i,j)
         end do
      end do

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

!compute cn0e,cn0i,cn0b,cn0e
      cn0e = 0.
      cn0i = 0.
      cn0c = 0.
      dum = 0.
      do i = 0,nr-1
         dum = dum+(jacoba(i)+jacoba(i+1))/2.0
      end do
      do i = 0,nr-1
         cn0e = cn0e+(xn0e(i)+xn0e(i+1))/2*(jacoba(i)+jacoba(i+1))/2.0
         cn0i = cn0i+(xn0i(i)+xn0i(i+1))/2*(jacoba(i)+jacoba(i+1))/2.0
         cn0c = cn0c+(xn0c(i)+xn0c(i+1))/2*(jacoba(i)+jacoba(i+1))/2.0
      end do
      cn0e = cn0e/dum
      cn0i = cn0i/dum
      cn0c = cn0c/dum
      do i = 0,nr
         xn0e(i) = xn0e(i)/cn0e
         xn0i(i) = xn0i(i)/cn0i
         xn0c(i) = xn0c(i)/cn0c
      end do

!assign value to xn0s...                                                                                                                                                                                                
      xn0s(1,:) = xn0i(:)
      xn0s(2,:) = xn0c(:)
!      xn0s(3,:) = xn0b(:)

      t0s(1,:) = t0i(:)
      t0s(2,:) = t0c(:)
!      t0s(3,:) = t0b(:)

      tgis(1) = t0s(1,1)
      tgis(2) = t0s(2,1)
      tgis(3) = t0s(3,1)
      tge = t0e(1)

      capts(1,:) = capti(:)
      capts(2,:) = captc(:)
!      capts(3,:) = captb(:)

      capns(1,:) = capni(:)
      capns(2,:) = capnc(:)
!      capns(3,:) = capnb(:)

      cn0s(1) = cn0i
      cn0s(2) = cn0c
      cn0s(3) = cn0b

      vpars(1,:) = 0. !vpari(:)
      vpars(2,:) = 0. !vparc(:)
!      vpars(3,:) = vparb(:)

      vparsp(1,:) = 0. !vparip(:)
      vparsp(2,:) = 0. !vparcp(:)
!      vparsp(3,:) = vparbp(:)


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
!      jfn = 1.

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
            bdcrvb(i,j) = f(i)/(bfld(i,j)**2*radius(i,j))*(psip2(i)*(gr(i,j))**2/radius(i,j)+curvbz(i,j))                       
         end do
      end do

      dipdr = 0.
      bdcrvb = 0.
      curvbz = 0.
      psip2 = 0.
  end subroutine new_equil

END MODULE gem_equil
