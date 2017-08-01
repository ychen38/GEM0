MODULE equil
  IMPLICIT NONE
  integer :: itube,iperi,iperidf,ibunit,icandy=0,isprime=0,ildu=1,eldu=0
  real(8) :: mimp=2,mcmp=12,chgi=1,chgc=6
  real(8) :: elon0=1.0,tria0=0.0,rmaj0=1000.0,r0,a=360.0,selon0=0.0,&
             stria0=0.0,rmaj0p=-0.0,q0p=0.006,q0=1.4, elonp0=0.,triap0=0.,erp=0.01,er0=0.,q0abs
  real(8) :: Rovera,shat0,teti,tcti,rhoia,Rovlni,Rovlti,Rovlne,Rovlte,Rovlnc,Rovltc,ncne,nuacs
  real(8) :: gamma_E,mach
  real(8) :: f0, f0p,bunit
  real(8) :: rin,rout,dr,dth,delz,jacmax,eadj
  real(8) :: cn0e,cn0i,cn0b,cn0c,n0emax,n0imax,n0bmax,n0cmax
  real(8) :: r0a,lxa,lymult,delra,delri,delre,delrn,rina,routa,betai, &
               tir0,xnir0,xu,frequ,vu,eru
  real(8) :: qoffset=0.  

  integer :: nr=300,nr2=150,ntheta=100,isgnq=-1,idiag=0
  real(8),dimension(:,:),allocatable :: bfld,qhat,radius,gr,gth,grdgt,grcgt, &
                                        gxdgy,dydr,dbdr,dbdth,dqhdr,jacob, &
                                        yfn,hght,thflx
  real(8),dimension(:),allocatable :: rmaj,rmajp,elon,selon,tria,stria, psi,&
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
                                      db2dl,db2drho,dbpsdl,dipdr
!for Miller local flux-tube
  real(8) :: candyf0p
  real(8),dimension(:),allocatable :: candyd0,candyd1,candyd2,candynus,candynu1,candydr

!for including bstar effects
  real(8),dimension(:),allocatable :: psip2
  real(8),dimension(:,:),allocatable :: curvbz,srbr,srbz,thbr,thbz,prsrbr,prsrbz,pthsrbr,pthsrbz,bdcrvb

  real(8),dimension(:,:),allocatable :: t0s,xn0s,capts,capns,vpars,vparsp
  real(8),dimension(:),allocatable :: cn0s,n0smax,tgis
  real(8) :: tge
  character(32) :: trflnm ! profile-data-file name
  real(8),external :: erf

contains
  subroutine new_equil()
      real(8) :: pi,pi2,r,th,s,ss,c1,c2,lti,s0,s1,delsi,lte,delse,ln,deln,lnh,delsh,sh
      parameter(c1=0.43236,c2=2.33528,lti=144.9,s1=0.5,delsi=0.6, &
                lte=44.9,delse=0.6,deln=0.2,ln=0.3,lnh=0.1,delsh=0.2,sh=0.5)
      integer :: i,j,k,j1,j2,ileft(0:nr)
      real(8) :: dum,x,denom,dum1,dum2,dum3,dum4,xp,elonp,triap,t,psitmax,dt,t0
      real(8) :: t0i0,t0e0,t0c0,ni0,nc0,ne0,t0i0p,t0e0p,t0c0p,ni0p,nc0p,ne0p
      real(8) :: dr24,cutoff,cutoffti,cutoffte,cutoffi,cutoffe,cutoffn
      real(8) :: e,proton,xu,vu,omegau,nu,Tu,Bu
      character(100) :: header,header2

      real(8) :: eps0=8.85e-12,me=0.9e-30,vte,pprime

      integer :: nraefl,nrgtc,nthgtc
      parameter(nraefl=100,nrgtc=100,nthgtc=200)
      real(8) :: Rgtc(0:nrgtc-1,0:nthgtc-1),Zgtc(0:nrgtc-1,0:nthgtc-1),Bgtc(0:nrgtc-1,0:nthgtc-1),rhogtc(0:nrgtc-1),thgtc(0:nthgtc-1),rminorgtc(0:nrgtc-1)
      real(8) :: thmid1(0:nrgtc-1),thmid2(0:nrgtc-1)
      real(8) :: ti_raefl(0:nraefl-1),te_raefl(0:nraefl-1),ni_raefl(0:nraefl-1),ne_raefl(0:nraefl-1),sf_raefl(0:nraefl-1),nh_raefl(0:nraefl-1),sf_ref(0:nraefl-1),tc_raefl(0:nraefl-1)      
      real(8) :: r1,r2,th1,th2,thgtcp,z1,z2,rdum,rdum1,thdum,thdum1,dthgtc,wx0,wx1,wy0,wy1,psitrin,psita,psitrout
      real(8) :: btor(0:nr,0:ntheta),rgemplot(0:ntheta-1),zgemplot(0:ntheta-1),rgtcplot(0:nthgtc-1),zgtcplot(0:nthgtc-1)
      real(8) :: sfgtc(0:nr),negtc(0:nr),nigtc(0:nr),ncgtc(0:nr),tegtc(0:nr),tigtc(0:nr),tcgtc(0:nr),psigtc(0:nr)   
      real(8) :: psit(0:nr),psitp(0:nr),rhogem(0:nr)

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
               upari(0:nr),dipdr(0:nr))

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

!Normalization
      e = 1.6e-19
      proton = 1.67e-27
      Bu = 1.98562
      Tu = 1000*e  !1KeV
      omegau = e*Bu/proton
      vu = sqrt(Tu/proton)
      xu = proton*vu/(e*Bu)
      nu = 2.5e19
      if(itube==0)then
         rmaj0 = 10./xu
         a = 1.0/xu
      end if
      betai = 4*3.14159*1e-7*nu*Tu/Bu**2
      open(51,file='gtcAE_1dProfiles.txt',status='old')
      do i = 0,nraefl-1
         read(51,*)sf_raefl(i),ne_raefl(i),nh_raefl(i),te_raefl(i),tc_raefl(i),ti_raefl(i)
      end do
      close(51)
!fix q(0 and q(a), adjust q(r)
      sf_ref = sf_raefl
      j = 60
      do i = 0,60
         sf_ref(i) = sf_raefl(0)+(sf_raefl(j)-sf_raefl(0))/float(j)*float(i)
      end do
      do i = 0,nraefl-1
         sf_raefl(i) = sf_ref(i)+(sf_raefl(i)-sf_ref(i))*1.-qoffset
      end do

      open(51,file='gtcAE_rho.txt',status='old')
      do i = 0,nrgtc-1
            read(51,*)rhogtc(i)
      end do
      close(51)

      open(51,file='gtcAE_rdata.txt',status='old')
      do i = 0,nrgtc-1
            read(51,*)(Rgtc(i,j),j = 0,nthgtc-1)
      end do
      close(51)
      Rgtc = Rgtc*1.72237/xu

      open(51,file='gtcAE_zdata.txt',status='old')
      do i = 0,nrgtc-1
         read(51,*)(zgtc(i,j),j = 0,nthgtc-1)
      end do
      close(51)
      Zgtc = Zgtc*1.72237/xu

      open(51,file='gtcAE_bdata_tor.txt',status='old')
!      open(51,file='gtcAE_bdata.txt',status='old')
      do i = 0,nrgtc-1
         read(51,*)(bgtc(i,j),j = 0,nthgtc-1)
      end do
      close(51)
      Bgtc = Bgtc*1.98562/Bu

!calculate the r grids corresponding to GTC rho grids, defined to be r=(R+-R-)/2 at Z=0
      pi = atan(1.0)*4
      pi2 = pi*2
      dthgtc = pi2/nthgtc
      rminorgtc = 0.
      do i = 1, nrgtc-1
         rdum = 0.
         rdum1 = 0.
!find the two major radii at Z=Zaxis
         k = 0
         do j = 0,nthgtc-1
            j1 = j
            j2 = j+1
            if(j2==nthgtc)j2 = 0
            th1 = j1*dthgtc
            th2 = th1+dthgtc
            z1 = zgtc(i,j1)
            z2 = zgtc(i,j2)
            dum = 0.
            if(z2==z1)dum = 1.e-10
            r1 = rgtc(i,j1)
            r2 = rgtc(i,j2)
            if(zgtc(i,j1)*zgtc(i,j2) .le. 0.0)then
               if(k==0)then
                  thdum = th1-dthgtc/(z2-z1+dum)*z1   !(thdum-th1) / (0-z1) = dthgtc/(z2-z1)
                  rdum = r1+(thdum-th1)/dthgtc*(r2-r1) !
                  k = 1
                  goto 200
               end if
               if(k==1)then
                  thdum1 = th1-dthgtc/(z2-z1+dum)*z1
                  rdum1 = r1+(thdum1-th1)/dthgtc*(r2-r1)
                  goto 201
               end if
            end if
 200        continue
         end do
 201     rminorgtc(i) = abs(rdum-rdum1)/2
         if(rdum<rdum1)then
            thmid1(i) = thdum1
            thmid2(i) = thdum
         end if
         if(rdum>rdum1)then
            thmid1(i) = thdum
            thmid2(i) = thdum1
         end if
         if(thmid1(i)>pi)thmid1(i) = thmid1(i)-pi2   !thdum1 around theta=0. Make all values on the same branch
      end do
      a = rminorgtc(nrgtc-1)
      rmaj0 = rgtc(0,0)

! if analytical equilibrium
      if(itube==1)then
         a = Rmaj0/Rovera
      endif

      r0=r0a*a
      rin=rina*a
      rout=routa*a
      pi = atan(1.0)*4.
      dr = (rout-rin)/nr
      dth = pi*2/ntheta
      if(idiag==1)write(*,*)pi,dr,dth,nr,ntheta

! My theta will be defined for Z>0 and Z<0 separately, so that up-down symmetry can be easily constructed
      k = 0
      do i = 0,nr
         r = rin+i*dr
         do j = k,nrgtc-2
            if((rminorgtc(j) .le. r) .and. (r .le. rminorgtc(j+1)))then
               ileft(i) = j
!               write(*,*)'in ileft loop', i,r,j
               goto 100
            end if
         end do
 100     k = ileft(i)
      end do


      do i = 0,nr
         r = rin+i*dr
         k = ileft(i)
         dum = rminorgtc(k+1)-rminorgtc(k)
         wx0 = (rminorgtc(k+1)-r)/dum
         wx1 = 1.-wx0
         psigtc(i) = wx0*rhogtc(k)+wx1*rhogtc(k+1)  !normalized toroidal flux, from GTC grids
         sfgtc(i) = wx0*sf_raefl(k)+wx1*sf_raefl(k+1)  
         negtc(i) = wx0*ne_raefl(k)+wx1*ne_raefl(k+1)  
         ncgtc(i) = wx0*nh_raefl(k)+wx1*nh_raefl(k+1)  
         tegtc(i) = wx0*te_raefl(k)+wx1*te_raefl(k+1)  
         tigtc(i) = wx0*ti_raefl(k)+wx1*ti_raefl(k+1)  
         tcgtc(i) = wx0*tc_raefl(k)+wx1*tc_raefl(k+1)  
         th1 = wx0*thmid1(k)+wx1*thmid1(k+1)
         th2 = wx0*thmid2(k)+wx1*thmid2(k+1)
         do j = 0,ntheta
            th = -pi+j*dth
            thgtcp = th/pi*(th2-th1)+th1
            if(th<0)thgtcp = (th+pi)/pi*(pi2+th1-th2)+th2
            thgtcp = mod(thgtcp+pi2,pi2)
            j1 = int(thgtcp/dthgtc)
            j1 = min(j1,nthgtc-1)
            j2 = j1+1
            if(j2==nthgtc)j2=0
            wy0 = ((j1+1)*dthgtc-thgtcp)/dthgtc
            wy1 = 1.-wy0
            radius(i,j) = wx0*wy0*rgtc(k,j1)+wx0*wy1*rgtc(k,j2) &
                         +wx1*wy0*rgtc(k+1,j1)+wx1*wy1*rgtc(k+1,j2) 
            hght(i,j) = wx0*wy0*zgtc(k,j1)+wx0*wy1*zgtc(k,j2) &
                         +wx1*wy0*zgtc(k+1,j1)+wx1*wy1*zgtc(k+1,j2) 
            btor(i,j) = wx0*wy0*bgtc(k,j1)+wx0*wy1*bgtc(k,j2) &
                         +wx1*wy0*bgtc(k+1,j1)+wx1*wy1*bgtc(k+1,j2) 
         end do
      end do
      negtc = negtc*0.32785e20/nu
      ncgtc = ncgtc*0.32785e20/nu
      tegtc = tegtc*1689.46*e/tu
      tigtc = tigtc*1689.46*e/tu
      tcgtc = tcgtc*1689.46*e/tu
      nigtc = negtc-ncgtc

!Make up-down symmetry
      do i = 0,nr
         do j = 1,ntheta/2-1
            radius(i,j) = radius(i,ntheta-j)
            hght(i,j) = -hght(i,ntheta-j)
            btor(i,j) = btor(i,ntheta-j)
         end do
      end do

! compute f
      do i = 0,nr
         dum = 0.
         do j = 0,ntheta-1
            dum = dum+btor(i,j)*radius(i,j)
         end do
         f(i) = dum/ntheta
      end do


! compute grad r,grad theta, grdgt,grcgt
      do i = 0,nr
         r = rin+i*dr
         x = asin(tria(i))
         xp = stria(i)/r
         elonp = selon(i)*elon(i)/r
         triap = stria(i)*sqrt(1-tria(i)**2)/r
         do j = 0,ntheta
            th = -pi+j*dth
            j1 = mod(j-1+ntheta,ntheta)
            j2 = mod(j+1+ntheta,ntheta)
            dum1 = (hght(i,j2)-hght(i,j1))/(2*dth)
            dum2 = (radius(i,j2)-radius(i,j1))/(2*dth)
            if(i==0)then
               dum3 = (hght(1,j)-hght(0,j))/dr
               dum4 = (radius(1,j)-radius(0,j))/dr
            elseif(i==nr)then
               dum3 = (hght(nr,j)-hght(nr-1,j))/dr
               dum4 = (radius(nr,j)-radius(nr-1,j))/dr
            else
               dum3 = (hght(i+1,j)-hght(i-1,j))/(2*dr)
               dum4 = (radius(i+1,j)-radius(i-1,j))/(2*dr)
            endif
!            dum1 = elon(i)*r*cos(th)                   !dZ d theta
!            dum2 = -r*sin(th+x*sin(th))*(1+x*cos(th))  !dR d theta
!            dum3 = elonp*r*sin(th)+elon(i)*sin(th)     !dZ d r 
!            dum4 = rmajp(i)+cos(th+x*sin(th))-r*sin(th+x*sin(th))* & 
!                    xp*sin(th)                         !dR d r
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


      do i = 0,nr
         r = rin+i*dr
         s = r/a
         sf(i) = 1.71+0.16*s**2
         sf(i) = sf(i)

         t0i(i) = 1.0
         xn0i(i) = 1.0

         xn0c(i) = 0.0375*exp(-deln/ln*tanh((s-s1)/deln))
         t0c(i) = 100.0

         xn0e(i) = xn0i(i)+xn0c(i)         
         t0e(i) = 1.0

         phincp(i) = 0.005*sin((r-rin)/(rout-rin)*2*pi)
         nue0(i) = 1.
         zeff(i) = 1.
      end do
      q0 = sf(nr2)
      q0abs = abs(q0)

! compute psitp(r), from f
      do i = 0,nr
         dum = 0.
         do j = 0,ntheta-1
            dum = dum+dth/radius(i,j)/grcgt(i,j)
         end do
         psitp(i) = f(i)/2/pi*dum
      end do

!compute psit(r)
      dum = 0.
      psit(0) = 0.
      do i = 1,nr
         psit(i) = psit(i-1)+dr*(psitp(i-1)+psitp(i))/2
      end do
!find psit(rin) and psit(a) so that the normalized psit matches gtc rho grids at rin and rout
      dum = psigtc(0)**2 !psit(rin)/psit(a)
      psita = psit(nr)/(psigtc(nr)**2-dum) !psit(a)
      psitrin = psigtc(0)**2*psita !psit(rin)
      psit = psit+psitrin

!assign q, T, n profiles using raefl.dat
      psitmax=psita !psit(nr)*(1./routa)**2 !extrapolate to t=1 using the scaling psit~r**2
      dt=1./nraefl
      t0 = sqrt(psit(0)/psitmax)
      do i = 0,nr
         t = sqrt(psit(i)/psitmax)
         rhogem(i) = t
         k = int(t/dt)
         k = min(k,nraefl-1)
         k = max(k,0)
         dum = ((k+1)*dt-t)/dt
         dum1 = 1.-dum
         t0i(i) = 1.0 !tigtc(i) !ti_raefl(k)*dum+ti_raefl(k+1)*dum1
         t0c(i) = tcgtc(i) !tc_raefl(k)*dum+tc_raefl(k+1)*dum1
         t0e(i) = 1.0 !tegtc(i) !te_raefl(k)*dum+te_raefl(k+1)*dum1
         xn0i(i) = nigtc(i) !(ni_raefl(k)*dum+ni_raefl(k+1)*dum1)/2.5 !n_u=2.5x10^13/cm^3

         xn0c(i) = ncgtc(i) !(nh_raefl(k)*dum+nh_raefl(k+1)*dum1)/2.5 !n_u=2.5x10^13/cm^3
!         xn0c(i) = exp(-delsh/lnh*tanh((t-sh)/delsh))
         xn0e(i) = 1.0 !negtc(i) !(ne_raefl(k)*dum+ne_raefl(k+1)*dum1)/2.5

         sf(i) = sfgtc(i)*isgnq !(sf_raefl(k)*dum+sf_raefl(k+1)*dum1)*isgnq

         phincp(i) = 0.005*sin((r-rin)/(rout-rin)*2*pi)
         nue0(i) = 1.
         zeff(i) = 1.
      end do
      q0 = sf(nr/2)
      q0abs = abs(q0)

      xn0c = xn0c*0.55
      xn0i = xn0e-xn0c

      capne = 0.
      capni = 0.
      capnc = 0.
      capte = 0.
      capti = 0.
      captc = 0.
      do i = 1,nr-1
         capne(i) = -(xn0e(i+1)-xn0e(i-1))/(2*dr)/xn0e(i)
         capte(i) = -(t0e(i+1)-t0e(i-1))/(2*dr)/t0e(i)

         capni(i) = -(xn0i(i+1)-xn0i(i-1))/(2*dr)/xn0i(i)
         capti(i) = -(t0i(i+1)-t0i(i-1))/(2*dr)/t0i(i)
         capnc(i) = -(xn0c(i+1)-xn0c(i-1))/(2*dr)/xn0c(i)
         captc(i) = -(t0c(i+1)-t0c(i-1))/(2*dr)/t0c(i)
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
!      jfn = 1.
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
!         lxa = 1./(q0p*lymult*a) 
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

  end subroutine new_equil

END MODULE equil
