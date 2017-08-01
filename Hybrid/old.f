!       3D Flux Tube Toroidal Electromagnetic GK Delta-f Code
!   global variables...

         use gem_com
         use equil
         use fft_wrapper
	implicit none
	integer :: n,i,j,k,ip

	call initialize
! use the following two lines for r-theta contour plot
        if(iCrs_Sec==1)then
           call pol2d
           goto 100
        end if
	if(iget.eq.0)call loader_wrapper
        if(iget.eq.1)then
           call restart(1,0)
        end if
        starttm=MPI_WTIME()

        do  timestep=ncurr,nm
           tcurr = tcurr+dt

	   call accumulate(timestep-1,0)
	   call poisson(timestep-1,0)
	   call ampere(timestep-1,0)
	   call field(timestep-1,0)
	   call split_weight(timestep-1,0)
	   call diagnose(timestep-1)
           call reporter(timestep-1)

	   call push_wrapper(timestep,1)

	   call accumulate(timestep,1)
	   call poisson(timestep,1)
	   call ampere(timestep,1)
	   call field(timestep,1)
	   call split_weight(timestep,1)

	   call push_wrapper(timestep,0)
           if(mod(timestep,10000)==0)then
              do i=0,last 
                 if(myid==i)write(*,*)myid,mmb
                 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
              end do
           end if

         end do
	 lasttm=MPI_WTIME()
	 tottm=lasttm-starttm
!	 write(*,*)'ps time=',pstm,'tot time=',tottm
         do i=0,last 
            if(myid==i)write(*,*)myid,mme
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         end do
 100     call MPI_FINALIZE(ierr)
         end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hybinit

      use gem_com
      use equil
      implicit none
      GCLR = int(MyId/ntube)
      GLST = numprocs/ntube-1
      TCLR = mod(MyId,ntube)
      TLST = ntube-1

!***MPI variables

!      if (GCLR.eq.GLST) 
!     %     mykm=km-GLST*mykm
      mykm = 1
      rngbr = GCLR+1
      if (GCLR.eq.GLST)rngbr = 0
      
      lngbr = GCLR-1
      if (GCLR.eq.0) lngbr = GLST

      idnxt = TCLR+1
      if (TCLR.eq.TLST)idnxt = 0
      
      idprv = TCLR-1
      if (TCLR.eq.0) idprv = TLST
!     
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine init
      
      use gem_com
      use equil
      implicit none
      character*(62) dumchar
      INTEGER :: i,j,k,n,ns,idum,i1,k1
      INTEGER :: mm1,lr1
      REAL(8) :: mims1,tets1,q1,kappan,kappat,r,qr,th,cost,dum,zdum
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp
      REAL(8) :: grp,gxdgyp,jacp,jfnp,gn0ep,gn0ip
      REAL(8) :: wx0,wx1,wz0,wz1,b

      IU=cmplx(0.,1.)
      pi=4.0*atan(1.0)
      pi2 = pi*2.

      open(115,file='gem.in')
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) imx,jmx,kmx,mmx,nmx,nsmx,modemx,ntube,icrs_sec,ipg,isphi
      read(115,*) dumchar
      read(115,*) im,jm,km,lxa,lymult,delra
      read(115,*) dumchar
      read(115,*) dt,nm,nsm,xshape,yshape,zshape,ipf,isham,ishgk
      read(115,*) dumchar
      read(115,*) kymode,iput,iget,ision,ish,peritr,llk,mlk,onemd,izon,kxz
      read(115,*) dumchar
      read(115,*) nplot,xnplt,modem,contu,pskip,wmax,imovie
      read(115,*) dumchar
      read(115,*) cut,amp,tor,ishift,fradi,kxcut,kycut,bcut
      read(115,*) dumchar
      read(115,*) dum,r0a,dum,dum,width,vpp,vt0,yd0
      read(115,*) dumchar
      read(115,*) c1,c2,c3,c4,ifluid,isg,amie,rneu
      read(115,*) dumchar
      read(115,*) beta,nonlin,nonline,ipara,vwidth,vwidthe,vcut,isuni,idg
!  read TRANSP profile-data-file name from gem.in
      read(115,*) dumchar
      read(115,"(a32)") trflnm

      call new_gem_com()
!      rin = r0-lx/2
!      rout = r0+lx/2
      rina = r0a-lxa/2.
      routa = r0a+lxa/2.

      call new_equil()
      lx=lxa*a
      ly=2.*pi*r0*lymult/q0
      br0 = rmaj0
      lr0 = r0
      qp = q0p
      lz = pi2*q0*rmaj0
      delz = lz/ntheta
!      write(*,*)'br0,lr0,q0,qp = ', br0,lr0,q0,qp
      if(iget.eq.0)then
         if(myid.eq.master)open(9, file='plot', status='unknown')
         if(myid.eq.master)open(11, file='flux', status='unknown')
         if(myid.eq.master)open(15, file='zonal_amp', status='unknown')
         if(myid.eq.master)open(16, file='indicator', status='unknown')
         if(myid.eq.master)open(17, file='yyre', status='unknown')
      end if
      if(iget.eq.1)then
         if(myid.eq.master)open(9, file='plot', & 
             status='unknown',position='append')
         if(myid.eq.master)open(11, file='flux', &
             status='unknown',position='append')
         if(myid.eq.master)open(15, file='zonal_amp', &
             status='unknown',position='append')
         if(myid.eq.master)open(16, file='indicator', &
             status='unknown',position='append')
         if(myid.eq.master)open(17, file='yyre', &
             status='unknown',position='append')
      end if

      if(isuni.eq.0)vwidthe=vwidth
      dte = dt
      iadi = 0
      if(isg.gt.0.)fradi = isg
      if(ifluid.eq.0)then
         iadi = 1
         fradi = 1.
      end if


!     begin reading species info, ns=1,nsm...
      if(nsm.le.0) write(*,*)'invalid nsm',nsm
      read(115,*) dumchar
      ptr(1)=1
      ns = 1
      read(115,*) dumchar
      read(115,*) mm1,mims(1),q(1),mims(2),q(2),lr1
      read(115,*) dumchar
      read(115,*) kappan,kapti,kapte
      read(115,*) dumchar
      read(115,*) nh,lh,ecpow,ehmax,ehmin,teth,iflr,iorb,iflrh
      tmm(1)=mm1
      mm(:)=int(mm1/numprocs)
      mme = int(mm1/numprocs)
      if (MyId.eq.Last) mm(ns)=mm1-Last*mm(ns)
!     write(*,*)'in init  ',Myid,mm(ns)
      tets(1)=1
      lr(1)=lr1
      lr(2)=lr1

      kapn(ns)=kappan
      kapt(ns)=kappat
      if (ns.lt.nsm) ptr(ns+1)=ptr(ns)+mm(ns)

      emass = 1./(amie-1.)
      qel = -1.

      mbeam = 2
      qbeam = 1

      if(iget.eq.1) amp=0.
!     totvol is the square for now...
      dx=lx/float(im)
      dy=ly/float(jm)
      dz=lz/float(km)
!      totvol=lx*ly*lz

      e0=lr0/q0/br0
!     
      do 10 i=0,nxpp
         xg(i)=i*dx !dx*(tclr*nxpp+i)
 10   continue
      do 12 j=0,jm
         yg(j)=dy*float(j)
 12   continue
      kcnt=1
      jcnt = jmx/ntube
      do 14 k=0,mykm
         n=GCLR*kcnt+k	
         zg(k)=dz*float(n)
 14   continue

!     set arrays for mode history plots...
      lmode(1)=0
      lmode(2)=0
      lmode(3)=0
      lmode(4)=0
      mmode(1)=1 !int(.33*ly/2/pi)-1
      mmode(2)=2 !int(.33*ly/2/pi)
      mmode(3)=3 !int(.33*ly/2/pi)+1
      mmode(4)=4 !int(.33*ly/2/pi)+2
      nmode(1)=km/2
      nmode(2)=km/2
      nmode(3)=km/2
      nmode(4)=km/2

!     initialize bfld
      zfnth(0) = 0.
      do j = 1,ntheta
         zfnth(j) = zfnth(j-1)+dth*q0*br0*(1./jfn(j-1)+1./jfn(j))/2
      end do

      thfnz(0) = -pi
      thfnz(ntheta) = pi
      k = 0
      do j = 1,ntheta-1
         zdum = j*lz/ntheta
         do i = k,ntheta-1
            if(zfnth(i)<=zdum.and.zfnth(i+1)>zdum)then
               k = i
               dum = (zdum-zfnth(i))*dth/(zfnth(i+1)-zfnth(i))
               thfnz(j) = i*dth-pi+dum
               go to 127
            end if
         end do
 127     continue
      end do

      do i1 = 0,nxpp
         r = xg(i1)-0.5*lx+lr0
         do k1 = 0,mykm
            k = int(zg(k1)/delz)
            k = min(k,ntheta-1)
            wz0 = ((k+1)*delz-zg(k1))/delz
            wz1 = 1-wz0
            th = wz0*thfnz(k)+wz1*thfnz(k+1)
            i = int((r-rin)/dr)
            i = min(i,nr-1)
            wx0 = (rin+(i+1)*dr-r)/dr
            wx1 = 1.-wx0
            k = int((th+pi)/dth)
            k = min(k,ntheta-1)
            wz0 = (-pi+(k+1)*dth-th)/dth
            wz1 = 1.-wz0
         dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
                 +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
         dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
                 +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
         grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
                 +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
         radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
                 +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
         dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
                 +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
         qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
                 +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
         grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
                 +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
         jacp = wx0*wz0*jacob(i,k)+wx0*wz1*jacob(i,k+1) &
                 +wx1*wz0*jacob(i+1,k)+wx1*wz1*jacob(i+1,k+1) 
         gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
                 +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 
         fp = wx0*f(i)+wx1*f(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         gn0ep = wx0*xn0e(i)+wx1*xn0e(i+1)        
         gn0ip = wx0*xn0i(i)+wx1*xn0i(i+1)        
         b=1.-tor+tor*bfldp
         cfx(i1,k1) = br0/b**3*fp/radiusp*dbdtp*grcgtp
         cfy(i1,k1) = br0/b**3*fp/radiusp* &
                      (dydrp*dbdtp-lr0/q0*qhatp*dbdrp)*grcgtp
         bdgxcgy(i1,k1) = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp             
         bmag(i1,k1) = b
         jac(i1,k1) = jacp*jfnp
         gn0i(i1) = gn0ip
         gn0e(i1) = gn0ep
         end do
      end do

      iseed = 1777+myid*13
      idum = ran2(iseed)
      phi = 0.
      apar = 0.
      dene = 0.
      upar = 0.
      upa0 = 0.

      if(myid.eq.master)then
         write(*,*)zfnth(ntheta),thfnz(ntheta/2),thfnz(ntheta/2+1)
         write(9,*)'dt,beta= ',dt, beta
         write(9,*)'amp,vpp,yd0 = ',amp,vpp,yd0
         write(9,*)'peritr,ifluid= ',peritr,ifluid
         write(9,*)'tor,nonlin= ',tor,nonlin
         write(9,*)'isuni= ',isuni, 'amie= ',amie
         write(9,*)'kxcut,kycut,bcut= ',kxcut,kycut,bcut
         write(9,*)'fradi,isg= ',fradi,isg
         write(9,*)'llk,mlk,onemd =',llk,mlk,onemd
         write(9,*)'vwidth (i,e), vcut= ', vwidth,vwidthe,vcut
         write(9,*)'rneu= ', rneu
         write(9,*)'mm1= ',mm1
         write(9,*)'ish,ishgk,isham = ',ish,ishgk,isham
         write(9,*)'nh,lh,th = ',nh,lh,teth
         write(9,*) 'lxa,lymult,delra,r0a,rina,routa=',lxa,lymult,delra,r0a,rina,routa
         write(9,*) 'a,r0,rmaj0,q0,lx,ly,lz=',a,r0,rmaj0,q0,lx,ly,lz
      end if

      if(myid.eq.master)then
         write(*,*)'dt,beta= ',dt, beta
         write(*,*)'amp,vpp,yd0 = ',amp,vpp,yd0
         write(*,*)'peritr,ifluid= ',peritr,ifluid
         write(*,*)'tor,nonlin= ',tor,nonlin
         write(*,*)'isuni= ',isuni, 'amie= ',amie
         write(*,*)'kxcut,kycut,bcut= ',kxcut,kycut,bcut
         write(*,*)'fradi,isg= ',fradi,isg
         write(*,*)'llk,mlk,onemd =',llk,mlk,onemd
         write(*,*)'vwidth (i,e), vcut= ', vwidth,vwidthe,vcut
         write(*,*)'rneu= ', rneu
         write(*,*)'mm1= ',mm1
         write(*,*)'ish,ishgk,isham = ',ish,ishgk,isham
         write(*,*)'nh,lh,th = ',nh,lh,teth
         write(*,*) 'lxa,lymult,delra,r0a,rina,routa=',lxa,lymult,delra,r0a,rina,routa
         write(*,*) 'a,r0,rmaj0,q0,lx,ly,lz=',a,r0,rmaj0,q0,lx,ly,lz
      end if
      close(115)
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ppush(n,ns)

      use gem_com
      use equil
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum,dum1
      INTEGER :: m,i,j,k,l,n,ns
      INTEGER :: np_old,np_new
      REAL(8) :: rhog,vfac,kap,vpar,pidum,kaptp,kapnp,xnp
      REAL(8) :: b,th,r,enerb,cost,sint,qr,laps,sz,ter
      REAL(8) :: xt,xs,yt,xdot,ydot,zdot,pzdot,edot,pzd0,vp0
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
      REAL(8) :: grp,gxdgyp,rhox(4),rhoy(4),psp,pzp

      pidum = 1./(pi*2)**1.5*vwidth**3
      do m=1,mm(ns)
         r=x2(ns,m)-0.5*lx+lr0
         k = int(z2(ns,m)/delz)
         wz0 = ((k+1)*delz-z2(ns,m))/delz
         wz1 = 1-wz0
         th = wz0*thfnz(k)+wz1*thfnz(k+1)

         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0
         dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
                 +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
         dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
                 +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
         grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
                 +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
         radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
                 +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
         dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
                 +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
         qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
                 +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
         grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
                 +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
         gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
                 +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 
         fp = wx0*f(i)+wx1*f(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         psp = wx0*psi(i)+wx1*psi(i+1)        
         ter = wx0*t0i(i)+wx1*t0i(i+1)        
         kaptp = wx0*capts(ns,i)+wx1*capts(ns,i+1)        
         kapnp = wx0*capns(ns,i)+wx1*capns(ns,i+1)        
         xnp = wx0*xn0s(ns,i)+wx1*xn0s(ns,i+1)        

         b=1.-tor+tor*bfldp
         pzp = mims(ns)*u2(ns,m)/b-psp/br0

         rhog=sqrt(2.*b*mu(ns,m)*mims(ns))/(q(ns)*b)*iflr

         rhox(1) = rhog*(1-tor)+rhog*grp*tor
         rhoy(1) = rhog*gxdgyp/grp*tor
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
         rhox(4) = 0
         rhoy(4) = -rhoy(3)
!    calculate avg. e-field...
!    do 1,2,4 point average, where lr is the no. of points...

         phip=0.
         exp1=0.
         eyp=0.
         ezp=0.
         delbxp=0.
         delbyp=0.
         dpdzp = 0.
         dadzp = 0.
         aparp = 0.

!  4 pt. avg. done explicitly for vectorization...
         do 200 l=1,lr(1)
!
            xs=x2(ns,m)+rhox(l) !rwx(1,l)*rhog
            yt=y2(ns,m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
!
!   particle can go out of bounds during gyroavg...
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            include "ppushngp.h"
 200     continue
         exp1 = exp1/4.
         eyp = eyp/4.
         ezp = ezp/4.
         delbxp = delbxp/4.
         delbyp = delbyp/4.
         dpdzp = dpdzp/4.
         dadzp = dadzp/4.
         aparp = aparp/4.
!
         vfac = 0.5*(mims(ns)*u2(ns,m)**2 + 2.*mu(ns,m)*b)
         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))
         kap = kapnp - (1.5-vfac/ter)*kaptp

         vpar = u2(ns,m)-q(ns)/mims(ns)*aparp*nonlin*0.
         enerb=(mu(ns,m)+mims(ns)*vpar*vpar/b)/q(ns)*tor
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
         xdot = vxdum*nonlin -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlin &
             +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0
         zdot =  vpar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         pzd0 = tor*(-mu(ns,m)/mims(ns)/radiusp/bfldp*psipp*dbdtp*grcgtp)
         pzdot = pzd0 + (q(ns)/mims(ns)*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
             +q(ns)/mims(ns)*(-xdot*delbyp+ydot*delbxp+zdot*dadzp))*ipara

         edot = xdot*exp1+(ydot-vp0)*eyp+zdot*ezp                      &
             +q(ns)*pzdot*aparp*tor     &
             +q(ns)*vpar*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)    &
             -q(ns)*vpar*delbxp*vp0

         x3(ns,m) = x2(ns,m) + 0.5*dt*xdot
         y3(ns,m) = y2(ns,m) + 0.5*dt*ydot
         z3(ns,m) = z2(ns,m) + 0.5*dt*zdot
         u3(ns,m) = u2(ns,m) + 0.5*dt*pzdot

         dum = 1-w2(ns,m)*nonlin*0.
         if(isuni.eq.1)dum = exp(-vfac)*pidum
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
!         vxdum = eyp+vpar/b*delbxp
         w3(ns,m)=w2(ns,m) + 0.5*dt*(vxdum*kap + edot/ter)*dum*xnp
         
!         if(x3(ns,m)>lx .or. x3(ns,m)<0.)w3(ns,m) = 0.


!         go to 333
         if(abs(pzp-pzi(ns,m))>3.0.or.abs(vfac-eki(ns,m))>0.2*eki(ns,m))then
            x3(ns,m) = xii(ns,m)
            z3(ns,m) = z0i(ns,m)
            r = x3(ns,m)-lx/2+lr0
            k = int(z3(ns,m)/delz)
            wz0 = ((k+1)*delz-z3(ns,m))/delz
            wz1 = 1-wz0
            th = wz0*thfnz(k)+wz1*thfnz(k+1)

            i = int((r-rin)/dr)
            wx0 = (rin+(i+1)*dr-r)/dr
            wx1 = 1.-wx0
            k = int((th+pi)/dth)
            wz0 = (-pi+(k+1)*dth-th)/dth
            wz1 = 1.-wz0
            b = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
            u3(ns,m) = sqrt(2/mims(ns)*abs(eki(ns,m)-mu(ns,m)*b))
            u2(ns,m) = u3(ns,m)
            w3(ns,m) = 0.
            w2(ns,m) = 0.
            x2(ns,m) = x3(ns,m)
            z2(ns,m) = z3(ns,m)
         end if

 333     continue
         laps=anint((z3(ns,m)/lz)-.5)*(1-peritr)
         r=x3(ns,m)-0.5*lx+lr0
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         qr = wx0*sf(i)+wx1*sf(i+1)
         y3(ns,m)=dmod(y3(ns,m)-laps*2*pi*qr*lr0/q0+800.*ly,ly)
         if(x3(ns,m)>lx)then
            x3(ns,m) = lx-1.e-8
            z3(ns,m)=lz-z3(ns,m)
            x2(ns,m) = x3(ns,m)
            z2(ns,m) = z3(ns,m)
            w2(ns,m) = 0.
            w3(ns,m) = 0.
         end if
         if(x3(ns,m)<0.)then
            x3(ns,m) = 1.e-8
            z3(ns,m)=lz-z3(ns,m)
            x2(ns,m) = x3(ns,m)
            z2(ns,m) = z3(ns,m)
            w2(ns,m) = 0.
            w3(ns,m) = 0.
         end if
         z3(ns,m)=dmod(z3(ns,m)+8.*lz,lz)
         x3(ns,m)=dmod(x3(ns,m)+8.*lx,lx)         
         x3(ns,m) = min(x3(ns,m),lx-1.0e-8)
         y3(ns,m) = min(y3(ns,m),ly-1.0e-8)
         z3(ns,m) = min(z3(ns,m),lz-1.0e-2)
         
      enddo

      np_old=mm(ns)
      call init_pmove(z3(ns,:),np_old,lz,ierr)

      call pmove(x2(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x3(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y2(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y3(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z2(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z3(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u2(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u3(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w2(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w3(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mu(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call pmove(xii(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0i(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(pzi(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(eki(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call end_pmove(ierr)
      mm(ns)=np_new

      return
      end

!-----------------------------------------------------------------------

      subroutine cpush(n,ns)

      use gem_com
      use equil
      implicit none
      INTEGER :: n
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,w3old
      INTEGER :: m,i,j,k,ns,l
      INTEGER :: np_old,np_new
      REAL(8) :: vfac,vpar,vxdum,dum,xdot,ydot,zdot,pzdot,edot,pzd0,vp0
      REAL(8) :: rhog,xt,yt,zt,kap,xs,pidum,dum1,kaptp,kapnp,xnp
      REAL(8) :: b,th,r,enerb,cost,sint,qr,laps,sz,ter
      REAL(8) :: myke,mypfl,myavewi
      REAL(8) :: myefl,mynos
      REAL(8) :: ketemp,pfltemp
      REAL(8) :: efltemp,nostemp
      REAL(8) :: sbuf(10),rbuf(10)
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
      REAL(8) :: grp,gxdgyp,rhox(4),rhoy(4),psp,pzp

      sbuf(1:10) = 0.
      rbuf(1:10) = 0.
      myavewi = 0.
      myke=0.    
      mypfl=0.    
      myefl=0. 
      mynos=0.   
      ketemp=0.
      pfltemp=0.
      efltemp=0.
      nostemp=0.
      pidum = 1./(pi*2)**1.5*vwidth**3

      do m=1,mm(ns)
         r=x3(ns,m)-0.5*lx+lr0

         k = int(z3(ns,m)/delz)
         wz0 = ((k+1)*delz-z3(ns,m))/delz
         wz1 = 1-wz0
         th = wz0*thfnz(k)+wz1*thfnz(k+1)

         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0
         dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
                 +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
         dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
                 +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
         grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
                 +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
         radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
                 +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
         dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
                 +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
         qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
                 +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
         grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
                 +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
         gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
                 +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 
         fp = wx0*f(i)+wx1*f(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         psp = wx0*psi(i)+wx1*psi(i+1)        
         ter = wx0*t0i(i)+wx1*t0i(i+1)        
         kaptp = wx0*capts(ns,i)+wx1*capts(ns,i+1)        
         kapnp = wx0*capns(ns,i)+wx1*capns(ns,i+1)        
         xnp = wx0*xn0s(ns,i)+wx1*xn0s(ns,i+1)        
         b=1.-tor+tor*bfldp
         pzp = mims(ns)*u3(ns,m)/b-psp/br0

         rhog=sqrt(2.*b*mu(ns,m)*mims(ns))/(q(ns)*b)*iflr

         rhox(1) = rhog*(1-tor)+rhog*grp*tor
         rhoy(1) = rhog*gxdgyp/grp*tor
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
         rhox(4) = 0
         rhoy(4) = -rhoy(3)
!    calculate avg. e-field...
!    do 1,2,4 point average, where lr is the no. of points...

         phip=0.
         exp1=0.
         eyp=0.
         ezp=0.
         delbxp = 0.
         delbyp = 0.
         dpdzp = 0.
         dadzp = 0.
         aparp = 0.

!  4 pt. avg. written out explicitly for vectorization...
         do 200 l=1,lr(1)
            xs=x3(ns,m)+rhox(l) !rwx(1,l)*rhog
            yt=y3(ns,m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
!   BOUNDARY
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            include "cpushngp.h"
 200     continue

         exp1=exp1/4.
         eyp=eyp/4.
         ezp=ezp/4.
         delbxp=delbxp/4.
         delbyp=delbyp/4.
         dpdzp = dpdzp/4.
         dadzp = dadzp/4.
         aparp = aparp/4.


         vfac = 0.5*(mims(ns)*u3(ns,m)**2 + 2.*mu(ns,m)*b)
         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))
         kap = kapnp -(1.5-vfac/ter)*kaptp

         vpar = u3(ns,m)-q(ns)/mims(ns)*aparp*nonlin*0.
         enerb=(mu(ns,m)+mims(ns)*vpar*vpar/b)/q(ns)*tor
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
         xdot = vxdum*nonlin -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlin     &
             +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0
         zdot =  vpar*(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         pzd0 = tor*(-mu(ns,m)/mims(ns)/radiusp/bfldp*psipp*dbdtp*grcgtp)
         pzdot = pzd0 + (q(ns)/mims(ns)*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
             +q(ns)/mims(ns)*(-xdot*delbyp+ydot*delbxp+zdot*dadzp))*ipara

         edot = xdot*exp1+(ydot-vp0)*eyp+zdot*ezp                      &
             +q(ns)*pzdot*aparp*tor     &
             +q(ns)*vpar*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)   &
             -q(ns)*vpar*delbxp*vp0

         x3(ns,m) = x2(ns,m) + dt*xdot
         y3(ns,m) = y2(ns,m) + dt*ydot
         z3(ns,m) = z2(ns,m) + dt*zdot
         u3(ns,m) = u2(ns,m) + dt*pzdot

         dum = 1-w3(ns,m)*nonlin*0.
         if(isuni.eq.1)dum = exp(-vfac)*pidum
!         vxdum = eyp+vpar/b*delbxp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
         w3old = w3(ns,m)
         w3(ns,m) = w2(ns,m) + dt*(vxdum*kap+edot/ter)*dum*xnp

         if(abs(w3(ns,m)).gt.1..and.nonlin==1)then
            w3(ns,m) = 0.
            w2(ns,m) = 0.
         end if


         laps=anint((z3(ns,m)/lz)-.5)*(1-peritr)
         r=x3(ns,m)-0.5*lx+lr0
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         qr = wx0*sf(i)+wx1*sf(i+1)
         y3(ns,m)=dmod(y3(ns,m)-laps*2*pi*qr*lr0/q0+800.*ly,ly)
         if(x3(ns,m)>lx)then
            x3(ns,m) = lx-1.e-8
            z3(ns,m)=lz-z3(ns,m)
            x2(ns,m) = x3(ns,m)
            z2(ns,m) = z3(ns,m)
            w2(ns,m) = 0.
            w3(ns,m) = 0.
         end if
         if(x3(ns,m)<0.)then
            x3(ns,m) = 1.e-8
            z3(ns,m)=lz-z3(ns,m)
            x2(ns,m) = x3(ns,m)
            z2(ns,m) = z3(ns,m)
            w2(ns,m) = 0.
            w3(ns,m) = 0.
         end if
         z3(ns,m)=dmod(z3(ns,m)+8.*lz,lz)
         x3(ns,m)=dmod(x3(ns,m)+8.*lx,lx)         
         x3(ns,m) = min(x3(ns,m),lx-1.0e-8)
         y3(ns,m) = min(y3(ns,m),ly-1.0e-8)
         z3(ns,m) = min(z3(ns,m),lz-1.0e-2)

!     particle diagnostics done here because info is available...
         mypfl=mypfl + w3old*(eyp+vpar*delbxp/b) 
         myefl=myefl + vfac*w3old*(eyp+vpar*delbxp/b)
         myke=myke + vfac*w3(ns,m)
         mynos=mynos + w3(ns,m)
         myavewi = myavewi+abs(w3(ns,m))

!     xn+1 becomes xn...
         u2(ns,m)=u3(ns,m)
         x2(ns,m)=x3(ns,m)
         y2(ns,m)=y3(ns,m)
         z2(ns,m)=z3(ns,m)
         w2(ns,m)=w3(ns,m)

!     100     continue
      enddo

      sbuf(1)=myke
      sbuf(2)=myefl
      sbuf(3)=mypfl
      sbuf(4)=mynos
      sbuf(5)=myavewi
      call MPI_ALLREDUCE(sbuf,rbuf,10,  &
          MPI_REAL8,MPI_SUM,           &
          MPI_COMM_WORLD,ierr)

      ketemp=rbuf(1)
      efltemp=rbuf(2)
      pfltemp=rbuf(3)
      nostemp=rbuf(4)
      avewi(n) = rbuf(5)/( float(tmm(1)) )
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      nos(1,n)=nostemp/( float(tmm(1)) )
      pfl(1,n)=pfltemp/( float(tmm(1)) )
      efl(1,n)=mims(ns)/tets(1)*efltemp/( float(tmm(1)) )
      ke(1,n)=ketemp/( 2.*float(tmm(1))*mims(ns) )
      np_old=mm(ns) 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call init_pmove(z3(ns,:),np_old,lz,ierr)
!     
      call pmove(x2(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x3(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y2(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y3(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z2(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z3(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u2(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u3(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w2(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w3(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mu(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xii(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0i(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(pzi(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(eki(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call end_pmove(ierr)
      mm(ns)=np_new
!     write(*,*)MyId,mm(ns)

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grad(ip)
  
!  currently set up for periodic in x,y,z

      use gem_com
      use equil
      implicit none
      INTEGER :: i,j,k,ip
      real(8) :: ux(0:imx,0:jmx,0:1),uy(0:imx,0:jmx,0:1)
      real(8) :: tmp(0:imx,0:jmx,0:1),uoverb(0:imx,0:jmx,0:1)

      if(idg==1)write(*,*)'enter grad'
      call gradu(phi(:,:,:),ux,uy)
      ex(:,:,:) = -ux(:,:,:)
      ey(:,:,:) = -uy(:,:,:)

      delbx = 0.
      delby = 0.
      if(ifluid.eq.1)then
         call gradu(apar(:,:,:),ux,uy)
         delbx(:,:,:) = uy(:,:,:)
         delby(:,:,:) = -ux(:,:,:)
      end if

      call gradu(dene(:,:,:),ux,uy)
      dnedx(:,:,:) = ux(:,:,:)
      dnedy(:,:,:) = uy(:,:,:)
      do i = 0,imx
         do j = 0,jmx
            do k = 0,1
               uoverb(i,j,k) = upar(i,j,k)  !/bfld(i,k)
            end do
         end do
      end do
      call gradu(uoverb(:,:,:),ux,uy)
      dupadx(:,:,:) = ux(:,:,:)
      dupady(:,:,:) = uy(:,:,:)

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grid1(ip,n)

!    source quantities are are calculated: n_i
!    right now only ion quantitities are calculated...

      use gem_com
      use equil
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: enerb,vxdum,dum,xdot,ydot
      INTEGER :: m,n,i,j,k,l,ns,ip
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,vte
      REAL(8) :: sz,wght,r,th,cost,sint,b,qr,dv
      REAL(8) :: xt,yt,rhog,pidum,vpar,xs,dely,vfac
      REAL(8) :: lbfs(0:imx,0:jmx)
      REAL(8) :: rbfs(0:imx,0:jmx)
      REAL(8) :: lbfr(0:imx,0:jmx)
      REAL(8) :: rbfr(0:imx,0:jmx)
      real(8) :: myden(0:imx,0:jmx,0:1),myjpar(0:imx,0:jmx,0:1)
      real(8) :: myhpar(0:imx,0:jmx,0:1),myhden(0:imx,0:jmx,0:1)
      real(8) :: mydene(0:imx,0:jmx,0:1),myupar(0:imx,0:jmx,0:1)
      real(8) :: mydti(0:imx,0:jmx,0:1),mydte(0:imx,0:jmx,0:1)
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
      REAL(8) :: grp,gxdgyp,rhox(4),rhoy(4)

      rho=0.
      jion = 0.
      myhpar = 0.
      myhden = 0.
      mydte = 0.
      ns=1
if(idg.eq.1)write(*,*)'enter ion grid1',mm(1)
      do ns = 1,nsm
      den(ns,:,:,:)=0.
      jpar(ns,:,:,:)=0.
      myden = 0.
      myjpar = 0.
      mydti = 0.
      do m=1,mm(ns)
         dv=float(lr(1))*(dx*dy*dz)
         r=x3(ns,m)-0.5*lx+lr0

         k = int(z3(ns,m)/delz)
         wz0 = ((k+1)*delz-z3(ns,m))/delz
         wz1 = 1-wz0
         th = wz0*thfnz(k)+wz1*thfnz(k+1)

         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0
         dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
                 +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
         dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
                 +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
         grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
                 +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
         radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
                 +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
         dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
                 +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
         qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
                 +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
         grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
                 +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
         gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
                 +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 
         fp = wx0*f(i)+wx1*f(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
!         b=1.-lr0/br0*cost
         b=1.-tor+tor*bfldp

         rhog=sqrt(2.*b*mu(ns,m)*mims(ns))/(q(ns)*b)*iflr

         rhox(1) = rhog*(1-tor)+rhog*grp*tor
         rhoy(1) = rhog*gxdgyp/grp*tor
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
         rhox(4) = 0
         rhoy(4) = -rhoy(3)

         vfac=0.5*(mims(ns)*u3(ns,m)**2 + 2.*mu(ns,m)*b )
         wght=w3(ns,m)/dv

         vpar = u3(ns,m) !linearly correct

!    now do 1,2,4 point average, where lr is the no. of points...
         do 200 l=1,lr(1)
            xs=x3(ns,m)+rhox(l) !rwx(1,l)*rhog
            yt=y3(ns,m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.e-8)
            yt = min(yt,ly-1.e-8)

            include "gridngp.h"
 200     continue
      enddo
if(idg.eq.1)write(*,*)myid,'pass ion grid1'
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!   enforce periodicity
      call enforce(myden(:,:,:))
      call enforce(myjpar)
      call enforce(mydti)
!      call filter(myden(:,:,:))
!      call filter(myjpar(:,:,:))

      do 410 i=0,im
         do 420 j=0,jm
            do 430 k=0,mykm
               den(ns,i,j,k)=q(ns)*myden(i,j,k)/n0/jac(i,k)*cn0s(ns)
               jpar(ns,i,j,k) = q(ns)*myjpar(i,j,k)/n0/jac(i,k)*cn0s(ns)
               mydti(i,j,k) = mydti(i,j,k)/n0/jac(i,k)*cn0s(ns)
 430        continue
 420     continue
 410  continue
      mydti(:,:,:) = mydti(:,:,:)-den(ns,:,:,:)
      call MPI_ALLREDUCE(mydti(0:im,0:jm,0:1),  &
     		dti(ns,0:im,0:jm,0:1),             &
     		(imx+1)*(jmx+1)*2,MPI_REAL8,       &
     		MPI_SUM,GRID_COMM,ierr)

      do 450 i=0,im
         do 460 j=0,jm
            do 470 k=0,mykm
               rho(i,j,k)=rho(i,j,k)+den(ns,i,j,k)*q(ns)
               jion(i,j,k) = jion(i,j,k)+jpar(ns,i,j,k)
 470        continue
 460     continue
 450  continue
end do

! electrons density and current
      vte = sqrt(amie)
      dene(:,:,:) = 0.
      upar(:,:,:) = 0.
      mydene = 0.
      myupar = 0.
if(idg.eq.1)write(*,*)'enter electron grid1'
      do m=1,mme
         dv=(dx*dy*dz)
         wght=w3e(m)/dv
         vpar = u3e(m) !linearly correct
         if(abs(vpar/vte).gt.vcut)wght = 0.

         xt=x3e(m)
         yt=y3e(m)
         include 'gridlie.h'
      enddo
if(idg.eq.1)write(*,*)'pass electron grid1'
!   enforce periodicity
      call enforce(mydene(:,:,:))
      call enforce(myupar(:,:,:))
!      call filter(mydene(:,:,:))
!      call filter(myupar(:,:,:))

      do  i=0,im
            do  j=0,jm
               do  k=0,mykm
                  dene(i,j,k)= mydene(i,j,k)/n0/jac(i,k)*cn0e*ifluid
                  upar(i,j,k) = myupar(i,j,k)/n0/jac(i,k)*cn0e*ifluid
               end do
            end do
         end do
 999  continue

      do i = 0,im
         do j = 0,jm
            do k = 0,mykm
               rho(i,j,k) = ision*rho(i,j,k) + dene(i,j,k)*qel
            enddo
         enddo
      enddo      
      if(ishgk==0)return
      do m=1,mmb
         dv=float(lr(1))*(dx*dy*dz)
         r=xh3(m)-0.5*lx+lr0

         k = int(zh3(m)/delz)
         wz0 = ((k+1)*delz-zh3(m))/delz
         wz1 = 1-wz0
         th = wz0*thfnz(k)+wz1*thfnz(k+1)

         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0
         dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
                 +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
         dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
                 +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
         grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
                 +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
         radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
                 +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
         dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
                 +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
         qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
                 +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
         grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
                 +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
         gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
                 +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 
         fp = wx0*f(i)+wx1*f(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
!         b=1.-lr0/br0*cost
         b=1.-tor+tor*bfldp

         rhog=sqrt(2.*b*muh(m)*mbeam)/(qbeam*b)*iflrh

         rhox(1) = rhog*(1-tor)+rhog*grp*tor
         rhoy(1) = rhog*gxdgyp/grp*tor
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
         rhox(4) = 0
         rhoy(4) = -rhoy(3)

         vfac=0.5*(mbeam*uh3(m)**2 + 2.*muh(m)*b )
         wght=wh3(m)/dv*nh

         vpar = uh3(m) !linearly correct

!    now do 1,2,4 point average, where lr is the no. of points...
         do 500 l=1,lr(1)
            xs=xh3(m)+rhox(l) !rwx(1,l)*rhog
            yt=yh3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.e-8)
            yt = min(yt,ly-1.e-8)

            i=int(xt/dx+0.5)
            j=int(yt/dy+0.5)
            k=int(zh3(m)/dz+0.5)-gclr*kcnt
     
            myhpar(i,j,k) = myhpar(i,j,k)+wght*vpar
            myhden(i,j,k) = myhden(i,j,k)+wght
 500     continue
      enddo
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call enforce(myhpar)
      call enforce(myhden)

      do 510 i=0,im
         do 520 j=0,jm
            do 530 k=0,mykm
               hpar(i,j,k) = qbeam*myhpar(i,j,k)/n0/jac(i,k)
               hden(i,j,k) = qbeam*myhden(i,j,k)/n0/jac(i,k)
 530        continue
 520     continue
 510  continue

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine modes2(u,modehis,n)

!     calculate mode histories. calculate modes for u at timestep n
!     and store in modehis(mode,n).

      use gem_com
      use equil
      use fft_wrapper
      implicit none
!     
      REAL(8) :: u(0:imx,0:jmx,0:1)
      COMPLEX(8) :: modehis(modemx,0:nmx)
      COMPLEX(8) :: modebuf
      INTEGER :: n,i,j,k,l,m,ii

      INTEGER :: mode,jj,thek,oproc,ifirst

!     

      if(n.eq.0) return

      do 100 mode=1,modem
         oproc=int(nmode(mode)/kcnt*ntube)

         if (MyId.eq.oproc) then
            thek=0
            do j=0,jm-1
               do i=0,im-1
                  tmpx(i)=u(i,j,thek)
               enddo

!     FT in x....
               call ccfft('x',1,imx,1.,tmpx,tmpx,coefx,workx,0)
               ii=lmode(mode) !+1
               if(lmode(mode).lt.0) write(*,*) 'lmode < 0, error'
               tmpy(j)=tmpx(ii)/float(im)
            enddo

!     FT in y....
            call ccfft('y',1,jmx,1.,tmpy,tmpy,coefy,worky,0)
            jj=mmode(mode)  !+1
            if(mmode(mode).lt.0) write(*,*) 'mmode < 0, error'
            modebuf=tmpy(jj)/float(jm)

         endif

         call MPI_BCAST(modebuf,1,MPI_COMPLEX16,oproc, &
             MPI_COMM_WORLD,ierr)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!         write(*,*)myid,modebuf
         modehis(mode,n)=modebuf
 100  continue
!     
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine restart(iflag,n)

      use gem_com
      use equil
      implicit none
      INTEGER :: m,ns,iflag,n,i,j,k,ip
      character*70 fname
      character*3 holdmyid

      write(holdmyid,'(I3.3)') MyId
      fname=directory//'dump_'//holdmyid//'.b'
      if(iflag.eq.1) then
         open(139+MyId,file=fname,form='unformatted',status='old')
         read(139+MyId)ncurr,tcurr,rmpp,rmaa
         do ns = 1,nsm
         read(139+MyId)mm(ns)
         do 110 m=1,mm(ns)
            read(139+MyId) mu(ns,m)
            read(139+MyId) x2(ns,m),y2(ns,m),z2(ns,m),u2(ns,m),w2(ns,m)
            read(139+MyId) xii(ns,m),z0i(ns,m),pzi(ns,m),eki(ns,m)
            w2(ns,m)=w2(ns,m)/cut
            x3(ns,m)=x2(ns,m)
            y3(ns,m)=y2(ns,m)
            z3(ns,m)=z2(ns,m)
            u3(ns,m)=u2(ns,m)
            w3(ns,m)=w2(ns,m)
 110     continue
         end do
         ns=3
         read(139+MyId)mm(ns)
         do 120 m=1,mm(ns)
            read(139+MyId) muh(m)
            read(139+MyId) xh2(m),yh2(m),zh2(m),uh2(m),wh2(m)
            read(139+MyId) xih(m),z0h(m),pzh(m),ekh(m)
            wh2(m)=wh2(m)/cut
            xh3(m)=xh2(m)
            yh3(m)=yh2(m)
            zh3(m)=zh2(m)
            uh3(m)=uh2(m)
            wh3(m)=wh2(m)
 120     continue

         read(139+MyId)mme
         do  m=1,mme
            read(139+MyId) x2e(m),y2e(m),z2e(m),u2e(m),w2e(m),mue2(m),ipass(m)
            read(139+MyId) xie(m),z0e(m),pze(m),eke(m),mue(m)
            w2e(m)=w2e(m)/cut
            x3e(m)=x2e(m)
            y3e(m)=y2e(m)
            z3e(m)=z2e(m)
            u3e(m)=u2e(m)
            w3e(m)=w2e(m)
            mue3(m)=mue2(m)
         end do

         read(139+myid)ke,fe,te,pfl,efl,nos,rmsphi,pmodehis
         close(139+MyId)
      endif

      if(iflag.eq.2) then
         open(139+MyId,file=fname,form='unformatted',status='unknown')
         write(139+MyId)n+1,tcurr-dt,rmpp,rmaa
         do ns = 1,nsm
         write(139+MyId)mm(ns)
         do 210 m=1,mm(ns)
            write(139+MyId) mu(ns,m)
            write(139+MyId) x2(ns,m),y2(ns,m),z2(ns,m),u2(ns,m),w2(ns,m)
            write(139+MyId) xii(ns,m),z0i(ns,m),pzi(ns,m),eki(ns,m)
 210     continue
         end do
         ns=3
         write(139+MyId)mm(ns)
         do 220 m=1,mm(ns)
            write(139+MyId) muh(m)
            write(139+MyId) xh2(m),yh2(m),zh2(m),uh2(m),wh2(m)
            write(139+MyId) xih(m),z0h(m),pzh(m),ekh(m)
 220     continue

         write(139+MyId)mme
         do  m=1,mme
            write(139+MyId) x2e(m),y2e(m),z2e(m),u2e(m),w2e(m),mue2(m),ipass(m)
            write(139+MyId) xie(m),z0e(m),pze(m),eke(m),mue(m)
         end do
         
         write(139+myid)ke,fe,te,pfl,efl,nos,rmsphi,pmodehis
         close(139+MyId)
      endif
!     
      return
      end
!-----------------------------------------------------------------------

!        Normal distribution random no. generator, stand. dev. = 1.
!        Version 2 does it Willy's way...

         subroutine parperp(vpar,vperp2,m,pi,cnt,MyId)

         REAL(8) :: vpar,vperp2,r1,r2,t,pi
         INTEGER :: m,iflag,cnt,MyId
         REAL(8) :: c0,c1,c2
         REAL(8) :: d1,d2,d3
         data c0,c1,c2/2.515517,0.802853,0.010328/
         data d1,d2,d3/1.432788,0.189269,0.001308/


          r1=revers(m+MyId*cnt,7)
          r2=revers(m+MyId*cnt,11)


!.....quiet start---see denavit pf '71(?) & abramowitz hand book
!.....fibonacci start---see denavit comm. pla. phy. & con. fus. '81
! warning: we have g1=1 in the x-direction. This surpresses all odd
!          modes in the x-direction!!!

         iflag=1
         if(r1.le.0.5) go to 110
         r1=1.-r1
         iflag=-1
  110    continue
         if(r1.ge.1.e-6) then
           t=sqrt(log(1.0/(r1*r1)))
         else
           t=5.0
           write(*,*)'parperp2 warning  m= ',m
         endif

         vpar=t-(c0+c1*t+c2*t**2)/(1.+d1*t+d2*t**2+d3*t**3)
         vpar=vpar*iflag

          vperp2=-2.0*dlog(r2)

        return
        end

!---------------------------------------------------------------------- 
!    calculate weights and delta j for proper periodicity 
!    in the toroidal direction, weights and delta j are
!    dependant only on minor radius, hence x, hence
!    weight is a vector in 0:imx

      subroutine weight

      use gem_com
      use equil
      implicit none
      INTEGER :: i,j
      REAL(8) :: dely,aweight,r,qr,wx0,wx1

!         peritr = 0
      if (GCLR.eq.Master.and.peritr.eq.0) then
	 do i=0,im
            r=float(i)*dx-0.5*lx+lr0
            j = int((r-rin)/dr)
            j = min(j,nr-1)
            wx0 = (rin+(j+1)*dr-r)/dr
            wx1 = 1.-wx0
            qr = wx0*sf(j)+wx1*sf(j+1)
            dely=dmod(2.*pi*lr0/q0*qr,ly)
            deljp(i)=int(dely/dy)
            deljm(i)=0
            aweight=dmod(dely,dy)
            weightp(i)=1.-aweight/dy
            weightm(i)=1.
	 enddo	
      elseif (GCLR.eq.GLST.and.peritr.eq.0) then
	 do i=0,im
            r=float(i)*dx-0.5*lx+lr0
            j = int((r-rin)/dr)
            j = min(j,nr-1)
            wx0 = (rin+(j+1)*dr-r)/dr
            wx1 = 1.-wx0
            qr = wx0*sf(j)+wx1*sf(j+1)
            dely=dmod(2.*pi*lr0/q0*qr,ly)
            deljm(i)=int(dely/dy)
            deljp(i)=0
            aweight=dmod(dely,dy)
            weightm(i)=1.-aweight/dy
            weightp(i)=1.
	 enddo	
      else
	 do i=0,im
            deljp(i)=0
            deljm(i)=0
            weightp(i)=1.
            weightm(i)=1.
	 enddo
      endif

      do j=0,jm
         do i=0,im
            jpl(i,j)=mod(j+deljp(i)+8*jm,jm)
            jpn(i,j)=mod(j+deljp(i)+1+8*jm,jm)
            jmi(i,j)=mod(j-deljm(i)+8*jm,jm)
            jmn(i,j)=mod(j-deljm(i)-1+8*jm,jm)
            weightpn(i)=1.-weightp(i)
            weightmn(i)=1.-weightm(i)
         enddo
      enddo

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gkps(nstep,ip)   
      use gem_com
      use equil
      use fft_wrapper
      implicit none
      REAL(8) :: b,b1,b2,gam0,gam1,delyz,th,bf,dum,r,qr,shat,ter
      REAL(8) :: kx1,kx2,ky
      REAL(8),dimension(:),allocatable :: akx,aky
      real(8),dimension(:,:,:,:),allocatable:: gamb1,gamb2
      complex(8),dimension(:,:,:,:),allocatable :: mx
      integer,dimension(:,:,:,:),allocatable :: ipiv
      REAL(8),dimension(:,:,:),allocatable :: formphi,formfe
      REAL(8) :: sgnx,sgny,sz,myfe
      INTEGER :: i,i1,j,j1,k,k1,l,m,n,ifirst,nstep,ip,INFO
      INTEGER :: l1,m1,myk,myj,izonal,ix,ikx
      COMPLEX(8) :: temp3dxy(0:imx-1,0:jmx-1,0:1),v(0:imx-1,0:jcnt-1,0:1)
      COMPLEX(8) :: sbuf(0:imx*jcnt*2-1),rbuf(0:imx*jmx*2-1)
      COMPLEX(8) :: sl(1:imx-1,0:jcnt-1,0:1)
      COMPLEX(8) :: aphik(0:nxpp-1),myaphik(0:nxpp-1)
      real(8) :: myaph(0:nxpp),aph(0:nxpp),u(0:imx,0:jmx,0:1)
      complex(8) :: cdum
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
      REAL(8) :: grp,gxdgyp,grdgtp,gthp
      REAL(8) :: wx0,wx1,wz0,wz1
      character*70 fname
      character*3 holdmyid

      save formphi,formfe,ifirst,akx,aky,mx,gamb1,gamb2,ipiv
      if(idg==1)write(*,*)'enter gkps'
      izonal = 1

      write(holdmyid,'(I3.3)') MyId
      fname='./matrix/'//'mx_phi_'//holdmyid
!     form factors....
      if (ifirst.ne.-99) then
         allocate(akx(0:imx-1),aky(0:jcnt-1), &
                  gamb1(0:imx-1,0:jcnt-1,0:imx-1,0:1), &
                  gamb2(0:imx-1,0:jcnt-1,0:imx-1,0:1))
         allocate(mx(imx-1,imx-1,0:jcnt-1,0:1),formphi(0:imx-1,0:jcnt-1,0:1))
         allocate(formfe(0:imx-1,0:jcnt-1,0:1),ipiv(imx-1,imx-1,0:jcnt-1,0:1))

         if(iget.eq.1) then
            open(10000+MyId,file=fname,form='unformatted',status='old')
            read(10000+MyId)mx,ipiv
            goto 200
         end if

         do 51 l=0,im-1
            do 52 m=0,jcnt-1
               j = tclr*jcnt+m
               do 54 i=0,im-1 
                  r = lr0-lx/2+i*dx
                  i1 = int((r-rin)/dr)
                  i1 = min(i1,nr-1)
                  wx0 = (rin+(i1+1)*dr-r)/dr
                  wx1 = 1.-wx0
                  ter = wx0*t0i(i1)+wx1*t0i(i1+1)
                  do 53 n = 0,1
                     k = gclr*kcnt+n
                     k1 = int(k*dz/delz)
                     k1 = min(k1,ntheta-1)
                     wz0 = ((k1+1)*delz-dz*k)/delz
                     wz1 = 1-wz0
                     th = wz0*thfnz(k1)+wz1*thfnz(k1+1)
                     k = int((th+pi)/dth)
                     k = min(k,ntheta-1)
                     wz0 = (-pi+(k+1)*dth-th)/dth
                     wz1 = 1.-wz0
                     bfldp = wx0*wz0*bfld(i1,k)+wx0*wz1*bfld(i1,k+1) &
                        +wx1*wz0*bfld(i1+1,k)+wx1*wz1*bfld(i1+1,k+1) 
                     dydrp = wx0*wz0*dydr(i1,k)+wx0*wz1*dydr(i1,k+1) &
                        +wx1*wz0*dydr(i1+1,k)+wx1*wz1*dydr(i1+1,k+1) 
                     qhatp = wx0*wz0*qhat(i1,k)+wx0*wz1*qhat(i1,k+1) &
                        +wx1*wz0*qhat(i1+1,k)+wx1*wz1*qhat(i1+1,k+1) 
                     grp = wx0*wz0*gr(i1,k)+wx0*wz1*gr(i1,k+1) &
                        +wx1*wz0*gr(i1+1,k)+wx1*wz1*gr(i1+1,k+1) 
                     gthp = wx0*wz0*gth(i1,k)+wx0*wz1*gth(i1,k+1) &
                        +wx1*wz0*gth(i1+1,k)+wx1*wz1*gth(i1+1,k+1) 
                     gxdgyp = wx0*wz0*gxdgy(i1,k)+wx0*wz1*gxdgy(i1,k+1) &
                        +wx1*wz0*gxdgy(i1+1,k)+wx1*wz1*gxdgy(i1+1,k+1) 
                     grdgtp = wx0*wz0*grdgt(i1,k)+wx0*wz1*grdgt(i1,k+1) &
                        +wx1*wz0*grdgt(i1+1,k)+wx1*wz1*grdgt(i1+1,k+1) 

                     if(j.ge.(jm/2+1)) then
                        m1=jm-j
                        sgny=-1.
                     else
                        m1=j
                        sgny=1.
                     endif
                     ky=sgny*2.*pi*float(m1)/ly
                     kx1=pi*float(l)/lx
                     kx2=-pi*float(l)/lx
                     bf=bfldp
                     b1=mims(1)*(kx1*kx1*grp**2 + &
                       ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                       +2*dydrp*lr0/q0*qhatp*grdgtp) &
                       +2*kx1*ky*gxdgyp)/(bf*bf)*ter

                     b2=mims(1)*(kx2*kx2*grp**2 + &
                       ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                       +2*dydrp*lr0/q0*qhatp*grdgtp) &
                       +2*kx2*ky*gxdgyp)/(bf*bf)*ter

                     call srcbes(b1,gam0,gam1)
                     gamb1(l,m,i,n)=gam0
                     call srcbes(b2,gam0,gam1)
                     gamb2(l,m,i,n)=gam0

!   formfactor in gkps
                     formfe(l,m,n) = 1.-gam0
                     formphi(l,m,n) = 1./jmx 
                     if(abs(ky)>kycut)formphi(l,m,n) = 0.
!                    if(abs(ky)==0.)formphi(l,m,n) = 0.
                     if(m1.ne.mlk.and.onemd==1)formphi(l,m,n) = 0.
 53               continue
 54            continue
 52         continue
 51      continue

         do k = 0,1
            do j = 0,jcnt-1
               do i = 1,imx-1
                  r = lr0-lx/2+i*dx
                  i1 = int((r-rin)/dr)
                  i1 = min(i1,nr-1)
                  wx0 = (rin+(i1+1)*dr-r)/dr
                  wx1 = 1.-wx0
                  ter = wx0*t0i(i1)+wx1*t0i(i1+1)
                  do ix = 1,imx-1
                     mx(i,ix,j,k) = 0.
                     if(i==ix)mx(i,ix,j,k) = fradi*cn0e
                     do ikx = 0,imx-1
                        mx(i,ix,j,k) = mx(i,ix,j,k)+sin(ix*ikx*pi/imx)* & 
                             ((1-gamb1(ikx,j,i,k))*exp(IU*ikx*i*pi/imx)- &
                              (1-gamb2(ikx,j,i,k))*exp(-IU*ikx*i*pi/imx)) &
                             /ter*cn0i*gn0i(i)/(IU*imx)
                     end do
                  end do
               end do
            end do
         end do
         do k = 0,1
            do j = 0,jcnt-1
               call ZGETRF(imx-1,imx-1,mx(:,:,j,k),imx-1,ipiv(:,:,j,k),INFO )
            end do
         end do

         if(iget.eq.0) then
            open(10000+MyId,file=fname,form='unformatted',status='unknown')
            write(10000+MyId)mx,ipiv
            goto 200
         end if

         if(gclr==kmx/2.and.tclr==0.and.nstep==0)then
 !           write(*,*)'setup in gkps'
            open(20,file="mx",status='unknown')
            j = 0
            k = 0
            do i = 1,imx-1
               do ix = 1,imx-1
                  if(abs(i-ix)<40)write(20,10)i,ix,mx(i,ix,j,k),mx(ix,i,j,k)
               end do
            end do
 10      format(1x,i5,1x,i5,2(2x,e12.5,2x,e12.5))
            close(20)
         end if
 200     ifirst=-99
      endif
      if(idg==1)write(*,*)'pass form factors'

!   now do field solve...

!      phi = 0.
!      return
!      temp3dxy = 0.
      aphik = 0.
      myaphik = 0.

!  find rho(kx,ky)
      do i = 0,nxpp
         do j = 0,jm
            do k = 0,1
               u(i,j,k) = rho(i,j,k)+hden(i,j,k)*ishgk+ &
                         fradi*(phi(i,j,k)/ntube*cn0e-den0(i,j,k))
            end do
         end do
      end do
      call dcmpy(u(0:imx-1,0:jmx-1,0:1),v)
      if(idg==1)write(*,*)'pass first fft'

      do k = 0,1
         do j = 0,jcnt-1
            myj=tclr*jcnt+j
            sl(1:imx-1,j,k) = v(1:imx-1,j,k)
            call ZGETRS('N',imx-1,1,mx(:,:,j,k),imx-1,ipiv(:,:,j,k), &
                         sl(:,j,k),imx-1,INFO) 
            temp3dxy(1:imx-1,myj,k) = sl(1:imx-1,j,k)
            temp3dxy(0,myj,k) = 0.
         end do
      end do

      do k = 0,1
         do j = 0,jcnt-1
            myj = tclr*jcnt+j
            do i = 0,imx-1
               m = j*2*imx+k*imx+i 
               sbuf(m) = temp3dxy(i,myj,k)
            end do
         end do
      end do

      call mpi_allgather(sbuf(0),jcnt*imx*2,mpi_complex16,rbuf,jcnt*imx*2, &
                         mpi_complex16,grid_comm,ierr)
                         
      do k = 0,1
         do j = 0,jmx-1
            do i = 0,imx-1
               m = j*2*imx+k*imx+i 
               temp3dxy(i,j,k) = rbuf(m)
            end do
         end do
      end do
      

!  from rho(kx,ky) to phi(kx,ky)
      do k = 0,1
         do j = 0,jmx-1
            do i = 0,imx-1
               temp3dxy(i,j,k)=temp3dxy(i,j,k)/jmx !*formphi(i,j,k)
            end do
         end do
      end do

!  from phi(kx,ky) to phi(x,y)
      do k=0,mykm
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3dxy(i,j,k)
            end do
            call ccfft('y',1,jmx,1.,tmpy,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3dxy(i,j,k) = tmpy(j)  ! phi(x,y)
            end do
         end do
      end do

 100  continue
      do i = 0,nxpp-1
         do j = 0,jm-1
            do k = 0,mykm
               phi(i,j,k) = temp3dxy(i,j,k)
            end do
         end do
      end do

!    x-y boundary points 
      call enfxy(phi(:,:,:))
      call enfz(phi(:,:,:))
      if(idg==1)write(*,*)'pass enfz', myid

      return
      end

!      End of gkps....
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine eqmo(ip)
      use gem_com
      use equil
      implicit none
      real(8) :: lbfr(0:imx,0:jmx)
      real(8) :: lbfs(0:imx,0:jmx)
      real(8) :: rbfr(0:imx,0:jmx)
      real(8) :: rbfs(0:imx,0:jmx)
      integer :: i,j,k,ip
      real(8) :: eta

      ez(:,:,:) = 0.
      dadz = 0.
      dpdz = 0.
!      return
      do i = 0,im
         do j = 0,jm
            do k = 1,mykm-1
               ez(i,j,k) = (phi(i,j,k-1)-phi(i,j,k+1))/(2.*dz)
               dadz(i,j,k) = (apar(i,j,k+1)-apar(i,j,k-1))/(2.*dz)
            end do
         end do
      end do   
      
      rbfs = phi(:,:,0)

      call MPI_SENDRECV(rbfs(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8,rngbr,204,                    &
          lbfr(0,0),(imx+1)*(jmx+1),              &
          MPI_REAL8,lngbr,204,                    &
          TUBE_COMM,stat,ierr)

      lbfs=phi(:,:,1)

      call MPI_SENDRECV(lbfs(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8,lngbr,205,                    &
          rbfr(0,0),(imx+1)*(jmx+1),              &
          MPI_REAL8,rngbr,205,                    &
          TUBE_COMM,stat,ierr)                    
      call MPI_BARRIER(TUBE_COMM,ierr)             
      do i=0,im
         do j=0,jm
            ez(i,j,0)=(weightp(i)*lbfr(i,jpl(i,j)) &
                +weightpn(i)*lbfr(i,jpn(i,j))        &
                -phi(i,j,1))/(2.*dz)

            ez(i,j,1)=( phi(i,j,0)-(weightm(i)*rbfr(i,jmi(i,j)) &
                +weightmn(i)*rbfr(i,jmn(i,j))))/(2.*dz)

         enddo
      enddo

      rbfs = apar(:,:,0)

      call MPI_SENDRECV(rbfs(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8,rngbr,206,                    &
          lbfr(0,0),(imx+1)*(jmx+1),              &
          MPI_REAL8,lngbr,206,                    &
          TUBE_COMM,stat,ierr)

      lbfs = apar(:,:,1)

      call MPI_SENDRECV(lbfs(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8,lngbr,207, &
          rbfr(0,0),(imx+1)*(jmx+1),   &
          MPI_REAL8,rngbr,207,         &
          TUBE_COMM,stat,ierr)
      call MPI_BARRIER(TUBE_COMM,ierr)
      do i=0,im
         do j=0,jm
            dadz(i,j,0)=-(weightp(i)*lbfr(i,jpl(i,j))  &
                +weightpn(i)*lbfr(i,jpn(i,j))  &
                -apar(i,j,1))/(2.*dz)

            dadz(i,j,1)=-( apar(i,j,0)-(weightm(i)*rbfr(i,jmi(i,j)) &
                +weightmn(i)*rbfr(i,jmn(i,j))))/(2.*dz)
         enddo
      enddo

      dpdz(:,:,:) = -ez(:,:,:)

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine spec(n)
      use gem_com
      use equil
      implicit none
      integer :: i,j,k,l,m,n

      i = tcurr-dt
      if(myid.eq.master)then
         if(izon.eq.0)then
            write(*,10)i,rmsphi(n),rmsapa(n),efl(1,n),efl(2,n),avewe(n)
!         write(*,11)pfl(1,n),pfl(2,n),pfl(3,n),efl(1,n),efl(2,n)
         end if
 10      format(1x,i6,5(2x,e10.3),2x,i7,2x,i3)
 11      format(6x,5(2x,e12.5))
 12      format(1x,i6,5(2x,e12.5))

         write(9,10)i,rmsphi(n),rmsapa(n),avewe(n),avewi(n),avewh(n)
         write(11,12)i,pfl(1,n),pfl(2,n),pfl(3,n),efl(1,n),efl(2,n)
         write(17,12)i,yyre(1,0),yyre(1,1),yyre(2,0),yyre(3,0)

!         write(*,*)ftrap
      end if   
      
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ezamp(nstep,ip)   

      use gem_com
      use equil
      use fft_wrapper
      implicit none
      REAL(8) :: b,b2,gam0,gam1,th,delyz,bf,rkper,r,qr,shat
      REAL(8) :: kx,ky
      REAL(8),dimension(:),allocatable :: akx,aky
      complex(8),dimension(:,:,:,:),allocatable:: nab1,nab2
      complex(8),dimension(:,:,:,:),allocatable :: mx
      REAL(8),dimension(:,:,:),allocatable :: formapa
      REAL(8) :: sgnx,sgny,sz,myfe
      integer,dimension(:,:,:,:),allocatable :: ipiv
      INTEGER :: i,i1,j,j1,k,k1,l,m,n,ifirst,nstep,ip,INFO
      INTEGER :: l1,m1,myk,myj,izonal,ix,ikx,iext
      COMPLEX(8) :: temp3dxy(0:imx-1,0:jmx-1,0:1),v(0:imx-1,0:jcnt-1,0:1)
      COMPLEX(8) :: sbuf(0:imx*jcnt*2-1),rbuf(0:imx*jmx*2-1)
      COMPLEX(8) :: sl(1:imx-1,0:jcnt-1,0:1)
      REAL(8) :: dum,u(0:imx,0:jmx,0:1)
      REAL(8) :: grp,gthp,gxdgyp,dydrp,qhatp,grdgtp,bfldp
      REAL(8) :: wx0,wx1,wz0,wz1
      character*70 fname
      character*3 holdmyid

      save formapa,ifirst,akx,aky,mx,IPIV

      write(holdmyid,'(I3.3)') MyId
      fname='./matrix/'//'mx_apa_'//holdmyid
      if (ifirst.ne.-99) then
         allocate(akx(0:imx-1),aky(0:jcnt-1),nab1(0:imx-1,0:jcnt-1,0:imx-1,0:1))
         allocate(nab2(0:imx-1,0:jcnt-1,0:imx-1,0:1),ipiv(imx-1,imx-1,0:jcnt-1,0:1))
         allocate(mx(imx-1,imx-1,0:jcnt-1,0:1),formapa(0:imx-1,0:jcnt-1,0:1))

         if(iget.eq.1) then
            open(20000+MyId,file=fname,form='unformatted',status='old')
            read(20000+MyId)mx,ipiv
            goto 200
         end if

         do 51 l=0,im-1
            do 52 m=0,jcnt-1
               j = tclr*jcnt+m
               do 54 i=0,im-1
                  r = lr0-lx/2+i*dx
                  do 53 n = 0,1
                     k = gclr*kcnt+n
                     k1 = int(dz*k/delz)
                     k1 = min(k1,ntheta-1)
                     wz0 = ((k1+1)*delz-dz*k)/delz
                     wz1 = 1-wz0
                     th = wz0*thfnz(k1)+wz1*thfnz(k1+1)
                     i1 = int((r-rin)/dr)
                     i1 = min(i1,nr-1)
                     wx0 = (rin+(i1+1)*dr-r)/dr
                     wx1 = 1.-wx0
                     k1 = int((th+pi)/dth)
                     k1 = min(k1,ntheta-1)
                     wz0 = (-pi+(k1+1)*dth-th)/dth
                     wz1 = 1.-wz0
                     bfldp = wx0*wz0*bfld(i1,k1)+wx0*wz1*bfld(i1,k1+1) &
                     +wx1*wz0*bfld(i1+1,k1)+wx1*wz1*bfld(i1+1,k1+1) 
                     dydrp = wx0*wz0*dydr(i1,k1)+wx0*wz1*dydr(i1,k1+1) &
                     +wx1*wz0*dydr(i1+1,k1)+wx1*wz1*dydr(i1+1,k1+1) 
                     qhatp = wx0*wz0*qhat(i1,k1)+wx0*wz1*qhat(i1,k1+1) &
                     +wx1*wz0*qhat(i1+1,k1)+wx1*wz1*qhat(i1+1,k1+1) 
                     grp = wx0*wz0*gr(i1,k1)+wx0*wz1*gr(i1,k1+1) &
                     +wx1*wz0*gr(i1+1,k1)+wx1*wz1*gr(i1+1,k1+1) 
                     gthp = wx0*wz0*gth(i1,k1)+wx0*wz1*gth(i1,k1+1) &
                     +wx1*wz0*gth(i1+1,k1)+wx1*wz1*gth(i1+1,k1+1) 
                     gxdgyp = wx0*wz0*gxdgy(i1,k1)+wx0*wz1*gxdgy(i1,k1+1) &
                     +wx1*wz0*gxdgy(i1+1,k1)+wx1*wz1*gxdgy(i1+1,k1+1) 
                     grdgtp = wx0*wz0*grdgt(i1,k1)+wx0*wz1*grdgt(i1,k1+1) &
                     +wx1*wz0*grdgt(i1+1,k1)+wx1*wz1*grdgt(i1+1,k1+1) 
                     
                     if(j.ge.(jm/2+1)) then
                        m1=jm-j
                        sgny=-1.
                     else
                        m1=j
                        sgny=1.
                     endif
                     ky=sgny*2.*pi*float(m1)/ly
                     kx=pi*float(l)/lx
                     akx(l) = kx
                     aky(m) = ky
                     bf=bfldp
                     b=mims(1)*(kx*kx*grp**2 + &
                        ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                        +2*dydrp*lr0/q0*qhatp*grdgtp) &
                        +2*kx*ky*gxdgyp)/bf/bf
                  
                     nab1(l,m,i,n) = kx**2*grp**2+ky**2*  &
                        (dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                        +2*dydrp*lr0/q0*qhatp*grdgtp)

                     nab2(l,m,i,n) = -IU*ky*kx*2*gxdgyp
!   formfactor in ezamp
                     formapa(l,m,n) = beta/jmx 
                     if(abs(ky)>kycut)formapa(l,m,n) = 0.
                     if(m1.ne.mlk.and.onemd==1)formapa(l,m,n) = 0.
 53               continue
 54            continue
 52         continue
 51      continue

         do k = 0,1
            do j = 0,jcnt-1
               do i = 1,imx-1
                  do ix = 1,imx-1
                     mx(i,ix,j,k) = 0.
                     if(i==ix)mx(i,ix,j,k) = beta*amie*gn0e(i)*cn0e
                     do ikx = 1,imx-1
                        mx(i,ix,j,k) = mx(i,ix,j,k) &
                             +nab1(ikx,j,i,k) &
                             *sin(ix*ikx*pi/imx)*sin(i*ikx*pi/imx)*2.0/imx &
                             +nab2(ikx,j,i,k) &
                              *sin(ix*ikx*pi/imx)*cos(i*ikx*pi/imx)*2.0/imx
                     end do
                  end do
               end do
            end do
         end do
         if(gclr==1.and.tclr==0.and.nstep==0)then
            open(20,file="mxamp",status='unknown')
            j = 0
            k = 0
            do i = 1,imx-1
               do ix = 1,imx-1
                  write(20,10)i,ix,mx(i,ix,j,k),mx(ix,i,j,k)
               end do
            end do
 10      format(1x,i5,1x,i5,2(2x,e12.5,2x,e12.5))
            close(20)
         end if
         do k = 0,1
            do j = 0,jcnt-1
               call ZGETRF( imx-1,imx-1,mx(:,:,j,k),imx-1,IPIV(:,:,j,k), INFO )
            end do
         end do

         if(iget.eq.0) then
            open(20000+MyId,file=fname,form='unformatted',status='unknown')
            write(20000+MyId)mx,ipiv
            goto 200
         end if

 200     ifirst=-99
      endif

!   now do field solve...

!  find jtot(kx,ky)
      iext = 1
      if(nstep>200)iext = 0
      do i = 0,imx
         do j = 0,jmx
            do k = 0,1
               u(i,j,k) = jion(i,j,k)*ision+hpar(i,j,k)*isham-upar(i,j,k) &
                 +amie*(apar(i,j,k)/ntube*gn0e(i)*cn0e-upa00(i,j,k)) &
                 +0.0*sin(mlk*yg(j)-0.01*tcurr)*exp(-(xg(i)-lx/2)**2/(lx/3)**2) &
                  *exp(-(zg(k)-lz/2)**2/(lz/2)**2)*iext
            end do
         end do
      end do
      call dcmpy(u(0:imx-1,0:jmx-1,0:1),v)

      do k = 0,1
         do j = 0,jcnt-1
            myj=tclr*jcnt+j
            sl(1:imx-1,j,k) = v(1:imx-1,j,k)
            call ZGETRS('N',imx-1,1,mx(:,:,j,k),imx-1,ipiv(:,:,j,k), &
                         sl(:,j,k),imx-1,INFO) 
            temp3dxy(1:imx-1,myj,k) = sl(1:imx-1,j,k)
            temp3dxy(0,myj,k) = 0.
         end do
      end do

      do k = 0,1
         do j = 0,jcnt-1
            myj = tclr*jcnt+j
            do i = 0,imx-1
               m = j*2*imx+k*imx+i 
               sbuf(m) = temp3dxy(i,myj,k)
            end do
         end do
      end do

      call mpi_allgather(sbuf(0),jcnt*imx*2,mpi_complex16,rbuf,jcnt*imx*2, &
                         mpi_complex16,grid_comm,ierr)
                         
      do k = 0,1
         do j = 0,jmx-1
            do i = 0,imx-1
               m = j*2*imx+k*imx+i 
               temp3dxy(i,j,k) = rbuf(m)
            end do
         end do
      end do
      

!  from rho(kx,ky) to phi(kx,ky)
      do k = 0,1
         do j = 0,jmx-1
            do i = 0,imx-1
               temp3dxy(i,j,k)=temp3dxy(i,j,k)*beta/jmx
            end do
         end do
      end do

!  from phi(kx,ky) to phi(x,y)
      do k=0,mykm
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3dxy(i,j,k)
            end do
            call ccfft('y',1,jmx,1.,tmpy,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3dxy(i,j,k) = tmpy(j)  ! phi(x,y)
            end do
         end do
      end do

      do i = 0,nxpp-1
         do j = 0,jm-1
            do k = 0,mykm
                  apar(i,j,k) = temp3dxy(i,j,k)
            end do
         end do
      end do

!    x-y boundary points 
      call enfxy(apar(:,:,:))
      call enfz(apar(:,:,:))
!      call filter(apar(:,:,:))

      return
      end

!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine filter(u)

      use gem_com
      use equil
      implicit none
      integer :: i,j,k,l,m,n,ip,NMDZ,NNI,ikz,ndum
      parameter(NMDZ=2)
      real(8) :: u(0:imx,0:jmx,0:1)
      real(8) :: temp(0:imx,0:jmx,0:1)
      real(8) :: cc0(0:imx,0:jmx,0:NMDZ-1),ss0(0:imx,0:jmx,0:NMDZ-1)
      real(8) :: cc1(0:imx,0:jmx,0:NMDZ-1),ss1(0:imx,0:jmx,0:NMDZ-1)
      real(8) ::  rkz,DW1,DW2
      parameter(NNI = 3)
      parameter(DW1 = 0.5)
      parameter(DW2 = -1./6.)
      REAL(8) :: ttemp3d(0:1,0:imx-1,0:jmx-1)
      REAL(8) :: htemp3d(0:1,0:imx-1,0:jmx-1)
      REAL(8) :: lbfs(0:imx-1,0:jmx-1) 
      REAL(8) :: lbfr(0:imx-1,0:jmx-1)
      REAL(8) :: rbfs(0:imx-1,0:jmx-1)
      REAL(8) :: rbfr(0:imx-1,0:jmx-1)
      real(8) :: tmpbc(0:imx,0:jmx)
      real(8) :: uz0(0:imx,0:jmx),uz1(0:imx,0:jmx)
      integer :: nflag

      nflag = 0

!      if(onemd.eq.1)goto 200

      do i = 0,im
         do j = 0,jm
            do k = 0,mykm
               temp(i,j,k) = u(i,j,k)
            end do
         end do
      end do

      do k=0,1
         do i = 0,im-1
            do j = 0,jm-1
               ttemp3d(k,i,j)=temp(i,j,k)
            end do
         end do
      enddo

!      include 'filtgt1p.h'
      include 'filt1p.h'

      do k=0,mykm
         do j=0,jm-1
            do i=0,im-1
	       temp(i,j,k)=ttemp3d(k,i,j)
            enddo
         enddo
      enddo

      call enfxy(temp)

 100  continue
      do i = 0,im
         do j = 0,jm
            do k = 0,mykm
               u(i,j,k) = temp(i,j,k)
            end do
         end do
      end do
      return

 200  continue

      if(GCLR.eq.0)uz0(:,:)=u(:,:,0)
      if(GCLR.eq.GLST)uz1(:,:)=u(:,:,mykm)

      call MPI_BCAST(uz0,(imx+1)*(jmx+1),MPI_REAL8,0, &
          TUBE_COMM,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      
      call MPI_BCAST(uz1,(imx+1)*(jmx+1),MPI_REAL8,GLST, &
          tube_comm,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do i = 0,imx
         do j = 0,jmx
            do k = 0,mykm
               temp(i,j,k) = u(i,j,k) !-(uz1(i,j)-uz0(i,j))/lz*zg(k)
!                   -uz0(i,j)
            end do
         end do
      end do

!   smooth for shearless slab, work on for peridical in z
      cc0 = 0.
      ss0 = 0.
      cc1 = 0.
      ss1 = 0.
      kzlook = 0
      do ikz = 0,NMDZ-1
         rkz = 2.*pi/lz*float(kzlook+ikz)
         do i = 0,im
            do j = 0,jm
               do k = 0,mykm-1
                  cc0(i,j,ikz) = cc0(i,j,ikz)+temp(i,j,k)*cos(zg(k)*rkz)
                  ss0(i,j,ikz) = ss0(i,j,ikz)+temp(i,j,k)*sin(zg(k)*rkz)
               end do
            end do
         end do
      end do

      ndum = (imx+1)*(jmx+1)*(NMDZ)
      call MPI_ALLREDUCE(cc0(0,0,0),cc1(0,0,0),ndum,MPI_REAL8, &
          MPI_SUM,tube_comm,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(ss0(0,0,0),ss1(0,0,0),ndum,MPI_REAL8, &
          MPI_SUM,tube_comm,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do i = 0,im
         do j = 0,jm
            cc1(i,j,0) = 1./float(km)*cc1(i,j,0)*0.
            ss1(i,j,0) = 1./float(km)*ss1(i,j,0)
         end do
      end do

      if(NMDZ.gt.1)then
         do i = 0,im
            do j = 0,jm
               do ikz = 1,NMDZ-1
                  cc1(i,j,ikz) = 2./float(km)*cc1(i,j,ikz)
                  ss1(i,j,ikz) = 2./float(km)*ss1(i,j,ikz)
               end do
            end do
         end do
      end if

      tmpbc = 0.
      do ikz = 0,NMDZ-1
         tmpbc(:,:) = tmpbc(:,:)+cc1(:,:,ikz)
      end do    

      do i = 0,im
         do j = 0,jm
            do k = 0,mykm
               temp(i,j,k) = 0.
               do ikz = 0,NMDZ-1
                  rkz = 2.*pi/lz*float(kzlook+ikz)
                  temp(i,j,k) = temp(i,j,k)+cc1(i,j,ikz)*cos(zg(k)*rkz) &
                   +ss1(i,j,ikz)*sin(zg(k)*rkz)
               end do
            end do
         end do
      end do

      do i = 0,im
         do j = 0,jm
            do k = 0,mykm
               u(i,j,k) = temp(i,j,k) !+(uz1(i,j)-uz0(i,j))/lz*zg(k)
!                   +uz0(i,j)
            end do
         end do
      end do

      call enfxy(u)

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine yveck(u,n)

!     calculate mode histories. calculate modes for u at timestep n
!     and store in modehis(mode,n).

      use gem_com
      use equil
      use fft_wrapper
      implicit none
!     
      REAL(8) :: u(0:imx,0:jmx,0:1)
      INTEGER :: n,i,j,k
      COMPLEX(8) :: tmp3d(0:imx-1,0:jmx-1,0:1-1)

      call ccfft('z',0,kmx,1.,tmpz,tmpz,coefz,workz,0)

!      if(n.eq.0) return

      do j=0,jm-1
         do i=0,im-1
            do k = 0,mykm-1
               tmp3d(i,j,k)=u(i,j,k)
            enddo
         end do
      end do

      do k = 0,mykm-1
         do j = 0,jmx-1
            do i = 0,imx-1
               tmpx(i) = tmp3d(i,j,k)
            end do
            call ccfft('x',1,imx,1.,tmpx,tmpx,coefx,workx,0)
            do i = 0,imx-1
               tmp3d(i,j,k) = tmpx(i)
            end do
         end do
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = tmp3d(i,j,k)
            end do
            call ccfft('y',1,jmx,1.,tmpy,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               tmp3d(i,j,k) = tmpy(j)
            end do
         end do
      end do   


      do j = 1,4
         do k = 0,mykm-1
            tmpz(k) = tmp3d(llk,j,k)
         end do
         
         if(GCLR.ne.master)then
            call MPI_SEND(tmpz(0),mykm,MPI_COMPLEX16,master, &
               gclr,tube_comm,ierr)
         end if

         if(gclr.eq.master)then
            do i = 1,GLST
               call MPI_RECV(tmpz(i*mykm),mykm,MPI_COMPLEX16,i, &
                   i,tube_comm,stat,ierr)
            end do
         end if

         if(GCLR.eq.master) then
            call ccfft('z',1,kmx,1.,tmpz,tmpz,coefz,workz,0)
         end if

         call MPI_BCAST(tmpz,kmx,MPI_COMPLEX16,master, &
            tube_comm,ierr)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         do i = 0,4
            yyamp(j,i) = abs(tmpz(i)) !cabs
            yyre(j,i) = real(tmpz(i))
            yyim(j,i) = aimag(tmpz(i))
         end do
      end do

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(8) function ran2(idum)
      parameter( IM1=2147483563,  &
                IM2=2147483399, &
                AM=1.0/IM1,&
                IMM1=IM1-1,&
                IA1=40014,&
                IA2=40692,&
                IQ1=53668,&
                IQ2=52774,&
                IR1=12211,&
                IR2=3791,&
                NTAB=32,&
                NDIV=1+IMM1/NTAB,&
                EPS=1.2e-7,&
                RNMX=1.0-EPS &
               )
      integer :: j,k,idum2=123456789,iy=0,iv(0:NTAB-1)
      real(8) :: temp

      save idum2, iy,iv
!      write(*,*)'idum2,iy  ',idum2,iy
      if(idum.le.0)then
         if(-idum.lt.1)then
            idum=1
         else
            idum = -idum
         end if
         idum2 = idum
         do j = NTAB+7,0,-1
            k = idum/IQ1
            idum = IA1*(idum-k*IQ1)-k*IR1
            if(idum.lt.0)idum = idum+IM1
            if(j.lt.NTAB)iv(j) = idum
         end do
         iy = iv(0)
      end if

      k = idum/IQ1
      idum = IA1*(idum-k*IQ1)-k*IR1
      if(idum.lt.0)idum = idum+IM1
      k = idum2/IQ2
      idum2 = IA2*(idum2-k*IQ2)-k*IR2
      if(idum2.lt.0)idum2 = idum2+IM2
      j = iy/NDIV
      iy = iv(j)-idum2
      iv(j) = idum
      if(iy<1)iy = iy+IMM1
      temp = AM*iy
      if(temp>RNMX)then
         ran2 = RNMX
      else
         ran2 = temp
      end if
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine loadi(ns)

      use gem_com
      use equil
      implicit none
      INTEGER :: i,j,k,m,idum,ns,m1
      INTEGER :: np_old,np_new
      REAL(8) :: vpar,vperp2,r,qr,th,b,cost,ter
      REAL(8) :: avgv,myavgv,avgw,myavgw
      real(8) :: dumx,dumy,dumz,jacp
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,psp
      REAL(8) :: grp,gxdgyp,zoldp
      REAL(8) :: wx0,wx1,wz0,wz1

      cnt=int(tmm(1)/numprocs)
!   write(*,*)'in loader cnt, mm(1)=',cnt,mm(1)

      myavgv=0.
      avgv=0.
      avgw = 0.
      myavgw = 0.

!      ns = 1
      m = 0
      do 160 j = 1,100000000

!     load a slab of ions...

         dumx=lx*ran2(iseed) !revers(MyId*cnt+j,2) !ran2(iseed)
         dumy=ly*ran2(iseed) !revers(MyId*cnt+j,3) !ran2(iseed)
         dumz=lz*ran2(iseed) !revers(MyId*cnt+j,5) !ran2(iseed)
         dumz = min(dumz,lz-1.e-8)
         r = lr0+dumx-0.5*lx
         th = (dumz-lz/2)/(q0*br0)
         i = int((r-rin)/dr)
         k = int((pi+th)/dth)
         jacp = jacob(i,k)
         if(ran2(iseed)<(0.5*jacp/jacmax))then
         m = m+1
         if(m>mm(ns))goto 170
         x2(ns,m)=min(dumx,lx-1.d-8)
         y2(ns,m)=min(dumy,ly-1.d-8)

         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1-wz0
         z2(ns,m) = wz0*zfnth(k)+wz1*zfnth(k+1)
         z2(ns,m)=min(z2(ns,m),lz-1.d-8)

         call parperp(vpar,vperp2,m,pi,cnt,MyId)
!   normalizations will be done in following loop...

         r=x2(ns,m)-0.5*lx+lr0
         cost=cos(th)
         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
         psp = wx0*psi(i)+wx1*psi(i+1)
         ter = wx0*t0s(ns,i)+wx1*t0s(ns,i+1)
         b=1.-tor+tor*bfldp

         u2(ns,m)=vpar/sqrt(mims(ns)/ter)
         mu(ns,m)=0.5*vperp2/b*ter
         eki(ns,m) = mu(ns,m)*b+0.5*mims(ns)*u2(ns,m)**2
         pzi(ns,m) = mims(ns)*u2(ns,m)/b-psp/br0
         z0i(ns,m) = z2(ns,m)
         xii(ns,m) = x2(ns,m)
         myavgv=myavgv+u2(ns,m)

!    LINEAR: perturb w(ns,m) to get linear growth...
         w2(ns,m)=2.*amp*(revers(MyId*cnt+m,13)-0.5) !(ran2(iseed) - 0.5 )
         if(izon.eq.1)w2(ns,m)=2.*amp*sin(x2(ns,m)/lx*2*pi)  
         myavgw=myavgw+w2(ns,m)
         end if
 160  continue
 170  continue
      myavgw = myavgw/mm(ns)
!      write(*,*)'myid ', myid,x2(10),y2(20),u2(20),w2(20)
!    subtract off avg. u...
      call MPI_ALLREDUCE(myavgv,avgv,1, &
          MPI_REAL8, &
          MPI_SUM,MPI_COMM_WORLD,ierr)
      if(idg.eq.1)write(*,*)'all reduce'
      avgv=avgv/float(tmm(1))
      do 180 m=1,mm(ns)
         u2(ns,m)=u2(ns,m)-avgv
         x3(ns,m)=x2(ns,m)
         y3(ns,m)=y2(ns,m)
         z3(ns,m)=z2(ns,m)
         u3(ns,m)=u2(ns,m)
!         w2(ns,m) = w2(ns,m)-myavgw
         w3(ns,m)=w2(ns,m)
 180  continue

      np_old=mm(ns)
      call init_pmove(z2(ns,:),np_old,lz,ierr)
      if(idg.eq.1)write(*,*)'pass init_pmove'
!     
      call pmove(x2(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x3(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y2(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y3(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z2(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z3(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u2(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u3(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w2(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w3(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mu(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call pmove(xii(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0i(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(pzi(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(eki(ns,:),np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
!     
      call end_pmove(ierr)
      mm(ns)=np_new

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine enforce(u)
      use gem_com
      use equil
      implicit none
      INTEGER :: m,i,j,k,ns,jj
      REAL(8) :: u(0:imx,0:jmx,0:1)      
      REAL(8) :: lbfs(0:imx,0:jmx)
      REAL(8) :: rbfs(0:imx,0:jmx)
      REAL(8) :: lbfr(0:imx,0:jmx)
      REAL(8) :: rbfr(0:imx,0:jmx)
      real(8) :: dum,dum1,dely,th,wy1,ydum

      do 300 j=0,jm-1
         do 310 k=0,mykm
            u(0,j,k) = u(0,j,k)+u(im,j,k)
 310     continue
 300  continue

      do 320 i=0,im-1 
         do 330 k=0,mykm
            u(i,0,k) = u(i,0,k)+u(i,jm,k)
            u(i,jm,k) = u(i,0,k)
 330     continue
 320  continue 

      rbfs=u(:,:,mykm)
      call MPI_SENDRECV(rbfs(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8, &
          rngbr,101, &
          lbfr(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8, &
          lngbr,101, &
          tube_comm,stat,ierr)	
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      lbfs=u(:,:,0)
      call MPI_SENDRECV(lbfs(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8, &
          lngbr,102, &
          rbfr(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8, &
          rngbr,102, &
          tube_comm,stat,ierr)	
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      do i=0,im
         do j=0,jm
            u(i,j,0)=u(i,j,0)  &
                +weightp(i)*lbfr(i,jpl(i,j))  &
                +weightpn(i)*lbfr(i,jpn(i,j)) 
         enddo
      enddo
      do i=0,im
         do j=0,jm
            u(i,j,mykm)=u(i,j,mykm)           &
                +weightm(i)*rbfr(i,jmi(i,j)) &
                +weightmn(i)*rbfr(i,jmn(i,j))
         enddo
      enddo

      call enfxy(u)
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine enfxy(u)
      use gem_com
      use equil
      implicit none
      real(8) :: u(0:imx,0:jmx,0:1)
      integer :: i,j,k,l,m,n,jj
      real(8) :: ydum,th,dely,wy1

!    periodic bc in y...
      do 611 k=0,mykm
         do 613 i=0,im-1
            u(i,jm,k)=u(i,0,k)
 613     continue
 611  continue

!   bc for x
      do 711 k=0,mykm
         do 712 j=0,jm
            u(im,j,k)=u(0,j,k)
 712     continue
 711  continue

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gradu(u,ux,uy)
      use gem_com
      use equil
      implicit none
      real(8) :: u(0:imx,0:jmx,0:1)
      real(8) :: ux(0:imx,0:jmx,0:1),uy(0:imx,0:jmx,0:1)
      integer :: i,j,k,l,m,n,jj,ju,jl
      real(8) :: ydum,th,dely,wy1,ul

      do j=0,jm-1
         ju = j+1
         jl = j-1
         if(j.eq.0)jl = jm-1
         do i=0,im-1
            do k=0,mykm
               uy(i,j,k)=(u(i,ju,k)-u(i,jl,k))/(2.*dy)
            enddo
         enddo
      enddo

      do i=1,im-1
         do j=0,jm-1
            do k=0,mykm
               ux(i,j,k)=(u(i+1,j,k)-u(i-1,j,k))/(2.*dx)
            enddo
         enddo
      enddo

! do boundary i=0
      do j=0,jm-1
         do k=0,mykm
            ul=u(im-1,j,k)
            ux(0,j,k)=(u(1,j,k)-ul)/(2.*dx)
         enddo
      enddo

      call enfxy(ux)
      call enfxy(uy)
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gradx(u,ux)
      use gem_com
      use equil
      implicit none
      real(8) :: u(0:imx,0:jmx,0:1)
      real(8) :: ux(0:imx,0:jmx,0:1)
      integer :: i,j,k,l,m,n,jj,ju,jl
      real(8) :: ydum,th,dely,wy1,ul

      do i=1,im-1
         do j=0,jm-1
            do k=0,mykm
               ux(i,j,k)=(u(i+1,j,k)-u(i-1,j,k))/(2.*dx)
            enddo
         enddo
      enddo

! do boundary i=0
      do j=0,jm-1
         do k=0,mykm
            ul=u(im-1,j,k)
            ux(0,j,k)=(u(1,j,k)-ul)/(2.*dx)
         enddo
      enddo

      call enfxy(ux)
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grady(u,uy)
      use gem_com
      use equil
      implicit none
      real(8) :: u(0:imx,0:jmx,0:1)
      real(8) :: uy(0:imx,0:jmx,0:1)
      integer :: i,j,k,l,m,n,jj,ju,jl
      real(8) :: ydum,th,dely,wy1,ul

      do j=0,jm-1
         ju = j+1
         jl = j-1
         if(j.eq.0)jl = jm-1
         do i=0,im-1
            do k=0,mykm
               uy(i,j,k)=(u(i,ju,k)-u(i,jl,k))/(2.*dy)
            enddo
         enddo
      enddo

      call enfxy(uy)
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine enfz(u)
      use gem_com
      use equil
      implicit none
      real(8) :: u(0:imx,0:jmx,0:1)
      real(8) :: lbfr(0:imx,0:jmx)
      real(8) :: lbfs(0:imx,0:jmx)
      real(8) :: rbfr(0:imx,0:jmx)
      real(8) :: rbfs(0:imx,0:jmx)
      integer :: i,j

      do i = 0,im
         do j = 0,jm
            rbfs(i,j)=u(i,j,mykm)
         end do   
      end do
      call MPI_SENDRECV(rbfs(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8,rngbr,204, &
          lbfr(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8,lngbr,204,       &
          tube_comm,stat,ierr)
      do i = 0,im
         do j = 0,jm
            lbfs(i,j)=u(i,j,0)
         end do
      end do   
      call MPI_SENDRECV(lbfs(0,0),(imx+1)*(jmx+1),  &
          MPI_REAL8,lngbr,205,  &
          rbfr(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8,rngbr,205,  &
          tube_comm,stat,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      do i=0,im
         do j=0,jm
            u(i,j,0)=(weightp(i)*lbfr(i,jpl(i,j))  &
                +weightpn(i)*lbfr(i,jpn(i,j))  &
                +u(i,j,0) )/2.  
            u(i,j,mykm)=(weightm(i)*rbfr(i,jmi(i,j)) &  
                +weightmn(i)*rbfr(i,jmn(i,j))  &
                +u(i,j,mykm) )/2.
         enddo
      enddo

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pint
      use gem_com
      use equil
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dgdtp,dpdzp,dadzp,aparp
      REAL(8) :: dum,vxdum,dum1,dum2,eps,x,h_x,h_coll
      INTEGER :: m,i,j,k,l,n,ipover,ieover
      INTEGER :: np_old,np_new
      REAL(8) :: vfac,kap,kapnp,kaptp,sz,vpar,ppar,vpdum,pidum,xnp
      REAL(8) :: b,th,r,enerb,cost,sint,qr,laps,ter
      REAL(8) :: xt,xs,yt,zt,xdot,ydot,zdot,pzdot,edot,pzd0,pzd1,vp0
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
      REAL(8) :: grp,gxdgyp,psp,pzp
      REAL(8) :: pzcrit,encrit,myavptch,myaven
      integer :: myopz,myoen
      REAL(8) :: x000,x001,x010,x011,x100,x101,x110,x111

      myopz = 0
      myoen = 0
      myavptch = 0.
      myaven = 0.
      pidum = 1./(pi*2)**1.5*(vwidthe)**3
      do m=1,mme
         r=x2e(m)-0.5*lx+lr0

         k = int(z2e(m)/delz)
         wz0 = ((k+1)*delz-z2e(m))/delz
         wz1 = 1-wz0
         th = wz0*thfnz(k)+wz1*thfnz(k+1)

         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0
         dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
                 +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
         dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
                 +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
         grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
                 +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
         radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
                 +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
         dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
                 +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
         qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
                 +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
         grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
                 +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
         gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
                 +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 
         fp = wx0*f(i)+wx1*f(i+1)       
         jfnp = wz0*jfn(k)+wz1*jfn(k+1) 
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         psp = wx0*psi(i)+wx1*psi(i+1)        
         ter = wx0*t0e(i)+wx1*t0e(i+1)        
         kaptp = wx0*capte(i)+wx1*capte(i+1)        
         kapnp = wx0*capne(i)+wx1*capne(i+1)        
         xnp = wx0*xn0e(i)+wx1*xn0e(i+1)        
         b=1.-tor+tor*bfldp
         pzp = emass*u2e(m)/b+psp/br0

         xt = x2e(m)
         yt = y2e(m)

         include 'ppushlie.h'

         vfac = 0.5*(emass*u2e(m)**2 + 2.*mue2(m)*b)
         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))
         kap = kapnp - (1.5-vfac/ter)*kaptp

         ppar = u2e(m)
         vpar = u2e(m)-qel/emass*aparp*ipara
         enerb=(mue2(m)+emass*vpar*vpar/b)/qel*tor
         dum2 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp+vpar/b*delbxp)*dum2
         xdot = vxdum*nonline  &
             -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1+vpar/b*delbyp)*dum2*nonline &
             +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0

         zdot =  vpar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         pzd0 = tor*(-mue2(m)/emass/radiusp/bfldp*psipp*dbdtp*grcgtp)
         pzd1 = qel/emass*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
             +qel/emass*(-xdot*delbyp+ydot*delbxp+zdot*dadzp) 
         pzdot = pzd0+pzd1*ipara
                      
         vpdum = vpar
         edot = qel*(xdot*exp1+(ydot-vp0)*eyp+zdot*ezp) &
                +qel*pzdot*aparp &
                +qel*vpdum*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)   &
                -qel*vpar*delbxp*vp0

         x3e(m) = x2e(m) + 0.5*dte*xdot
         y3e(m) = y2e(m) + 0.5*dte*ydot
         z3e(m) = z2e(m) + 0.5*dte*zdot
         u3e(m) = u2e(m) + 0.5*dte*pzdot
         mue3(m) = mue2(m)

         eps = b*mue2(m)+0.5*emass*u2e(m)*u2e(m)
         x = sqrt(eps)
         h_x    = 4.0*x*x/(3.0*sqrt(pi))
         h_coll = h_x/sqrt(1.0+h_x**2)
!         h_coll = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*DERF(x)
         dum1 = 0.5*rneu/eps**1.5*(1+h_coll)
         if(x<0.136)dum1=0.0

         dum = 1-w2e(m)*nonline*0.
         if(isuni.eq.1)dum = exp(-vfac)*pidum
         vxdum = (eyp+vpdum/b*delbxp)*dum2
         w3e(m)=w2e(m) + 0.5*dte*(  &
             (vxdum*kap + edot/ter-2*dum1*ppar*aparp)*xnp     &
             +isg*(-dgdtp-zdot*dpdzp+xdot*exp1+ydot*eyp))*dum

!         if(x3e(m)>lx .or. x3e(m)<0.)w3e(m) = 0. 

!         go to 333

         pzcrit = 3.
         encrit = 5.
         ieover = 0
         ipover = 0
         if(abs(pzp-pze(m))>pzcrit)then
            myopz = myopz+1
            ipover = 1
         end if
         if(abs(vfac-eke(m))>encrit.and.abs(x2e(m)-lx/2)<(lx/2-1))then
            myoen = myoen+1
            ieover = 1
            myaven = myaven+eke(m)
            myavptch = myavptch+abs(vpar)/sqrt(2/emass*vfac)
         end if
!         goto 333
         if(ieover==1.or.ipover==1)then
            x3e(m) = xie(m)
            z3e(m) = z0e(m)
            r = x3e(m)-lx/2+lr0
            k = int(z3e(m)/delz)
            wz0 = ((k+1)*delz-z3e(m))/delz
            wz1 = 1-wz0
            th = wz0*thfnz(k)+wz1*thfnz(k+1)

            i = int((r-rin)/dr)
            wx0 = (rin+(i+1)*dr-r)/dr
            wx1 = 1.-wx0
            k = int((th+pi)/dth)
            wz0 = (-pi+(k+1)*dth-th)/dth
            wz1 = 1.-wz0
            b = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
            u3e(m) = sqrt(2/emass*abs(eke(m)-mue(m)*b))
            u2e(m) = u3e(m)
            w3e(m) = 0.
            w2e(m) = 0.
            x2e(m) = x3e(m)
            z2e(m) = z3e(m)
            mue3(m) = mue(m)
            mue2(m) = mue(m)
         end if
 333     continue
         laps=anint((z3e(m)/lz)-.5)*(1-peritr)
         r=x3e(m)-0.5*lx+lr0
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         qr = wx0*sf(i)+wx1*sf(i+1)
         y3e(m)=dmod(y3e(m)-laps*2*pi*qr*lr0/q0+800.*ly,ly)
         if(x3e(m)>lx)then
            x3e(m) = lx-1.e-8
            z3e(m)=lz-z3e(m)
            x2e(m) = x3e(m)
            z2e(m) = z3e(m)
            w2e(m) = 0.
            w3e(m) = 0.
         end if
         if(x3e(m)<0.)then
            x3e(m) = 1.e-8
            z3e(m)=lz-z3e(m)
            x2e(m) = x3e(m)
            z2e(m) = z3e(m)
            w2e(m) = 0.
            w3e(m) = 0.
         end if
         z3e(m)=dmod(z3e(m)+8.*lz,lz)
         x3e(m)=dmod(x3e(m)+8.*lx,lx)         
         x3e(m) = min(x3e(m),lx-1.0e-8)
         y3e(m) = min(y3e(m),ly-1.0e-8)
         z3e(m) = min(z3e(m),lz-1.0e-2)

      end do
      call MPI_ALLREDUCE(myopz,nopz,1,MPI_integer, &
          MPI_SUM, MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(myoen,noen,1,MPI_integer, &
          MPI_SUM, MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(myaven,aven,1,MPI_real8, &
          MPI_SUM, MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(myavptch,avptch,1,MPI_real8, &
          MPI_SUM, MPI_COMM_WORLD,ierr)
      aven = aven/(noen+0.1)
      avptch = avptch/(noen+0.1)

      np_old=mme
      call init_pmove(z3e,np_old,lz,ierr)
      call pmove(x2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mue2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mue3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call pmove(index,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ipass,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call pmove(mue,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xie,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(eke,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(pze,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call end_pmove(ierr)
      mme=np_new

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cint(n)
      use gem_com
      use equil
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dgdtp,dpdzp,dadzp,aparp
      REAL(8) :: dum,vxdum,dum1,dum2,eps,x,h_x,h_coll
      INTEGER :: m,i,j,k,l,n,mynowe
      INTEGER :: np_old,np_new
      REAL(8) :: vfac,kap,kapnp,kaptp,sz,vpar,ppar,vpdum,pidum,xnp
      REAL(8) :: b,th,r,enerb,cost,sint,qr,laps,ter
      REAL(8) :: xt,xs,yt,zt,xdot,ydot,zdot,pzdot,edot,pzd0,pzd1,vp0
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,w3old
      REAL(8) :: myke,mypfl,myptrp
      REAL(8) :: myefl,mynos,myavewe
      REAL(8) :: ketemp,pfltemp,ptrptmp
      REAL(8) :: efltemp,nostemp
      REAL(8) :: sbuf(10),rbuf(10)
      REAL(8) :: mytotn,mytrap,totn,ttrap
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
      REAL(8) :: grp,gxdgyp,psp,pzp
      REAL(8) :: x000,x001,x010,x011,x100,x101,x110,x111

      sbuf(1:10) = 0.
      rbuf(1:10) = 0.
      myavewe = 0.
      myke=0.    
      mypfl=0.    
      myptrp=0.
      myefl=0. 
      mynos=0.   
      ketemp=0.
      pfltemp=0.
      efltemp=0.
      nostemp=0.
      mytotn = 0.
      mytrap = 0.
      mynowe = 0

      pidum = 1./(pi*2)**1.5*(vwidthe)**3
      do m=1,mme
         r=x3e(m)-0.5*lx+lr0

         k = int(z3e(m)/delz)
         wz0 = ((k+1)*delz-z3e(m))/delz
         wz1 = 1-wz0
         th = wz0*thfnz(k)+wz1*thfnz(k+1)

         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0
         dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
                 +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
         dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
                 +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
         grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
                 +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
         radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
                 +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
         dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
                 +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
         qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
                 +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
         grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
                 +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
         gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
                 +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 
         fp = wx0*f(i)+wx1*f(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         psp = wx0*psi(i)+wx1*psi(i+1)        
         ter = wx0*t0e(i)+wx1*t0e(i+1)        
         kaptp = wx0*capte(i)+wx1*capte(i+1)        
         kapnp = wx0*capne(i)+wx1*capne(i+1)        
         xnp = wx0*xn0e(i)+wx1*xn0e(i+1)        
         b=1.-tor+tor*bfldp
         pzp = emass*u3e(m)/b+psp/br0

         xt = x3e(m)
         yt = y3e(m)

         include 'cpushlie.h'

         vfac = 0.5*(emass*u3e(m)**2 + 2.*mue3(m)*b)
         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))
         kap = kapnp - (1.5-vfac/ter)*kaptp

         ppar = u3e(m)
         vpar = u3e(m)-qel/emass*aparp*ipara
         dum2 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         enerb=(mue3(m)+emass*vpar*vpar/b)/qel*tor
         vxdum = (eyp+vpar/b*delbxp)*dum2
         xdot = vxdum*nonline &
             -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1+vpar/b*delbyp)*dum2*nonline  &
             +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0

         zdot =  vpar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         pzd0 = tor*(-mue3(m)/emass/radiusp/bfldp*psipp*dbdtp*grcgtp)
         pzd1 = qel/emass*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
             +qel/emass*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)
         pzdot = pzd0+pzd1*ipara

         vpdum = vpar
         edot = qel*(xdot*exp1+(ydot-vp0)*eyp+zdot*ezp) &
                +qel*pzdot*aparp  &
                +qel*vpdum*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)   &
                -qel*vpar*delbxp*vp0

         x3e(m) = x2e(m) + dte*xdot
         y3e(m) = y2e(m) + dte*ydot
         z3e(m) = z2e(m) + dte*zdot
         u3e(m) = u2e(m) + dte*pzdot
         mue3(m) = mue2(m)

         eps = b*mue3(m)+0.5*emass*u3e(m)*u3e(m)
         x = sqrt(eps)
         h_x    = 4.0*x*x/(3.0*sqrt(pi))
         h_coll = h_x/sqrt(1.0+h_x**2)
!         h_coll = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*DERF(x)
         dum1 = 0.5*rneu/eps**1.5*(1+h_coll)
         if(x<0.136)dum1=0.0

         dum = 1-w3e(m)*nonline*0.
         if(isuni.eq.1)dum = exp(-vfac)*pidum
         vxdum = (eyp+vpdum/b*delbxp)*dum2
!             -(2*u3e(m)*aparp+qel/emass*aparp*aparp)/(b*b)/br0*sint
         w3old = w3e(m)
         w3e(m)=w2e(m) + dte*(  &
             (vxdum*kap + edot/ter  -2*dum1*ppar*aparp)*xnp    & 
             +isg*(-dgdtp-zdot*dpdzp+xdot*exp1+ydot*eyp))*dum

         if(abs(w3e(m)).gt.0.6.and.nonline==1)then
            w3e(m) = 0.
            w2e(m) = 0.
            mynowe = mynowe+1
         end if


         mytotn = mytotn+1
         if(isuni.eq.1)mytotn = mytotn+exp(-vfac)
         mytrap = mytrap+1-ipass(m)
         if(isuni.eq.1)mytrap = mytrap+exp(-vfac)*(1-ipass(m))

         laps=anint((z3e(m)/lz)-.5)*(1-peritr)
         r=x3e(m)-0.5*lx+lr0
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         qr = wx0*sf(i)+wx1*sf(i+1)
         y3e(m)=dmod(y3e(m)-laps*2*pi*qr*lr0/q0+800.*ly,ly)
         if(x3e(m)>lx)then
            x3e(m) = lx-1.e-8
            z3e(m)=lz-z3e(m)
            x2e(m) = x3e(m)
            z2e(m) = z3e(m)
            w2e(m) = 0.
            w3e(m) = 0.
         end if
         if(x3e(m)<0.)then
            x3e(m) = 1.e-8
            z3e(m)=lz-z3e(m)
            x2e(m) = x3e(m)
            z2e(m) = z3e(m)
            w2e(m) = 0.
            w3e(m) = 0.
         end if
         z3e(m)=dmod(z3e(m)+8.*lz,lz)
         x3e(m)=dmod(x3e(m)+8.*lx,lx)         
         x3e(m) = min(x3e(m),lx-1.0e-8)
         y3e(m) = min(y3e(m),ly-1.0e-8)
         z3e(m) = min(z3e(m),lz-1.0e-2)

!     particle diagnostics done here because info is available...
         mypfl=mypfl + w3old*(eyp+vpar*delbxp/b) 
         myptrp=myptrp + w3old*eyp*(1-ipass(m)) 
         myefl=myefl + vfac*w3old*(eyp+vpar*delbxp/b)
         myke=myke + vfac*w3e(m)
         mynos=mynos + w3e(m)
         myavewe = myavewe+abs(w3e(m))

         if(abs(z3e(m)-z2e(m)).gt.lz/2)ipass(m)=1

!    xn+1 becomes xn...
         u2e(m)=u3e(m)
         x2e(m)=x3e(m)
         y2e(m)=y3e(m)
         z2e(m)=z3e(m)
         w2e(m)=w3e(m)

      enddo
      call MPI_ALLREDUCE(mynowe,nowe,1,MPI_integer, &
          MPI_SUM, MPI_COMM_WORLD,ierr)

      sbuf(1)=myke
      sbuf(2)=myefl
      sbuf(3)=mypfl
      sbuf(4)=mynos
      sbuf(5)=mytotn
      sbuf(6)=mytrap
      sbuf(7)=myptrp
      sbuf(8)=myavewe
      call MPI_ALLREDUCE(sbuf,rbuf,10,  &
          MPI_REAL8,MPI_SUM,  &
          MPI_COMM_WORLD,ierr)

      ketemp=rbuf(1)
      efltemp=rbuf(2)
      pfltemp=rbuf(3)
      nostemp=rbuf(4)
      totn=rbuf(5)
      ttrap=rbuf(6)
      ptrptmp=rbuf(7)
      avewe(n) = rbuf(8)/( float(tmm(1)) )
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      nos(2,n)=nostemp/( float(tmm(1)) )
      pfl(2,n)=pfltemp/( float(tmm(1)) )
      pfl(3,n)=ptrptmp/( float(tmm(1)) )
      efl(2,n)=efltemp/( float(tmm(1)) )
!      ke(2,n)=ketemp/( 2.*float(tmm(1))*mims(ns) )
      ftrap = ttrap/totn

      np_old=mme
      call init_pmove(z3e,np_old,lz,ierr)

      call pmove(x2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      if (ierr.ne.0) call ppexit
      call pmove(y2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mue2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mue3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call pmove(index,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ipass,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call pmove(mue,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xie,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(eke,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(pze,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call end_pmove(ierr)
      mme=np_new

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine drdt(ip)
      use gem_com
      use equil
      implicit none
      integer :: i,j,k,ip,isdndt
      real(8) :: lbfr(0:imx,0:jmx)
      real(8) :: lbfs(0:imx,0:jmx)
      real(8) :: rbfr(0:imx,0:jmx)
      real(8) :: rbfs(0:imx,0:jmx)
      real(8) :: dum
      real(8) :: djdx(0:imx,0:jmx,0:1),djdy(0:imx,0:jmx,0:1)

      call gradx(jionx*ision+hpex*ishgk-upex,djdx)
      call grady(jiony*ision+hpey*ishgk-upey,djdy)
!      djdx = 0.
!      djdy = 0.
      isdndt = 1
!      upa0(:,:,:) = apar(ip,:,:,:)
      do i = 0,im
         do j = 0,jm
            rbfs(i,j)=jion(i,j,0)*ision+hpar(i,j,0)*ishgk-(upazd(i,j,0)+amie*upa0(i,j,0))
         end do   
      end do
      call MPI_SENDRECV(rbfs,(imx+1)*(jmx+1),  &
          MPI_REAL8,rngbr,404, &
          lbfr,(imx+1)*(jmx+1), &
          MPI_REAL8,lngbr,404, &
          tube_comm,stat,ierr)
      do i = 0,im
         do j = 0,jm
            lbfs(i,j)=jion(i,j,1)*ision+hpar(i,j,1)*ishgk-(upazd(i,j,1)+amie*upa0(i,j,1))
         end do   
      end do
      call MPI_SENDRECV(lbfs,(imx+1)*(jmx+1),  &
          MPI_REAL8,lngbr,405,  &
          rbfr,(imx+1)*(jmx+1),&
          MPI_REAL8,rngbr,405, &
          tube_comm,stat,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do i = 0,im-1
         do j = 0,jm-1
            dum =  weightp(i)*lbfr(i,jpl(i,j)) &
       	+weightpn(i)*lbfr(i,jpn(i,j))
            drhodt(i,j,0) = -(jion(i,j,1)*ision+hpar(i,j,1)*ishgk-upazd(i,j,1) &
                -amie*upa0(i,j,1)-dum)/(2.*dz) &
                -djdx(i,j,0)-djdy(i,j,0) &
                +(drhoidt(i,j,0)*ision+dnhdt(i,j,0)*ishgk-dnedt(i,j,0))

            dum = weightm(i)*rbfr(i,jmi(i,j)) &
                +weightmn(i)*rbfr(i,jmn(i,j))
            drhodt(i,j,1) = -(dum-(jion(i,j,0)*ision+hpar(i,j,0)*ishgk-upazd(i,j,0)  &
                -amie*upa0(i,j,0)))/(2.*dz)  &
                -djdx(i,j,1)-djdy(i,j,1) &
                +(drhoidt(i,j,1)*ision+dnhdt(i,j,1)*ishgk-dnedt(i,j,1))
         end do
      end do
      call enfxy(drhodt(:,:,:))
      call enfz(drhodt(:,:,:))
!      call filter(drhodt(:,:,:))
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dpdt(ip)   
      use gem_com
      use equil
      use fft_wrapper
      implicit none
      REAL(8) :: b,b1,b2,gam0,gam1,delyz,th,bf,dum,r,qr,shat,ter
      REAL(8) :: kx1,kx2,ky
      REAL(8),dimension(:),allocatable :: akx,aky
      real(8),dimension(:,:,:,:),allocatable:: gamb1,gamb2
      complex(8),dimension(:,:,:,:),allocatable :: mx
      integer,dimension(:,:,:,:),allocatable :: ipiv
      REAL(8),dimension(:,:,:),allocatable :: formdpt
      REAL(8) :: sgnx,sgny,sz,myfe,u(0:imx,0:jmx,0:1)
      INTEGER :: i,i1,j,j1,k,k1,l,m,n,ifirst,nstep,ip,INFO
      INTEGER :: l1,m1,myk,myj,izonal,ix,ikx
      COMPLEX(8) :: temp3dxy(0:imx-1,0:jmx-1,0:1),v(0:imx-1,0:jcnt-1,0:1)
      COMPLEX(8) :: sbuf(0:imx*jcnt*2-1),rbuf(0:imx*jmx*2-1)
      COMPLEX(8) :: sl(1:imx-1,0:jcnt-1,0:1)
      complex(8) :: cdum
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
      REAL(8) :: grp,gxdgyp,grdgtp,gthp
      REAL(8) :: wx0,wx1,wz0,wz1
      character*70 fname
      character*3 holdmyid

      save formdpt,ifirst,akx,aky,mx,gamb1,gamb2,ipiv
      if(idg==1)write(*,*)'enter gkps'
      izonal = 1

!     form factors....
      write(holdmyid,'(I3.3)') MyId
      fname='./matrix/'//'mx_dpt_'//holdmyid
      if (ifirst.ne.-99) then
         allocate(akx(0:imx-1),aky(0:jcnt-1),&
                  gamb1(0:imx-1,0:jcnt-1,0:imx-1,0:1),&
                  gamb2(0:imx-1,0:jcnt-1,0:imx-1,0:1))
         allocate(mx(imx-1,imx-1,0:jcnt-1,0:1),formdpt(0:imx-1,0:jcnt-1,0:1))
         allocate(ipiv(imx-1,imx-1,0:jcnt-1,0:1))

         if(iget.eq.1) then
            open(30000+MyId,file=fname,form='unformatted',status='old')
            read(30000+MyId)mx,ipiv
            goto 200
         end if

         do 51 l=0,im-1
            do 52 m=0,jcnt-1
               j = tclr*jcnt+m
               do 54 i=0,im-1 
                  r = lr0-lx/2+i*dx
                  i1 = int((r-rin)/dr)
                  i1 = min(i1,nr-1)
                  wx0 = (rin+(i1+1)*dr-r)/dr
                  wx1 = 1.-wx0
                  ter = wx0*t0i(i1)+wx1*t0i(i1+1)
                  do 53 n = 0,1
                     k = gclr*kcnt+n
                     k1 = int(k*dz/delz)
                     k1 = min(k1,ntheta-1)
                     wz0 = ((k1+1)*delz-dz*k)/delz
                     wz1 = 1-wz0
                     th = wz0*thfnz(k1)+wz1*thfnz(k1+1)
                     k = int((th+pi)/dth)
                     k = min(k,ntheta-1)
                     wz0 = (-pi+(k+1)*dth-th)/dth
                     wz1 = 1.-wz0
                     bfldp = wx0*wz0*bfld(i1,k)+wx0*wz1*bfld(i1,k+1) &
                        +wx1*wz0*bfld(i1+1,k)+wx1*wz1*bfld(i1+1,k+1) 
                     dydrp = wx0*wz0*dydr(i1,k)+wx0*wz1*dydr(i1,k+1) &
                        +wx1*wz0*dydr(i1+1,k)+wx1*wz1*dydr(i1+1,k+1) 
                     qhatp = wx0*wz0*qhat(i1,k)+wx0*wz1*qhat(i1,k+1) &
                        +wx1*wz0*qhat(i1+1,k)+wx1*wz1*qhat(i1+1,k+1) 
                     grp = wx0*wz0*gr(i1,k)+wx0*wz1*gr(i1,k+1) &
                        +wx1*wz0*gr(i1+1,k)+wx1*wz1*gr(i1+1,k+1) 
                     gthp = wx0*wz0*gth(i1,k)+wx0*wz1*gth(i1,k+1) &
                        +wx1*wz0*gth(i1+1,k)+wx1*wz1*gth(i1+1,k+1) 
                     gxdgyp = wx0*wz0*gxdgy(i1,k)+wx0*wz1*gxdgy(i1,k+1) &
                        +wx1*wz0*gxdgy(i1+1,k)+wx1*wz1*gxdgy(i1+1,k+1) 
                     grdgtp = wx0*wz0*grdgt(i1,k)+wx0*wz1*grdgt(i1,k+1) &
                        +wx1*wz0*grdgt(i1+1,k)+wx1*wz1*grdgt(i1+1,k+1) 

                     if(j.ge.(jm/2+1)) then
                        m1=jm-j
                        sgny=-1.
                     else
                        m1=j
                        sgny=1.
                     endif
                     ky=sgny*2.*pi*float(m1)/ly
                     kx1=pi*float(l)/lx
                     kx2=-pi*float(l)/lx
                     bf=bfldp
                     b1=mims(1)*(kx1*kx1*grp**2 + &
                       ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                       +2*dydrp*lr0/q0*qhatp*grdgtp) &
                       +2*kx1*ky*gxdgyp)/(bf*bf)*ter

                     b2=mims(1)*(kx2*kx2*grp**2 + &
                       ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                       +2*dydrp*lr0/q0*qhatp*grdgtp) &
                       +2*kx2*ky*gxdgyp)/(bf*bf)*ter
                  
                     call srcbes(b1,gam0,gam1)
                     gamb1(l,m,i,n)=gam0
                     call srcbes(b2,gam0,gam1)
                     gamb2(l,m,i,n)=gam0

!   formfactor in gkps
                     formdpt(l,m,n) = 1./jmx 
                     if(abs(ky)>kycut)formdpt(l,m,n) = 0.
!                    if(abs(ky)==0.)formdpt(l,m,n) = 0.
                     if(m1.ne.mlk.and.onemd==1)formdpt(l,m,n) = 0.
 53               continue
 54            continue
 52         continue
 51      continue

         do k = 0,1
            do j = 0,jcnt-1
               do i = 1,imx-1
                  r = lr0-lx/2+i*dx
                  i1 = int((r-rin)/dr)
                  i1 = min(i1,nr-1)
                  wx0 = (rin+(i1+1)*dr-r)/dr
                  wx1 = 1.-wx0
                  ter = wx0*t0i(i1)+wx1*t0i(i1+1)
                  do ix = 1,imx-1
                     mx(i,ix,j,k) = 0.
                     do ikx = 0,imx-1
                        mx(i,ix,j,k) = mx(i,ix,j,k)+sin(ix*ikx*pi/imx)* & 
                             ((1-gamb1(ikx,j,i,k))*exp(IU*ikx*i*pi/imx)- &
                              (1-gamb2(ikx,j,i,k))*exp(-IU*ikx*i*pi/imx)) &
                             /ter*cn0i*gn0i(i)/(IU*imx)
                     end do
                  end do
               end do
            end do
         end do
         do k = 0,1
            do j = 0,jcnt-1
               call ZGETRF(imx-1,imx-1,mx(:,:,j,k),imx-1,ipiv(:,:,j,k),INFO )
            end do
         end do

         if(iget.eq.0) then
            open(30000+MyId,file=fname,form='unformatted',status='unknown')
            write(30000+MyId)mx,ipiv
            goto 200
         end if

         if(gclr==kmx/2.and.tclr==0.and.nstep==0)then
 !           write(*,*)'setup in gkps'
            open(20,file="mx",status='unknown')
            j = 0
            k = 0
            do i = 1,imx-1
               do ix = 1,imx-1
                  if(abs(i-ix)<40)write(20,10)i,ix,mx(i,ix,j,k),mx(ix,i,j,k)
               end do
            end do
 10      format(1x,i5,1x,i5,2(2x,e12.5,2x,e12.5))
            close(20)
         end if
 200     ifirst=-99
      endif
      if(idg==1)write(*,*)'pass form factors'

      call dcmpy(drhodt(0:imx-1,0:jmx-1,0:1),v)

      if(tclr==0)then
         do k = 0,1
            cdum = 0.
            do i = 0,imx-1
               cdum = cdum+v(i,0,k)
            end do
            cdum = cdum/imx
            do i = 0,imx-1
               v(i,0,k) = v(i,0,k)-cdum
            end do
         end do
      end if

      do k = 0,1
         do j = 0,jcnt-1
            myj=tclr*jcnt+j
            sl(1:imx-1,j,k) = v(1:imx-1,j,k)
            call ZGETRS('N',imx-1,1,mx(:,:,j,k),imx-1,ipiv(:,:,j,k), &
                         sl(:,j,k),imx-1,INFO) 
            temp3dxy(1:imx-1,myj,k) = sl(1:imx-1,j,k)
            temp3dxy(0,myj,k) = 0.
         end do
      end do

      do k = 0,1
         do j = 0,jcnt-1
            myj = tclr*jcnt+j
            do i = 0,imx-1
               m = j*2*imx+k*imx+i 
               sbuf(m) = temp3dxy(i,myj,k)
            end do
         end do
      end do

      call mpi_allgather(sbuf(0),jcnt*imx*2,mpi_complex16,rbuf,jcnt*imx*2, &
                         mpi_complex16,grid_comm,ierr)
                         
      do k = 0,1
         do j = 0,jmx-1
            do i = 0,imx-1
               m = j*2*imx+k*imx+i 
               temp3dxy(i,j,k) = rbuf(m)
            end do
         end do
      end do
      

!  from rho(kx,ky) to phi(kx,ky)
      do k = 0,1
         do j = 0,jmx-1
            do i = 0,imx-1
               temp3dxy(i,j,k)=temp3dxy(i,j,k)/jmx !*formdpt(i,j,k)
            end do
         end do
      end do

!  from phi(kx,ky) to phi(x,y)
      do k=0,mykm
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3dxy(i,j,k)
            end do
            call ccfft('y',1,jmx,1.,tmpy,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3dxy(i,j,k) = tmpy(j)  ! phi(x,y)
            end do
         end do
      end do

 100  continue
      do i = 0,nxpp-1
         do j = 0,jm-1
            do k = 0,mykm
               dphidt(i,j,k) = temp3dxy(i,j,k)
            end do
         end do
      end do

!    x-y boundary points 
      call enfxy(dphidt(:,:,:))
      call enfz(dphidt(:,:,:))
      call fltx(dphidt(:,:,:),0)
      call filter(dphidt(:,:,:))
      if(idg==1)write(*,*)'pass enfz', myid

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine jie(ip,n)

!    source quantities are are calculated: n_i
!    right now only ion quantitities are calculated...

      use gem_com
      use equil
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: enerb,vxdum,dum,xdot,ydot,zdot,pidum,dum1,dum2
      INTEGER :: m,n,i,j,k,l,ns,ip,nonfi,nonfe
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,vte
      REAL(8) :: sz,wght,wght0,wght1,wght2,r,th,cost,sint,b,qr,dv,kap,ter
      REAL(8) :: kapnp,kaptp,xnp
      REAL(8) :: xt,yt,rhog,vpar,xs,dely,vfac,vp0
      REAL(8) :: lbfs(0:imx,0:jmx)
      REAL(8) :: rbfs(0:imx,0:jmx)
      REAL(8) :: lbfr(0:imx,0:jmx)
      REAL(8) :: rbfr(0:imx,0:jmx)
      real(8) :: myjpar(0:imx,0:jmx,0:1),myjpex(0:imx,0:jmx,0:1)
      real(8) :: myjpey(0:imx,0:jmx,0:1),myupar(0:imx,0:jmx,0:1)

      real(8) :: myhpar(0:imx,0:jmx,0:1),myhpex(0:imx,0:jmx,0:1)
      real(8) :: myhpey(0:imx,0:jmx,0:1),mydnhdt(0:imx,0:jmx,0:1)

      real(8) :: myupex(0:imx,0:jmx,0:1),myupey(0:imx,0:jmx,0:1)
      real(8) :: myupazd(0:imx,0:jmx,0:1)
      real(8) :: mydnidt(0:imx,0:jmx,0:1),mydnedt(0:imx,0:jmx,0:1)
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
      REAL(8) :: grp,gxdgyp,rhox(4),rhoy(4)

      nonfi = 1 
      nonfe = 1 
      jion = 0.
      jionx = 0.
      jiony = 0.
      drhoidt = 0.
      myhpar = 0.
      myhpex = 0.
      myhpey = 0.
      mydnhdt = 0.
      ns=1

      pidum = 1./(pi*2)**1.5*(vwidth)**3
      if(isuni.eq.0)pidum = 1.

      do ns = 1,nsm
      myjpar = 0.
      myjpex = 0.
      myjpey = 0.
      mydnidt = 0.
      do m=1,mm(ns)
         dv=float(lr(ns))*(dx*dy*dz)

         r=x3(ns,m)-0.5*lx+lr0

         k = int(z3(ns,m)/delz)
         wz0 = ((k+1)*delz-z3(ns,m))/delz
         wz1 = 1-wz0
         th = wz0*thfnz(k)+wz1*thfnz(k+1)

         cost=cos(th)
         sint=sin(th)
         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0
         dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
                 +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
         dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
                 +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
         grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
                 +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
         radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
                 +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
         dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
                 +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
         qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
                 +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
         grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
                 +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
         gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
                 +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 
         fp = wx0*f(i)+wx1*f(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         ter = wx0*t0i(i)+wx1*t0i(i+1)        
         kaptp = wx0*capts(ns,i)+wx1*capts(ns,i+1)        
         kapnp = wx0*capns(ns,i)+wx1*capns(ns,i+1)        
         xnp = wx0*xn0s(ns,i)+wx1*xn0s(ns,i+1)        

!         b=1.-lr0/br0*cost
         b=1.-tor+tor*bfldp

         rhog=sqrt(2.*b*mu(ns,m)*mims(ns))/(q(ns)*b)*iflr
         vfac = 0.5*(mims(ns)*u3(ns,m)**2 + 2.*mu(ns,m)*b)
         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))
         kap = kapnp - (1.5-vfac/ter)*kaptp
         vpar = u3(ns,m)
         wght=w3(ns,m)/dv
         wght0 = exp(-vfac)/dv
         if(isuni.eq.0)wght0 = 1./dv

         rhox(1) = rhog*(1-tor)+rhog*grp*tor
         rhoy(1) = rhog*gxdgyp/grp*tor
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
         rhox(4) = 0
         rhoy(4) = -rhoy(3)

         exp1=0.
         eyp=0.
         delbxp = 0.
         delbyp = 0.
         aparp = 0.
         do 100 l=1,lr(1)
            xs=x3(ns,m)+rhox(l) !rwx(1,l)*rhog
            yt=y3(ns,m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            i=int(xt/dx+0.5)
            j=int(yt/dy+0.5)
            k=int(z3(ns,m)/dz+0.5)-gclr*kcnt
            exp1=exp1 + ex(i,j,k)
            eyp=eyp + ey(i,j,k)
            delbxp = delbxp+delbx(i,j,k)
            delbyp = delbyp+delby(i,j,k)
            aparp = aparp+apar(i,j,k)
 100     continue

         exp1=exp1/4.
         eyp=eyp/4.
         delbxp=delbxp/4.
         delbyp=delbyp/4.
         aparp = aparp/4.

         vpar = u3(ns,m)-q(ns)/mims(ns)*aparp*0.
         enerb=(mu(ns,m)+mims(ns)*vpar*vpar/b)/q(ns)*tor
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp+vpar/b*delbxp)*dum1

         xdot = vxdum*nonlin  &
              -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1+vpar/b*delbyp)*dum1*nonlin  &
             +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0
         zdot =  vpar*(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

!    now do 1,2,4 point average, where lr is the no. of points...
         do 200 l=1,lr(1)
            xs=x3(ns,m)+rhox(l) !rwx(1,l)*rhog
            yt=y3(ns,m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            i=int(xt/dx+0.5)
            j=int(yt/dy+0.5)
            k=int(z3(ns,m)/dz+0.5)-gclr*kcnt

            wght1 = wght0*(vxdum*kap+xdot*exp1/ter+(ydot-vp0)*eyp/ter)*xnp
            myjpar(i,j,k) = myjpar(i,j,k)+wght*zdot
            myjpex(i,j,k) = myjpex(i,j,k)+wght*xdot
            myjpey(i,j,k) = myjpey(i,j,k)+wght*ydot
            mydnidt(i,j,k) = mydnidt(i,j,k)+wght1
 200     continue

      enddo

!   enforce periodicity
      call enforce(myjpar)
      call enforce(myjpex)
      call enforce(myjpey)
      call enforce(mydnidt)
!      call filter(myjpar)
!      call filter(mydnidt)

      do 410 i=0,im
         do 420 j=0,jm
            do 430 k=0,mykm
               jpar(ns,i,j,k) = q(ns)*myjpar(i,j,k)/n0/jac(i,k)*cn0s(ns)
               jpex(ns,i,j,k) = q(ns)*myjpex(i,j,k)/n0/jac(i,k)*cn0s(ns)
               jpey(ns,i,j,k) = q(ns)*myjpey(i,j,k)/n0/jac(i,k)*cn0s(ns)
               dnidt(ns,i,j,k) = q(ns)*mydnidt(i,j,k)/n0/jac(i,k)*pidum*cn0s(ns)
 430        continue
 420     continue
 410  continue

      do i = 0,im
         do j = 0,jm
            do k = 0,mykm
               jion(i,j,k) = jion(i,j,k)+jpar(ns,i,j,k)
               jionx(i,j,k) = jionx(i,j,k)+jpex(ns,i,j,k)
               jiony(i,j,k) = jiony(i,j,k)+jpey(ns,i,j,k)
               drhoidt(i,j,k) = drhoidt(i,j,k)+dnidt(ns,i,j,k)
            end do
         end do
      end do
      end do

! electrons current
      pidum = 1./(pi*2)**1.5*(vwidthe)**3
      if(isuni.eq.0)pidum = 1.

      vte = sqrt(amie)
      myupex = 0.
      myupey = 0.
      myupazd = 0.
      mydnedt = 0.
      do m=1,mme
         dv=(dx*dy*dz)

         r=x3e(m)-0.5*lx+lr0

         k = int(z3e(m)/delz)
         wz0 = ((k+1)*delz-z3e(m))/delz
         wz1 = 1-wz0
         th = wz0*thfnz(k)+wz1*thfnz(k+1)

         cost=cos(th)
         sint=sin(th)
         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0
         dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
                 +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
         dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
                 +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
         grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
                 +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
         radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
                 +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
         dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
                 +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
         qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
                 +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
         grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
                 +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
         gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
                 +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 
         fp = wx0*f(i)+wx1*f(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         ter = wx0*t0e(i)+wx1*t0e(i+1)        
         kaptp = wx0*capte(i)+wx1*capte(i+1)        
         kapnp = wx0*capne(i)+wx1*capne(i+1)        
         xnp = wx0*xn0e(i)+wx1*xn0e(i+1)        

!         b=1.-lr0/br0*cost
         b=1.-tor+tor*bfldp

         vfac = 0.5*(emass*u3e(m)**2 + 2.*mue3(m)*b)
         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))
         kap = kapnp - (1.5-vfac/ter)*kaptp

         wght=w3e(m)/dv
         wght0=exp(-vfac)/dv
         if(isuni.eq.0)wght0=1./dv
         wght0 = wght0
         vpar = u3e(m)
         if(abs(vpar/vte).gt.vcut)wght = 0.
         if(abs(vpar/vte).gt.vcut)wght0 = 0.

         enerb=(mue3(m)+emass*vpar*vpar/b)/qel*tor
         xdot = -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0

         zdot =  u3e(m)*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         xt = x3e(m)
         yt = y3e(m)
         i=int(xt/dx)
         j=int(yt/dy)
         k=0 !int(z3e(m)/dz)-gclr*kcnt

         eyp = w000(m)*ey(i,j,k)  &
             + w100(m)*ey(i+1,j,k) &
             + w010(m)*ey(i,j+1,k) &
             + w110(m)*ey(i+1,j+1,k) &
             + w001(m)*ey(i,j,k+1) &
             + w101(m)*ey(i+1,j,k+1) &
             + w011(m)*ey(i,j+1,k+1) &
             + w111(m)*ey(i+1,j+1,k+1)

         exp1 = w000(m)*ex(i,j,k)  &
             + w100(m)*ex(i+1,j,k) &
             + w010(m)*ex(i,j+1,k) &
             + w110(m)*ex(i+1,j+1,k) &
             + w001(m)*ex(i,j,k+1) &
             + w101(m)*ex(i+1,j,k+1) &
             + w011(m)*ex(i,j+1,k+1) &
             + w111(m)*ex(i+1,j+1,k+1)

         ezp = w000(m)*ez(i,j,k)  &
             + w100(m)*ez(i+1,j,k) &
             + w010(m)*ez(i,j+1,k) &
             + w110(m)*ez(i+1,j+1,k) &
             + w001(m)*ez(i,j,k+1) &
             + w101(m)*ez(i+1,j,k+1) &
             + w011(m)*ez(i,j+1,k+1) &
             + w111(m)*ez(i+1,j+1,k+1)

         delbxp = w000(m)*delbx(i,j,k)  &
             + w100(m)*delbx(i+1,j,k) &
             + w010(m)*delbx(i,j+1,k) &
             + w110(m)*delbx(i+1,j+1,k) &
             + w001(m)*delbx(i,j,k+1) &
             + w101(m)*delbx(i+1,j,k+1) &
             + w011(m)*delbx(i,j+1,k+1) &
             + w111(m)*delbx(i+1,j+1,k+1)

         delbyp = w000(m)*delby(i,j,k)  &
             + w100(m)*delby(i+1,j,k) &
             + w010(m)*delby(i,j+1,k) &
             + w110(m)*delby(i+1,j+1,k) &
             + w001(m)*delby(i,j,k+1) &
             + w101(m)*delby(i+1,j,k+1) &
             + w011(m)*delby(i,j+1,k+1) &
             + w111(m)*delby(i+1,j+1,k+1)

         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp+vpar/b*delbxp)*dum1
         xdot = xdot+vxdum*nonline
         ydot = ydot+(-exp1+vpar/b*delbyp)*dum1*nonline

         phip = w000(m)*phi(i,j,k)  &
             + w100(m)*phi(i+1,j,k) &
             + w010(m)*phi(i,j+1,k) &
             + w110(m)*phi(i+1,j+1,k) &
             + w001(m)*phi(i,j,k+1) &
             + w101(m)*phi(i+1,j,k+1) &
             + w011(m)*phi(i,j+1,k+1) &
             + w111(m)*phi(i+1,j+1,k+1)

         wght1 = (wght+1./dv*phip*isg)
         wght2 = wght*zdot

         myupex(i,j,k)      =myupex(i,j,k)+wght1*w000(m)*xdot
         myupex(i+1,j,k)    =myupex(i+1,j,k)+wght1*w100(m)*xdot
         myupex(i,j+1,k)    =myupex(i,j+1,k)+wght1*w010(m)*xdot
         myupex(i+1,j+1,k)  =myupex(i+1,j+1,k)+wght1*w110(m)*xdot
         myupex(i,j,k+1)    =myupex(i,j,k+1)+wght1*w001(m)*xdot
         myupex(i+1,j,k+1)  =myupex(i+1,j,k+1)+wght1*w101(m)*xdot
         myupex(i,j+1,k+1)  =myupex(i,j+1,k+1)+wght1*w011(m)*xdot
         myupex(i+1,j+1,k+1)=myupex(i+1,j+1,k+1)+wght1*w111(m)*xdot
         
         myupey(i,j,k)      =myupey(i,j,k)+wght1*w000(m)*ydot
         myupey(i+1,j,k)    =myupey(i+1,j,k)+wght1*w100(m)*ydot
         myupey(i,j+1,k)    =myupey(i,j+1,k)+wght1*w010(m)*ydot
         myupey(i+1,j+1,k)  =myupey(i+1,j+1,k)+wght1*w110(m)*ydot
         myupey(i,j,k+1)    =myupey(i,j,k+1)+wght1*w001(m)*ydot
         myupey(i+1,j,k+1)  =myupey(i+1,j,k+1)+wght1*w101(m)*ydot
         myupey(i,j+1,k+1)  =myupey(i,j+1,k+1)+wght1*w011(m)*ydot
         myupey(i+1,j+1,k+1)=myupey(i+1,j+1,k+1)+wght1*w111(m)*ydot

         myupazd(i,j,k)      =myupazd(i,j,k)+wght2*w000(m)
         myupazd(i+1,j,k)    =myupazd(i+1,j,k)+wght2*w100(m)
         myupazd(i,j+1,k)    =myupazd(i,j+1,k)+wght2*w010(m)
         myupazd(i+1,j+1,k)  =myupazd(i+1,j+1,k)+wght2*w110(m)
         myupazd(i,j,k+1)    =myupazd(i,j,k+1)+wght2*w001(m)
         myupazd(i+1,j,k+1)  =myupazd(i+1,j,k+1)+wght2*w101(m)
         myupazd(i,j+1,k+1)  =myupazd(i,j+1,k+1)+wght2*w011(m)
         myupazd(i+1,j+1,k+1)=myupazd(i+1,j+1,k+1)+wght2*w111(m)

         dum2 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp+vpar/b*delbxp)*dum2 
         wght0 = wght0*(vxdum*kap-xdot*exp1/ter-(ydot-vp0)*eyp/ter)*xnp

         mydnedt(i,j,k)      =mydnedt(i,j,k)+wght0*w000(m)
         mydnedt(i+1,j,k)    =mydnedt(i+1,j,k)+wght0*w100(m)
         mydnedt(i,j+1,k)    =mydnedt(i,j+1,k)+wght0*w010(m)
         mydnedt(i+1,j+1,k)  =mydnedt(i+1,j+1,k)+wght0*w110(m)
         mydnedt(i,j,k+1)    =mydnedt(i,j,k+1)+wght0*w001(m)
         mydnedt(i+1,j,k+1)  =mydnedt(i+1,j,k+1)+wght0*w101(m)
         mydnedt(i,j+1,k+1)  =mydnedt(i,j+1,k+1)+wght0*w011(m)
         mydnedt(i+1,j+1,k+1)=mydnedt(i+1,j+1,k+1)+wght0*w111(m)

      enddo

!   enforce periodicity
      call enforce(myupex(:,:,:))
      call enforce(myupey(:,:,:))
      call enforce(myupazd(:,:,:))
      call enforce(mydnedt(:,:,:))
!      call filter(mydnedt(:,:,:))

      do  i=0,im
         do  j=0,jm
            do  k=0,mykm
               upex(i,j,k) = myupex(i,j,k)/n0/jac(i,k)*cn0e
               upey(i,j,k) = myupey(i,j,k)/n0/jac(i,k)*cn0e
               upazd(i,j,k) = myupazd(i,j,k)/n0/jac(i,k)*cn0e
               dnedt(i,j,k) = mydnedt(i,j,k)/n0/jac(i,k)*pidum*cn0e
            end do
         end do
      end do

! energetic particles
      if(ishgk==0)return
      do m=1,mmb
         dv=float(lr(1))*(dx*dy*dz)

         r=xh3(m)-0.5*lx+lr0

         k = int(zh3(m)/delz)
         wz0 = ((k+1)*delz-zh3(m))/delz
         wz1 = 1-wz0
         th = wz0*thfnz(k)+wz1*thfnz(k+1)

         cost=cos(th)
         sint=sin(th)
         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0
         dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
                 +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
         dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
                 +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
         grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
                 +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
         radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
                 +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
         dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
                 +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
         qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
                 +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
         grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
                 +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
         gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
                 +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 
         fp = wx0*f(i)+wx1*f(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
!         ter = wx0*t0i(i)+wx1*t0i(i+1)        
!         kaptp = wx0*capti(i)+wx1*capti(i+1)        

!         b=1.-lr0/br0*cost
         b=1.-tor+tor*bfldp

         rhog=sqrt(2.*b*muh(m)*mbeam)/(qbeam*b)*iflrh
         vfac = 0.5*(mbeam*uh3(m)**2 + 2.*muh(m)*b)
         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))
         vpar = uh3(m)
         wght=wh3(m)/dv*nh
         wght0 = exp(-vfac)/dv*nh
         if(isuni.eq.0)wght0 = 1./dv*nh

         rhox(1) = rhog*(1-tor)+rhog*grp*tor
         rhoy(1) = rhog*gxdgyp/grp*tor
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
         rhox(4) = 0
         rhoy(4) = -rhoy(3)

         exp1=0.
         eyp=0.
         delbxp = 0.
         delbyp = 0.
         aparp = 0.
         do 300 l=1,lr(1)
            xs=xh3(m)+rhox(l) !rwx(1,l)*rhog
            yt=yh3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            i=int(xt/dx+0.5)
            j=int(yt/dy+0.5)
            k=int(zh3(m)/dz+0.5)-gclr*kcnt
            exp1=exp1 + ex(i,j,k)
            eyp=eyp + ey(i,j,k)
            delbxp = delbxp+delbx(i,j,k)
            delbyp = delbyp+delby(i,j,k)
            aparp = aparp+apar(i,j,k)
 300     continue

         exp1=exp1/4.
         eyp=eyp/4.
         delbxp=delbxp/4.
         delbyp=delbyp/4.
         aparp = aparp/4.

         vpar = uh3(m)
         enerb=(muh(m)+mbeam*vpar*vpar/b)/qbeam*tor
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp+vpar/b*delbxp)*dum1

         xdot = vxdum*nonlin  &
              -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1+vpar/b*delbyp)*dum1*nonlin  &
             +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0
         zdot =  vpar*(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

!    now do 1,2,4 point average, where lr is the no. of points...
         do 400 l=1,lr(1)
            xs=xh3(m)+rhox(l) !rwx(1,l)*rhog
            yt=yh3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            i=int(xt/dx+0.5)
            j=int(yt/dy+0.5)
            k=int(zh3(m)/dz+0.5)-gclr*kcnt

            kap = 1/lh
            ter = teth
            wght1 = wght0*(vxdum*kap+xdot*exp1/ter+(ydot-vp0)*eyp/ter)
            myhpar(i,j,k) = myhpar(i,j,k)+wght*zdot
            myhpex(i,j,k) = myhpex(i,j,k)+wght*xdot
            myhpey(i,j,k) = myhpey(i,j,k)+wght*ydot
            mydnhdt(i,j,k) = mydnhdt(i,j,k)+wght1
 400     continue

      enddo

!   enforce periodicity
      call enforce(myhpar)
      call enforce(myhpex)
      call enforce(myhpey)
      call enforce(mydnhdt)
!      call filter(myjpar)
!      call filter(mydnidt)

      do 510 i=0,im
         do 520 j=0,jm
            do 530 k=0,mykm
               hpar(i,j,k) = qbeam*myhpar(i,j,k)/n0/jac(i,k)
               hpex(i,j,k) = qbeam*myhpex(i,j,k)/n0/jac(i,k)
               hpey(i,j,k) = qbeam*myhpey(i,j,k)/n0/jac(i,k)
               dnhdt(i,j,k) = qbeam*mydnhdt(i,j,k)/n0/jac(i,k)*pidum
 530        continue
 520     continue
 510  continue

 999  continue
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lorentz(ip,n)
      use gem_com
      use equil
      implicit none
      integer :: i,ip,k,m,n,ncol,icol
      real(8) :: edum,vdum,dum,dum1,ptch,vte,r,qr,th,cost,b
      real(8) :: h_x,h_coll,x,eps,dtcol,uold
      real(8) :: wx0,wx1,wz0,wz1

      ncol = 1
      if(ip.eq.1)dtcol = dt/ncol*2
      if(ip.eq.0)dtcol = dt/ncol
      vte = sqrt(amie)
      if(rneu==0.0)return
      do k = 1,mme
         r=x3e(k)-0.5*lx+lr0

         m = int(z3e(k)/delz)
         wz0 = ((m+1)*delz-z3e(k))/delz
         wz1 = 1-wz0
         th = wz0*thfnz(m)+wz1*thfnz(m+1)
         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         m = int((th+pi)/dth)
         wz0 = (-pi+(m+1)*dth-th)/dth
         wz1 = 1.-wz0
         b = wx0*wz0*bfld(i,m)+wx0*wz1*bfld(i,m+1) &
                 +wx1*wz0*bfld(i+1,m)+wx1*wz1*bfld(i+1,m+1) 
         uold = u3e(k)
         edum = b*mue3(k)+0.5*emass*u3e(k)*u3e(k)
         vdum = sqrt(2.*edum/emass)
         ptch = u3e(k)/vdum
         eps = edum
         x = sqrt(eps)
         h_x    = 4.0*x*x/(3.0*sqrt(pi))
         h_coll = h_x/sqrt(1.0+h_x**2)
!         h_coll = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*DERF(x)
         dum1 = 1/eps**1.5*(1+h_coll)
         dum = dtcol*rneu*dum1
         if(x<0.136)dum=0.0
         do icol = 1,ncol
            ptch =ptch-dum*ptch+sqrt(dum*(1.-ptch*ptch)) &
                *dsign(1.d0,ran2(iseed)-0.5)
            ptch = dmin1(ptch,0.999d0)
            ptch = dmax1(ptch,-0.999d0)
         end do
         u3e(k) = vdum*ptch
         mue3(k) = 0.5*emass*vdum*vdum*(1.-ptch*ptch)/b
!         u1e(k) = u1e(k)+u3e(k)-uold
         if(ip.eq.0)then
            mue2(k) = mue3(k)
            u2e(k) = u3e(k)
         end if 
!         if(myid.eq.master)write(*,*)mue(k),vdum,k,ptch
      end do
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine initialize
         use gem_com
      use equil
      use fft_wrapper
	implicit none
        real(8) :: dum,dum1,dum2,jacp,xndum,r,wx0,wx1
        complex(8),dimension(0:1) :: x,y
	integer :: n,i,j,k,ip
	
        call ppinit(MyId,numprocs,ntube,TUBE_COMM,GRID_COMM)
!	write(*,*)'ppinit  ',myid,numprocs,ntube,TUBE_COMM,GRID_COMM

!     program begins....

!     reset timestep counter.
         Last=numprocs-1
         timestep=0
         tcurr = 0.

!     read input data and initialize...

         call hybinit
         do i=0,Last
            if (MyId.eq.i) call init
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         enddo
         if(idg.eq.1)write(*,*)'past init'

      dum = 0.
      dum1 = 0.
      do i = 0,im-1
         r = xg(i)-0.5*lx+lr0
         j = int((r-rin)/dr)
         j = min(j,nr-1)
         wx0 = (rin+(j+1)*dr-r)/dr
         wx1 = 1.-wx0
         xndum = wx0*xn0e(j)+wx1*xn0e(j+1)
         dum = dum+(jac(i,0)+jac(i+1,0)+jac(i,1)+jac(i+1,1))/4
         dum1 = dum1+xndum*(jac(i,0)+jac(i+1,0)+jac(i,1)+jac(i+1,1))/4
      end do
      call MPI_ALLREDUCE(dum,jacp,1,  &
          MPI_REAL8,MPI_SUM,           &
          tube_comm,ierr)

      call MPI_ALLREDUCE(dum1,dum2,1,  &
          MPI_REAL8,MPI_SUM,           &
          tube_comm,ierr)

      totvol = dx*ly*dz*jacp    
      n0=float(tmm(1))/totvol
      if(myid==0)then
         write(*,*)'totvol,jacp,dum2=',totvol,jacp,dum2
      end if
	call weight
!     initialize particle quantities...
         if( cut.eq.0.) cut=1.
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         call ccfft('x',0,imx,0.,tmpx,tmpx,coefx,workx,0)
         call ccfft('y',0,jmx,0.,tmpy,tmpy,coefy,worky,0)
         call ccfft('z',0,kmx,0.,tmpz,tmpz,coefz,workz,0)
         call dsinf(1,x,1,0,y,1,0,imx*2,1,1.d0,aux1,50000,aux2,20000)

         ncurr = 1
         call blendf
	end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine loader_wrapper
         use gem_com
      use equil
	implicit none

	integer :: n,i,j,k,ip,ns

        do ns = 1,nsm
           if(isuni.eq.0)call loadi(ns)
        end do
        if(ifluid.eq.1)call ldel
        if(ishgk==1)call loadh
        if(izon.eq.1)then
           do i = 1,200
              call accumulate(0,0)
              call poisson(0,0)
              if(mod(i,20).eq.0)call outd(0)
              if(ifluid.eq.1)call ldel
           end do
        end if
        if(idg.eq.1)write(*,*)'past loader'
end subroutine loader_wrapper
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine accumulate(n,ip)
         use gem_com
         use equil
	implicit none

	integer :: n,i,j,k,ip
        if(ifluid==1)call setw(ip,n)
	call grid1(ip,n)
	if(idg.eq.1)write(*,*)'pass grid1'
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine accumulate
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine poisson(n,ip)
         use gem_com
         use equil
	implicit none
	integer :: n,i,i1,j,k,ip,it,iter=5
        real(8) :: myrmsphi,rmp(20)

        do it = 1,iter
           call den0_phi(ip,n,it)
           call gkps(n,ip)
           myrmsphi=0.
           rmp(it)=0.
           do k=0,mykm-1
              do j=0,jm-1
                 do i1=0,im-1
                    myrmsphi=myrmsphi+phi(i1,j,k)*phi(i1,j,k)
                 enddo
              enddo
           enddo
           call MPI_ALLREDUCE(myrmsphi,rmp(it),1, &
              MPI_REAL8,                               &
              MPI_SUM,TUBE_COMM,ierr)
           rmp(it)=sqrt(rmp(it)/(im*jm*km))
        end do
        rmsphi(n)=rmp(iter)        
        if(ip==1)ipred = 1
        if(ip==0)icorr = 1
        if(iter==1)goto 100
        if(rmp(iter)/rmp(iter-1)>1.1)then
           phi(:,:,:) = 0.
           if(ip==1)ipred = 0
           if(ip==0)icorr = 0
           rmsphi(n)=0
        end if
        if(n>100.and.rmp(iter)/rmpp>1.5)then
           phi(:,:,:) = 0.
           if(ip==1)ipred = -1
           if(ip==0)icorr = -1
           rmsphi(n)=0           
        end if
        rmpp = rmp(iter)
 100    continue
        call fltx(phi(:,:,:),0)
        call filter(phi(:,:,:))
	if(idg.eq.1)write(*,*)'pass gkps'
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine poisson
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine ampere(n,ip)
         use gem_com
         use equil
	implicit none
        integer :: iter=10
	integer :: n,i,i1,j,k,ip
        real(8) :: myrmsapa,rma(20)

        if(ifluid==1.and.beta.gt.1.e-8)then
           do i = 1,iter
              call jpar0(ip,n,i,0)
              if(idg.eq.1)write(*,*)'pass jpar0'
              call ezamp(n,ip)
              if(idg.eq.1)write(*,*)'pass ezamp'
              if(i==iter)call jpar0(ip,n,i,1)
              myrmsapa=0.
              rma(i)=0.
              do k=0,mykm-1
                 do j=0,jm-1
                    do i1=0,im-1
                       myrmsapa=myrmsapa+apar(i1,j,k)*apar(i1,j,k)
                    enddo
                 enddo
              enddo
              call MPI_ALLREDUCE(myrmsapa,rma(i),1, &
                 MPI_REAL8,                               &
                 MPI_SUM,TUBE_COMM,ierr)
              rma(i)=sqrt(rma(i)/(im*jm*km))
           end do
        end if
        rmsapa(n) = rma(iter)
        if(ip==1)jpred = 1
        if(ip==0)jcorr = 1
        if(rma(iter)/rma(iter-1)>1.1)then
           apar(:,:,:) = 0.
           if(ip==1)jpred = 0
           if(ip==0)jcorr = 0
           rmsapa(n) = 0.
        end if
        if(n>100.and.rma(iter)/rmaa>1.5)then
           apar(:,:,:) = 0.
           if(ip==1)jpred = -1
           if(ip==0)jcorr = -1
           rmsapa(n) = 0.
        end if
        rmaa=rma(iter)
        call fltx(apar(:,:,:),0)
        call filter(apar(:,:,:))
	if(idg.eq.1)write(*,*)'pass filter(apar)'
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine ampere
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine split_weight(n,ip)
         use gem_com
         use equil
	implicit none

	integer :: n,i,j,k,ip
            if(isg.gt.0..and.ifluid.eq.1)then
               call jie(ip,n)
	if(idg.eq.1)write(*,*)'pass jie'
               call drdt(ip)
	if(idg.eq.1)write(*,*)'pass drdt'
               call dpdt(ip)
	if(idg.eq.1)write(*,*)'pass dpdt'
            end if
	if(idg.eq.1)write(*,*)'pass split_weight'
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine split_weight
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine field(n,ip)
         use gem_com
         use equil
	implicit none
	integer :: n,i,j,k,ip
	call grad(ip)
	if(idg.eq.1)write(*,*)'pass grad'
	call eqmo(ip)
	if(idg.eq.1)write(*,*)'pass eqmo'
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine field
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine push_wrapper(n,ip)
         use gem_com
         use equil
	implicit none
	integer :: n,i,j,k,ip,ns

        do ns = 1,nsm
            if(ip.eq.1.and.ision==1)call ppush(n,ns)
            if(ip.eq.0.and.ision==1)call cpush(n,ns)
         end do
         if(idg.eq.1)write(*,*)'pass ppush'

         if(ip.eq.1.and.ishgk==1)call hpush(n)
         if(ip.eq.0.and.ishgk==1)call hcush(n)
         if(idg.eq.1)write(*,*)'pass hcush'

         if(ip.eq.1.and.ifluid==1)call pint
         if(ip.eq.0.and.ifluid==1)call cint(n)
         if(idg.eq.1)write(*,*)'pass pint'
         if(ifluid==1.and.ip==0)call lorentz(ip,n)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)        		
end subroutine push_wrapper
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine diagnose(n)
         use gem_com
         use equil
	implicit none
	integer :: n,i,j,k,ip

        call modes2(phi,pmodehis,n)
        if(idg.eq.1)write(*,*)'pass modes2'  
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)          
        call yveck(phi(:,:,:),n)
        if(idg.eq.1)write(*,*)'pass yvec'    
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)        

end subroutine diagnose
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine reporter(n)
         use gem_com
         use equil
	implicit none
	integer :: n,i,j,k,ip

        if(mod(n,xnplt).eq.0) then
           call spec(n)
        endif
 13     format(1x,i6,7(2x,i7))
        if(myid==0)write(16,13)n,ipred,icorr,jpred,jcorr,nopz,noen,nowe
!        if(myid==0)write(*,13)n,ipred,icorr,jpred,jcorr,nopz,noen
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!     save particle arrays for restart if iput=1...
!     do this before the code crashes due to graphics problems
        if((iput.eq.1).and.mod(n+1,100).eq.0)call restart(2,n)

!     periodically make output for plots
        call outd(n)
        if(idg.eq.1)write(*,*)'pass outd'

end subroutine reporter
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine jpar0(ip,n,it,itp)

!    source quantities are are calculated: n_i
!    right now only ion quantitities are calculated...

      use gem_com
      use equil
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: enerb,vxdum,dum,xdot,ydot,avede0
      INTEGER :: m,n,i,j,k,l,ns,ip,it,itp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,vte
      REAL(8) :: sz,wght,wght0,wght1,r,th,cost,sint,b,qr,dv,xnp
      REAL(8) :: xt,yt,zt,rhog,pidum,vpar,xs,dely,vfac,ter
      REAL(8) :: lbfs(0:imx,0:jmx)
      REAL(8) :: rbfs(0:imx,0:jmx)
      REAL(8) :: lbfr(0:imx,0:jmx)
      REAL(8) :: rbfr(0:imx,0:jmx)
      real(8) :: myupa(0:imx,0:jmx,0:1),myupa0(0:imx,0:jmx,0:1)
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
      REAL(8) :: grp,gxdgyp
      REAL(8) :: zdot

      pidum = 1./(pi*2)**1.5*(vwidthe)**3
      if(isuni.eq.0)pidum = 1.
      ns=1

!    electrons current due to f_M (p_para)
      vte = sqrt(amie)
      myupa = 0.
      myupa0 = 0.
!      upa0(:,:,:) = apar(ip,:,:,:)
!      return
      if(it.eq.1)then
         apar(:,:,:) = 0.
         upa0(:,:,:) = 0.
         upa00(:,:,:) = 0.
         return
      end if
      do m=1,mme
         r=x3e(m)-0.5*lx+lr0
         vpar = u3e(m)

         k = int(z3e(m)/delz)
         wz0 = ((k+1)*delz-z3e(m))/delz
         wz1 = 1-wz0
         th = wz0*thfnz(k)+wz1*thfnz(k+1)
            
         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0

         ter = wx0*t0e(i)+wx1*t0e(i+1)        
         xnp = wx0*xn0e(i)+wx1*xn0e(i+1)        
         if(itp==1)then
            grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
                 +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
            radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
                 +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
            bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
            dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
                 +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
            jfnp = wz0*jfn(k)+wz1*jfn(k+1)
            psipp = wx0*psip(i)+wx1*psip(i+1)        
            fp = wx0*f(i)+wx1*f(i+1)
            b=1.-tor+tor*bfldp
            enerb=(mue3(m)+emass*vpar*vpar/b)/qel*tor
            zdot =  vpar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp
         end if

         dv=(dx*dy*dz)
         wght0 = 1/dv
         wght1 = 1/dv
         if(abs(vpar/vte).gt.vcut)then
            wght0 = 0.
            wght1 = 0.
         end if

         xt=x3e(m)
         yt=y3e(m)
         zt=z3e(m)
         i=int(xt/dx)
         j=int(yt/dy)
         k=0 !int(z3e(m)/dz)-gclr*kcnt
         aparp = w000(m)*apar(i,j,k)  &
             + w100(m)*apar(i+1,j,k) &
             + w010(m)*apar(i,j+1,k) &
             + w110(m)*apar(i+1,j+1,k) &
             + w001(m)*apar(i,j,k+1) &
             + w101(m)*apar(i+1,j,k+1) &
             + w011(m)*apar(i,j+1,k+1) &
             + w111(m)*apar(i+1,j+1,k+1)

         if(itp==0)then
            wght1 = wght1*aparp*vpar*vpar/amie/ter*xnp 
            myupa0(i,j,k)      =myupa0(i,j,k)+wght1*w000(m)
            myupa0(i+1,j,k)    =myupa0(i+1,j,k)+wght1*w100(m)
            myupa0(i,j+1,k)    =myupa0(i,j+1,k)+wght1*w010(m)
            myupa0(i+1,j+1,k)  =myupa0(i+1,j+1,k)+wght1*w110(m)
            myupa0(i,j,k+1)    =myupa0(i,j,k+1)+wght1*w001(m)
            myupa0(i+1,j,k+1)  =myupa0(i+1,j,k+1)+wght1*w101(m)
            myupa0(i,j+1,k+1)  =myupa0(i,j+1,k+1)+wght1*w011(m)
            myupa0(i+1,j+1,k+1)=myupa0(i+1,j+1,k+1)+wght1*w111(m)
         end if
         if(itp==1)then
            wght0 = wght0*aparp*vpar*zdot/amie/ter*xnp 
            myupa(i,j,k)      =myupa(i,j,k)+wght0*w000(m)
            myupa(i+1,j,k)    =myupa(i+1,j,k)+wght0*w100(m)
            myupa(i,j+1,k)    =myupa(i,j+1,k)+wght0*w010(m)
            myupa(i+1,j+1,k)  =myupa(i+1,j+1,k)+wght0*w110(m)
            myupa(i,j,k+1)    =myupa(i,j,k+1)+wght0*w001(m)
            myupa(i+1,j,k+1)  =myupa(i+1,j,k+1)+wght0*w101(m)
            myupa(i,j+1,k+1)  =myupa(i,j+1,k+1)+wght0*w011(m)
            myupa(i+1,j+1,k+1)=myupa(i+1,j+1,k+1)+wght0*w111(m)
         end if
      enddo

!   enforce periodicity
      if(itp==1)call enforce(myupa(:,:,:))
      if(itp==0)call enforce(myupa0(:,:,:))
!      call filter(myupa(:,:,:))

      if(itp==1)then
         do  i=0,im
            do  j=0,jm
               do  k=0,mykm
                  upa0(i,j,k)= myupa(i,j,k)/n0/jac(i,k)*pidum*cn0e
               end do
            end do
         end do
      end if
      if(itp==0)then
         do  i=0,im
            do  j=0,jm
               do  k=0,mykm
                  upa00(i,j,k)= myupa0(i,j,k)/n0/jac(i,k)*pidum*cn0e
               end do
            end do
         end do
      end if

      upa0(:,:,:) = upa0(:,:,:)+ &
                   (isg*phi(:,:,:)+dene(:,:,:))*apar(:,:,:)*nonline*0.

      upa00(:,:,:) = upa00(:,:,:)+ &
                   (isg*phi(:,:,:)+dene(:,:,:))*apar(:,:,:)*nonline*0.
 999  continue

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine den0_phi(ip,n,it)

!    source quantities are are calculated: n_i
!    right now only ion quantitities are calculated...

      use gem_com
      use equil
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: enerb,vxdum,dum,xdot,ydot,avede0
      INTEGER :: m,n,i,j,k,l,ns,ip,it
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,vte
      REAL(8) :: sz,wght,wght0,r,th,cost,sint,b,qr,dv
      REAL(8) :: xt,yt,zt,rhog,pidum,vpar,xs,dely,vfac
      REAL(8) :: lbfs(0:imx,0:jmx)
      REAL(8) :: rbfs(0:imx,0:jmx)
      REAL(8) :: lbfr(0:imx,0:jmx)
      REAL(8) :: rbfr(0:imx,0:jmx)
      real(8) :: myden0(0:imx,0:jmx,0:1)

      pidum = 1./(pi*2)**1.5*(vwidthe)**3
      if(isuni.eq.0)pidum = 1.
      ns=1

! electrons density due to phi*f_M (p_para)
      vte = sqrt(amie)
      myden0 = 0.
      den0 = 0.
      if(it.eq.1)then
         phi(:,:,:) = 0.
         den0(:,:,:) = 0.
         return
      end if
      do m=1,mme
         dv=(dx*dy*dz)
         wght0 = 1./dv
         xt=x3e(m)
         yt=y3e(m)
         zt=z3e(m)
         i=int(xt/dx)
         j=int(yt/dy)
         k=0 !int(z3e(m)/dz)-gclr*kcnt

         phip = 0.
         phip = w000(m)*phi(i,j,k)  &
             + w100(m)*phi(i+1,j,k) &
             + w010(m)*phi(i,j+1,k) &
             + w110(m)*phi(i+1,j+1,k) &
             + w001(m)*phi(i,j,k+1) &
             + w101(m)*phi(i+1,j,k+1) &
             + w011(m)*phi(i,j+1,k+1) &
             + w111(m)*phi(i+1,j+1,k+1)

         wght0 = wght0*phip
         myden0(i,j,k)      =myden0(i,j,k)+wght0*w000(m)
         myden0(i+1,j,k)    =myden0(i+1,j,k)+wght0*w100(m)
         myden0(i,j+1,k)    =myden0(i,j+1,k)+wght0*w010(m)
         myden0(i+1,j+1,k)  =myden0(i+1,j+1,k)+wght0*w110(m)
         myden0(i,j,k+1)    =myden0(i,j,k+1)+wght0*w001(m)
         myden0(i+1,j,k+1)  =myden0(i+1,j,k+1)+wght0*w101(m)
         myden0(i,j+1,k+1)  =myden0(i,j+1,k+1)+wght0*w011(m)
         myden0(i+1,j+1,k+1)=myden0(i+1,j+1,k+1)+wght0*w111(m)
      enddo

!   enforce periodicity
      call enforce(myden0(:,:,:))

      do  i=0,im
         do  j=0,jm
            do  k=0,mykm
               den0(i,j,k)= myden0(i,j,k)/n0/jac(i,k)*pidum*cn0e
            end do
         end do
      end do

      return
      end
!!!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine countw(n)
      use gem_com
      use equil
      implicit none

      INTEGER :: n,nbin,m,i,j,k
      parameter(nbin=200)
      REAL(8) :: myavwe,avwe,mymaxw,myminw,maxw,minw
      REAL(8) :: sbuf(10),rbuf(10)
      REAL(8),dimension(:),allocatable :: wghtmax,wghtmin
      REAL(8) :: dw
      integer :: npbin(0:nbin),mynpbin(0:nbin)
      real(8) :: wpbin(0:nbin),mywpbin(0:nbin)

      allocate (wghtmax(0:numprocs-1),wghtmin(0:numprocs-1))

      sbuf(1:10) = 0.
      rbuf(1:10) = 0.
      myavwe = 0.
      avwe = 0.
      mymaxw = 0.
      myminw = 0.
      maxw = 0.
      minw = 0.
      do m=1,mme
         if((w3e(m)).gt.mymaxw)mymaxw = (w3e(m))
         if((w3e(m)).lt.myminw)myminw = (w3e(m))
         myavwe = myavwe+abs(w3e(m))
      end do
      call MPI_GATHER(mymaxw,1, &
               MPI_REAL8, &
               wghtmax,1,MPI_REAL8, &
               Last,mpi_comm_world,ierr)
      call MPI_GATHER(myminw,1, &
               MPI_REAL8, &
               wghtmin,1,MPI_REAL8, &
               Last,mpi_comm_world,ierr)
      if(myid.eq.Last)then
         do i = 0,Last
            if(wghtmax(i).gt.maxw)maxw = wghtmax(i)
            if(wghtmin(i).lt.minw)minw = wghtmin(i)
         end do
      end if
      call MPI_BCAST(maxw,1,MPI_REAL8,Last,Mpi_comm_world,ierr)
      call MPI_BCAST(minw,1,MPI_REAL8,Last,Mpi_comm_world,ierr)

      dw = (maxw-minw)/nbin
      mynpbin = 0
      npbin = 0
      wpbin = 0.
      mywpbin = 0.
      do m = 1,mme
         i = int(((w3e(m))-minw)/dw)
         mynpbin(i) = mynpbin(i)+1
         mywpbin(i) = mywpbin(i)+abs(w3e(m))
         if(abs(w3e(m)).gt.4.)then
            w3e(m) = 0.
            w2e(m) = 0.
         end if
      end do
      call MPI_ALLREDUCE(mynpbin(0:nbin),npbin(0:nbin),nbin+1, &
          MPI_integer, &
          MPI_SUM, &
          MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(mywpbin(0:nbin),wpbin(0:nbin),nbin+1, &
          MPI_real8, &
          MPI_SUM, &
          MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(myavwe,avwe,1, &
          MPI_real8, &
          MPI_SUM, &
          MPI_COMM_WORLD,ierr)
      avwe = avwe/tmm(1)
      wpbin = wpbin/tmm(1)
      if(myid.eq.0.and.mod(n,xnplt).eq.0)then
!         write(*,*)'maxw,minw,avwe= ', maxw, minw,avwe
         do i = 0,nbin
!            write(*,*)n,(minw+i*dw)/avwe,npbin(i),wpbin(i)
         end do
      end if

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine setw(ip,n)
      use gem_com
      implicit none
      INTEGER :: m,n,ip,i,j,k
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,vte
      REAL(8) :: wght
      REAL(8) :: xt,yt,zt

      vte = sqrt(amie)
      do m=1,mme
         xt=x3e(m)
         yt=y3e(m)
         zt=z3e(m)
         i=int(xt/dx)
         j=int(yt/dy)
         k=0 !int(z3e(m)/dz)-gclr*kcnt
         wx0=float(i+1)-xt/dx 
         wx1=1.-wx0
         wy0=float(j+1)-yt/dy
         wy1=1.-wy0
         wz0=float(gclr*kcnt+k+1)-z3e(m)/dz
         wz1=1.-wz0

         wght = 1.
         if(abs(u3e(m)/vte).gt.vcut)wght = 0.

         w000(m)=wx0*wy0*wz0*wght
         w001(m)=wx0*wy0*wz1*wght
         w010(m)=wx0*wy1*wz0*wght
         w011(m)=wx0*wy1*wz1*wght
         w100(m)=wx1*wy0*wz0*wght
         w101(m)=wx1*wy0*wz1*wght
         w110(m)=wx1*wy1*wz0*wght
         w111(m)=wx1*wy1*wz1*wght
      end do
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      subroutine fltx(u,isbl)   
      use gem_com
      use fft_wrapper
      implicit none
      INTEGER :: i,j,k,k1,l,m,n,ifirst,nstep,ip,INFO,isbl
      INTEGER :: l1,m1,myk,myj,izonal,ix,ikx
      real(8) :: u(0:imx,0:jmx,0:1)
      COMPLEX(8) :: temp3dxy(0:imx-1,0:jmx-1,0:1)
      real(8) :: filterx(0:imx-1,0:jmx-1),gr(0:imx-1),gi(0:imx-1)
      real :: kx,ky,kx0,th,shat,sgny

      do j = 0,jmx-1
         do i = 0,imx-1
            if(j.ge.(jm/2+1)) then
               m1=jm-j
               sgny=-1.
            else
               m1=j
               sgny=1.
            endif
            ky=sgny*2.*pi*float(m1)/ly
            kx = i*pi/lx

            filterx(i,j) = 2.0/imx
            if(abs(ky)>kycut)filterx(i,j) = 0.
            if(kx>kxcut)filterx(i,j) = 0.
            if(onemd==1.and.m1.ne.mlk)filterx(i,j) = 0.
            if(j==0.and.i<=1)filterx(i,j) = 0.
         end do
      end do

      do k=0,mykm
         do j=0,jm-1
            do i=0,imx-1
               temp3dxy(i,j,k)=u(i,j,k)
            enddo
         enddo
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3dxy(i,j,k)
            end do
            call ccfft('y',-1,jmx,1.,tmpy,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3dxy(i,j,k) = tmpy(j)   !rho(ky,x)
            end do
         end do
      enddo

      do k = 0,1
         do j = 0,jmx-1
            gr(:) = real(temp3dxy(:,j,k))
            gi(:) = aimag(temp3dxy(:,j,k))
       call dsinf(0,gr(0),1,0,gr(0),1,0,imx*2,1, 1.d0,aux1,50000,aux2,20000)  
       call dsinf(0,gi(0),1,0,gi(0),1,0,imx*2,1, 1.d0,aux1,50000,aux2,20000)  
            gr(:) = gr(:)*filterx(:,j)
            gi(:) = gi(:)*filterx(:,j)
       call dsinf(0,gr(0),1,0,gr(0),1,0,imx*2,1, 1.d0,aux1,50000,aux2,20000)  
       call dsinf(0,gi(0),1,0,gi(0),1,0,imx*2,1, 1.d0,aux1,50000,aux2,20000)  
            temp3dxy(:,j,k) = cmplx(gr(:),gi(:))
         end do
      end do

      temp3dxy(:,:,:)=temp3dxy(:,:,:)/jmx
      if(isbl==1)call filtbl(temp3dxy(0:imx-1,0:jmx-1,0:1))

      do k=0,mykm
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3dxy(i,j,k)
            end do
            call ccfft('y',1,jmx,1.,tmpy,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3dxy(i,j,k) = tmpy(j)  ! phi(x,y)
            end do
         end do
      end do

 100  continue
      do i = 0,imx-1
         do j = 0,jm-1
            do k = 0,mykm
               u(i,j,k) = temp3dxy(i,j,k)
            end do
         end do
      end do

      call enfxy(u(:,:,:))
      call enfz(u)
!      call filter(dphidt(:,:,:))

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dcmpy(u,v)   
      use gem_com
      use fft_wrapper
      implicit none
      INTEGER :: i,j,k,k1,l,m,n,ifirst,nstep,ip,INFO
      INTEGER :: l1,m1,myk,myj,izonal,ix,ikx,id
      INTEGER :: recvcnt(0:ntube-1)      
      real(8) :: u(0:imx-1,0:jmx-1,0:1)
      complex(8) :: v(0:imx-1,0:jcnt-1,0:1)
      complex(8) :: sbuf(0:imx*jmx*2-1),rbuf(0:imx*jcnt*2-1)
      COMPLEX(8) :: temp3d(0:imx-1,0:jmx-1,0:1)
      real :: kx,ky,kx0,th,shat,sgny

      do k=0,mykm
         do j=0,jm-1
            do i=0,imx-1
               temp3d(i,j,k)=u(i,j,k)
            enddo
         enddo
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3d(i,j,k)
            end do
            call ccfft('y',-1,jmx,1.,tmpy,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3d(i,j,k) = tmpy(j)
            end do
         end do
      enddo

      cnt = 2*jcnt*imx
      recvcnt(0:ntube-1) = cnt
      do j = 0,jmx-1
         do k =0,1
            do i = 0,imx-1
               sbuf(j*2*imx+k*imx+i) = temp3d(i,j,k)
            end do
         end do
      end do
      call mpi_reduce_scatter(sbuf,rbuf,recvcnt,mpi_complex16,mpi_sum, &
                              grid_comm,ierr)

      do j = 0,jcnt-1
         do k =0,1
            do i = 0,imx-1
               v(i,j,k) = rbuf(j*2*imx+k*imx+i)
            end do
         end do
      end do

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hpush(n)

      use gem_com
      use equil
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum,dum1
      INTEGER :: m,i,j,k,l,n
      INTEGER :: np_old,np_new
      REAL(8) :: rhog,vfac,sqvfac,kap,vpar,pidum,kaptp,kapnp
      REAL(8) :: b,th,r,enerb,cost,sint,qr,laps,sz,ter
      REAL(8) :: xt,xs,yt,xdot,ydot,zdot,pzdot,edot,pzd0,vp0
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
      REAL(8) :: grp,gxdgyp,rhox(4),rhoy(4),psp,pzp,cvbzp

      pidum = 1./(pi*2)**1.5*vwidth**3
      do m=1,mmb
         r=xh2(m)-0.5*lx+lr0
         k = int(zh2(m)/delz)
         wz0 = ((k+1)*delz-zh2(m))/delz
         wz1 = 1-wz0
         th = wz0*thfnz(k)+wz1*thfnz(k+1)

         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0
         dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
                 +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
         dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
                 +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
         grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
                 +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
         radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
                 +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
         dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
                 +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
         qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
                 +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
         grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
                 +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
         gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
                 +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 
         cvbzp = wx0*wz0*curvbz(i,k)+wx0*wz1*curvbz(i,k+1) &
                 +wx1*wz0*curvbz(i+1,k)+wx1*wz1*curvbz(i+1,k+1) 
         fp = wx0*f(i)+wx1*f(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         psp = wx0*psi(i)+wx1*psi(i+1)        

         b=1.-tor+tor*bfldp
         pzp = mbeam*uh2(m)/b-psp/br0

         rhog=sqrt(2.*b*muh(m)*mbeam)/(qbeam*b)*iflrh

         rhox(1) = rhog*(1-tor)+rhog*grp*tor
         rhoy(1) = rhog*gxdgyp/grp*tor
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
         rhox(4) = 0
         rhoy(4) = -rhoy(3)
!    calculate avg. e-field...
!    do 1,2,4 point average, where lr is the no. of points...

         phip=0.
         exp1=0.
         eyp=0.
         ezp=0.
         delbxp=0.
         delbyp=0.
         dpdzp = 0.
         dadzp = 0.
         aparp = 0.

!  4 pt. avg. done explicitly for vectorization...
         do 200 l=1,lr(1)
!
            xs=xh2(m)+rhox(l) !rwx(1,l)*rhog
            yt=yh2(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
!
!   particle can go out of bounds during gyroavg...
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            include "hpushngp.h"
 200     continue
         exp1 = exp1/4.
         eyp = eyp/4.
         ezp = ezp/4.
         delbxp = delbxp/4.
         delbyp = delbyp/4.
         dpdzp = dpdzp/4.
         dadzp = dadzp/4.
         aparp = aparp/4.
!
         vfac = 0.5*(mbeam*uh2(m)**2 + 2.*muh(m)*b)
         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))

         vpar = uh2(m)
         enerb=(muh(m)+mbeam*vpar*vpar/b)/qbeam*tor
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
         dum = 1/(b*b)*fp/radiusp*cvbzp
         xdot = vxdum*nonlin -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlin &
             +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0
         zdot =  (vpar) &
                 *(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         pzd0 = tor*(-muh(m)/mbeam/radiusp/bfldp*psipp*dbdtp*grcgtp)
         pzdot = pzd0 + (qbeam/mbeam*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
             +qbeam/mbeam*(-xdot*delbyp+ydot*delbxp+zdot*dadzp))*ipara

         edot = xdot*exp1+(ydot-vp0)*eyp+zdot*ezp                      &
             +qbeam*pzdot*aparp*tor     &
             +qbeam*vpar*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)    &
             -qbeam*vpar*delbxp*vp0

         xh3(m) = xh2(m) + 0.5*dt*xdot
         yh3(m) = yh2(m) + 0.5*dt*ydot
         zh3(m) = zh2(m) + 0.5*dt*zdot
         uh3(m) = uh2(m) + 0.5*dt*pzdot

         dum = 1-wh2(m)*nonlin*0.
         if(isuni.eq.1)dum = exp(-vfac)*pidum
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
!         vxdum = eyp+vpar/b*delbxp

         sqvfac = sqrt(vfac)
         kap = 1/lh
         wh3(m)=wh2(m) + 0.5*dt*(vxdum*kap + edot/teth)*dum
!         wh3(m)=wh2(m) + 0.5*dt*(vxdum/lh  &
!                  + edot*3/2*sqvfac/(vfac*sqvfac+ecpow))*dum
         
!         if(xh3(m)>lx .or. xh3(m)<0.)wh3(m) = 0.

!         go to 333
         if(abs(pzp-pzh(m))>3.0.or.abs(vfac-ekh(m))>0.2*ekh(m))then
            xh3(m) = xih(m)
            zh3(m) = z0h(m)
            r = xh3(m)-lx/2+lr0
            k = int(zh3(m)/delz)
            wz0 = ((k+1)*delz-zh3(m))/delz
            wz1 = 1-wz0
            th = wz0*thfnz(k)+wz1*thfnz(k+1)

            i = int((r-rin)/dr)
            wx0 = (rin+(i+1)*dr-r)/dr
            wx1 = 1.-wx0
            k = int((th+pi)/dth)
            wz0 = (-pi+(k+1)*dth-th)/dth
            wz1 = 1.-wz0
            b = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
            uh3(m) = sqrt(2/mbeam*abs(ekh(m)-muh(m)*b))
            uh2(m) = uh3(m)
            wh3(m) = 0.
            wh2(m) = 0.
            xh2(m) = xh3(m)
            zh2(m) = zh3(m)
         end if

 333     continue
         laps=anint((zh3(m)/lz)-.5)*(1-peritr)
         r=xh3(m)-0.5*lx+lr0
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         qr = wx0*sf(i)+wx1*sf(i+1)
         yh3(m)=dmod(yh3(m)-laps*2*pi*qr*lr0/q0+800.*ly,ly)
         if(xh3(m)>lx)then
            xh3(m) = lx-1.e-8
            zh3(m)=lz-zh3(m)
            zh2(m)=zh3(m)
            xh2(m) = xh3(m)
!            wh2(m) = 0.
!            wh3(m) = 0.
         end if
         if(xh3(m)<0.)then
            xh3(m) = 1.e-8
            zh3(m)=lz-zh3(m)
            zh2(m)=zh3(m)
            xh2(m) = xh3(m)
!            wh2(m) = 0.
!            wh3(m) = 0.
         end if
         zh3(m)=dmod(zh3(m)+8.*lz,lz)
         xh3(m)=dmod(xh3(m)+8.*lx,lx)         
         xh3(m) = min(xh3(m),lx-1.0e-8)
         yh3(m) = min(yh3(m),ly-1.0e-8)
         zh3(m) = min(zh3(m),lz-1.0e-2)
         
      enddo

      np_old=mmb
      call init_pmove(zh3,np_old,lz,ierr)

      call pmove(xh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xh3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yh3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zh3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(uh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(uh3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wh3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(muh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call pmove(xih,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0h,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(pzh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ekh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call end_pmove(ierr)
      mmb=np_new

      return
      end

!-----------------------------------------------------------------------

      subroutine hcush(n)

      use gem_com
      use equil
      implicit none
      INTEGER :: n
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,w3old
      INTEGER :: m,i,j,k,ns,l
      INTEGER :: np_old,np_new
      REAL(8) :: vfac,vpar,vxdum,dum,xdot,ydot,zdot,pzdot,edot,pzd0,vp0
      REAL(8) :: rhog,xt,yt,zt,kap,xs,pidum,dum1,sqvfac,kaptp,kapnp
      REAL(8) :: b,th,r,enerb,cost,sint,qr,laps,sz,ter
      REAL(8) :: myke,mypfl,myavewh
      REAL(8) :: myefl,mynos
      REAL(8) :: ketemp,pfltemp
      REAL(8) :: efltemp,nostemp
      REAL(8) :: sbuf(10),rbuf(10)
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
      REAL(8) :: grp,gxdgyp,rhox(4),rhoy(4),psp,pzp,cvbzp

      sbuf(1:10) = 0.
      rbuf(1:10) = 0.
      myavewh = 0.
      myke=0.    
      mypfl=0.    
      myefl=0. 
      mynos=0.   
      ketemp=0.
      pfltemp=0.
      efltemp=0.
      nostemp=0.
      pidum = 1./(pi*2)**1.5*vwidth**3

      do m=1,mmb
         r=xh3(m)-0.5*lx+lr0

         k = int(zh3(m)/delz)
         wz0 = ((k+1)*delz-zh3(m))/delz
         wz1 = 1-wz0
         th = wz0*thfnz(k)+wz1*thfnz(k+1)

         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0
         dbdrp = wx0*wz0*dbdr(i,k)+wx0*wz1*dbdr(i,k+1) &
                 +wx1*wz0*dbdr(i+1,k)+wx1*wz1*dbdr(i+1,k+1) 
         dbdtp = wx0*wz0*dbdth(i,k)+wx0*wz1*dbdth(i,k+1) &
                 +wx1*wz0*dbdth(i+1,k)+wx1*wz1*dbdth(i+1,k+1) 
         grcgtp = wx0*wz0*grcgt(i,k)+wx0*wz1*grcgt(i,k+1) &
                 +wx1*wz0*grcgt(i+1,k)+wx1*wz1*grcgt(i+1,k+1) 
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
         radiusp = wx0*wz0*radius(i,k)+wx0*wz1*radius(i,k+1) &
                 +wx1*wz0*radius(i+1,k)+wx1*wz1*radius(i+1,k+1) 
         dydrp = wx0*wz0*dydr(i,k)+wx0*wz1*dydr(i,k+1) &
                 +wx1*wz0*dydr(i+1,k)+wx1*wz1*dydr(i+1,k+1) 
         qhatp = wx0*wz0*qhat(i,k)+wx0*wz1*qhat(i,k+1) &
                 +wx1*wz0*qhat(i+1,k)+wx1*wz1*qhat(i+1,k+1) 
         grp = wx0*wz0*gr(i,k)+wx0*wz1*gr(i,k+1) &
                 +wx1*wz0*gr(i+1,k)+wx1*wz1*gr(i+1,k+1) 
         gxdgyp = wx0*wz0*gxdgy(i,k)+wx0*wz1*gxdgy(i,k+1) &
                 +wx1*wz0*gxdgy(i+1,k)+wx1*wz1*gxdgy(i+1,k+1) 
         cvbzp = wx0*wz0*curvbz(i,k)+wx0*wz1*curvbz(i,k+1) &
                 +wx1*wz0*curvbz(i+1,k)+wx1*wz1*curvbz(i+1,k+1) 
         fp = wx0*f(i)+wx1*f(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         psp = wx0*psi(i)+wx1*psi(i+1)        
!         ter = wx0*t0i(i)+wx1*t0i(i+1)        
!         kaptp = wx0*capti(i)+wx1*capti(i+1)        
         b=1.-tor+tor*bfldp
         pzp = mbeam*uh3(m)/b-psp/br0

         rhog=sqrt(2.*b*muh(m)*mbeam)/(qbeam*b)*iflrh

         rhox(1) = rhog*(1-tor)+rhog*grp*tor
         rhoy(1) = rhog*gxdgyp/grp*tor
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
         rhox(4) = 0
         rhoy(4) = -rhoy(3)
!    calculate avg. e-field...
!    do 1,2,4 point average, where lr is the no. of points...

         phip=0.
         exp1=0.
         eyp=0.
         ezp=0.
         delbxp = 0.
         delbyp = 0.
         dpdzp = 0.
         dadzp = 0.
         aparp = 0.

!  4 pt. avg. written out explicitly for vectorization...
         do 200 l=1,lr(1)
            xs=xh3(m)+rhox(l) !rwx(1,l)*rhog
            yt=yh3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
!   BOUNDARY
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            include "hcushngp.h"
 200     continue

         exp1=exp1/4.
         eyp=eyp/4.
         ezp=ezp/4.
         delbxp=delbxp/4.
         delbyp=delbyp/4.
         dpdzp = dpdzp/4.
         dadzp = dadzp/4.
         aparp = aparp/4.


         vfac = 0.5*(mbeam*uh3(m)**2 + 2.*muh(m)*b)
         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))

         vpar = uh3(m)
         enerb=(muh(m)+mbeam*vpar*vpar/b)/qbeam*tor
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
         dum = 1/(b*b)*fp/radiusp*cvbzp
         xdot = vxdum*nonlin -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlin     &
             +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0
         zdot =  (vpar) &
                 *(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         pzd0 = tor*(-muh(m)/mbeam/radiusp/bfldp*psipp*dbdtp*grcgtp)
         pzdot = pzd0 + (qbeam/mbeam*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
             +qbeam/mbeam*(-xdot*delbyp+ydot*delbxp+zdot*dadzp))*ipara

         edot = xdot*exp1+(ydot-vp0)*eyp+zdot*ezp                      &
             +qbeam*pzdot*aparp*tor     &
             +qbeam*vpar*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)   &
             -qbeam*vpar*delbxp*vp0

         xh3(m) = xh2(m) + dt*xdot
         yh3(m) = yh2(m) + dt*ydot
         zh3(m) = zh2(m) + dt*zdot
         uh3(m) = uh2(m) + dt*pzdot

         dum = 1-wh3(m)*nonlin*0.
         if(isuni.eq.1)dum = exp(-vfac)*pidum
!         vxdum = eyp+vpar/b*delbxp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
         w3old = wh3(m)
         kap = 1/lh
         wh3(m)=wh2(m) + dt*(vxdum*kap + edot/teth)*dum
!         wh3(m)=wh2(m) + dt*(vxdum/lh  &
!                  + edot*3/2*sqvfac/(vfac*sqvfac+ecpow))*dum

         if(abs(wh3(m)).gt.1..and.nonlin==1)then
            wh3(m) = 0.
            wh2(m) = 0.
         end if


         laps=anint((zh3(m)/lz)-.5)*(1-peritr)
         r=xh3(m)-0.5*lx+lr0
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         qr = wx0*sf(i)+wx1*sf(i+1)
         yh3(m)=dmod(yh3(m)-laps*2*pi*qr*lr0/q0+800.*ly,ly)
         if(xh3(m)>lx)then
            xh3(m) = lx-1.e-8
            zh3(m)=lz-zh3(m)
            zh2(m)=zh3(m)
            xh2(m) = xh3(m)
!            wh2(m) = 0.
!            wh3(m) = 0.
         end if
         if(xh3(m)<0.)then
            xh3(m) = 1.e-8
            zh3(m)=lz-zh3(m)
            zh2(m)=zh3(m)
            xh2(m) = xh3(m)
!            wh2(m) = 0.
!            wh3(m) = 0.
         end if
         zh3(m)=dmod(zh3(m)+8.*lz,lz)
         xh3(m)=dmod(xh3(m)+8.*lx,lx)         
         xh3(m) = min(xh3(m),lx-1.0e-8)
         yh3(m) = min(yh3(m),ly-1.0e-8)
         zh3(m) = min(zh3(m),lz-1.0e-2)

!     particle diagnostics done here because info is available...
         mypfl=mypfl + w3old*(eyp+vpar*delbxp/b) 
         myefl=myefl + vfac*w3old*(eyp+vpar*delbxp/b)
         myke=myke + vfac*wh3(m)
         mynos=mynos + wh3(m)
         myavewh = myavewh+abs(wh3(m))

!     xn+1 becomes xn...
         uh2(m)=uh3(m)
         xh2(m)=xh3(m)
         yh2(m)=yh3(m)
         zh2(m)=zh3(m)
         wh2(m)=wh3(m)

!     100     continue
      enddo

      sbuf(1)=myke
      sbuf(2)=myefl
      sbuf(3)=mypfl
      sbuf(4)=mynos
      sbuf(5)=myavewh
      call MPI_ALLREDUCE(sbuf,rbuf,10,  &
          MPI_REAL8,MPI_SUM,           &
          MPI_COMM_WORLD,ierr)

      ketemp=rbuf(1)
      efltemp=rbuf(2)
      pfltemp=rbuf(3)
      nostemp=rbuf(4)
      avewh(n) = rbuf(5)/( float(tmm(1)) )
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      nos(3,n)=nostemp/( float(tmm(1)) )
      pfl(3,n)=pfltemp/( float(tmm(1)) )
      efl(3,n)=mims(ns)/tets(1)*efltemp/( float(tmm(1)) )
      ke(3,n)=ketemp/( 2.*float(tmm(1))*mims(ns) )
      np_old=mmb 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call init_pmove(zh3,np_old,lz,ierr)

      call pmove(xh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xh3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yh3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zh3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(uh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(uh3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wh3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(muh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xih,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0h,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(pzh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ekh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call end_pmove(ierr)
      mmb=np_new
!     write(*,*)MyId,mm(1)

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine loadh

      use gem_com
      use equil
      implicit none
      INTEGER :: i,j,k,m,idum,ns,m1
      INTEGER :: np_old,np_new
      REAL(8) :: vpar,vperp2,r,qr,th,b,cost,ter
      REAL(8) :: avgv,myavgv,avgw,myavgw
      real(8) :: dumx,dumy,dumz,jacp
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,psp
      REAL(8) :: grp,gxdgyp,zoldp
      REAL(8) :: wx0,wx1,wz0,wz1

      cnt=int(tmm(1)/numprocs)
!   write(*,*)'in loader cnt, mm(1)=',cnt,mm(1)

      myavgv=0.
      avgv=0.
      avgw = 0.
      myavgw = 0.

      mmb=int(tmm(1)/numprocs)
      m = 0
      do 160 j = 1,100000000

!     load a slab of ions...

         dumx=lx*ran2(iseed) !revers(MyId*cnt+j,2) !ran2(iseed)
         dumy=ly*ran2(iseed) !revers(MyId*cnt+j,3) !ran2(iseed)
         dumz=lz*ran2(iseed) !revers(MyId*cnt+j,5) !ran2(iseed)
         dumz = min(dumz,lz-1.e-8)
         r = lr0+dumx-0.5*lx
         th = (dumz-lz/2)/(q0*br0)
         i = int((r-rin)/dr)
         k = int((pi+th)/dth)
         jacp = jacob(i,k)
         if(ran2(iseed)<(0.5*jacp/jacmax))then
         m = m+1
         if(m>mmb)goto 170
         xh2(m)=min(dumx,lx-1.d-8)
         yh2(m)=min(dumy,ly-1.d-8)

         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1-wz0
         zh2(m) = wz0*zfnth(k)+wz1*zfnth(k+1)
         zh2(m)=min(zh2(m),lz-1.d-8)

         call parperp(vpar,vperp2,m,pi,cnt,MyId)
!   normalizations will be done in following loop...

         r=xh2(m)-0.5*lx+lr0
         cost=cos(th)
         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
         psp = wx0*psi(i)+wx1*psi(i+1)
         ter = teth !wx0*t0i(i)+wx1*t0i(i+1)
         b=1.-tor+tor*bfldp

         uh2(m)=vpar/sqrt(mbeam/ter)
         muh(m)=0.5*vperp2/b*ter
         ekh(m) = muh(m)*b+0.5*mbeam*uh2(m)**2
         pzh(m) = mbeam*uh2(m)/b-psp/br0
         z0h(m) = zh2(m)
         xih(m) = xh2(m)
         myavgv=myavgv+uh2(m)

!    LINEAR: perturb w(m) to get linear growth...
         wh2(m)=2.*amp*(revers(MyId*cnt+m,13)-0.5) !(ran2(iseed) - 0.5 )
         myavgw=myavgw+wh2(m)
         end if
 160  continue
 170  continue
      myavgw = myavgw/mmb
!      write(*,*)'myid ', myid,x2(10),y2(20),u2(20),w2(20)
!    subtract off avg. u...
      call MPI_ALLREDUCE(myavgv,avgv,1, &
          MPI_REAL8, &
          MPI_SUM,MPI_COMM_WORLD,ierr)
      if(idg.eq.1)write(*,*)'all reduce'
      avgv=avgv/float(tmm(1))
      do 180 m=1,mmb
         uh2(m)=uh2(m)-avgv
         xh3(m)=xh2(m)
         yh3(m)=yh2(m)
         zh3(m)=zh2(m)
         uh3(m)=uh2(m)
         wh3(m)=wh2(m)
 180  continue

      np_old=mmb
      call init_pmove(zh3,np_old,lz,ierr)
      if(idg.eq.1)write(*,*)'pass init_pmove'
!     
      call pmove(xh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xh3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yh3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zh3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(uh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(uh3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wh3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(muh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call pmove(xih,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0h,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(pzh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ekh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
!     
      call end_pmove(ierr)
      mmb=np_new
     write(*,*)MyId,j,mmb

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real(8) function en3(s)
      real(8) :: s

      if(0 .le. s .and. s .le. 1.0)then
         en3 = s*s/2
      else if(1.0 .le. s .and. s .le. 2.0) then
         en3 = -3./2.+3*s-s*s
      else if(2.0 .le. s .and. s .le. 3.0) then
         en3 = (3-s)*(3-s)/2.0
      else
         en3 = 0.
      end if

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine blendf
      use gem_com
      use equil
      implicit none

      INTEGER :: m,i,j,k,m1,m2
      complex(8) :: dum,dum1
      complex(8) :: tmp(nb,nb),work(100)
      integer :: IPIV(10),INFO
      REAL(8) :: r,qr,s1,s2,s3,dth1,wx0,wx1
      real(8) :: aky(0:jmx-1),dely(0:imx-1)

      dth1 = pi2/nb
      do 51 j = 0,im-1
         r = rin+xg(j)
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         qr = wx0*sf(i)+wx1*sf(i+1)
         dely(j) = dmod(-pi2*lr0/q0*qr+80.*ly,ly)*tor
 51   continue
      do 52 j = 0,jm-1
         if(j.ge.(jm/2+1)) then
            aky(j) = -2.*pi*float(jm-j)/ly
         else
            aky(j) = 2.*pi*float(j)/ly
         end if
 52   continue
      do 54 j = 0,jm-1
         do i = 0,im-1
            pfac(i,j) = exp(IU*aky(j)*dely(i))
         end do
 54   continue

      do m = 1,nb
         do i = 0,imx-1
            do j = 0,jmx-1
               do k = 0,kmx
                  s1 = float(k)*pi2/kmx/dth1-m
                  s2 = float(k)*pi2/kmx/dth1-m-nb
                  s3 = float(k)*pi2/kmx/dth1-m+nb
                  pol(m,i,j,k) = en3(s1)+en3(s2)*pfac(i,j)+en3(s3)/pfac(i,j)
               end do 
            end do 
         end do
      end do

      do i = 0,imx-1
         do j = 0,jmx-1
            do m1 = 1,nb
               do m2 = 1,nb
                  dum = 0.
                  do k = 0,km-1
                     dum = dum+dconjg(pol(m1,i,j,k))*pol(m2,i,j,k)
                  end do
                  tmp(m1,m2) = dum
                  pmtrx(i,j,m1,m2) = dum
               end do
            end do
            call F07ARF(nb,nb,tmp,nb,IPIV,INFO)
            call F07AWF(nb,tmp,nb,IPIV,work,100,INFO)
            do m1 = 1,nb
               do m2 = 1,nb
                  pmtrxi(i,j,m1,m2) = tmp(m1,m2)
               end do
            end do
         end do
      end do

      return
      if(myid==0)then
         i = 5
         j = 1
         do k = 0,kmx
            write(*,10)(pol(m,i,j,k),m=1,nb)
         end do

         write(*,*)'pfac= ', pfac(i,j)
         do m = 1,nb
            write(*,10)pol(m,i,j,km)-pol(m,i,j,0)*pfac(i,j)
         end do
 10   format(12(1x,e10.3))

         write(*,*)' '
         write(*,*)' '

         do m1=1,nb
            write(*,11)(pmtrx(i,j,m1,k),k=1,nb)
         end do
         write(*,*)' '
         write(*,*)' '
         do m1=1,nb
            write(*,11)(pmtrxi(i,j,m1,k),k=1,nb)
         end do

         do m1=1,nb
            do m2=1,nb
               dum = 0.
               do k = 1,nb
                  dum = dum+pmtrx(i,j,m1,k)*pmtrxi(i,j,k,m2)
               end do
               tmp(m1,m2) = dum
            end do
         end do

         write(*,*)' '
         write(*,*)' '
         do m1=1,nb
            write(*,11)(tmp(m1,k),k=1,nb)
         end do
      end if

 11   format(12(1x,e10.3))

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine filtbl(u)   
      use gem_com
      use equil
      use fft_wrapper
      implicit none
      complex(8) :: lbfs(0:imx-1,0:jmx-1)
      complex(8) :: rbfr(0:imx-1,0:jmx-1)
      complex(8) :: u(0:imx-1,0:jmx-1,0:1)
      complex(8) :: bm(1:nb),cm(1:nb),dum
      REAL(8) :: qr
      real(8) :: aky(0:jmx-1),dely(0:imx-1),dklz(0:imx-1,0:jmx-1)
      INTEGER :: i,j,k,l,m,n,ind,m1,m2
      INTEGER :: myk,mynum,id
      complex(8),dimension(:),allocatable :: holdu,rbuf,sbuf

!      write(*,*)'enter filtor'
      mynum = imx*jmx/numprocs      
      allocate(holdu(0:imx*jmx*kmx/numprocs-1),  &
               rbuf(0:mynum-1), &
               sbuf(0:mynum-1))

!      return
!pack data to send

       do id = 0,last
          if(id.ne.myid)then
             do ind = 0,mynum-1 !id*mynum,(id+1)*mynum
                j = (id*mynum+ind)/im
                i = id*mynum+ind-j*im
                sbuf(ind) = u(i,j,0)
             end do
!             send sbuf to id
               call MPI_SENDRECV(sbuf(0),mynum, &
                         mpi_complex16,id,10, &
                         rbuf(0),mynum, &
                         mpi_complex16,id,10, &
                         mpi_comm_world,stat,ierr)
!unpack and put into holdu
             do ind = 0,mynum-1
                holdu(ind*km+int(id/ntube)) = rbuf(ind)
             end do
          end if
       end do

!put own u(:,:,) in holdu
       do ind = 0,mynum-1
          j = (myid*mynum+ind)/im
          i = myid*mynum+ind-j*im
          holdu(ind*km+int(myid/ntube)) = u(i,j,0)
       end do

       do ind = 0,mynum-1
          j = (myid*mynum+ind)/im
          i = myid*mynum+ind-j*im
          do k = 0,km-1
             tmpz(k) = holdu(ind*km+k)
          end do

!do blending function filtering
          do m1 = 1,nb
             dum = 0.
             do k = 0,km-1
                dum = dum+dconjg(pol(m1,i,j,k))*tmpz(k)
             end do
             bm(m1) = dum
          end do

          do m1 = 1,nb
             dum = 0.
             do m2 = 1,nb
                dum = dum+pmtrxi(i,j,m1,m2)*bm(m2)
             end do
             cm(m1) = dum
          end do
          
          do k = 0,km-1
             dum = 0.
             do m1 = 1,nb
                dum = dum+cm(m1)*pol(m1,i,j,k)
             end do
             tmpz(k) = dum
          end do

          do k = 0,km-1
             holdu(ind*km+k) = tmpz(k)  
          end do
       end do

       do id = 0,last
          if(id.ne.myid)then
             do ind = 0,mynum-1
                sbuf(ind) = holdu(ind*km+int(id/ntube))
             end do
!             send sbuf to id
               call MPI_SENDRECV(sbuf(0),mynum, &
                         mpi_complex16,id,20, &
                         rbuf(0),mynum, &
                         mpi_complex16,id,20, &
                         mpi_comm_world,stat,ierr)
!unpack and put into u
             do ind = 0,mynum-1
                j = (id*mynum+ind)/im
                i = id*mynum+ind-j*im
                u(i,j,0) = rbuf(ind)
             end do
          end if
       end do
                
!put own holdu in u
       do ind = 0,mynum-1
          j = (myid*mynum+ind)/im
          i = myid*mynum+ind-j*im
          u(i,j,0) = holdu(ind*km+int(myid/ntube))
       end do

!       write(*,*)'before assign',myid
!assign u(:,:,1)

       do i = 0,im-1
          do j = 0,jm-1
             lbfs(i,j) = u(i,j,0)
          end do
       end do

       call MPI_SENDRECV(lbfs(0:imx-1,0:jmx-1),imx*jmx, &
                         MPI_complex16,lngbr,10, &
                         rbfr(0:imx-1,0:jmx-1),imx*jmx, &
                         MPI_complex16,rngbr,10, &
                         TUBE_COMM,stat,ierr)
!       write(*,*)'after send_recv',myid
       if(gclr.ne.glst)then                  
          do i = 0,im-1
             do j = 0,jm-1
                u(i,j,1) = rbfr(i,j)
             end do
          end do
       end if
       if(gclr==glst)then
          do i = 0,im-1
             do j = 0,jm-1
                u(i,j,1) = rbfr(i,j)*pfac(i,j)
             end do
          end do
       end if
          
!       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       return
       end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ldel

      use gem_com
      use equil
      implicit none
      INTEGER :: i,j,k,m,idum,ns,m1
      INTEGER :: np_old,np_new
      REAL(8) :: vpar,vperp2,r,qr,th,b,cost,ter
      REAL(8) :: avgv,myavgv,avgw,myavgw
      real(8) :: dumx,dumy,dumz,jacp
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,psp
      REAL(8) :: grp,gxdgyp,zoldp
      REAL(8) :: wx0,wx1,wz0,wz1

      cnt=int(tmm(1)/numprocs)
!   write(*,*)'in loader cnt, mm(1)=',cnt,mm(1)

      myavgv=0.
      avgv=0.
      avgw = 0.
      myavgw = 0.

      m = 0
      do 160 j = 1,100000000

!     load a slab of electrons according to J*xn0e...

         dumx=lx*ran2(iseed) !revers(MyId*cnt+j,2) !ran2(iseed)
         dumy=ly*ran2(iseed) !revers(MyId*cnt+j,3) !ran2(iseed)
         dumz=lz*ran2(iseed) !revers(MyId*cnt+j,5) !ran2(iseed)
         dumz = min(dumz,lz-1.e-8)
         r = lr0+dumx-0.5*lx
         th = (dumz-lz/2)/(q0*br0)
         i = int((r-rin)/dr)
         k = int((pi+th)/dth)
         jacp = jacob(i,k)
         if(ran2(iseed)<(0.5*jacp/jacmax))then
         m = m+1
         if(m>mme)goto 170
         x2e(m)=min(dumx,lx-1.d-8)
         y2e(m)=min(dumy,ly-1.d-8)

         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1-wz0
         z2e(m) = wz0*zfnth(k)+wz1*zfnth(k+1)
         z2e(m)=min(z2e(m),lz-1.d-8)

         call parperp(vpar,vperp2,m,pi,cnt,MyId)
!   normalizations will be done in following loop...

         r=x2e(m)-0.5*lx+lr0
         cost=cos(th)
         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
         psp = wx0*psi(i)+wx1*psi(i+1)
         ter = wx0*t0e(i)+wx1*t0e(i+1)
         b=1.-tor+tor*bfldp

         u2e(m)=vpar*sqrt(amie*ter)
         mue(m)=0.5*vperp2/b*ter
         eke(m) = mue(m)*b+0.5*emass*u2e(m)**2
         pze(m) = emass*u2e(m)/b+psp/br0
         z0e(m) = z2e(m)
         xie(m) = x2e(m)
         myavgv=myavgv+u2e(m)

!    LINEAR: perturb w(m) to get linear growth...
         w2e(m)=2.*amp*(revers(MyId*cnt+m,13)-0.5) !(ran2(iseed) - 0.5 )
         myavgw=myavgw+w2e(m)
         end if
 160  continue
 170  continue
      myavgw = myavgw/mme
!      write(*,*)'myid ', myid,x2(10),y2(20),u2(20),w2(20)
!    subtract off avg. u...
      call MPI_ALLREDUCE(myavgv,avgv,1, &
          MPI_REAL8, &
          MPI_SUM,MPI_COMM_WORLD,ierr)
      if(idg.eq.1)write(*,*)'all reduce'
      avgv=avgv/float(tmm(1))
      do 180 m=1,mme
         u2e(m)=u2e(m)-avgv
         x3e(m)=x2e(m)
         y3e(m)=y2e(m)
         z3e(m)=z2e(m)
         u3e(m)=u2e(m)
         w3e(m)=w2e(m)
         mue2(m) = mue(m)
         mue3(m) = mue(m)
 180  continue

      np_old=mme
      call init_pmove(z2e,np_old,lz,ierr)
      if(idg.eq.1)write(*,*)'pass init_pmove'
!     
      call pmove(x2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mue2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mue3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call pmove(mue,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xie,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(pze,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(eke,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
!     
      call end_pmove(ierr)
      mme=np_new

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
