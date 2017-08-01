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
        if(iCrs_Sec==2)then
           call constrdte
           goto 100
        end if
        if(isft==1)then
           call ftcamp
           goto 100
        end if
        starttm=MPI_WTIME()

        do  timestep=ncurr,nm
           tcurr = tcurr+dt

	   call accumulate(timestep-1,0)
	   call ampere(timestep-1,0) !need upar prior to drdt()
           if(isgkm==1)then
              call jie(timestep-1)
              call drdt(timestep-1)
              call dpdt(timestep-1)
           end if
	   if(isgkm==0)call poisson(timestep-1,0)
	   call field(timestep-1,0)
	   call diagnose(timestep-1)
           call reporter(timestep-1)

	   call push_wrapper(timestep,1)

	   call accumulate(timestep,1)
	   call ampere(timestep,1)     
           if(isgkm==1)then
              call jie(timestep)
              call drdt(timestep)
              call dpdt(timestep)
           end if
	   if(isgkm==0)call poisson(timestep,1)
	   call field(timestep,1)

	   call push_wrapper(timestep,0)
           if(mod(timestep,10000)==0)then
              do i=0,last 
!                 if(myid==i)write(*,*)myid,mm(3)
                 call MPI_BARRIER(MPI_COMM_WORLD,ierr)
              end do
           end if
         end do
         call sweep
	 lasttm=MPI_WTIME()
	 tottm=lasttm-starttm
!	 write(*,*)'ps time=',pstm,'tot time=',tottm
!         do i=0,last 
!            if(myid==i)write(*,*)myid,mm(3)
!            call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!         end do
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
      INTEGER :: i,j,k,n,ns,idum,i1,k1,m,j1
      INTEGER :: mm1,lr1
      REAL(8) :: mims1,tets1,q1,kappan,kappat,r,qr,th,cost,dum,zdum
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,upae0p,dnuobdrp,dnuobdtp
      REAL(8) :: grp,gxdgyp,g2xp,g2yp,jacp,jfnp,gn0ip,gn0ep,gt0ip,gt0ep,capnpp,captep
      REAL(8) :: wx0,wx1,wz0,wz1,b

      IU=cmplx(0.,1.)
      pi=4.0*atan(1.0)
      pi2 = pi*2.

      open(115,file='gem.in')
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) imx,jmx,kmx,mmx,mmxe,nmx,nsmx,modemx,ntube,icrs_sec,ipg,isphi,itrace,isft,irlk,qoffset
      read(115,*) dumchar
      read(115,*) lxa,lyfra,nrst,eprs,epsnx,epsny,epsnz,epsax,epsay,epsaz
      read(115,*) dumchar
      read(115,*) dt,nm,nsm,xshape,yshape,zshape,isham,ishgk,icmprs,iez,isgkm,icncl,ispre,isdte
      read(115,*) dumchar
      read(115,*) kymode,iput,iget,ision,ish,peritr,llk,mlk,onemd,izon,kxz,oml,omu,nom,delom
      read(115,*) dumchar
      read(115,*) nplot,xnplt,modem,nzcrt,npzh,nenh,dtref,islow,mynf,frmax,ifskp,jcnt
      read(115,*) dumchar
      read(115,*) cut,amp,tor,ishift,fradi,kxcut,kycut,bcut,epspx,epspz
      read(115,*) dumchar
      read(115,*) dum,r0a,dum,dum,width,vpp,vt0,yd0,etaohm
      read(115,*) dumchar
      read(115,*) c1,c2,c3,c4,ifluid,isg,amie,rneu,gamion,gamele,gambm,gamdne,gamphi,gamapa,gamte,nugam
      read(115,*) dumchar
      read(115,*) beta,nonlin,nonline,nonlinh,ipara,iparah,nonlint,iparat,vwidth,vwidthe,vcut,isuni,idg
!     begin reading species info, ns=1,nsm...
      if(nsm.le.0) write(*,*)'invalid nsm',nsm
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) mm1,tmme,tets1,mims1,q1,lr1
      read(115,*) dumchar
      read(115,*) kappan,kapti,kapte
      read(115,*) dumchar
      read(115,*) betah,lh,delth,vc,vi,vmin,iflr,iorb,iflrh,iske,ifrzt,idenek,iflut,ieqmo814
      isbam = 0
      isbgk = 0
      iscam = 1
      iscgk = 1
      im=imx
      jm=jmx
      km=kmx

!set up omlk(0:nom-1)
      iomlk = 0.
      do i = 0,nom-1
         omlk(i) = -omu+2*omu/nom*i
         if(abs(omlk(i)-oml)<delom .or. abs(omlk(i)+oml)<delom) iomlk(i)=1.
      end do
      nfreq = kmx*mynf
      call new_gem_com()
      ns = 1
      ptr(1) = 1
      tmm(ns)=mm1
      mm(ns)=int(mm1/numprocs)
      mme = int(tmme/numprocs)
      mtrace=tmm(1) !/8
      mmt = mtrace/numprocs
      if (MyId.eq.Last) mm(ns)=mm1-Last*mm(ns)
!     write(*,*)'in init  ',Myid,mm(ns)
      tets(ns)=tets1
      mims(ns)=mims1
      q(ns)=q1
      lr(ns)=lr1
      kapn(ns)=kappan
      kapt(ns)=kappat
      tets(2) = 1.
      mims(2) = mims(1)
      lr(2) = 1
      q(2) = 1.
      kapn(2) = kapn(1)

      mims(3) = 2.
      q(3) = 1

      mims(4) = 2
      q(4) = 1

      emass = 1./amie
      qel = -1

      rina = r0a-lxa/2
      routa = r0a+lxa/2
      call new_equil(myid)
      lx = lxa*a
      ly = 2.*pi*r0/q0abs/lyfra
      br0 = rmaj0
      lr0 = r0
!      qp = q0p                !useless
      lz = pi2*q0abs*rmaj0
      delz = lz/ntheta
      tload = t0i(10)
      tloadc = t0c(10)

      pzcrit(3) = q(3)*abs(psi(nr)-psi(0))/br0/npzh
      encrit(3) = 0.5*mims(3)*vc**2/nenh

!      write(*,*)'br0,lr0,q0,qp = ', br0,lr0,q0,qp
      if(iget.eq.0)then
         if(myid.eq.master)open(9, file='plot', status='unknown')
         if(myid.eq.master)open(11, file='flux', status='unknown')
         if(myid.eq.master)open(15, file='freq', status='unknown')
         if(gclr==kmx/2 .and. tclr==0)open(16, file='indicator', status='unknown')
         if(gclr==kmx/2 .and. tclr==0)open(17, file='yyre', status='unknown')
         if(myid.eq.master)open(18, file='yyre1', status='unknown')
         if(gclr==kmx/2 .and. tclr==0)open(22, file='yyre2', status='unknown')
         if(myid.eq.master)open(21, file='sweep', status='unknown')
         if(gclr==kmx/2 .and. tclr==0)open(23, file='mdhis', status='unknown')
         if(gclr==kmx/2 .and. tclr==0)open(24, file='mdhisa', status='unknown')
         if(gclr==kmx/2 .and. tclr==0)open(25, file='stress', status='unknown')
      end if
      if(iget.eq.1)then
         if(myid.eq.master)open(9, file='plot', & 
             status='unknown',position='append')
         if(myid.eq.master)open(11, file='flux', &
             status='unknown',position='append')
         if(myid.eq.master)open(15, file='freq', &
             status='unknown',position='append')
         if(gclr==kmx/2 .and. tclr==0)open(16, file='indicator', &
             status='unknown',position='rewind')
         if(gclr==kmx/2 .and. tclr==0)open(17, file='yyre', &
             status='unknown',position='append')
         if(myid.eq.master)open(18, file='yyre1', &
             status='unknown',position='append')
         if(gclr==kmx/2 .and. tclr==0)open(22, file='yyre2', &
             status='unknown',position='append')
         if(myid.eq.master)open(21, file='sweep', &
             status='unknown',position='rewind')
         if(gclr==kmx/2 .and. tclr==0)open(23, file='mdhis', status='unknown',position='append')
         if(gclr==kmx/2 .and. tclr==0)open(24, file='mdhisa', status='unknown',position='append')
         if(gclr==kmx/2 .and. tclr==0)open(25, file='stress', status='unknown',position='append')
      end if

      if(isuni.eq.0)vwidthe=vwidth
      dte = dt
      iadi = 0
      if(isg.gt.0.)fradi = isg
      if(ifluid.eq.0)then
         iadi = 1
         fradi = 1.
      end if



      if(iget.eq.1) amp=0.
!     totvol is the square for now...
      dx=lx/float(im)
      dy=ly/float(jm)
      dz=lz/float(km)
!      totvol=lx*ly*lz

!      e0=lr0/q0abs/br0   !useless
!     
      do 10 i=0,nxpp
         xg(i)=i*dx !dx*(tclr*nxpp+i)
 10   continue
      do 12 j=0,jm
         yg(j)=dy*float(j)
 12   continue
      kcnt=1

!assign mode sign and its index in fft calls
!      jcnt = 3  !jmx/ntube
      mstart = 0
      ntor0 = mstart+1
      do m = 0,jcnt-1
         if(m==0)ntor(m) = 0
         if(m>0)ntor(m) = ntor0+int((m-1)/2)

         isgnft(m) = 1
         j1 = mstart+int((float(m)+1.0)/2)
         jft(m) = j1
         if(m==0)then
            isgnft(m) = 1
            jft(m) = 0
         end if
         if(m>0.and.mod(m,2)==0)then
            isgnft(m) = -1
            jft(m) = jmx-j1
         end if
      end do

      do 14 k=0,mykm
         n=GCLR*kcnt+k	
         zg(k)=dz*float(n)
 14   continue

!     set ring weighting array, lr is the no. of points...
      do 300 ns=1,nsm
         if (lr(ns).eq.1) then
            rwx(ns,1)=0.
            rwy(ns,1)=0.
         elseif(lr(ns).eq.2) then
            rwx(ns,1)=1.
            rwx(ns,2)=-1.
            rwy(ns,1)=0.
            rwy(ns,2)=0.
         elseif(lr(ns).eq.4) then
            rwx(ns,1)=1.
            rwx(ns,2)=-1.
            rwx(ns,3)=0.
            rwx(ns,4)=0.
            rwy(ns,1)=0.
            rwy(ns,2)=0.
            rwy(ns,3)=1.
            rwy(ns,4)=-1.
         else
            write(*,*)'incorrect value of lr',lr(ns)
         endif
 300  continue
          
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
      if(q0<0)zfnth = zfnth+lz

      thfnz(0) = -pi
      thfnz(ntheta) = pi
      if(q0<0.)then
         thfnz(0) = pi
         thfnz(ntheta) = -pi
      end if

      if(q0>0)then
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
 127        continue
         end do
      end if

      if(q0<0)then
         k = 0
         do j = 1,ntheta-1
            zdum = lz-j*lz/ntheta
            do i = k,ntheta-1
               if(zfnth(i)>=zdum.and.zfnth(i+1)<zdum)then
                  k = i
                  dum = (zdum-zfnth(i))*dth/(zfnth(i+1)-zfnth(i))
                  thfnz(ntheta-j) = i*dth-pi+dum
                  go to 128
               end if
            end do
 128        continue
         end do
      end if

      if(myid.eq.master)then
         open(19, file='eqdat', status='unknown')
         do i = 0,nr
            write(19, 100)i, psi(i)
         end do
 100     format(1x, i4, 5(2x,e10.3))

         do j = 0,ntheta
            write(19,99)j,bdcrvb(nr/2,j),thfnz(j),zfnth(j)
         end do
         close(19)
      end if
 99      format(1x,i3,2x,8(1x,e10.3))

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

         dnuobdrp = wx0*wz0*dnuobdr(i,k)+wx0*wz1*dnuobdr(i,k+1) &
                 +wx1*wz0*dnuobdr(i+1,k)+wx1*wz1*dnuobdr(i+1,k+1)
         dnuobdtp = wx0*wz0*dnuobdt(i,k)+wx0*wz1*dnuobdt(i,k+1) &
                 +wx1*wz0*dnuobdt(i+1,k)+wx1*wz1*dnuobdt(i+1,k+1)
         upae0p = wx0*wz0*upae0(i,k)+wx0*wz1*upae0(i,k+1) &
                 +wx1*wz0*upae0(i+1,k)+wx1*wz1*upae0(i+1,k+1)

         fp = wx0*f(i)+wx1*f(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         gt0ip = wx0*t0i(i)+wx1*t0i(i+1)
         gt0ep = wx0*t0e(i)+wx1*t0e(i+1)
         gn0ep = wx0*xn0e(i)+wx1*xn0e(i+1)
         gn0ip = wx0*xn0i(i)+wx1*xn0i(i+1)
         capnpp = wx0*capne(i)+wx1*capne(i+1)
         captep = wx0*capte(i)+wx1*capte(i+1)
         b=1.-tor+tor*bfldp
         cfx(i1,k1) = br0/b**3*fp/radiusp*dbdtp*grcgtp
         cfy(i1,k1) = br0/b**3*fp/radiusp* &
                      (dydrp*dbdtp-lr0/q0*qhatp*dbdrp)*grcgtp
         bdgxcgy(i1,k1) = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp             
         ggr(i1,k1) = grp
         bmag(i1,k1) = b
         jac(i1,k1) = jacp*jfnp
         bdgrzn(i1,k1) = q0*br0/radiusp/b*psipp*grcgtp/jfnp
         gt0i(i1) = gt0ip
         gt0e(i1) = gt0ep
         gn0e(i1) = gn0ep
         gn0i(i1) = gn0ip
         gcpne(i1) = capnpp
         gcpte(i1) = captep
         gupae0(i1,k1) = upae0p
         gnuoby(i1,k1) = (-dydrp*dnuobdtp+r0/q0*qhatp*dnuobdrp)*fp/radiusp*grcgtp
         gnuobx(i1,k1) = dnuobdtp*fp/radiusp*grcgtp
         end do
      end do

!set up dfltz for digital filtering
      do k = 0,1
         n = gclr+k
         dfltz(k) = float(abs(n-kmx/2))/float(kmx/2)
      end do
      dfltz = dfltz*dt/dtref

      iseed = -(1777+myid*13)
      idum = ran2(iseed)
      phi = 0.
      apar = 0.
      dene = 0.
      upar = 0.
      camp = 0.
      phiom = 0.
      dteom = 0.
      delte = 0.

      do i = 0,im
         do j = 0,jm
            do k = 0,mykm
               phi(i,j,k) = amp*(ran2(idum)-0.5)*ifluid*1e-8  !amp*sin(nzcrt*xg(i)*pi/lx) !
               dene(i,j,k) = amp*(ran2(idum)-0.5)*ifluid *1e-8
               apar(i,j,k) = amp*(ran2(idum)-0.5)*ifluid *1e-10 
            end do
         end do
      end do

      if(myid.eq.master)then
         write(9,*)'dt,beta= ',dt, beta
         write(9,*)'kapt,kapn,kpte= ', kapti, kappan,kapte
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
         write(9,*)'nh,lh,psimax = ',nh,lh,psimax
         write(9,*)'vc,vi,vmin = ',vc,vi,vmin
         write(9,*)'a,rmaj0,lx,ly,lz = ',a,rmaj0,lx,ly,lz
      end if

      if(myid.eq.master)then
         write(*,*)'dt,beta= ',dt, beta
         write(*,*)'kapt,kapn,kpte= ', kapti, kappan,kapte
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
         write(*,*)'nh,lh,psimax = ',nh,lh,psimax
         write(*,*)'vc,vi,vmin = ',vc,vi,vmin
         write(*,*)'a,rmaj0,lx,ly,lz = ',a,rmaj0,lx,ly,lz,erf(0.7)
      end if
      close(115)
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ppush(n)

      use gem_com
      use equil
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum,dum1
      INTEGER :: m,i,j,k,l,n
      INTEGER :: np_old,np_new
      REAL(8) :: rhog,vfac,kap,vpar,pidum,kaptp,kapnp,xnp
      REAL(8) :: b,th,r,enerb,cost,sint,qr,laps,sz,ter
      REAL(8) :: xt,xs,yt,xdot,ydot,zdot,xdt,ydt,pzdot,edot,pzd0,vp0
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
      REAL(8) :: grp,gxdgyp,rhox(4),rhoy(4),psp,pzp

      pidum = 1./(pi*2)**1.5*vwidth**3
      do m=1,mm(1)
         r=x2(m)-0.5*lx+lr0
         k = int(z2(m)/delz)
         wz0 = ((k+1)*delz-z2(m))/delz
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
         kaptp = wx0*capti(i)+wx1*capti(i+1)
         kapnp = wx0*capni(i)+wx1*capni(i+1)
         xnp = wx0*xn0i(i)+wx1*xn0i(i+1)

         b=1.-tor+tor*bfldp
         pzp = mims(1)*u2(m)/b*fp/br0-q(1)*psp/br0

         rhog=sqrt(2.*b*mu(m)*mims(1))/(q(1)*b)*iflr

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
            xs=x2(m)+rhox(l) !rwx(1,l)*rhog
            yt=y2(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
!
!   particle can go out of bounds during gyroavg...
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            include "ppushli.h"
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
         vfac = 0.5*(mims(1)*u2(m)**2 + 2.*mu(m)*b)
         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))
         kap = kapnp - (1.5-vfac/ter)*kaptp

         vpar = u2(m)-q(1)/mims(1)*aparp*nonlin*0.
         enerb=(mu(m)+mims(1)*vpar*vpar/b)/q(1)*tor
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
         xdot = vxdum*nonlin -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         xdt = -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlin &
             +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0
         ydt = iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0
         zdot =  vpar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         pzd0 = tor*(-mu(m)/mims(1)/radiusp/bfldp*psipp*dbdtp*grcgtp)
         pzdot = pzd0 + (q(1)/mims(1)*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
             +q(1)/mims(1)*(-xdot*delbyp+ydot*delbxp+zdot*dadzp))*ipara

         edot = q(1)*(xdt*exp1+ydt*eyp+vpar*ezp)

         x3(m) = x2(m) + 0.5*dt*xdot
         y3(m) = y2(m) + 0.5*dt*ydot
         z3(m) = z2(m) + 0.5*dt*zdot
         u3(m) = u2(m) + 0.5*dt*pzdot

         dum = 1-w2(m)*nonlin*0.
         if(isuni.eq.1)dum = exp(-vfac)*pidum
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
!         vxdum = eyp+vpar/b*delbxp
         w3(m)=w2(m) + 0.5*dt*(vxdum*kap + edot/ter)*dum*xnp*(tload/ter)**1.5*exp(vfac*(1/tload-1./ter))-0.5*dt*gamion*w2(m)
         
!         if(x3(m)>lx .or. x3(m)<0.)w3(m) = 0.

!         go to 333
         if(abs(pzp-pzi(m))>3.0.or.abs(vfac-eki(m))>0.2*eki(m))then
            x3(m) = xii(m)
            z3(m) = z0i(m)
            r = x3(m)-lx/2+lr0
            k = int(z3(m)/delz)
            wz0 = ((k+1)*delz-z3(m))/delz
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
            u3(m) = sqrt(2/mims(1)*abs(eki(m)-mu(m)*b))
            u2(m) = u3(m)
            w3(m) = 0.
            w2(m) = 0.
            x2(m) = x3(m)
            z2(m) = z3(m)
         end if

 333     continue
         laps=anint((z3(m)/lz)-.5)*(1-peritr)
         r=x3(m)-0.5*lx+lr0
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         i = max(i,0)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         qr = wx0*sf(i)+wx1*sf(i+1)
         y3(m)=dmod(y3(m)-laps*2*pi*qr*lr0/q0*dsign(1.d0,q0)+800.*ly,ly)
         if(x3(m)>lx)then
            x3(m) = lx-1.e-8
            z3(m)=lz-z3(m)
            x2(m) = x3(m)
            z2(m) = z3(m)
            w2(m) = 0.
            w3(m) = 0.
         end if
         if(x3(m)<0.)then
            x3(m) = 1.e-8
            z3(m)=lz-z3(m)
            x2(m) = x3(m)
            z2(m) = z3(m)
            w2(m) = 0.
            w3(m) = 0.
         end if
         z3(m)=dmod(z3(m)+8.*lz,lz)
         x3(m)=dmod(x3(m)+8.*lx,lx)         
         x3(m) = min(x3(m),lx-1.0e-8)
         y3(m) = min(y3(m),ly-1.0e-8)
         z3(m) = min(z3(m),lz-1.0e-2)
         
      enddo

      np_old=mm(1)
      call init_pmove(z3,np_old,lz,ierr)

      call pmove(x2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mu,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call pmove(xii,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0i,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(pzi,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(eki,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call end_pmove(ierr)
      mm(1)=np_new

      return
      end

!-----------------------------------------------------------------------

      subroutine cpush(n)

      use gem_com
      use equil
      implicit none
      INTEGER :: n
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,w3old
      INTEGER :: m,i,j,k,ns,l
      INTEGER :: np_old,np_new
      REAL(8) :: vfac,vpar,vxdum,dum,xdot,ydot,xdt,ydt,zdot,pzdot,edot,pzd0,vp0
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

      do m=1,mm(1)
         r=x3(m)-0.5*lx+lr0

         k = int(z3(m)/delz)
         wz0 = ((k+1)*delz-z3(m))/delz
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
         kaptp = wx0*capti(i)+wx1*capti(i+1)
         kapnp = wx0*capni(i)+wx1*capni(i+1)
         xnp = wx0*xn0i(i)+wx1*xn0i(i+1)
         b=1.-tor+tor*bfldp
         pzp = mims(1)*u3(m)/b*fp/br0-q(1)*psp/br0

         rhog=sqrt(2.*b*mu(m)*mims(1))/(q(1)*b)*iflr

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
            xs=x3(m)+rhox(l) !rwx(1,l)*rhog
            yt=y3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
!   BOUNDARY
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            include "cpushli.h"
 200     continue

         exp1=exp1/4.
         eyp=eyp/4.
         ezp=ezp/4.
         delbxp=delbxp/4.
         delbyp=delbyp/4.
         dpdzp = dpdzp/4.
         dadzp = dadzp/4.
         aparp = aparp/4.


         vfac = 0.5*(mims(1)*u3(m)**2 + 2.*mu(m)*b)
         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))
         kap = kapnp -(1.5-vfac/ter)*kaptp

         vpar = u3(m)-q(1)/mims(1)*aparp*nonlin*0.
         enerb=(mu(m)+mims(1)*vpar*vpar/b)/q(1)*tor
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
         xdot = vxdum*nonlin -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         xdt = -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlin     &
             +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0
         ydt = iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0
         zdot =  vpar*(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         pzd0 = tor*(-mu(m)/mims(1)/radiusp/bfldp*psipp*dbdtp*grcgtp)
         pzdot = pzd0 + (q(1)/mims(1)*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
             +q(1)/mims(1)*(-xdot*delbyp+ydot*delbxp+zdot*dadzp))*ipara

         edot = q(1)*(xdt*exp1+ydt*eyp+vpar*ezp)

         x3(m) = x2(m) + dt*xdot
         y3(m) = y2(m) + dt*ydot
         z3(m) = z2(m) + dt*zdot
         u3(m) = u2(m) + dt*pzdot

         dum = 1-w3(m)*nonlin*0.
         if(isuni.eq.1)dum = exp(-vfac)*pidum
!         vxdum = eyp+vpar/b*delbxp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
         w3old = w3(m)
         w3(m) = w2(m) + dt*(vxdum*kap+edot/ter)*dum*xnp*(tload/ter)**1.5*exp(vfac*(1/tload-1./ter))-dt*gamion*w2(m)

         if(abs(w3(m)).gt.1..and.nonlin==1)then
            w3(m) = 0.
            w2(m) = 0.
         end if


         laps=anint((z3(m)/lz)-.5)*(1-peritr)
         r=x3(m)-0.5*lx+lr0
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         qr = wx0*sf(i)+wx1*sf(i+1)
         y3(m)=dmod(y3(m)-laps*2*pi*qr*lr0/q0*dsign(1.d0,q0)+800.*ly,ly)
         if(x3(m)>lx)then
            x3(m) = lx-1.e-8
            z3(m)=lz-z3(m)
            x2(m) = x3(m)
            z2(m) = z3(m)
            w2(m) = 0.
            w3(m) = 0.
         end if
         if(x3(m)<0.)then
            x3(m) = 1.e-8
            z3(m)=lz-z3(m)
            x2(m) = x3(m)
            z2(m) = z3(m)
            w2(m) = 0.
            w3(m) = 0.
         end if
         z3(m)=dmod(z3(m)+8.*lz,lz)
         x3(m)=dmod(x3(m)+8.*lx,lx)         
         x3(m) = min(x3(m),lx-1.0e-8)
         y3(m) = min(y3(m),ly-1.0e-8)
         z3(m) = min(z3(m),lz-1.0e-2)

!     particle diagnostics done here because info is available...
         mypfl=mypfl + w3old*(eyp+vpar*delbxp/b) 
         myefl=myefl + vfac*w3old*(eyp+vpar*delbxp/b)
         myke=myke + vfac*w3(m)
         mynos=mynos + w3(m)
         myavewi = myavewi+abs(w3(m))

!     xn+1 becomes xn...
         u2(m)=u3(m)
         x2(m)=x3(m)
         y2(m)=y3(m)
         z2(m)=z3(m)
         w2(m)=w3(m)

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
      efl(1,n)=efltemp/( float(tmm(1)) )
      ke(1,n)=ketemp/(float(tmm(1)))
      np_old=mm(1) 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call init_pmove(z3,np_old,lz,ierr)
!     
      call pmove(x2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mu,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xii,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0i,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(pzi,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(eki,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call end_pmove(ierr)
      mm(1)=np_new
!     write(*,*)MyId,mm(1)

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
      real(8) :: mydbr(0:imx-1),myjac(0:imx-1),v(0:imx-1),dum,dum1

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

!compute deltbr
      do i = 0,nxpp-1
         dum = 0.
         dum1 = 0.
         do j = 0,jm-1
            dum = dum+(delbx(i,j,0)*bdgxcgy(i,0)/ggr(i,0))**2*jac(i,0)
            dum1 = dum1+jac(i,0)
         end do
         mydbr(i) = dum
         myjac(i) = dum1
      end do
      call MPI_ALLREDUCE(mydbr,dbr,nxpp,  &
          MPI_REAL8,MPI_SUM,           &
          tube_comm,ierr)
      call MPI_ALLREDUCE(myjac,v,nxpp,  &
          MPI_REAL8,MPI_SUM,           &
          tube_comm,ierr)
      dbr(0:imx-1) = sqrt(dbr(0:imx-1)/v(0:imx-1))


      do k = 0,1
         do i = 0,imx
            do j = 0,jmx
               if(iske==0)tmp(i,j,k) = (gt0e(i)*dene(i,j,k)+gn0e(i)*delte(i,j,k)*isdte)/gt0e(i)
               if(iske==1)tmp(i,j,k) = (dltpe(i,j,k)+dltpa(i,j,k))/2/gt0e(i)
            end do
         end do
      end do

      call gradu(tmp(:,:,:),ux,uy)
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
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp,nhp,nbp
      REAL(8) :: enerb,vxdum,dum,xdot,ydot
      INTEGER :: m,n,i,j,k,l,ns,ip
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,vte,ter
      REAL(8) :: sz,wght,wght0,wght1,wght2,wght3,wght4,r,th,cost,sint,b,qr,dv
      REAL(8) :: xt,yt,rhog,pidum,vpar,xs,dely,vfac
      REAL(8) :: lbfs(0:imx,0:jmx)
      REAL(8) :: rbfs(0:imx,0:jmx)
      REAL(8) :: lbfr(0:imx,0:jmx)
      REAL(8) :: rbfr(0:imx,0:jmx)
      real(8) :: myden(0:imx,0:jmx,0:1),myjpar(0:imx,0:jmx,0:1)
      real(8) :: myhpar(0:imx,0:jmx,0:1),myhden(0:imx,0:jmx,0:1),myhden0(0:imx,0:jmx,0:1)
      real(8) :: mybpar(0:imx,0:jmx,0:1),mybden(0:imx,0:jmx,0:1)
      real(8) :: mydene(0:imx,0:jmx,0:1),myupar(0:imx,0:jmx,0:1)
      real(8) :: mydltpa(0:imx,0:jmx,0:1),mydltpe(0:imx,0:jmx,0:1)
      real(8) :: mydti(0:imx,0:jmx,0:1),mydte(0:imx,0:jmx,0:1)
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
      REAL(8) :: grp,gxdgyp,rhox(4),rhoy(4)
      real(8) :: mydbr(0:imx-1),myjac(0:imx-1),v(0:imx-1),dum1

      rho=0.
      den=0.
      jpar = 0.
      myden = 0.
      myjpar = 0.
      mydti = 0.
      mydte = 0.
      ns=1
if(idg.eq.1)write(*,*)'enter ion grid1',mm(1)
      do m=1,mm(1)
         dv=float(lr(1))*(dx*dy*dz)
         r=x3(m)-0.5*lx+lr0

         k = int(z3(m)/delz)
         wz0 = ((k+1)*delz-z3(m))/delz
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

         rhog=sqrt(2.*b*mu(m)*mims(1))/(q(1)*b)*iflr

         rhox(1) = rhog*(1-tor)+rhog*grp*tor
         rhoy(1) = rhog*gxdgyp/grp*tor
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
         rhox(4) = 0
         rhoy(4) = -rhoy(3)

         vfac=0.5*(mims(1)*u3(m)**2 + 2.*mu(m)*b )
         wght=w3(m)/dv

         vpar = u3(m) !linearly correct

!    now do 1,2,4 point average, where lr is the no. of points...
         do 100 l=1,lr(1)
            xs=x3(m)+rhox(l) !rwx(1,l)*rhog
            yt=y3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.e-8)
            yt = min(yt,ly-1.e-8)

            include "gridli.h"
 100     continue
      enddo
if(idg.eq.1)write(*,*)myid,'pass ion grid1'
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!   enforce periodicity
      call enforce(myden(:,:,:))
      call enforce(myjpar)
      call enforce(mydti)
!      call filter(myden(:,:,:))
!      call filter(myjpar(:,:,:))

      do 110 i=0,im
         do 120 j=0,jm
            do 130 k=0,mykm
               den(1,i,j,k)=q(ns)*myden(i,j,k)/n0/jac(i,k)*cn0i
               jpar(i,j,k) = q(ns)*myjpar(i,j,k)/n0/jac(i,k)*ifluid*cn0i
               mydti(i,j,k) = mydti(i,j,k)/n0/jac(i,k)*cn0i
 130        continue
 120     continue
 110  continue

      mydti(:,:,:) = den(1,:,:,:)
      call MPI_ALLREDUCE(mydti(0:im,0:jm,0:1),  &
     		dti(0:im,0:jm,0:1),             &
     		(imx+1)*(jmx+1)*2,MPI_REAL8,       &
     		MPI_SUM,GRID_COMM,ierr)

      do 150 i=0,im
         do 160 j=0,jm
            do 170 k=0,mykm
               rho(i,j,k)=rho(i,j,k)+den(1,i,j,k)
 170        continue
 160     continue
 150  continue

      if(ishgk==0)goto 299
      myhpar = 0.
      myhden = 0.
      myhden0 = 0.
      do m=1,mm(3)
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
         nhp = wx0*nhi(i)+wx1*nhi(i+1)
!         b=1.-lr0/br0*cost
         b=1.-tor+tor*bfldp

         rhog=sqrt(2.*b*muh(m)*mims(3))/(q(3)*b)*iflrh

         rhox(1) = rhog*(1-tor)+rhog*grp*tor
         rhoy(1) = rhog*gxdgyp/grp*tor
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
         rhox(4) = 0
         rhoy(4) = -rhoy(3)

         vfac=0.5*(mims(3)*uh3(m)**2 + 2.*muh(m)*b )
         wght=wh3(m)/dv
         wght0 = nh*nhp/dv

         vpar = uh3(m) !linearly correct

!    now do 1,2,4 point average, where lr is the no. of points...
         do 200 l=1,lr(1)
            xs=xh3(m)+rhox(l) !rwx(1,l)*rhog
            yt=yh3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.e-8)
            yt = min(yt,ly-1.e-8)

!            i=int(xt/dx+0.5)
!            j=int(yt/dy+0.5)
!            k=int(zh3(m)/dz+0.5)-gclr*kcnt
     
!            myhpar(i,j,k) = myhpar(i,j,k)+wght*vpar
!            myhden(i,j,k) = myhden(i,j,k)+wght
            include 'gridlih.h'
 200     continue
      enddo
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call enforce(myhpar)
      call enforce(myhden)
      call enforce(myhden0)

      do 210 i=0,im
         do 220 j=0,jm
            do 230 k=0,mykm
               hpar(i,j,k) = q(3)*myhpar(i,j,k)/n0/jac(i,k)*cv
               hden(i,j,k) = q(3)*myhden(i,j,k)/n0/jac(i,k)*cv
               myhden0(i,j,k) = q(3)*myhden0(i,j,k)/n0/jac(i,k)*cv
 230        continue
 220     continue
 210  continue

      call MPI_ALLREDUCE(hden(0:im,0:jm,0:1),  &
     		denh(0:im,0:jm,0:1),             &
     		(imx+1)*(jmx+1)*2,MPI_REAL8,       &
     		MPI_SUM,GRID_COMM,ierr)
      call MPI_ALLREDUCE(myhden0(0:im,0:jm,0:1),  &
     		hden0(0:im,0:jm,0:1),             &
     		(imx+1)*(jmx+1)*2,MPI_REAL8,       &
     		MPI_SUM,GRID_COMM,ierr)

!compute denh
      do i = 0,nxpp-1
         dum = 0.
         dum1 = 0.
         do j = 0,jm-1
            dum = dum+(denh(i,j,0))**2*jac(i,0)
            dum1 = dum1+jac(i,0)
         end do
         mydbr(i) = dum
         myjac(i) = dum1
      end do
      call MPI_ALLREDUCE(mydbr,denhr,nxpp,  &
          MPI_REAL8,MPI_SUM,           &
          tube_comm,ierr)
      call MPI_ALLREDUCE(myjac,v,nxpp,  &
          MPI_REAL8,MPI_SUM,           &
          tube_comm,ierr)
      denhr(0:imx-1) = sqrt(denhr(0:imx-1)/v(0:imx-1))

!compute hden0
      do i = 0,nxpp-1
         dum = 0.
         dum1 = 0.
         do j = 0,jm-1
            dum = dum+(hden0(i,j,0))**2*jac(i,0)
            dum1 = dum1+jac(i,0)
         end do
         mydbr(i) = dum
         myjac(i) = dum1
      end do
      call MPI_ALLREDUCE(mydbr,hden0r,nxpp,  &
          MPI_REAL8,MPI_SUM,           &
          tube_comm,ierr)
      call MPI_ALLREDUCE(myjac,v,nxpp,  &
          MPI_REAL8,MPI_SUM,           &
          tube_comm,ierr)
      hden0r(0:imx-1) = sqrt(hden0r(0:imx-1)/v(0:imx-1))
 299  continue

!Beam ions
      if(isbgk==0)goto 399
      mybpar = 0.
      mybden = 0.
      do m=1,mm(4)
         dv=float(lr(1))*(dx*dy*dz)
         r=xb3(m)-0.5*lx+lr0

         k = int(zb3(m)/delz)
         wz0 = ((k+1)*delz-zb3(m))/delz
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
         nbp = wx0*nbi(i)+wx1*nbi(i+1)
!         b=1.-lr0/br0*cost
         b=1.-tor+tor*bfldp

         rhog=sqrt(2.*b*mub(m)*mims(4))/(q(4)*b)

         rhox(1) = rhog*(1-tor)+rhog*grp*tor
         rhoy(1) = rhog*gxdgyp/grp*tor
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
         rhox(4) = 0
         rhoy(4) = -rhoy(3)

         vfac=0.5*(mims(4)*ub3(m)**2 + 2.*mub(m)*b )
         wght=wb3(m)/dv
         wght0 = nbeam*nbp/dv

         vpar = ub3(m) !linearly correct

!    now do 1,2,4 point average, where lr is the no. of points...
         do 300 l=1,lr(1)
            xs=xb3(m)+rhox(l) !rwx(1,l)*rhog
            yt=yb3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.e-8)
            yt = min(yt,ly-1.e-8)

!            i=int(xt/dx+0.5)
!            j=int(yt/dy+0.5)
!            k=int(zh3(m)/dz+0.5)-gclr*kcnt
     
!            myhpar(i,j,k) = myhpar(i,j,k)+wght*vpar
!            myhden(i,j,k) = myhden(i,j,k)+wght
            include 'gridlib.h'
 300     continue
      enddo
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call enforce(mybpar)
      call enforce(mybden)

      do 310 i=0,im
         do 320 j=0,jm
            do 330 k=0,mykm
               bpar(i,j,k) = q(4)*mybpar(i,j,k)/n0/jac(i,k)*cvbeam
               bden(i,j,k) = q(4)*mybden(i,j,k)/n0/jac(i,k)*cvbeam
 330        continue
 320     continue
 310  continue

      call MPI_ALLREDUCE(bden(0:im,0:jm,0:1),  &
     		denb(0:im,0:jm,0:1),             &
     		(imx+1)*(jmx+1)*2,MPI_REAL8,       &
     		MPI_SUM,GRID_COMM,ierr)

 399  continue

      if(iscgk==0)goto 499
      cpar = 0.
      myden = 0.
      myjpar = 0.
if(idg.eq.1)write(*,*)'enter ion grid1',mm(1)
      do m=1,mm(2)
         dv=float(lr(1))*(dx*dy*dz)
         r=xc3(m)-0.5*lx+lr0

         k = int(zc3(m)/delz)
         wz0 = ((k+1)*delz-zc3(m))/delz
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

         rhog=sqrt(2.*b*muc(m)*mims(2))/(q(2)*b)*iflr

         rhox(1) = rhog*(1-tor)+rhog*grp*tor
         rhoy(1) = rhog*gxdgyp/grp*tor
         rhox(2) = -rhox(1)
         rhoy(2) = -rhoy(1)
         rhox(3) = 0
         rhoy(3) = rhog*(1-tor)+rhog/b/grp*fp/radiusp*qhatp*lr0/q0*grcgtp*tor
         rhox(4) = 0
         rhoy(4) = -rhoy(3)

         vfac=0.5*(mims(2)*uc3(m)**2 + 2.*muc(m)*b )
         wght=wc3(m)/dv

         vpar = uc3(m) !linearly correct

!    now do 1,2,4 point average, where lr is the no. of points...
         do 400 l=1,lr(1)
            xs=xc3(m)+rhox(l) !rwx(1,l)*rhog
            yt=yc3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.e-8)
            yt = min(yt,ly-1.e-8)

            include "gridlic.h"
 400     continue
      enddo
if(idg.eq.1)write(*,*)myid,'pass ion grid1'
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!   enforce periodicity
      call enforce(myden(:,:,:))
      call enforce(myjpar)
!      call enforce(mydti)
!      call filter(myden(:,:,:))
!      call filter(myjpar(:,:,:))

      do 410 i=0,im
         do 420 j=0,jm
            do 430 k=0,mykm
               den(2,i,j,k)=q(2)*myden(i,j,k)/n0/jac(i,k)*cn0c
               cpar(i,j,k) = q(2)*myjpar(i,j,k)/n0/jac(i,k)*ifluid*cn0c
 430        continue
 420     continue
 410  continue

      mydti(:,:,:) = den(2,:,:,:)
      call MPI_ALLREDUCE(mydti(0:im,0:jm,0:1),  &
     		dtc(0:im,0:jm,0:1),             &
     		(imx+1)*(jmx+1)*2,MPI_REAL8,       &
     		MPI_SUM,GRID_COMM,ierr)


 499  continue
      do i = 0,im
         do j = 0,jm
            do k = 0,mykm
               rho(i,j,k) = ision*rho(i,j,k) + dene(i,j,k)*qel/ntube
            enddo
         enddo
      enddo      

      if(iske==0)return
! electrons density and current
      vte = sqrt(amie*t0e(nr/2))
      dltpe(:,:,:) = 0.
      dltpa(:,:,:) = 0.
      mydltpe(:,:,:) = 0.
      mydltpa(:,:,:) = 0.
      mydene = 0.
      myupar = 0.
if(idg.eq.1)write(*,*)'enter electron grid1'
      do m=1,mme
         dv=(dx*dy*dz)
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
         b = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
         ter = wx0*t0e(i)+wx1*t0e(i+1)        
         vpar = u3e(m) !linearly correct
         vfac = 0.5*(emass*u3e(m)**2 + 2.*mue3(m)*b)

         xt = x3e(m)
         yt = y3e(m)
         i=int(xt/dx)
         j=int(yt/dy)
         k=0 !int(z3e(m)/dz)-gclr*kcnt
         phip = w000(m)*phi(i,j,k)  &
             + w100(m)*phi(i+1,j,k) &
             + w010(m)*phi(i,j+1,k) &
             + w110(m)*phi(i+1,j+1,k) &
             + w001(m)*phi(i,j,k+1) &
             + w101(m)*phi(i+1,j,k+1) &
             + w011(m)*phi(i,j+1,k+1) &
             + w111(m)*phi(i+1,j+1,k+1)

         wght=w3e(m)/dv*mue3(m)*b
         wght1=w3e(m)/dv*emass*vpar**2
         wght2=w3e(m)/dv
         wght3=w3e(m)/dv*vpar

         if(abs(vfac/ter).gt.vcut)then
            wght = 0.
            wght1 = 0.
            wght2 = 0.
            wght3 = 0.
         end if

         xt=x3e(m)
         yt=y3e(m)
         include 'gridlie.h'
      enddo
if(idg.eq.1)write(*,*)'pass electron grid1'
!   enforce periodicity
      call enforce(mydene(:,:,:))
      call enforce(myupar(:,:,:))
      call enforce(mydltpe(:,:,:))
      call enforce(mydltpa(:,:,:))

      do  i=0,im
         do  j=0,jm
            do  k=0,mykm
               mydene(i,j,k)= mydene(i,j,k)/n0e/jac(i,k)*cn0e*ifluid
               myupar(i,j,k) = myupar(i,j,k)/n0e/jac(i,k)*cn0e*ifluid
               mydltpe(i,j,k)= mydltpe(i,j,k)/n0e/jac(i,k)*cn0e*ifluid
               mydltpa(i,j,k) = mydltpa(i,j,k)/n0e/jac(i,k)*cn0e*ifluid
            end do
         end do
      end do

      call MPI_ALLREDUCE(mydene(0:im,0:jm,0:1),  &
     		denek(0:im,0:jm,0:1),             &
     		(imx+1)*(jmx+1)*2,MPI_REAL8,       &
     		MPI_SUM,GRID_COMM,ierr)

      call MPI_ALLREDUCE(myupar(0:im,0:jm,0:1),  &
     		upark(0:im,0:jm,0:1),             &
     		(imx+1)*(jmx+1)*2,MPI_REAL8,       &
     		MPI_SUM,GRID_COMM,ierr)

      call MPI_ALLREDUCE(mydltpe(0:im,0:jm,0:1),  &
     		dltpe(0:im,0:jm,0:1),             &
     		(imx+1)*(jmx+1)*2,MPI_REAL8,       &
     		MPI_SUM,GRID_COMM,ierr)

      call MPI_ALLREDUCE(mydltpa(0:im,0:jm,0:1),  &
     		dltpa(0:im,0:jm,0:1),             &
     		(imx+1)*(jmx+1)*2,MPI_REAL8,       &
     		MPI_SUM,GRID_COMM,ierr)
      
!      call flty(denek,0)
!      call flty(upark,0)
      call filter(denek)
      call filter(upark)

      call fltx(dltpe,0,0,1,0)
      call fltx(dltpa,0,0,1,0)
      call filter(dltpe)
      call filter(dltpa)
 999  continue

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
               call ccfft('x',1,imx,1.d0,tmpx,coefx,workx,0)
               ii=lmode(mode) !+1
               if(lmode(mode).lt.0) write(*,*) 'lmode < 0, error'
               tmpy(j)=tmpx(ii)/float(im)
            enddo

!     FT in y....
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
            jj=mmode(mode)  !+1
            if(mmode(mode).lt.0) write(*,*) 'mmode < 0, error'
            modebuf=tmpy(jj)/float(jm)

         endif

         call MPI_BCAST(modebuf,1,MPI_DOUBLE_COMPLEX,oproc, &
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
      character*4 holdmyid

      write(holdmyid,'(I4.4)') MyId
      fname=directory//'dump_'//holdmyid//'.b'
      if(iflag.eq.1) then
         open(139+MyId,file=fname,form='unformatted',status='old')
         ns=1
         read(139+MyId)mm(ns),ncurr,tcurr,rmpp,rmaa
         do 110 m=ptr(ns),ptr(ns)+mm(ns)-1
            read(139+MyId) mu(m)
            read(139+MyId) x2(m),y2(m),z2(m),u2(m),w2(m)
            read(139+MyId) xii(m),z0i(m),pzi(m),eki(m)
            w2(m)=w2(m)/cut
            x3(m)=x2(m)
            y3(m)=y2(m)
            z3(m)=z2(m)
            u3(m)=u2(m)
            w3(m)=w2(m)
 110     continue

         if(ishgk==1)then
         ns=3
         read(139+MyId)mm(ns)
         do 120 m=1,mm(ns)
            read(139+MyId) muh(m),muh2(m),mui(m)
            read(139+MyId) xh2(m),yh2(m),zh2(m),uh2(m),wh2(m)
            read(139+MyId) xih(m),z0h(m),pzh(m),ekh(m),vih(m),index(m)
            wh2(m)=wh2(m)/cut
            xh3(m)=xh2(m)
            yh3(m)=yh2(m)
            zh3(m)=zh2(m)
            uh3(m)=uh2(m)
            wh3(m)=wh2(m)
 120     continue
         end if

         if(isbgk==1)then
         ns=4
         read(139+MyId)mm(ns)
         do 140 m=1,mm(ns)
            read(139+MyId) mub(m),mub2(m)
            read(139+MyId) xb2(m),yb2(m),zb2(m),ub2(m),wb2(m)
            wb2(m)=wb2(m)/cut
            xb3(m)=xb2(m)
            yb3(m)=yb2(m)
            zb3(m)=zb2(m)
            ub3(m)=ub2(m)
            wb3(m)=wb2(m)
 140     continue
         end if

         if(iscgk==1)then
         ns=2
         read(139+MyId)mm(ns)
         do 150 m=1,mm(ns)
            read(139+MyId) muc(m)
            read(139+MyId) xc2(m),yc2(m),zc2(m),uc2(m),wc2(m)
            wc2(m)=wc2(m)/cut
            xc3(m)=xc2(m)
            yc3(m)=yc2(m)
            zc3(m)=zc2(m)
            uc3(m)=uc2(m)
            wc3(m)=wc2(m)
 150     continue
         end if

 160     continue
         if(iske==1)then
         read(139+MyId)mme
         do  m=1,mme
            read(139+MyId) x2e(m),y2e(m),z2e(m),u2e(m),w2e(m),mue2(m)
            read(139+MyId) xie(m),z0e(m),pze(m),eke(m),mue(m),u0e(m)
            w2e(m)=w2e(m)/cut
            x3e(m)=x2e(m)
            y3e(m)=y2e(m)
            z3e(m)=z2e(m)
            u3e(m)=u2e(m)
            w3e(m)=w2e(m)
            mue3(m)=mue2(m)
         end do
         end if

         read(139+myid)apar,dene,phi,phiom,delte !,dteom
         read(139+myid)ke,fe,te,pfl,efl,nos,rmsphi,pmodehis,camp
         close(139+MyId)
      endif

      if(iflag.eq.2) then
         open(139+MyId,file=fname,form='unformatted',status='unknown')
         ns=1
         write(139+MyId)mm(ns),n+1,tcurr-dt,rmpp,rmaa
         do 210 m=ptr(ns),ptr(ns)+mm(ns)-1
            write(139+MyId) mu(m)
            write(139+MyId) x2(m),y2(m),z2(m),u2(m),w2(m)
            write(139+MyId) xii(m),z0i(m),pzi(m),eki(m)
 210     continue

         if(ishgk==1)then
         ns=3
         write(139+MyId)mm(ns)
         do 220 m=1,mm(ns)
            write(139+MyId) muh(m),muh2(m),mui(m)
            write(139+MyId) xh2(m),yh2(m),zh2(m),uh2(m),wh2(m)
            write(139+MyId) xih(m),z0h(m),pzh(m),ekh(m),vih(m),index(m)
 220     continue
         end if

         if(isbgk==1)then
         ns=4
         write(139+MyId)mm(ns)
         do 240 m=1,mm(ns)
            write(139+MyId) mub(m),mub2(m)
            write(139+MyId) xb2(m),yb2(m),zb2(m),ub2(m),wb2(m)
 240     continue
         end if

         if(iscgk==1)then
         ns=2
         write(139+MyId)mm(ns)
         do 250 m=1,mm(ns)
            write(139+MyId) muc(m)
            write(139+MyId) xc2(m),yc2(m),zc2(m),uc2(m),wc2(m)
 250     continue
         end if

 260     continue
         if(iske==1)then
         write(139+MyId)mme
         do  m=1,mme
            write(139+MyId) x2e(m),y2e(m),z2e(m),u2e(m),w2e(m),mue2(m)
            write(139+MyId) xie(m),z0e(m),pze(m),eke(m),mue(m),u0e(m)
         end do
         end if

         write(139+myid)apar,dene,phi,phiom,delte !,dteom
         write(139+myid)ke,fe,te,pfl,efl,nos,rmsphi,pmodehis,camp
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
            dely=dmod(2.*pi*lr0/q0*qr*dsign(1.d0,q0)+800.*ly,ly)
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
            dely=dmod(2.*pi*lr0/q0*qr*dsign(1.d0,q0)+800.*ly,ly)
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
      character*4 holdmyid

      save formphi,formfe,ifirst,akx,aky,mx,gamb1,gamb2,ipiv
      if(idg==1)write(*,*)'enter gkps'
      izonal = 1

      write(holdmyid,'(I4.4)') MyId
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

                     m1 = mstart+int((float(m)+1.0)/2)
                     if(m==0)m1=0
                     sgny = isgnft(m)

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
                     if(i==ix)mx(i,ix,j,k) = fradi
                     do ikx = 0,imx-1
                        mx(i,ix,j,k) = mx(i,ix,j,k)+q(1)*sin(ix*ikx*pi/imx)* & 
                             ((1-gamb1(ikx,j,i,k))*exp(IU*ikx*i*pi/imx)- &
                              (1-gamb2(ikx,j,i,k))*exp(-IU*ikx*i*pi/imx)) &
                             /ter*gn0e(i)/(IU*imx)
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
      u = rho+hden*ishgk+bden*isbgk+den(2,:,:,:)*iscgk
      call dcmpy(u(0:imx-1,0:jmx-1,0:1),v)
      if(idg==1)write(*,*)'pass first fft'

      temp3dxy = 0.
      do k = 0,1
         do j = 0,jcnt-1   !j=0,jcnt-1 to include n=0 mode
            myj=jft(j)
            sl(1:imx-1,j,k) = v(1:imx-1,j,k)
            call ZGETRS('N',imx-1,1,mx(:,:,j,k),imx-1,ipiv(:,:,j,k), &
                         sl(:,j,k),imx-1,INFO) 
            temp3dxy(1:imx-1,myj,k) = sl(1:imx-1,j,k)
            temp3dxy(0,myj,k) = 0.
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
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
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
      if(iez==0)return

      if(iske==1)dene = denek
      rbfs = dene(:,:,0)

      call MPI_SENDRECV(rbfs(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8,rngbr,204,                    &
          lbfr(0,0),(imx+1)*(jmx+1),              &
          MPI_REAL8,lngbr,204,                    &
          TUBE_COMM,stat,ierr)

      lbfs=dene(:,:,1)

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
                -dene(i,j,1))/(2.*dz)*bdgrzn(i,0)*gt0e(i)/gn0e(i) &
                +gt0e(i)*delbx(i,j,0)*(gcpne(i)+gcpte(i)*isdte)*bdgxcgy(i,0)+etaohm*djpa(i,j,0)/gn0e(i)

            ez(i,j,1)=( dene(i,j,0)-(weightm(i)*rbfr(i,jmi(i,j)) &
                +weightmn(i)*rbfr(i,jmn(i,j))))/(2.*dz)*bdgrzn(i,1)*gt0e(i)/gn0e(i) &
                +gt0e(i)*delbx(i,j,1)*(gcpne(i)+gcpte(i)*isdte)*bdgxcgy(i,1)+etaohm*djpa(i,j,1)/gn0e(i)

         enddo
      enddo

      call filter(ez)
      return
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
      real(8) :: tdum

      m = imx/8
      k = imx/16
      tdum = tcurr-dt
      if(myid.eq.master)then
         write(*,10)tdum,rmsphi(n),rmsapa(n),rmsez(n),efl(1,n),efl(2,n),avewi(n),avewc(n),yyre(mlk,0),yyre(mlk,1)
!            write(*,12)i,dbr(m),dbr(2*m),dbr(3*m),dbr(4*m),dbr(5*m),dbr(6*m),dbr(7*m),dtr(m),dtr(2*m),dtr(3*m),dtr(4*m),dtr(5*m),dtr(6*m),dtr(7*m)
 10      format(1x,f10.1,15(2x,e10.3),2x,i8,2x,i8)
 11      format(6x,5(2x,e12.5))
 12      format(1x,f10.1,18(2x,e10.3))

         write(9,10)tdum,rmsphi(n),rmsapa(n),rmsez(n),efl(1,n),efl(2,n),avewi(n),avewc(n),yyre(mlk,0),yyre(mlk,1)
         write(11,12)tdum,pfle(4,n),pfl(1,n),pfl(2,n),pfl(4,n),efle(4,n),efl(1,n),efl(2,n),efl(4,n)
!         write(17,12)i,yyre(mlk,0),yyre(mlk,1),yyre(mlk,2),yyre(mlk,3),yyre(mlk,4)!, &
!             yyre(6,0),yyre(7,0),yyre(8,0),yyre(9,0),yyre(10,0),yyre(11,0),yyre(12,0)
         write(18,12)tdum,dbr(m),dbr(2*m),dbr(3*m),dbr(4*m),dbr(5*m),dbr(6*m),dbr(7*m),dtr(m),dtr(2*m),dtr(3*m),dtr(4*m),dtr(5*m),dtr(6*m),dtr(7*m)
!         write(18,12)i,dtr(2*k),dtr(3*k),dtr(4*k),dtr(5*k),dtr(6*k),dtr(7*k),dtr(8*k),dtr(9*k),dtr(10*k),dtr(11*k),dtr(12*k),dtr(13*k),dtr(14*k),dtr(15*k)
      end if   
      if(gclr==kmx/2 .and. tclr==0)then
         write(23,12)tdum,mdhis(0),mdhis(1),mdhis(2),mdhisa(0),mdhisa(1),mdhisa(2)

         write(24,12)tdum,mdhisb(0),mdhisb(1),mdhisc(0),mdhisc(1),mdhisd(0),mdhisd(1)
         
         do  i = 0,6             
            write(22,13)tdum,i,real(phihis(i,0)),(real(phihis(i,j)), aimag(phihis(i,j)), j = 1,jcnt-2,2)
            write(22,13)tdum,i,real(aparhis(i,0)),(real(aparhis(i,j)), aimag(aparhis(i,j)), j = 1,jcnt-2,2)
         end do
 13      format(1x,f10.1,1x,i2,33(2x,e10.3))

         write(25,14)tdum,(reyn00(i),i = 0,imx-1)
         write(25,14)tdum,(reynh00(i),i = 0,imx-1)
         write(25,14)tdum,(maxw00(i),i = 0,imx-1)
         write(25,14)tdum,(maxwe00(i),i = 0,imx-1)
!         write(25,14)tdum,(maxwh00(i),i = 0,imx-1)
         write(25,14)tdum,(drdt00(i),i = 0,imx-1)

 14      format(1x,f10.1,128(2x,e12.5))
      end if
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ezamp(nstep,ip)   

      use gem_com
      use equil
      use fft_wrapper
      implicit none
      REAL(8) :: b2,gam0,gam1,th,delyz,bf,rkper,r,qr,shat,ter !,b
      REAL(8) :: kx,ky
      complex(8) :: lbfr(0:imx-1,0:jcnt-1)
      complex(8) :: lbfs(0:imx-1,0:jcnt-1)
      complex(8) :: rbfr(0:imx-1,0:jcnt-1)
      complex(8) :: rbfs(0:imx-1,0:jcnt-1)
      REAL(8),dimension(:),allocatable :: dely,aky
      complex(8),dimension(:,:,:,:),allocatable:: nab1,nab2
      complex(8),dimension(:,:,:,:),allocatable :: mx
      REAL(8),dimension(:,:,:),allocatable :: formapa
      REAL(8) :: sgnx,sgny,sz,myfe
      INTEGER :: i,i1,j,j1,k,k1,l,m,n,ifirst,nstep,ip,INFO
      INTEGER :: l1,m1,myk,myj,izonal,ix,ikx,iext
      COMPLEX(8) :: temp3dxy(0:imx-1,0:jmx-1,0:1),va(0:imx-1,0:jmx-1,0:1),vb(0:imx-1,0:jmx-1,0:1),&
            v(0:imx-1,0:jcnt-1,0:1)
      COMPLEX(8) :: sbuf(0:imx*jcnt*2-1),rbuf(0:imx*jmx*2-1)
      COMPLEX(8) :: sl(1:imx-1,0:jcnt-1,0:1)
      REAL(8) :: dum,u(0:imx,0:jmx,0:1)
      complex(8) :: ua(0:imx-1,0:jcnt-1,0:1),calph,cbeta
      REAL(8) :: grp,gthp,gxdgyp,dydrp,qhatp,grdgtp,bfldp,g2xp,g2yp, &
                 radiusp,g2zp,gzp,gxdgzp,gydgzp
      REAL(8) :: wx0,wx1,wz0,wz1
      character*70 fname
      character*4 holdmyid

      save formapa,ifirst,dely,aky,mx
      calph = 1.
      cbeta = 0.
      write(holdmyid,'(I4.4)') MyId
      fname='./matrix/'//'mx_apa_'//holdmyid
      if (ifirst.ne.-99) then
         allocate(dely(0:imx-1),aky(0:jmx-1),nab1(0:imx-1,0:jcnt-1,0:imx-1,0:1))
         allocate(nab2(0:imx-1,0:jcnt-1,0:imx-1,0:1))
         allocate(mx(imx-1,imx-1,0:jcnt-1,0:1),formapa(0:imx-1,0:jcnt-1,0:1))

         if(iget.eq.1) then
            open(20000+MyId,file=fname,form='unformatted',status='old')
            read(20000+MyId)mx
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
                     ter = wx0*t0i(i1)+wx1*t0i(i1+1)
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
                     radiusp = wx0*wz0*radius(i1,k1)+wx0*wz1*radius(i1,k1+1) &
                     +wx1*wz0*radius(i1+1,k1)+wx1*wz1*radius(i1+1,k1+1)

                     m1 = mstart+int((float(m)+1.0)/2)
                     if(m==0)m1=0
                     sgny = isgnft(m)

                     ky=sgny*2.*pi*float(m1)/ly
                     kx=pi*float(l)/lx
                     bf=bfldp
!                     b=mims(1)*(kx*kx*grp**2 + &
!                        ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
!                        +2*dydrp*lr0/q0*qhatp*grdgtp) &
!                        +2*kx*ky*gxdgyp)/bf/bf
                  
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
!                     if(i==ix)mx(i,ix,j,k) = beta*amie*(1+nh*ish)
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
!         do k = 0,1
!            do j = 0,jcnt-1
!               call ZGETRF( imx-1,imx-1,mx(:,:,j,k),imx-1,IPIV(:,:,j,k), INFO )
!            end do
!         end do

         if(iget.eq.0) then
            open(20000+MyId,file=fname,form='unformatted',status='unknown')
            write(20000+MyId)mx
            goto 200
         end if

 200     ifirst=-99
      endif

!   now do field solve...

!  find jtot(kx,ky)
      iext = 0
      if(nstep>200)iext = 0
      do i = 0,imx
         do j = 0,jmx
            do k = 0,1
               u(i,j,k) = jpar(i,j,k)*ision+hpar(i,j,k)*isham+bpar(i,j,k)*isbam+cpar(i,j,k)*iscam &
                 +0.1*sin(mlk*yg(j)-0.04*tcurr)*exp(-(xg(i)-lx/2)**2/(lx/3)**2) &
                  *exp(-(zg(k)-lz/2)**2/(lz/2)**2)*iext
            end do
         end do
      end do
      call dcmpy(u(0:imx-1,0:jmx-1,0:1),v)
      call dcmpy(apar(0:imx-1,0:jmx-1,0:1)/ntube,ua)

      temp3dxy = 0.
      do k = 0,1
         do j = 0,jcnt-1  !no n=0 !!!
            myj=jft(j)
            call ZGEMV('N',imx-1,imx-1,calph,mx(:,:,j,k),imx-1, &
                       ua(1:imx-1,j,k), 1,cbeta,sl(1:imx-1,j,k), 1)
            temp3dxy(1:imx-1,myj,k) = sl(1:imx-1,j,k)
            temp3dxy(0,myj,k) = 0.
         end do
      end do

!  from rho(kx,ky) to phi(kx,ky)
      do k = 0,1
         do j = 0,jmx-1
            do i = 0,imx-1
               temp3dxy(i,j,k)=temp3dxy(i,j,k)/jmx
            end do
         end do
      end do

      v = v/jmx
      va = 0.
      do k = 0,1
         do j = 0,jcnt-1
            myj = jft(j)
            va(0:imx-1,myj,k) = v(0:imx-1,j,k)
         end do
      end do

!djpa for Ohm's law
      vb = 0.
      do i = 0,nxpp-1
         do j = 0,jm-1
            do k = 0,mykm
                  vb(i,j,k) = 1/beta*temp3dxy(i,j,k)
            end do
         end do
      end do

      do i = 0,nxpp-1
         do j = 0,jm-1
            do k = 0,mykm
                  temp3dxy(i,j,k) = (va(i,j,k)*1.-1/beta*temp3dxy(i,j,k))/gn0e(i)
            end do
         end do
      end do


!  from phi(kx,ky) to phi(x,y)
      do k=0,mykm
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3dxy(i,j,k)
            end do
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3dxy(i,j,k) = tmpy(j)  ! phi(x,y)
            end do
         end do
      end do

      do i = 0,nxpp-1
         do j = 0,jm-1
            do k = 0,mykm
                  upar(i,j,k) = temp3dxy(i,j,k)
            end do
         end do
      end do

!    x-y boundary points 
      call enfxy(upar(:,:,:))
      call enfz(upar(:,:,:))
!      call filter(upar(:,:,:))

!  djpa
      do k=0,mykm
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = vb(i,j,k)
            end do
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               vb(i,j,k) = tmpy(j)  ! phi(x,y)
            end do
         end do
      end do

      do i = 0,nxpp-1
         do j = 0,jm-1
            do k = 0,mykm
                  djpa(i,j,k) = vb(i,j,k)
            end do
         end do
      end do

!    x-y boundary points 
      call enfxy(djpa(:,:,:))
      call enfz(djpa(:,:,:))

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
      parameter(NNI = 5)
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

      call ccfft('z',0,kmx,1.d0,tmpz,coefz,workz,0)

!      if(n.eq.0) return

      do j=0,jm-1
         do i=0,im-1
            do k = 0,mykm-1
               tmp3d(i,j,k)=u(i,j,k)
            enddo
         end do
      end do

      do k = 0,mykm-1
!         do j = 0,jmx-1
!            do i = 0,imx-1
!               tmpx(i) = tmp3d(i,j,k)
!            end do
!            call ccfft('x',1,imx,1.d0,tmpx,coefx,workx,0)
!            do i = 0,imx-1
!               tmp3d(i,j,k) = tmpx(i)
!            end do
!         end do
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = tmp3d(i,j,k)
            end do
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               tmp3d(i,j,k) = tmpy(j)
            end do
         end do
      end do   

      if(gclr==kmx/2 .and. tclr==0 .and. mod(n,1000)==0)then
         write(17,*)n
         do i = 0,imx-1
            write(17,10)i,(rin+i*dx)/a,real(tmp3d(i,mlk,0)),aimag(tmp3d(i,mlk,0)),abs(tmp3d(i,mlk,0))
         end do
      end if
 10   format(1x,i5,4(2x,e10.3))

!      do n=1,12
      j=mlk
!      do j = 1,4
         do k = 0,mykm-1
            tmpz(k) = tmp3d(llk,j,k)
         end do
         
         if(GCLR.ne.master)then
            call MPI_SEND(tmpz(0),mykm,MPI_DOUBLE_COMPLEX,master, &
               gclr,tube_comm,ierr)
         end if

         if(gclr.eq.master)then
            do i = 1,GLST
               call MPI_RECV(tmpz(i*mykm),mykm,MPI_DOUBLE_COMPLEX,i, &
                   i,tube_comm,stat,ierr)
            end do
         end if

         if(GCLR.eq.master) then
            call ccfft('z',1,kmx,1.d0,tmpz,coefz,workz,0)
         end if

         call MPI_BCAST(tmpz,kmx,MPI_DOUBLE_COMPLEX,master, &
            tube_comm,ierr)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         do i = 0,4
            yyamp(mlk,i) = abs(tmpz(i)) !cabs
            yyre(mlk,i) = real(tmpz(i))
            yyim(mlk,i) = aimag(tmpz(i))
         end do
!      end do
!      end do

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine yveck1(u,n)

!     calculate mode histories. calculate modes for u at timestep n
!     and store in modehis(mode,n).

      use gem_com
      use equil
      use fft_wrapper
      implicit none
!     
      REAL(8) :: u(0:imx,0:jmx,0:1)
      INTEGER :: n,i,j,k,m
      COMPLEX(8) :: tmp3d(0:imx,0:jmx,0:1)

      call ccfft('z',0,kmx,1.d0,tmpz,coefz,workz,0)

!      if(n.eq.0) return

      tmp3d = 0.
      do j=0,jm-1
         do i=0,im-1
            do k = 0,mykm
               tmp3d(i,j,k)=u(i,j,k)
            enddo
         end do
      end do

      do k = 0,mykm
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = tmp3d(i,j,k)
            end do
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               tmp3d(i,j,k) = tmpy(j)
            end do
         end do
      end do   
      do i = 0,im
         do k = 0,1
            do m = 0,nom-1
               phiom(i,0,k,m) = phiom(i,0,k,m)+tmp3d(i,mlk,k)*exp(-IU*omlk(m)*n*dt)
               phiom(i,1,k,m) = phiom(i,1,k,m)+tmp3d(i,jmx-mlk,k)*exp(-IU*omlk(m)*n*dt)
            end do
         end do
      end do

      j = int(n/ifskp)
      if(mod(n,ifskp)==0)then
         camp(:,j)=0.
         if(j>0)camp(:,j-1)=camp(:,j-1)/ifskp
      end if
      if(gclr==kmx/2)then
         do i = 0,6
            camp(i,j) = camp(i,j)+u(imx/8*(i+1),jmx/2,0)
         end do
      end if
      call MPI_BCAST(camp,7*50000,MPI_DOUBLE_COMPLEX,kmx/2, &
            tube_comm,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine yveck2(u,n)

!     calculate mode histories. calculate modes for u at timestep n
!     and store in modehis(mode,n).

      use gem_com
      use equil
      use fft_wrapper
      implicit none
!     
      REAL(8) :: u(0:imx,0:jmx,0:1)
      INTEGER :: n,i,j,k,m
      COMPLEX(8) :: tmp3d(0:imx,0:jmx,0:1)

      call ccfft('z',0,kmx,1.d0,tmpz,coefz,workz,0)

!      if(n.eq.0) return

      tmp3d = 0.
      do j=0,jm-1
         do i=0,im-1
            do k = 0,mykm
               tmp3d(i,j,k)=u(i,j,k)
            enddo
         end do
      end do

      do k = 0,mykm
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = tmp3d(i,j,k)
            end do
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               tmp3d(i,j,k) = tmpy(j)
            end do
         end do
      end do   
      do i = 0,im
         do k = 0,1
            do m = 0,nom-1
               dteom(i,0,k,m) = dteom(i,0,k,m)+tmp3d(i,1,k)*exp(-IU*omlk(m)*n*dt)
               dteom(i,1,k,m) = dteom(i,1,k,m)+tmp3d(i,jmx-1,k)*exp(-IU*omlk(m)*n*dt)
            end do
         end do
      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

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
      subroutine loadi

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

!     load a slab of ions...

         dumx=lx*ran2(iseed) !revers(MyId*cnt+j,2) !ran2(iseed)
         dumy=ly*ran2(iseed) !revers(MyId*cnt+j,3) !ran2(iseed)
         dumz=lz*ran2(iseed) !revers(MyId*cnt+j,5) !ran2(iseed)
         dumz = min(dumz,lz-1.e-8)
         r = lr0+dumx-0.5*lx
         th = (dumz-lz/2)/(q0*br0) !should be ok for q0<0
         i = int((r-rin)/dr)
         k = int((pi+th)/dth)
         jacp = jacob(i,k)
         if(ran2(iseed)<(0.5*jacp/jacmax))then
         m = m+1
         if(m>mm(1))goto 170
         x2(m)=min(dumx,lx-1.d-8)
         y2(m)=min(dumy,ly-1.d-8)

         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1-wz0
         z2(m) = wz0*zfnth(k)+wz1*zfnth(k+1)
         z2(m)=min(z2(m),lz-1.d-8)

         call parperp(vpar,vperp2,m,pi,cnt,MyId)
!   normalizations will be done in following loop...

         r=x2(m)-0.5*lx+lr0
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
         fp = wx0*f(i)+wx1*f(i+1)
         ter = wx0*t0i(i)+wx1*t0i(i+1)
         b=1.-tor+tor*bfldp

         u2(m)=vpar/sqrt(mims(1)/tload)
         mu(m)=0.5*vperp2/b*tload
         eki(m) = mu(m)*b+0.5*mims(1)*u2(m)**2
         pzi(m) = mims(1)*u2(m)/b*fp/br0-q(1)*psp/br0
         z0i(m) = z2(m)
         xii(m) = x2(m)
         myavgv=myavgv+u2(m)

!    LINEAR: perturb w(m) to get linear growth...
         w2(m)=2.*amp*sin(pi2/ly*y2(m))*exp(-(z2(m)-lz/2)**2/(lz/8)**2)*exp(-(x2(m)-0.4*lx)**2/(lx/8)**2)
!         w2(m)=2.*amp*(revers(MyId*cnt+m,13)-0.5)*0. !(ran2(iseed) - 0.5 )
!         if(izon.eq.1)w2(m)=2.*amp*sin(x2(m)/lx*2*pi)  
         myavgw=myavgw+w2(m)
         end if
 160  continue
 170  continue
      myavgw = myavgw/mm(1)
!      write(*,*)'myid ', myid,x2(10),y2(20),u2(20),w2(20)
!    subtract off avg. u...
      call MPI_ALLREDUCE(myavgv,avgv,1, &
          MPI_REAL8, &
          MPI_SUM,MPI_COMM_WORLD,ierr)
      if(idg.eq.1)write(*,*)'all reduce'
      avgv=avgv/float(tmm(1))
      do 180 m=1,mm(1)
         u2(m)=u2(m)-avgv
         x3(m)=x2(m)
         y3(m)=y2(m)
         z3(m)=z2(m)
         u3(m)=u2(m)
!         w2(m) = w2(m)-myavgw
         w3(m)=w2(m)
 180  continue

      np_old=mm(1)
      call init_pmove(z3,np_old,lz,ierr)
      if(idg.eq.1)write(*,*)'pass init_pmove'
!     
      call pmove(x2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mu,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call pmove(xii,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0i,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(pzi,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(eki,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
!     
      call end_pmove(ierr)
      mm(1)=np_new
      mm(2)=mm(1)
!     write(*,*)MyId,j,mm(1)
!     
!    species two loaded on top of ions, for now species one
!    is ions and species two is electrons, additional species will
!    require modification of loader to initialize particle array...

      do 260 ns=2,nsm
         do 250 m=1,mm(1)
            m1=ptr(ns)+m-1
            x2(m1)=x2(m)
            y2(m1)=y2(m)
            z2(m1)=z2(m)
            u2(m1)=sqrt(mims(ns)*tets(ns))*u2(m)
            mu(m1)=mims(ns)*tets(ns)*mu(m)

 250     continue
 260  continue

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine loadc

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

      mm(2)=int(tmm(1)/numprocs)
      m = 0
      do 160 j = 1,100000000

!     load a slab of ions...

         dumx=lx*ran2(iseed) !revers(MyId*cnt+j,2) !ran2(iseed)
         dumy=ly*ran2(iseed) !revers(MyId*cnt+j,3) !ran2(iseed)
         dumz=lz*ran2(iseed) !revers(MyId*cnt+j,5) !ran2(iseed)
         dumz = min(dumz,lz-1.e-8)
         r = lr0+dumx-0.5*lx
         th = (dumz-lz/2)/(q0*br0) !should be ok for q0<0
         i = int((r-rin)/dr)
         k = int((pi+th)/dth)
         jacp = jacob(i,k)
         if(ran2(iseed)<(0.5*jacp/jacmax))then
         m = m+1
         if(m>mm(1))goto 170
         xc2(m)=min(dumx,lx-1.d-8)
         yc2(m)=min(dumy,ly-1.d-8)

         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1-wz0
         zc2(m) = wz0*zfnth(k)+wz1*zfnth(k+1)
         zc2(m)=min(zc2(m),lz-1.d-8)

         call parperp(vpar,vperp2,m,pi,cnt,MyId)
!   normalizations will be done in following loop...

         r=xc2(m)-0.5*lx+lr0
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
         fp = wx0*f(i)+wx1*f(i+1)
         ter = wx0*t0c(i)+wx1*t0c(i+1)
         b=1.-tor+tor*bfldp

         uc2(m)=vpar/sqrt(mims(2)/tloadc)
         muc(m)=0.5*vperp2/b*tloadc
         myavgv=myavgv+uc2(m)

!    LINEAR: perturb w(m) to get linear growth...
         wc2(m)=0. ! 2.*amp*(revers(MyId*cnt+m,13)-0.5) !(ran2(iseed) - 0.5 )
!         if(izon.eq.1)w2(m)=2.*amp*sin(x2(m)/lx*2*pi)  
         myavgw=myavgw+wc2(m)
         end if
 160  continue
 170  continue
      myavgw = myavgw/mm(2)
!      write(*,*)'myid ', myid,x2(10),y2(20),u2(20),w2(20)
!    subtract off avg. u...
      call MPI_ALLREDUCE(myavgv,avgv,1, &
          MPI_REAL8, &
          MPI_SUM,MPI_COMM_WORLD,ierr)
      if(idg.eq.1)write(*,*)'all reduce'
      avgv=avgv/float(tmm(1))
      do 180 m=1,mm(2)
         uc2(m)=uc2(m)-avgv
         xc3(m)=xc2(m)
         yc3(m)=yc2(m)
         zc3(m)=zc2(m)
         uc3(m)=uc2(m)
!         w2(m) = w2(m)-myavgw
         wc3(m)=wc2(m)
 180  continue

      np_old=mm(2)
      call init_pmove(zc3,np_old,lz,ierr)
      if(idg.eq.1)write(*,*)'pass init_pmove'
!     
      call pmove(xc2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xc3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yc2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yc3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zc2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zc3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(uc2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(uc3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wc2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wc3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(muc,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

!     
      call end_pmove(ierr)
      mm(2)=np_new
!     write(*,*)MyId,j,mm(1)
!     
!    species two loaded on top of ions, for now species one
!    is ions and species two is electrons, additional species will
!    require modification of loader to initialize particle array...

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
	subroutine initialize
         use gem_com
      use equil
      use fft_wrapper
	implicit none
        real(8) :: dum,dum1,dum2,jacp,xndum,r,wx0,wx1
!        complex(8),dimension(0:1) :: x,y
        real(8),dimension(0:1) :: x,y
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
         do i = 0,im-1
            dum = dum+(jac(i,0)+jac(i+1,0)+jac(i,1)+jac(i+1,1))/4
         end do
         call MPI_ALLREDUCE(dum,jacp,1,  &
             MPI_REAL8,MPI_SUM,           &
             tube_comm,ierr)
         totvol = dx*ly*dz*jacp    
         n0=float(tmm(1))/totvol
         n0e=float(tmme)/totvol
         cv = 4.*pi/3.*log((vc**3+vi**3)/(vmin**3+vi**3))

         cvbeam = 4.*pi/3.*log((vcbeam**3+vibeam**3)/(vmin**3+vibeam**3))

!  calculate the volume of each radial subdomain
         do k = 1,nsubd
            dum = 0.
            do i = (k-1)*im/nsubd,k*im/nsubd-1
               r = xg(i)-0.5*lx+lr0
               j = int((r-rin)/dr)
               j = min(j,nr-1)
               wx0 = (rin+(j+1)*dr-r)/dr
               wx1 = 1.-wx0
               dum = dum+(jac(i,0)+jac(i+1,0)+jac(i,1)+jac(i+1,1))/4
            end do
            call MPI_ALLREDUCE(dum,jacp,1,  &
                MPI_REAL8,MPI_SUM,           &
                tube_comm,ierr)
            vol(k) = dx*ly*dz*jacp    
         end do

	call weight
!     initialize particle quantities...
         if( cut.eq.0.) cut=1.
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         call ccfft('x',0,imx,0.d0,tmpx,coefx,workx,0)
         call ccfft('y',0,jmx,0.d0,tmpy,coefy,worky,0)
         call ccfft('z',0,kmx,0.d0,tmpz,coefz,workz,0)
         call dsinf(1,x,1,0,1,0,imx*2,1,1.d0,aux1,50000,aux2,20000)
         
         ncurr = 1
         call blendf
	end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine loader_wrapper
         use gem_com
      use equil
	implicit none

	integer :: n,i,j,k,ip

         if(isuni.eq.0)call loadi
         if(ishgk==1)call loadh
         if(isbgk==1)call loadb
         if(iscgk==1)call loadc
         if(iske.eq.1)call ldel
         if(idg.eq.1)write(*,*)'past loader'
end subroutine loader_wrapper
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine accumulate(n,ip)
         use gem_com
         use equil
	implicit none

	integer :: n,i,j,k,ip
        if(iske==1)call setw(ip,n)
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

        call gkps(n,ip)
        call fltx(phi(:,:,:),0,0,0,0)
        call filter(phi)

        myrmsphi=0.
        rmsphi(n)=0.
        do k=0,mykm-1
           do j=0,jm-1
              do i=0,im-1
                 myrmsphi=myrmsphi+phi(i,j,k)*phi(i,j,k)
              enddo
           enddo
        enddo
        call MPI_ALLREDUCE(myrmsphi,rmsphi(n),1, &
            MPI_REAL8, &
            MPI_SUM,tube_comm,ierr)
        rmsphi(n)=sqrt(rmsphi(n)/(im*jm*km))

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

        call ezamp(n,ip)
!        upar = upark
        if(idg.eq.1)write(*,*)'pass ezamp'
        call fltx(upar(:,:,:),0,0,0,0)
        call filter(upar(:,:,:))

	if(idg.eq.1)write(*,*)'pass filter(apar)'
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine ampere
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine field(n,ip)
         use gem_com
         use equil
	implicit none
        integer :: n,i,j,k,ip,i1
        real(8) :: lbfr(0:imx,0:jmx)
        real(8) :: lbfs(0:imx,0:jmx)
        real(8) :: rbfr(0:imx,0:jmx)
        real(8) :: rbfs(0:imx,0:jmx)
        real(8) :: dum
        real(8) :: myrmsphi,rmp(20),myavap(0:imx-1)
        real(8) :: myjaca(0:imx-1),jaca(0:imx-1)

	if(idg.eq.2)write(*,*)'enter field'
	call grad(ip)
	if(idg.eq.2)write(*,*)'pass grad in field'
        if(iez==0)then
           ez = 0
           goto 100
        end if

        do i = 0,im
           do j = 0,jm
              rbfs(i,j)=phi(i,j,0)
           end do
        end do
        call MPI_SENDRECV(rbfs,(imx+1)*(jmx+1),  &
          MPI_REAL8,rngbr,404, &
          lbfr,(imx+1)*(jmx+1), &
          MPI_REAL8,lngbr,404, &
          tube_comm,stat,ierr)
        do i = 0,im
           do j = 0,jm
              lbfs(i,j)=phi(i,j,1)
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
              dpdz(i,j,0) = (phi(i,j,1)-dum)/(2.*dz)

              dum = weightm(i)*rbfr(i,jmi(i,j)) &
                +weightmn(i)*rbfr(i,jmn(i,j))
              dpdz(i,j,1) = (dum-phi(i,j,0))/(2.*dz)
           end do
        end do
        call enfxy(dpdz)
        call enfz(dpdz)
        do i = 0,imx
           do j = 0,jmx
              do k = 0,1
                 ezs(i,j,k) = -dpdz(i,j,k)*bdgrzn(i,k)
              end do
           end do
        end do
	if(idg.eq.2)write(*,*)'pass ezs in field'
!        call filter(ez)                                                                                                                              
!        goto 100                                                                                                                                     

        if(iske==1 .or. ieqmo814==1)then
           call jie814(ip,n)
           call drdt814(ip)
           call eqmo814(ip)
        end if
	if(idg.eq.2)write(*,*)'pass eqmo814 in field'
        if(iske==0 .and. ieqmo814==0)call eqmo(ip)
	if(idg.eq.2)write(*,*)'pass eqmo in field'

        call fltx(ez(:,:,:),1,0,1,0)
	if(idg.eq.2)write(*,*)'pass fltx in field'
!        call filter(ez(:,:,:))                                                                                                                       

 100         continue
        myrmsphi=0.
        do k=0,mykm-1
           do j=0,jm-1
              do i1=0,im-1
                 myrmsphi=myrmsphi+ez(i1,j,k)*ez(i1,j,k)
              enddo
           enddo
        enddo
        call MPI_ALLREDUCE(myrmsphi,rmsez(n),1, &
              MPI_REAL8,                               &
              MPI_SUM,TUBE_COMM,ierr)
        rmsez(n)=sqrt(rmsez(n)/(im*jm*km))

	if(idg.eq.2)write(*,*)'pass rmsez in field'
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine field
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine push_wrapper(n,ip)
         use gem_com
         use equil
	implicit none
	integer :: n,i,j,k,ip

            if(ip.eq.1.and.ision==1)call ppush(n)
            if(ip.eq.0.and.ision==1)call cpush(n)
            if(eprs>0.and.ip.eq.0.and.mod(n+1,nrst).eq.0)call rstpi
            if(idg.eq.1)write(*,*)'pass ppush'

            if(ip.eq.1.and.ishgk==1)call hpush(n)
            if(ip.eq.0.and.ishgk==1)call hcush(n)
            if(ip.eq.0.and.ishgk==1)call collh(n)
            if(eprs>0.and.ip.eq.0.and.mod(n+1,nrst).eq.0)call rstph

            if(ip.eq.1.and.isbgk==1)call bpush(n)
            if(ip.eq.0.and.isbgk==1)call bcush(n)
            if(ip.eq.0.and.isbgk==1)call collb(n)

            if(ip.eq.1.and.iscgk==1)call mpush(n)
            if(ip.eq.0.and.iscgk==1)call mcush(n)
!            if(ip.eq.0.and.iscgk==1)call collc(n)

            if(ip.eq.1.and.iske==1)call pint
            if(ip.eq.0.and.iske==1)call cint(n)

            if(ip.eq.1.and.ifluid==1)call pintef
            if(ip.eq.0.and.ifluid==1)call cintef(n)
            if(idg.eq.1)write(*,*)'pass pint'
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)        		
end subroutine push_wrapper
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine diagnose(n)
         use gem_com
         use equil
	implicit none
	integer :: n,i,j,k,m,ip

        if(idg.eq.1)write(*,*)'before modes2 in diagnose'  
        call modes2(phi,pmodehis,n)
        if(idg.eq.1)write(*,*)'pass modes2'  
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)          
        call yveck(phi(:,:,:),n)
        call yveck1(phi(:,:,:),n)
        call yveck2(delte(:,:,:),n)

        call mdampd(phi(:,:,:),mdhis)
        call mdampd(apar(:,:,:),mdhisa)
        call fltd(dti(:,:,:))
        call fltd(dtc(:,:,:))
        dne = dene
        call fltd(dne(:,:,:))
        call mdampd(dti(:,:,:),mdhisb)
        call mdampd(dtc(:,:,:),mdhisc)
        call mdampd(dne(:,:,:),mdhisd)

        if(idg.eq.1)write(*,*)'pass yvec'    
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)        

end subroutine diagnose
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine reporter(n)
         use gem_com
         use equil
	implicit none
	integer :: n,i,j,k,m,ip

        if(mod(n,xnplt).eq.0) then
           call spec(n)
        endif

      if(n==(nm-1) .and. gclr==kmx/2 .and. tclr==0)then
         do i = 0,7
            j = imx/8*i+imx/16
            do m = 0,nom-1
               write(16,13)i,m,omlk(m),real(phiom(j,0,0,m)),aimag(phiom(j,0,0,m)),abs(phiom(j,0,0,m)) 
            end do
         end do
      end if
 13            format(1x,i3,1x,i4,4(2x,e10.3))

!        write(16,13)n,nopz,noen
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!     save particle arrays for restart if iput=1...
!     do this before the code crashes due to graphics problems
        if((iput.eq.1).and.mod(n+1,500).eq.0)call restart(2,n)

!     periodically make output for plots
        call outd(n)
        if(idg.eq.1)write(*,*)'pass outd'

end subroutine reporter
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fltx(u,isbl,idigit,ism,isdamp)   
      use gem_com
      use fft_wrapper
      implicit none
      INTEGER :: i,j,k,k1,l,m,n,ifirst,nstep,ip,INFO,isbl,idigit,NNI1
      INTEGER :: l1,m1,myk,myj,izonal,ix,ikx,ism,isdamp
      real(8) :: u(0:imx,0:jmx,0:1)
      COMPLEX(8) :: v(0:imx-1),u1(0:imx-1)
      COMPLEX(8),dimension(:,:,:),allocatable :: temp3dxy,tmp3d,t3d
      real(8) :: filterx(0:imx-1,0:jmx-1),gr(0:imx-1),gi(0:imx-1)
      real :: kx,ky,kx0,th,shat,sgny,dum,dum1
      real(8) ::  DW1,DW2
      parameter(NNI1 = 3)
      parameter(DW1 = 0.5)
      parameter(DW2 = -1./6.)

      allocate(temp3dxy(0:imx-1,0:jmx-1,0:1),tmp3d(0:imx,0:jmx-1,0:1),t3d(0:imx,0:jmx-1,0:1))

      if(idg.eq.3)write(*,*)'enter fltx'

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
            if(1 .le. j .and. j .le. mstart)filterx(i,j) = 0.
            if(((jmx-mstart) .le. j) .and. (j .le. (jmx-1)))filterx(i,j) = 0.
            if(abs(ky)>kycut)filterx(i,j) = 0.
            if(kx>kxcut)filterx(i,j) = 0.
            if(onemd==1.and.m1.ne.mlk)filterx(i,j) = 0.
!            if(m1>(jcnt-1)/2)filterx(i,j) = 0.
            if(izon==0 .and. j==0)filterx(i,j) = 0.
            if(izon==1 .and. j==0 .and. i<nzcrt)filterx(i,j) = 0.
         end do
      end do
      
      if(idg.eq.3)write(*,*)'pass filter fltx'         
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
            call ccfft('y',-1,jmx,1.d0,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3dxy(i,j,k) = tmpy(j)   !rho(ky,x)
            end do
         end do
      enddo

      if(idigit==1)then
         tmp3d = 0.
         t3d(0:imx-1,:,:) = temp3dxy(0:imx-1,:,:)
         t3d(imx,:,:) = 0.
         do l = 1,NNI1
            do k = 0,1
               do j = 0,jmx-1
                  do i = 1,imx-1
                     tmp3d(i,j,k) = (DW1*(t3d(i-1,j,k)+t3d(i+1,j,k))+t3d(i,j,k))/(1.+2*DW1)
                  end do
               end do
            end do
            t3d = tmp3d
            do k = 0,1
               do j = 0,jmx-1
                  do i = 1,imx-1
                     tmp3d(i,j,k) = (DW2*(t3d(i-1,j,k)+t3d(i+1,j,k))+t3d(i,j,k))/(1.+2*DW2)
                  end do
               end do
            end do
            t3d = tmp3d
         end do
         temp3dxy(0:imx-1,:,:) = t3d(0:imx-1,:,:)
         goto 200
      end if
      if(idg.eq.3)write(*,*)'pass idigit fltx'

      do k = 0,1
         do j = 0,jmx-1
            gr(:) = real(temp3dxy(:,j,k))
            gi(:) = aimag(temp3dxy(:,j,k))
       call dsinf(0,gr,1,0,1,0,imx*2,1,1.d0,aux1,50000,aux2,20000)
       call dsinf(0,gi,1,0,1,0,imx*2,1,1.d0,aux1,50000,aux2,20000)
            if(j==0.and.ism==0)then
               u1(:) = cmplx(gr(:),gi(:))
               call gam(u1(:),v(:),isdamp)
               gr(:) = real(v(:))
               gi(:) = aimag(v(:))
            end if
            gr(:) = gr(:)*filterx(:,j)
            gi(:) = gi(:)*filterx(:,j)
       call dsinf(0,gr,1,0,1,0,imx*2,1,1.d0,aux1,50000,aux2,20000)
       call dsinf(0,gi,1,0,1,0,imx*2,1,1.d0,aux1,50000,aux2,20000)
            temp3dxy(:,j,k) = cmplx(gr(:),gi(:))
         end do
      end do
      if(idg.eq.3)write(*,*)'pass dsin fltx'

200  continue
      temp3dxy(:,:,:)=temp3dxy(:,:,:)/jmx
      if(isbl==1)call filtbl(temp3dxy(0:imx-1,0:jmx-1,0:1))

      do k=0,mykm
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3dxy(i,j,k)
            end do
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
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
      return

!      call filter(dphidt(:,:,:))
      do k = 0,1
         dum = 0.
         dum1 = 0.
         do j=0,jm-1
            do i = 0,im-1
               dum = dum+u(i,j,k)*jac(i,k)
               dum1 = dum1+jac(i,k)
            end do
         end do
         dum = dum/dum1
         u(:,:,k) = u(:,:,k)-dum
      end do
      if(idg.eq.3)write(*,*)'end of fltx'

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine fltd(u)   
        use gem_com
      use fft_wrapper
      implicit none
      INTEGER :: i,j,k,k1,l,m,n,ifirst,nstep,ip,INFO,isbl,ism
      INTEGER :: l1,m1,myk,myj,ix,ikx,ierror
      real(8) :: u(0:imx,0:jmx,0:1)
      COMPLEX(8) :: temp3dxy(0:imx-1,0:jmx-1,0:1),v(0:imx-1),u1(0:imx-1)
      real(8) :: filterx(0:imx-1,0:jmx-1),gr(0:imx-1),gi(0:imx-1)
      real(8) :: kx,ky,sgny,dum,dum1,b2,myuav(0:imx-1)
      real(8) :: myjaca(0:imx-1),jaca(0:imx-1),uav(0:imx-1)

      xshape=1.0
      yshape=1.0

      do j = 0,jmx-1
         do i = 0,imx-1
            if(j.ge.(jm/2+1)) then
               m1=jm-j
               sgny=-1.d0
            else
               m1=j
               sgny=1.d0
            endif
            ky=sgny*2.d0*pi*float(m1)/ly
            kx = i*pi/lx
            b2 = xshape**2*kx**2+yshape**2*ky**2

            filterx(i,j) = 2.d0/imx*exp(-b2**2)   !be careful with always filtering with hyper-Gaussian

            if(abs(ky)>kycut)filterx(i,j) = 0.d0
            if(kx>kxcut)filterx(i,j) = 0.d0
            if(j==0.and.i<=nzcrt)filterx(i,j) = 0.d0
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
            call ccfft('y',-1,jmx,1.d0,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3dxy(i,j,k) = tmpy(j)   !rho(ky,x)
            end do
         end do
      enddo

      do k = 0,1
         do j = 0,jmx-1
            gr(:) = real(temp3dxy(:,j,k))
            gi(:) = aimag(temp3dxy(:,j,k))
       call dsinf(0,gr,1,0,1,0,imx*2,1,1.d0,aux1,50000,aux2,20000)
       call dsinf(0,gi,1,0,1,0,imx*2,1,1.d0,aux1,50000,aux2,20000)
       if(j==0)then
          u1(:) = cmplx(gr(:),gi(:))
          call gam1(u1(:),v(:))
          gr(:) = real(v(:))
          gi(:) = aimag(v(:))
       end if
       gr(:) = gr(:)*filterx(:,j)
       gi(:) = gi(:)*filterx(:,j)
       call dsinf(0,gr,1,0,1,0,imx*2,1,1.d0,aux1,50000,aux2,20000)
       call dsinf(0,gi,1,0,1,0,imx*2,1,1.d0,aux1,50000,aux2,20000)
            temp3dxy(:,j,k) = cmplx(gr(:),gi(:))
         end do
      end do

      temp3dxy(:,:,:)=temp3dxy(:,:,:)/jmx

      do k=0,mykm
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3dxy(i,j,k)
            end do
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
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
      subroutine flty(u)   
      use gem_com
      use fft_wrapper
      implicit none
      INTEGER :: i,j,k,l,NNI1
      real(8) :: u(0:imx,0:jmx,0:1)
      real(8) :: tmp3d(0:imx,0:jmx,0:1),t3d(0:imx,0:jmx,0:1)
      real :: dum,dum1
      real(8) ::  DW1,DW2
      parameter(NNI1 = 3)
      parameter(DW1 = 0.5)
      parameter(DW2 = -1./6.)

      t3d = u
      do l = 1,NNI1
         do k = 0,1
            do i = 0,imx
               do j = 0,jmx
                  if(j>0)dum = t3d(i,j-1,k)
                  if(j==0)dum = t3d(i,jmx-1,k)
                  if(j<jmx)dum1 = t3d(i,j+1,k)
                  if(j==jmx)dum1 = t3d(i,1,k)
                  tmp3d(i,j,k) = (DW1*(dum+dum1)+t3d(i,j,k))/(1.+2*DW1)
               end do
            end do
         end do
         t3d = tmp3d
         do k = 0,1
            do i = 0,imx
               do j = 0,jmx
                  if(j>0)dum = t3d(i,j-1,k)
                  if(j==0)dum = t3d(i,jmx-1,k)
                  if(j<jmx)dum1 = t3d(i,j+1,k)
                  if(j==jmx)dum1 = t3d(i,1,k)
                  tmp3d(i,j,k) = (DW2*(dum+dum1)+t3d(i,j,k))/(1.+2*DW2)
               end do
            end do
         end do
         t3d = tmp3d
      end do
      u = t3d

      call enfxy(u(:,:,:))
      call enfz(u)
!      call filter(dphidt(:,:,:))
      do k = 0,1
         dum = 0.
         dum1 = 0.
         do j=0,jm-1
            do i = 0,im-1
               dum = dum+u(i,j,k)*jac(i,k)
               dum1 = dum1+jac(i,k)
            end do
         end do
         dum = dum/dum1
         u(:,:,k) = u(:,:,k)-dum
      end do

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
      complex(8) :: v(0:imx-1,0:jcnt-1,0:1),myv(0:imx-1,0:jcnt-1,0:1)
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
            call ccfft('y',-1,jmx,1.d0,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3d(i,j,k) = tmpy(j)
            end do
         end do
      enddo

      do m = 0,jcnt-1
         myv(:,m,:) = temp3d(:,jft(m),:)
      end do

      cnt = 2*jcnt*imx
      call mpi_allreduce(myv,v,cnt,MPI_DOUBLE_COMPLEX,mpi_sum, &
                              grid_comm,ierr)

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hpush(n)

      use gem_com
      use equil
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum,dum1
      INTEGER :: m,i,j,k,l,n,mynopz
      INTEGER :: np_old,np_new
      REAL(8) :: rhog,vfac,kap,vpar,pidum,kaphp
      REAL(8) :: b,th,r,enerb,cost,sint,qr,laps,sz,ter,bstar
      REAL(8) :: xt,xs,yt,xdot,ydot,xdt,ydt,zdot,pzdot,edot,pzd0,vp0
      REAL(8) :: pzd1,dpzdt,dthdt,dpsidt,grdgtp,psip2p
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,kaphip,nhip
      REAL(8) :: grp,gxdgyp,rhox(4),rhoy(4),psp,pzp,cvbzp
      REAL(8) :: v,dist,s,s0,rpsi,drdpsi,dfdpsi

      pidum = 1./(pi*2)**1.5*vwidth**3
      mynopz = 0
      do m=1,mm(3)
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
         grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
                 +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1) 
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
         psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         psp = wx0*psi(i)+wx1*psi(i+1)        
         kaphp = wx0*kaphi(i)+wx1*kaphi(i+1)
         nhip = wx0*nhi(i)+wx1*nhi(i+1)        

         b=1.-tor+tor*bfldp
         pzp = mims(3)*uh2(m)/b*fp/br0-q(3)*psp/br0

         rhog=sqrt(2.*b*muh(m)*mims(3))/(q(3)*b)*iflrh

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

            include "hpushli.h"
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
         vfac = 0.5*(mims(3)*uh2(m)**2 + 2.*muh(m)*b)
         v = sqrt(2.*vfac/mims(3))

         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))

         vpar = uh2(m)
         bstar = b*(1+mims(3)*vpar/(b*q(3))*fp/(b**2*radiusp)*(psip2p*grp**2/radiusp+cvbzp))
         enerb=(muh(m)+mims(3)*vpar*vpar/b)/q(3)*b/bstar*tor
         
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1*b/bstar
         dum = 1/(b*b)*fp/radiusp*cvbzp
         xdot = vxdum*nonlinh -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         xdt = -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1/b+vpar/b*delbyp)*dum1*b/bstar*nonlinh &
             +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0  &
             -mims(3)*vpar**2/(q(3)*bstar*b)*(psip2p*grp**2/radiusp+cvbzp)*lr0/(radiusp*q0) 
         ydt = tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0  &
             -mims(3)*vpar**2/(q(3)*bstar*b)*(psip2p*grp**2/radiusp+cvbzp)*lr0/(radiusp*q0) 
         zdot =  (vpar)*b/bstar &
                 *(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         pzd0 = tor*(-muh(m)/mims(3)/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar
         pzd1 = q(3)/mims(3)*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp 
         pzdot = pzd0 + pzd1*iparah

         edot = q(3)*(xdt*exp1+ydt*eyp+vpar*ezp )

         xh3(m) = xh2(m) + 0.5*dt*xdot
         yh3(m) = yh2(m) + 0.5*dt*ydot
         zh3(m) = zh2(m) + 0.5*dt*zdot
         uh3(m) = uh2(m) + 0.5*dt*pzdot

         dist = nh/(v**3+vi**3)*nhip
         dum = v**3+vi**3
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1*b/bstar

         wh3(m)=wh2(m) + 0.5*dt*(kaphp*vxdum*dist &
                       + edot*3*v/mims(3)*dist/(v**3+vi**3))*dum-0.5*dt*gambm*wh2(m)
         
!         if(xh3(m)>lx .or. xh3(m)<0.)wh3(m) = 0.

         go to 333
         if(abs(pzp-pzh(m))>pzcrit(3).or.abs(vfac-ekh(m))>encrit(3))then
            mynopz = mynopz+1
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
            uh3(m) = vih(m) !sqrt(2/mims(3)*abs(ekh(m)-muh(m)*b))
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
         i = max(i,0)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         qr = wx0*sf(i)+wx1*sf(i+1)
         yh3(m)=dmod(yh3(m)-laps*2*pi*qr*lr0/q0*dsign(1.d0,q0)+800.*ly,ly)
         if(xh3(m)>lx)then
            xh3(m) = lx-1.e-8
            zh3(m)=lz-zh3(m)
            zh2(m)=zh3(m)
            xh2(m) = xh3(m)
            yh3(m) = ly*ran2(iseed)
            yh2(m) = yh3(m)
!            wh2(m) = 0.
!            wh3(m) = 0.
         end if
         if(xh3(m)<0.)then
            xh3(m) = 1.e-8
            zh3(m)=lz-zh3(m)
            zh2(m)=zh3(m)
            xh2(m) = xh3(m)
            yh3(m) = ly*ran2(iseed)
            yh2(m) = yh3(m)
!            wh2(m) = 0.
!            wh3(m) = 0.
         end if
         zh3(m)=dmod(zh3(m)+8.*lz,lz)
         xh3(m)=dmod(xh3(m)+8.*lx,lx)         
         xh3(m) = min(xh3(m),lx-1.0e-8)
         yh3(m) = min(yh3(m),ly-1.0e-8)
         zh3(m) = min(zh3(m),lz-1.0e-2)
         
      enddo

!      call MPI_ALLREDUCE(mynopz,nopz,1,MPI_integer, &
!          MPI_SUM, MPI_COMM_WORLD,ierr)

      np_old=mm(3)
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
      call pmove(muh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mui,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call pmove(xih,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0h,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(pzh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ekh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(vih,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(index,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call end_pmove(ierr)
      mm(3)=np_new

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
      REAL(8) :: vfac,vpar,vxdum,dum,xdot,ydot,xdt,ydt,zdot,pzdot,edot,pzd0,vp0
      REAL(8) :: pzd1,dpzdt,dthdt,dpsidt,grdgtp,psip2p
      REAL(8) :: rhog,xt,yt,zt,kap,xs,pidum,dum1,kaphp
      REAL(8) :: b,th,r,enerb,cost,sint,qr,laps,sz,ter,bstar
      REAL(8) :: myke,mypfl,myavewh
      REAL(8) :: myefl,mynos
      REAL(8) :: ketemp,pfltemp
      REAL(8) :: efltemp,nostemp
      REAL(8) :: sbuf(10),rbuf(10)
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,kaphip,nhip
      REAL(8) :: grp,gxdgyp,rhox(4),rhoy(4),psp,pzp,cvbzp
      REAL(8) :: v,dist,s,s0,rpsi,drdpsi,dfdpsi

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

      do m=1,mm(3)
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
         grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
                 +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1) 
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
         psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         psp = wx0*psi(i)+wx1*psi(i+1)        
         ter = wx0*t0i(i)+wx1*t0i(i+1)        
         kapt(1) = wx0*capti(i)+wx1*capti(i+1)        
         b=1.-tor+tor*bfldp
         pzp = mims(3)*uh3(m)/b*fp/br0-q(3)*psp/br0
         kaphp = wx0*kaphi(i)+wx1*kaphi(i+1)
         nhip = wx0*nhi(i)+wx1*nhi(i+1)        

         rhog=sqrt(2.*b*muh(m)*mims(3))/(q(3)*b)*iflrh

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

            include "hcushli.h"
 200     continue

         exp1=exp1/4.
         eyp=eyp/4.
         ezp=ezp/4.
         delbxp=delbxp/4.
         delbyp=delbyp/4.
         dpdzp = dpdzp/4.
         dadzp = dadzp/4.
         aparp = aparp/4.


         vfac = 0.5*(mims(3)*uh3(m)**2 + 2.*muh(m)*b)
         v = sqrt(2.*vfac/mims(3))

         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))

         vpar = uh3(m)
         bstar = b*(1+mims(3)*vpar/(b*q(3))*fp/(b**2*radiusp)*(psip2p*grp**2/radiusp+cvbzp))
         enerb=(muh(m)+mims(3)*vpar*vpar/b)/q(3)*b/bstar*tor

         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1*b/bstar
         dum = 1/(b*b)*fp/radiusp*cvbzp
         xdot = vxdum*nonlinh -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         xdt = -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1/b+vpar/b*delbyp)*dum1*b/bstar*nonlinh     &
             +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
             -mims(3)*vpar**2/(q(3)*bstar*b)*(psip2p*grp**2/radiusp+cvbzp)*lr0/(radiusp*q0) 
         ydt = tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
             -mims(3)*vpar**2/(q(3)*bstar*b)*(psip2p*grp**2/radiusp+cvbzp)*lr0/(radiusp*q0) 
         zdot =  (vpar)*b/bstar &
                 *(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         pzd0 = tor*(-muh(m)/mims(3)/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar
         pzd1 = q(3)/mims(3)*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp 
         pzdot = pzd0 + pzd1*iparah

         edot = q(3)*(xdot*exp1+ydot*eyp+vpar*ezp )

         xh3(m) = xh2(m) + dt*xdot
         yh3(m) = yh2(m) + dt*ydot
         zh3(m) = zh2(m) + dt*zdot
         uh3(m) = uh2(m) + dt*pzdot

         if(uh2(m)<0..and.uh3(m)>0.)index(m) = 0.
         if(uh2(m)>0..and.uh3(m)<0.)index(m) = 0.
         if(zh3(m)>lz)index(m) = 1.
         if(zh3(m)<0.)index(m) = 1.

         dist = nh/(v**3+vi**3)*nhip
         dum = v**3+vi**3
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1*b/bstar

         w3old = wh3(m)
         wh3(m)=wh2(m) + dt*(kaphp*vxdum*dist &
                       + edot*3*v/mims(3)*dist/(v**3+vi**3))*dum-dt*gambm*wh2(m)

         if(abs(wh3(m)).gt.nh.and.nonlinh==1)then
!            wh3(m) = 0.
!            wh2(m) = 0.
         end if


         laps=anint((zh3(m)/lz)-.5)*(1-peritr)
         r=xh3(m)-0.5*lx+lr0
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         i = max(i,0)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         qr = wx0*sf(i)+wx1*sf(i+1)
         yh3(m)=dmod(yh3(m)-laps*2*pi*qr*lr0/q0*dsign(1.d0,q0)+800.*ly,ly)
         if(xh3(m)>lx)then
            xh3(m) = lx-1.e-8
            zh3(m)=lz-zh3(m)
            zh2(m)=zh3(m)
            xh2(m) = xh3(m)
            yh3(m) = ly*ran2(iseed)
            yh2(m) = yh3(m)
!            wh2(m) = 0.
!            wh3(m) = 0.
         end if
         if(xh3(m)<0.)then
            xh3(m) = 1.e-8
            zh3(m)=lz-zh3(m)
            zh2(m)=zh3(m)
            xh2(m) = xh3(m)
            yh3(m) = ly*ran2(iseed)
            yh2(m) = yh3(m)
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
         myke=myke + vfac !*wh3(m)
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
      efl(3,n)=efltemp/( float(tmm(1)) )
      ke(3,n)=ketemp/(float(tmm(1)))
      np_old=mm(3) 
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
      call pmove(muh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mui,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xih,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0h,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(pzh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ekh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(vih,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(index,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call end_pmove(ierr)
      mm(3)=np_new
!     write(*,*)MyId,mm(1)

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine loadh

      use gem_com
      use equil
      implicit none
      INTEGER :: i,j,k,m,idum,ns,m1,nvgrd=100
      INTEGER :: np_old,np_new
      REAL(8) :: vpar,vperp2,r,qr,th,b,cost,ter,x
      REAL(8) :: avgv,myavgv,avgw,myavgw
      real(8) :: dumx,dumy,dumz,jacp
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,psp
      REAL(8) :: grp,gxdgyp,zoldp
      REAL(8) :: wx0,wx1,wz0,wz1
      REAL(8) :: vx,vy,v,vg(0:10000),dv,f1,f2

      cnt=int(tmm(1)/numprocs)
!   write(*,*)'in loader cnt, mm(1)=',cnt,mm(1)

      myavgv=0.
      avgv=0.
      avgw = 0.
      myavgw = 0.

!find a set of v grids for slowing-down distribution
      dv = (vc-vmin)/nvgrd
      vg(0) = vmin
      vg(1) = vmin+dv
      do i=1,1000000
         f1 = 1./(vg(i-1)**3+vi**3)
         f2 = 1./(vg(i)**3+vi**3)
         vg(i+1) = vg(i)+(vg(i)-vg(i-1))*vg(i-1)**2*f1/(vg(i)**2*f2)
         if(vg(i+1)>vc)then
            nvgrd = i
            goto 150
         end if
      end do
 150  continue
!compute x=int v^4/(v^3+vi^3)
      x = 0.
      dv = (vc-vmin)/nvgrd
      do i = 0, nvgrd-1
         v = vmin+(i+0.5)*dv
         x = x+v**4/(v**3+vi**3)*dv
      end do
      if(myid==0)write(*,*)'nvgrd, x= ', nvgrd,x*pi*mims(3)*nh

      mm(3)=int(tmm(1)/numprocs)
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
         v = vg(int(ran2(iseed)*nvgrd))+(ran2(iseed)-0.5)*dv
         vpar = 2.*v*(ran2(iseed)-0.5)

         if((ran2(iseed)<(0.5*jacp/jacmax)).and.(v<vc).and.(v>vmin))then
         m = m+1
         if(m>mm(3))goto 170
         xh2(m)=min(dumx,lx-1.d-8)
         yh2(m)=min(dumy,ly-1.d-8)

         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1-wz0
         zh2(m) = wz0*zfnth(k)+wz1*zfnth(k+1)
         zh2(m)=min(zh2(m),lz-1.d-8)

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
         fp = wx0*f(i)+wx1*f(i+1)        
!         ter = tets(3) !wx0*t0i(i)+wx1*t0i(i+1)
         b=1.-tor+tor*bfldp

         uh2(m)=vpar
         muh(m)=0.5*mims(3)*(v**2-vpar**2)/b
         muh2(m) = muh(m)
         ekh(m) = muh(m)*b+0.5*mims(3)*uh2(m)**2
         pzh(m) = mims(3)*uh2(m)/b*fp/br0-q(3)*psp/br0
         z0h(m) = zh2(m)
         xih(m) = xh2(m)
         vih(m) = uh2(m)
         mui(m) = muh(m)
         myavgv=myavgv+uh2(m)

!    LINEAR: perturb w(m) to get linear growth...
         wh2(m)=2.*amp*nh/cv*(revers(MyId*cnt+m,13)-0.5) !(ran2(iseed) - 0.5 )
         myavgw=myavgw+wh2(m)
         index(m) = 1
         end if
 160  continue
 170  continue
      myavgw = myavgw/mm(3)
!      write(*,*)'myid ', myid,x2(10),y2(20),u2(20),w2(20)
!    subtract off avg. u...
      call MPI_ALLREDUCE(myavgv,avgv,1, &
          MPI_REAL8, &
          MPI_SUM,MPI_COMM_WORLD,ierr)
      if(idg.eq.1)write(*,*)'all reduce'
      avgv=avgv/float(tmm(1))
      do 180 m=1,mm(3)
         uh2(m)=uh2(m)-avgv
         xh3(m)=xh2(m)
         yh3(m)=yh2(m)
         zh3(m)=zh2(m)
         uh3(m)=uh2(m)
         wh3(m)=wh2(m)
 180  continue

      np_old=mm(3)
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
      call pmove(muh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mui,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xih,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0h,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(pzh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ekh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(vih,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(index,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
!     
      call end_pmove(ierr)
      mm(3)=np_new
!     write(*,*)MyId,j,mm(3)

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pintef
      use gem_com
      use equil
      implicit none

      integer :: i,j,k,ip
      real*8 :: lbfr(0:imx,0:jmx,2)
      real*8 :: lbfs(0:imx,0:jmx,2)
      real*8 :: rbfr(0:imx,0:jmx,2)
      real*8 :: rbfs(0:imx,0:jmx,2)
      real*8 :: dum,dumphi,dum1,dum2,denel,dener,aparl,aparr
      real*8 :: dmnl1(0:imx,0:jmx,0:1),dmnl2(0:imx,0:jmx,0:1),dmnl3(0:imx,0:jmx,0:1),dmnl4(0:imx,0:jmx,0:1)
      real*8 :: ux(0:imx,0:jmx,0:1),uy(0:imx,0:jmx,0:1)

      phis = phi
      denes = dene
      apars = apar
      deltes = delte

      do i = 0,im
         do j = 0,jm
            rbfs(i,j,1)=upar(i,j,0)/bmag(i,0)
            rbfs(i,j,2) = dene(i,j,0)
         end do   
      end do
      call MPI_SENDRECV(rbfs(0,0,1),(imx+1)*(jmx+1)*2, &
          MPI_REAL8,rngbr,404, &
          lbfr(0,0,1),(imx+1)*(jmx+1)*2, &
          MPI_REAL8,lngbr,404, &
          tube_comm,stat,ierr)
      do i = 0,im
         do j = 0,jm
            lbfs(i,j,1)=upar(i,j,1)/bmag(i,1)
            lbfs(i,j,2) = dene(i,j,1)
         end do   
      end do
      call MPI_SENDRECV(lbfs(0,0,1),(imx+1)*(jmx+1)*2, &
          MPI_REAL8,lngbr,405, &
          rbfr(0,0,1),(imx+1)*(jmx+1)*2, &
          MPI_REAL8,rngbr,405, &
          tube_comm,stat,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(iske==0)then
         call gradu(dene*ex,ux,uy)
         dmnl1 = uy
         call gradu(dene*ey,ux,uy)
         dmnl2 = ux
         call gradu(delby*upar,ux,uy)
         dmnl3 = uy
         call gradu(delbx*upar,ux,uy)
         dmnl4 = ux
      end if

      if(iske==1)then
         call gradu(denek*ex,ux,uy)
         dmnl1 = uy
         call gradu(denek*ey,ux,uy)
         dmnl2 = ux
         call gradu(delby*upark,ux,uy)
         dmnl3 = uy
         call gradu(delbx*upark,ux,uy)
         dmnl4 = ux
      end if

!dene_p for PC, denes for RK
      do i = 0,im-1
         do j = 0,jm-1
            dum =  weightp(i)*lbfr(i,jpl(i,j),1) &
       	+weightpn(i)*lbfr(i,jpn(i,j),1)
            ddedt(i,j,0) = -(upar(i,j,1)/bmag(i,1)-dum) &
                /(2.*dz)*bmag(i,0)*bdgrzn(i,0)*gn0e(i) &
                +gn0e(i)*gcpne(i)*ey(i,j,0)/bmag(i,0)*bdgxcgy(i,0)*icmprs &
                +2.*(cfx(i,0)*(-dnedx(i,j,0)*gt0e(i)/gn0e(i)*ispre-ex(i,j,0)) &
                +cfy(i,0)*(-dnedy(i,j,0)*gt0e(i)/gn0e(i)*ispre-ey(i,j,0)))/br0*icmprs*gn0e(i) &
                -(delby(i,j,0)*gnuobx(i,0)+delbx(i,j,0)*gnuoby(i,0))    &
                -bdgxcgy(i,0)*(-dmnl1(i,j,0)+dmnl2(i,j,0))/bmag(i,0)*nonline &
                -bdgxcgy(i,0)*(dmnl3(i,j,0)+dmnl4(i,j,0))*gn0e(i)/bmag(i,0)*nonline*iflut

            dene(i,j,0) = denes(i,j,0)+0.5*dt*ddedt(i,j,0)-0.5*dt*gamdne*denes(i,j,0)

            dum = weightm(i)*rbfr(i,jmi(i,j),1) &
                +weightmn(i)*rbfr(i,jmn(i,j),1)
            ddedt(i,j,1) = -(dum-upar(i,j,0)/bmag(i,0)) &
                /(2.*dz)*bmag(i,1)*bdgrzn(i,1)*gn0e(i) &
                +gn0e(i)*gcpne(i)*ey(i,j,1)/bmag(i,1)*bdgxcgy(i,1)*icmprs &
                +2.*(cfx(i,1)*(-dnedx(i,j,1)*gt0e(i)/gn0e(i)*ispre-ex(i,j,1)) &
                +cfy(i,1)*(-dnedy(i,j,1)*gt0e(i)/gn0e(i)*ispre-ey(i,j,1)))/br0*icmprs*gn0e(i) &
                -(delby(i,j,1)*gnuobx(i,1)+delbx(i,j,1)*gnuoby(i,1))    &
                -bdgxcgy(i,1)*(-dmnl1(i,j,1)+dmnl2(i,j,1))/bmag(i,1)*nonline &
                -bdgxcgy(i,1)*(dmnl3(i,j,1)+dmnl4(i,j,1))*gn0e(i)/bmag(i,1)*nonline*iflut

            dene(i,j,1) = denes(i,j,1)+0.5*dt*ddedt(i,j,1)-0.5*dt*gamdne*denes(i,j,1)
         end do
      end do

      do i = 0,im
         do j = 0,jm
            rbfs(i,j,1)=phi(i,j,0)
            rbfs(i,j,2) = apar(i,j,0)            
         end do   
      end do
      call MPI_SENDRECV(rbfs(0,0,1),(imx+1)*(jmx+1)*2, &
          MPI_REAL8,rngbr,406, &
          lbfr(0,0,1),(imx+1)*(jmx+1)*2, &
          MPI_REAL8,lngbr,406, &
          tube_comm,stat,ierr)
      do i = 0,im
         do j = 0,jm
            lbfs(i,j,1)=phi(i,j,1)
            lbfs(i,j,2) = apar(i,j,1)
         end do   
      end do
      call MPI_SENDRECV(lbfs(0,0,1),(imx+1)*(jmx+1)*2, &
          MPI_REAL8,lngbr,407, &
          rbfr(0,0,1),(imx+1)*(jmx+1)*2, &
          MPI_REAL8,rngbr,407, &
          tube_comm,stat,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!apar_p for PC, apars for RK
      do i = 0,im-1
         do j = 0,jm-1
            dum =  weightp(i)*lbfr(i,jpl(i,j),1) &
       	+weightpn(i)*lbfr(i,jpn(i,j),1)
            dumphi = (phi(i,j,1)-dum)/(2.*dz)*bdgrzn(i,0)
            apar(i,j,0) = apars(i,j,0)-0.5*dt*(dumphi+ez(i,j,0))-0.5*dt*gamapa*apars(i,j,0)
      
            dum = weightm(i)*rbfr(i,jmi(i,j),1) &
       	+weightmn(i)*rbfr(i,jmn(i,j),1)
            dumphi = (dum-phi(i,j,0))/(2.*dz)*bdgrzn(i,1)
            apar(i,j,1) = apars(i,j,1) &
                -0.5*dt*(dumphi+ez(i,j,1))-0.5*dt*gamapa*apars(i,j,1)
         end do
      end do

      call enfxy(apar(:,:,:))
      call enfxy(dene(:,:,:))
      call enfz(apar(:,:,:))
      call enfz(dene(:,:,:))

      do i = 0,im-1
         do j = 0,jm-1
            delte(i,j,0) = deltes(i,j,0)+0.5*dt*gt0e(i)*gcpte(i)*ey(i,j,0)*bdgxcgy(i,0)/bmag(i,0) &
                           -0.5*dt*gamte*deltes(i,j,0)
            delte(i,j,1) = deltes(i,j,1)+0.5*dt*gt0e(i)*gcpte(i)*ey(i,j,1)*bdgxcgy(i,1)/bmag(i,1) &
                           -0.5*dt*gamte*deltes(i,j,1)
         end do
      end do
      call enfxy(delte(:,:,:))
      call enfz(delte(:,:,:))


      if(isgkm==0)return
      do i = 0,im-1
         do j = 0,jm-1
            phi(i,j,0) = phis(i,j,0)+0.5*dt*dphidt(i,j,0)-0.5*dt*gamphi*phis(i,j,0)
            phi(i,j,1) = phis(i,j,1)+0.5*dt*dphidt(i,j,1)-0.5*dt*gamphi*phis(i,j,1)
         end do
      end do
      call enfxy(phi(:,:,:))
      call enfz(phi(:,:,:))

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cintef(n)
      use gem_com
      use equil
      implicit none

      integer :: i,j,k,n,ip
      real*8 :: tmp(0:imx,0:jmx,0:1)
      real*8 :: lbfr(0:imx,0:jmx,2)
      real*8 :: lbfs(0:imx,0:jmx,2)
      real*8 :: rbfr(0:imx,0:jmx,2)
      real*8 :: rbfs(0:imx,0:jmx,2)
      real*8 :: tmpa(0:imx,0:jmx,0:1),tmpd(0:imx,0:jmx,0:1)
      real*8 :: dum,dumphi,dum1,dum2,denel,dener,aparl,aparr
      REAL*8 :: myrmsapa
      real*8 :: dmnl1(0:imx,0:jmx,0:1),dmnl2(0:imx,0:jmx,0:1),dmnl3(0:imx,0:jmx,0:1),dmnl4(0:imx,0:jmx,0:1)
      real*8 :: ux(0:imx,0:jmx,0:1),uy(0:imx,0:jmx,0:1)
      real(8) :: mydbr(0:imx-1),myjac(0:imx-1),v(0:imx-1)

      do i = 0,im
         do j = 0,jm
            rbfs(i,j,1)=upar(i,j,0)/bmag(i,0)
            rbfs(i,j,2) = dene(i,j,0)
         end do   
      end do
      call MPI_SENDRECV(rbfs(0,0,1),(imx+1)*(jmx+1)*2, &
          MPI_REAL8,rngbr,404, &
          lbfr(0,0,1),(imx+1)*(jmx+1)*2, &
          MPI_REAL8,lngbr,404, &
          tube_comm,stat,ierr)
      do i = 0,im
         do j = 0,jm
            lbfs(i,j,1)=upar(i,j,1)/bmag(i,1)
            lbfs(i,j,2) = dene(i,j,1)
         end do   
      end do
      call MPI_SENDRECV(lbfs(0,0,1),(imx+1)*(jmx+1)*2, &
          MPI_REAL8,lngbr,405, &
          rbfr(0,0,1),(imx+1)*(jmx+1)*2, &
          MPI_REAL8,rngbr,405, &
          tube_comm,stat,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(iske==0)then
         call gradu(dene*ex,ux,uy)
         dmnl1 = uy
         call gradu(dene*ey,ux,uy)
         dmnl2 = ux
         call gradu(delby*upar,ux,uy)
         dmnl3 = uy
         call gradu(delbx*upar,ux,uy)
         dmnl4 = ux
      end if

      if(iske==1)then
         call gradu(denek*ex,ux,uy)
         dmnl1 = uy
         call gradu(denek*ey,ux,uy)
         dmnl2 = ux
         call gradu(delby*upark,ux,uy)
         dmnl3 = uy
         call gradu(delbx*upark,ux,uy)
         dmnl4 = ux
      end if

      do i = 0,im-1
         do j = 0,jm-1
            dum =  weightp(i)*lbfr(i,jpl(i,j),1) &
       	+weightpn(i)*lbfr(i,jpn(i,j),1)
            tmp(i,j,0) = -(upar(i,j,1)/bmag(i,1)-dum) &
                /(2.*dz)*bmag(i,0)*bdgrzn(i,0)*gn0e(i) &
                +gn0e(i)*gcpne(i)*ey(i,j,0)/bmag(i,0)*bdgxcgy(i,0)*icmprs &
                +2.*(cfx(i,0)*(-dnedx(i,j,0)*gt0e(i)/gn0e(i)*ispre-ex(i,j,0)) &
                +cfy(i,0)*(-dnedy(i,j,0)*gt0e(i)/gn0e(i)*ispre-ey(i,j,0)))/br0*icmprs*gn0e(i) &
                -(delby(i,j,0)*gnuobx(i,0)+delbx(i,j,0)*gnuoby(i,0))    &
                -bdgxcgy(i,0)*(-dmnl1(i,j,0)+dmnl2(i,j,0))/bmag(i,0)*nonline &
                -bdgxcgy(i,0)*(dmnl3(i,j,0)+dmnl4(i,j,0))*gn0e(i)/bmag(i,0)*nonline*iflut

            tmpd(i,j,0) = denes(i,j,0)+dt*(ddedt(i,j,0)*0.+tmp(i,j,0))-dt*gamdne*denes(i,j,0)

            dum = weightm(i)*rbfr(i,jmi(i,j),1) &
                +weightmn(i)*rbfr(i,jmn(i,j),1)
            tmp(i,j,1) = -(dum-upar(i,j,0)/bmag(i,0)) &
                /(2.*dz)*bmag(i,1)*bdgrzn(i,1)*gn0e(i) &
                +gn0e(i)*gcpne(i)*ey(i,j,1)/bmag(i,1)*bdgxcgy(i,1)*icmprs &
                +2.*(cfx(i,1)*(-dnedx(i,j,1)*gt0e(i)/gn0e(i)*ispre-ex(i,j,1)) &
                +cfy(i,1)*(-dnedy(i,j,1)*gt0e(i)/gn0e(i)*ispre-ey(i,j,1)))/br0*icmprs*gn0e(i) &
                -(delby(i,j,1)*gnuobx(i,1)+delbx(i,j,1)*gnuoby(i,1))    &
                -bdgxcgy(i,1)*(-dmnl1(i,j,1)+dmnl2(i,j,1))/bmag(i,1)*nonline &
                -bdgxcgy(i,1)*(dmnl3(i,j,1)+dmnl4(i,j,1))*gn0e(i)/bmag(i,1)*nonline*iflut

            tmpd(i,j,1) = denes(i,j,1)+dt*(ddedt(i,j,1)*0.+tmp(i,j,1))-dt*gamdne*denes(i,j,1)
         end do
      end do

!phis for PC, phi for RK
      do i = 0,im
         do j = 0,jm
            rbfs(i,j,1)=0.5*(phi(i,j,0)+phi(i,j,0))
            rbfs(i,j,2) = apar(i,j,0)            
         end do   
      end do
      call MPI_SENDRECV(rbfs(0,0,1),(imx+1)*(jmx+1)*2, &
          MPI_REAL8,rngbr,406, &
          lbfr(0,0,1),(imx+1)*(jmx+1)*2, &
          MPI_REAL8,lngbr,406, &
          tube_comm,stat,ierr)
      do i = 0,im
         do j = 0,jm
            lbfs(i,j,1)=0.5*(phi(i,j,1)+phi(i,j,1))
            lbfs(i,j,2) = apar(i,j,1)
         end do   
      end do
      call MPI_SENDRECV(lbfs(0,0,1),(imx+1)*(jmx+1)*2, &
          MPI_REAL8,lngbr,407, &
          rbfr(0,0,1),(imx+1)*(jmx+1)*2, &
          MPI_REAL8,rngbr,407, &
          tube_comm,stat,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do i = 0,im-1
         do j = 0,jm-1
            dum =  weightp(i)*lbfr(i,jpl(i,j),1) &
       	+weightpn(i)*lbfr(i,jpn(i,j),1)
            dumphi = (0.5*(phi(i,j,1)+phi(i,j,1))-dum)/(2.*dz)*bdgrzn(i,0)
            tmpa(i,j,0) = apars(i,j,0)-dt*(dumphi+ez(i,j,0))-dt*gamapa*apars(i,j,0)
      
            dum = weightm(i)*rbfr(i,jmi(i,j),1) &
       	+weightmn(i)*rbfr(i,jmn(i,j),1)
            dumphi = (dum-0.5*(phi(i,j,0)+phi(i,j,0)))/(2.*dz)*bdgrzn(i,1)
            tmpa(i,j,1) = apars(i,j,1) &
                -dt*(dumphi+ez(i,j,1))-dt*gamapa*apars(i,j,1)
         end do
      end do
 100          continue

!            remove k_perp = 0 component from dene
      do k = 0,1
         dum = 0.
         do i = 0,im-1
            do j = 0,jm-1
               dum = dum+tmpd(i,j,k)
            end do
         end do
         dum = dum/float(im*jm)
         do i = 0,im-1
            do j = 0,jm-1
               tmpd(i,j,k) = tmpd(i,j,k)-dum
            end do
         end do
      end do

      apar(:,:,:) = tmpa(:,:,:)
      dene(:,:,:) = tmpd(:,:,:)
      
      call enfxy(apar(:,:,:))
      call enfxy(dene(:,:,:))
      call enfz(apar(:,:,:))
      call enfz(dene(:,:,:))
      call fltx(apar,0,0,0,1)
!      call fltx(dene,0,0,0,1)

!      dene_p = dene
!      call flty(dene_p)
!      dene = dene*(1-epsny)+dene_p*epsny

!      dene_p = dene
!      call fltx(dene_p,0,1,0,0)
!      dene = dene*(1-epsnx)+dene_p*epsnx
      dene_p = dene
      call filter(dene_p)
      dene = dene*(1-epsnz)+dene_p*epsnz

!      apar_p = apar
!      call flty(apar_p)
!      apar = apar*(1-epsay)+apar_p*epsay
      apar_p = apar
      call fltx(apar_p,0,1,0,0)
      apar = apar*(1-epsax)+apar_p*epsax
      apar_p = apar
      call filter(apar_p)
      apar = apar*(1-epsaz)+apar_p*epsaz

      myrmsapa = 0.
      rmsapa(n) = 0.
      do k=0,mykm-1
         do j=0,jm-1
            do i=0,im-1
               myrmsapa = myrmsapa+apar(i,j,k)*apar(i,j,k)
            enddo
         enddo
      enddo
      call MPI_ALLREDUCE(myrmsapa,rmsapa(n),1, &
          MPI_REAL8, &
          MPI_SUM,tube_comm,ierr)
      rmsapa(n) = sqrt(rmsapa(n)/(im*jm*km))

      do i = 0,im-1
         do j = 0,jm-1
            delte(i,j,0) = deltes(i,j,0)+dt*gt0e(i)*gcpte(i)*ey(i,j,0)*bdgxcgy(i,0)/bmag(i,0) &
                           -dt*gamte*deltes(i,j,0)
            delte(i,j,1) = deltes(i,j,1)+dt*gt0e(i)*gcpte(i)*ey(i,j,1)*bdgxcgy(i,1)/bmag(i,1) &
                           -dt*gamte*deltes(i,j,1)
         end do
      end do

      call enfxy(delte(:,:,:))
      call enfz(delte(:,:,:))
      call fltx(delte,0,0,0,1)

      phi_p = delte
      call fltx(phi_p,0,1,0,0)
      delte = delte*(1-epspx)+phi_p*epspx
      phi_p = delte
      call filter(phi_p)
      delte = delte*(1-epspz)+phi_p*epspz

!compute delter
      do i = 0,nxpp-1
         dum = 0.
         dum1 = 0.
         do j = 0,jm-1
            dum = dum+(delte(i,j,0))**2*jac(i,0)
            dum1 = dum1+jac(i,0)
         end do
         mydbr(i) = dum
         myjac(i) = dum1
      end do
      call MPI_ALLREDUCE(mydbr,dtr,nxpp,  &
          MPI_REAL8,MPI_SUM,           &
          tube_comm,ierr)
      call MPI_ALLREDUCE(myjac,v,nxpp,  &
          MPI_REAL8,MPI_SUM,           &
          tube_comm,ierr)
      dtr(0:imx-1) = sqrt(dtr(0:imx-1)/v(0:imx-1))

      if(isgkm==0)return
      do i = 0,im-1
         do j = 0,jm-1
            phi(i,j,0) = phis(i,j,0)+dt*dphidt(i,j,0)-dt*gamphi*phis(i,j,0)
            phi(i,j,1) = phis(i,j,1)+dt*dphidt(i,j,1)-dt*gamphi*phis(i,j,1)
         end do
      end do
      call enfxy(phi(:,:,:))
      call enfz(phi(:,:,:))
      call fltx(phi,0,0,0,1)

      phi_p = phi
      call fltx(phi_p,0,1,0,0)
      phi = phi*(1-epspx)+phi_p*epspx
      phi_p = phi
      call filter(phi_p)
      phi = phi*(1-epspz)+phi_p*epspz

      myrmsapa = 0.
      rmsphi(n) = 0.
      do k=0,mykm-1
         do j=0,jm-1
            do i=0,im-1
               myrmsapa = myrmsapa+phi(i,j,k)*phi(i,j,k)
            enddo
         enddo
      enddo
      call MPI_ALLREDUCE(myrmsapa,rmsphi(n),1, &
          MPI_REAL8, &
          MPI_SUM,tube_comm,ierr)
      rmsphi(n) = sqrt(rmsphi(n)/(im*jm*km))

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
      integer :: nbin=50,nplty=9

      dth1 = pi2/nb
      do 51 j = 0,im-1
         r = rin+xg(j)
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         qr = wx0*sf(i)+wx1*sf(i+1)
         dely(j) = dmod(-pi2*lr0/q0*qr*dsign(1.d0,q0)+800.*ly,ly)*tor
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
!            call F07ARF(nb,nb,tmp,nb,IPIV,INFO)
!            call F07AWF(nb,tmp,nb,IPIV,work,100,INFO)
!  call by lapack name instead
            call ZGETRF(nb,nb,tmp,nb,IPIV,INFO)
            call ZGETRI(nb,tmp,nb,IPIV,work,100,INFO)
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
                         MPI_DOUBLE_COMPLEX,id,10, &
                         rbuf(0),mynum, &
                         MPI_DOUBLE_COMPLEX,id,10, &
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
                         MPI_DOUBLE_COMPLEX,id,20, &
                         rbuf(0),mynum, &
                         MPI_DOUBLE_COMPLEX,id,20, &
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
                         MPI_DOUBLE_COMPLEX,lngbr,10, &
                         rbfr(0:imx-1,0:jmx-1),imx*jmx, &
                         MPI_DOUBLE_COMPLEX,rngbr,10, &
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
      subroutine collh(n)
      use gem_com
      use equil
      implicit none
      integer :: i,ip,k,m,n,ncol,icol,isdum,np_old,np_new
      real(8) :: edum,vdum,dum,dum1,ptch,vte,r,qr,th,cost,b,psp
      real(8) :: h_x,h_coll,x,eps,dtcol,uold,hee,nue,ter
      real(8) :: wx0,wx1,wz0,wz1

      ncol = 1
      dtcol = dt/ncol
      if(rneu==0.0)return
      do k = 1,mm(3)
         r=xh3(k)-0.5*lx+lr0

         m = int(zh3(k)/delz)
         wz0 = ((m+1)*delz-zh3(k))/delz
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
         psp = wx0*psi(i)+wx1*psi(i+1)
         uold = uh3(k)
         edum = b*muh(k)+0.5*mims(3)*uh3(k)*uh3(k)
         vdum = sqrt(2.*edum/mims(3))
         ptch = uh3(k)/vdum

! collision frequency for experimental profiles

         nue=rneu*vi**3/(vdum**3+1.0)
         dum = mims(1)/mims(3)*dtcol*nue

         do icol = 1,ncol
            ptch =ptch-dum*ptch+sqrt(dum*(1.-ptch*ptch)) &
                *dsign(1.d0,ran2(iseed)-0.5)
            ptch = dmin1(ptch,0.999d0)
            ptch = dmax1(ptch,-0.999d0)

            vdum = vdum-rneu*(vdum+vi**3/vdum**2)*dtcol*islow
         end do
!source and sink
!         if(vdum>vc)vdum = vmin
         isdum = 0
         if(vdum<vmin)then
            vdum = vc
            ptch = 2.*(ran2(iseed)-0.5)
            isdum=1
         end if

         uh3(k) = vdum*ptch
         muh(k) = 0.5*mims(3)*vdum*vdum*(1.-ptch*ptch)/b
         muh2(k) = muh(k)
         uh2(k) = uh3(k)
         if(isdum==1)then
!with this recycling, ek and pzet as criteria are not needed
!            ekh(k) = muh(k)*b+0.5*mims(3)*uh2(k)**2
!            pzh(k) = mims(3)*uh2(k)/b-q(3)*psp/br0  
!            vih(k) = uh2(k)
!            mui(k) = muh(k)
            yh2(k) = ly*ran2(iseed)
            yh3(k) = yh2(k)
            xh2(k) = xih(k)
            xh3(k) = xh2(k)
            zh2(k) = z0h(k)
            zh3(k) = zh2(k)

            r=xh3(k)-0.5*lx+lr0
            m = int(zh3(k)/delz)
            wz0 = ((m+1)*delz-zh3(k))/delz
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
            muh(k) = 0.5*mims(3)*vdum*vdum*(1.-ptch*ptch)/b
            muh2(k) = muh(k)
            wh2(k) = 0.
            wh3(k) = 0.
         end if
      end do

      np_old=mm(3)
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
      call pmove(muh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
!      call pmove(mui,np_old,np_new,ierr)
!      if (ierr.ne.0) call ppexit
      call pmove(xih,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0h,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
!      call pmove(pzh,np_old,np_new,ierr)
!      if (ierr.ne.0) call ppexit
!      call pmove(ekh,np_old,np_new,ierr)
!      if (ierr.ne.0) call ppexit
!      call pmove(vih,np_old,np_new,ierr)
!      if (ierr.ne.0) call ppexit
      call pmove(index,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
!     
      call end_pmove(ierr)
      mm(3)=np_new

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rstph
      use gem_com
      use equil
      implicit none
      INTEGER :: m,n,ip,i,j,k,j1,j2,j3,j4,j5,nobge
      real(8) :: wx0,wx1,wz0,wz1,bfldp,ter,gdum
      real(8) :: favx=1,favy=1,favz=1,fave=0,favl=1
      REAL(8) :: wx(0:1),wy(0:1),wz(0:1),we(0:1),wp(0:1),vte
      REAL(8) :: r,qr,th,jacp,b,vfac,js,jv,dum,dum1,dum2,pidum,vel,jfnp
      REAL(8) :: xt,yt,zt,vmac,dvfac,dlamb,lamb,dely,g,myavewe,avwh
      integer :: isign,il,ie,bfcnt,mynobge
      real(8) :: w(0:1,0:1,0:1,0:1,0:1)
      real(8) :: mytotw(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1), &
                 totw(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1)
      real(8) :: wght(1:mmx),y(1:mmx)
      real(8) :: myg(0:imx,0:jmx,0:1,0:negrd,0:nlgrd), &
                 tg(0:imx,0:jmx,0:1,0:negrd,0:nlgrd)
      integer :: ng(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1), &
                 myng(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1)
      real(8) ::    h(0:imx,0:jmx,0:1,0:negrd,0:nlgrd), &
                    h1(0:imx,0:jmx,0:1,0:negrd,0:nlgrd), &
                    lbfs(0:imx,0:jmx,0:negrd,0:nlgrd), &
                    rbfs(0:imx,0:jmx,0:negrd,0:nlgrd), &
                    lbfr(0:imx,0:jmx,0:negrd,0:nlgrd), &
                    rbfr(0:imx,0:jmx,0:negrd,0:nlgrd)
      real(8) :: myncell(0:imx-1,0:jmx-1),ncell(0:imx-1,0:jmx-1)  
      real(8) :: mydwcell(0:imx-1,0:jmx-1),dwcell(0:imx-1,0:jmx-1)  
      real(8) :: mydwpcell(0:imx-1,0:jmx-1),dwpcell(0:imx-1,0:jmx-1)  
      real(8) :: mydwecell(0:imx-1,0:jmx-1),dwecell(0:imx-1,0:jmx-1)  
      real(8) :: myavwcell(0:imx-1,0:jmx-1),avwcell(0:imx-1,0:jmx-1)  
      real(8) :: dwn(0:imx-1,0:jmx-1),dwe(0:imx-1,0:jmx-1)  
      real(8) :: dwp(0:imx-1,0:jmx-1)

      vmac = vc
      dvfac = vmac/negrd
      dlamb = 2.0/nlgrd
      pidum = sqrt(pi)*pi2

      if(idg==1)write(*,*)'before 1 loop',myid
      h = 0.
      mytotw = 0.
      mynobge = 0
      myg = 0.
      myng = 0
      myncell = 0.
      myavwcell = 0.

      do m=1,mm(3)
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
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         jacp = wx0*wz0*jacob(i,k)+wx0*wz1*jacob(i,k+1) &
                 +wx1*wz0*jacob(i+1,k)+wx1*wz1*jacob(i+1,k+1)
         jacp = jacp*jfnp
         b=1.-tor+tor*bfldp
         vfac = sqrt(uh3(m)**2 + 2.*muh(m)*b/mims(3))
         vel = vfac
         lamb = uh3(m)/(vel+1.e-6)

         xt = xh3(m)
         i=int(xt/dx)
         i = min(i,im-1)
         yt = yh3(m)
         j = int(yh3(m)/dy)
         j = min(j,jmx-1)
         myncell(i,j) = myncell(i,j)+1
         myavwcell(i,j) = myavwcell(i,j)+abs(wh3(m))
         k=0
         wx(0)=1.-favx+(float(i+1)-xt/dx)*favx
         wx(1)=1.-wx(0)
         wy(0)=1.-favy+(float(j+1)-yt/dy)*favy
         wy(1)=1.-wy(0)
         wz(0)=1.-favz+(float(gclr*kcnt+k+1)-zh3(m)/dz)*favz
         wz(1)=1.-wz(0)
         ie = int(vfac/dvfac)
         il = int((lamb+1.0)/dlamb)
         il = min(il,nlgrd-1)
         if(ie.gt.(negrd-1))then
            mynobge = mynobge+1
            goto 100
         end if
         we(0) = 1.-fave+(float(ie+1)-vfac/dvfac)*fave
         we(1) = 1.-we(0)
         wp(0) = 1-favl+(float(il+1)-(lamb+1.)/dlamb)*favl
         wp(1) = 1.-wp(0)

         do j1 = 0,1
            do j2 = 0,1
               do j3 = 0,1
                  do j4 = 0,1
                     do j5 = 0,1
                        w(j1,j2,j3,j4,j5) = wx(j1)*wy(j2)*wz(j3)*we(j4)*wp(j5)
                     end do
                  end do
               end do
            end do
         end do


         js = jacp
         jv = 1.
         dum = totvol/tmm(1)/(dx*dy*dz*dvfac*dlamb*js*jv)

         myng(i,j,ie,il) = myng(i,j,ie,il)+1
         mytotw(i,j,ie,il) = mytotw(i,j,ie,il)+wh3(m)

         do j1 = 0,1
            do j2 = 0,1
               do j3 = 0,1
                  do j4 = 0,1
                     do j5 = 0,1
                        h(i+j1,j+j2,j3,ie+j4,il+j5) =  &
                           h(i+j1,j+j2,j3,ie+j4,il+j5) + &
                           dum*w(j1,j2,j3,j4,j5)*wh3(m)

                        myg(i+j1,j+j2,j3,ie+j4,il+j5) =  &
                           myg(i+j1,j+j2,j3,ie+j4,il+j5) + &
                           dum*w(j1,j2,j3,j4,j5)
                     end do
                  end do
               end do
            end do
         end do
 100               continue
      end do

      if(idg==1)write(*,*)'before enforce'
      bfcnt = (imx+1)*(jmx+1)*(negrd+1)*(nlgrd+1)
      rbfs(:,:,:,:)=h(:,:,1,:,:)
      call MPI_SENDRECV(rbfs(0,0,0,0),bfcnt, &
          MPI_REAL8, &
          rngbr,101, &
          lbfr(0,0,0,0),bfcnt, &
          MPI_REAL8, &
          lngbr,101, &
          tube_comm,stat,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      lbfs(:,:,:,:)=h(:,:,0,:,:)
      call MPI_SENDRECV(lbfs(0,0,0,0),bfcnt, &
          MPI_REAL8, &
          lngbr,102, &
          rbfr(0,0,0,0),bfcnt, &
          MPI_REAL8, &
          rngbr,102, &
          tube_comm,stat,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do i = 0,imx
         do j = 0,jmx
            do ie = 0,negrd
               do il = 0,nlgrd
                  h(i,j,0,ie,il) = h(i,j,0,ie,il) &
                +weightp(i)*lbfr(i,jpl(i,j),ie,il)  &
                +weightpn(i)*lbfr(i,jpn(i,j),ie,il)
               end do
            end do
         end do
      enddo

      do i = 0,imx
         do j = 0,jmx
            do ie = 0,negrd
               do il = 0,nlgrd
                  h(i,j,1,ie,il) = h(i,j,1,ie,il) &
                +weightm(i)*rbfr(i,jmi(i,j),ie,il) &
                +weightmn(i)*rbfr(i,jmn(i,j),ie,il)
               end do
            end do
         end do
      enddo

      rbfs(:,:,:,:)=myg(:,:,1,:,:)
      call MPI_SENDRECV(rbfs(0,0,0,0),bfcnt, &
          MPI_REAL8, &
          rngbr,103, &
          lbfr(0,0,0,0),bfcnt, &
          MPI_REAL8, &
          lngbr,103, &
          tube_comm,stat,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      lbfs(:,:,:,:)=myg(:,:,0,:,:)
      call MPI_SENDRECV(lbfs(0,0,0,0),bfcnt, &
          MPI_REAL8, &
          lngbr,104, &
          rbfr(0,0,0,0),bfcnt, &
          MPI_REAL8, &
          rngbr,104, &
          tube_comm,stat,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do i = 0,imx
         do j = 0,jmx
            do ie = 0,negrd
               do il = 0,nlgrd
                  myg(i,j,0,ie,il) = myg(i,j,0,ie,il) &
                +weightp(i)*lbfr(i,jpl(i,j),ie,il)  &
                +weightpn(i)*lbfr(i,jpn(i,j),ie,il)
               end do
            end do
         end do
      enddo

      do i = 0,imx
         do j = 0,jmx
            do ie = 0,negrd
               do il = 0,nlgrd
                  myg(i,j,1,ie,il) = myg(i,j,1,ie,il) &
                +weightm(i)*rbfr(i,jmi(i,j),ie,il) &
                +weightmn(i)*rbfr(i,jmn(i,j),ie,il)
               end do
            end do
         end do
      enddo

      if(idg==1)write(*,*)'before reduce'
      call MPI_ALLREDUCE(h,  &
                h1,             &
                2*bfcnt,MPI_REAL8,       &
                MPI_SUM,GRID_COMM,ierr)

      call MPI_ALLREDUCE(mytotw, totw, &
              imx*jmx*negrd*nlgrd,MPI_REAL8,       &
                MPI_SUM,GRID_COMM,ierr)

      call MPI_ALLREDUCE(myg, tg, &
                2*bfcnt,MPI_REAL8,       &
                MPI_SUM,GRID_COMM,ierr)

      call MPI_ALLREDUCE(myng, ng, &
              imx*jmx*negrd*nlgrd,MPI_INTEGER,       &
                MPI_SUM,GRID_COMM,ierr)

      call MPI_ALLREDUCE(mynobge, nobge, &
              1,MPI_integer,       &
                MPI_SUM,MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(myavwcell, avwcell, &
              imx*jmx,MPI_REAL8,       &
                MPI_SUM,GRID_COMM,ierr)

      call MPI_ALLREDUCE(myncell, ncell, &
              imx*jmx,MPI_REAL8,       &
                MPI_SUM,GRID_COMM,ierr)

      if(idg==1)write(*,*)'before 2 loop',myid
      myavewe = 0.
      mydwcell = 0.
      mydwecell = 0.
      mydwpcell = 0.
      do m=1,mm(3)
         wght(m) = wh3(m)
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
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
         jacp = wx0*wz0*jac(i,k)+wx0*wz1*jac(i,k+1) &
                 +wx1*wz0*jac(i+1,k)+wx1*wz1*jac(i+1,k+1)
         b=1.-tor+tor*bfldp

         vfac = sqrt(uh3(m)**2 + 2.*muh(m)*b/mims(3))
         vel = vfac
         lamb = uh3(m)/(vel+1.e-6)

         xt = xh3(m)
         i=int(xt/dx)
         i = min(i,im-1)
         yt = yh3(m)
         j = int(yh3(m)/dy)
         j = min(j,jmx-1)
         k=0
         wx(0)=1.-favx+(float(i+1)-xt/dx)*favx
         wx(1)=1.-wx(0)
         wy(0)=1.-favy+(float(j+1)-yt/dy)*favy
         wy(1)=1.-wy(0)
         wz(0)=1.-favz+(float(gclr*kcnt+k+1)-zh3(m)/dz)*favz
         wz(1)=1.-wz(0)
         ie = int(vfac/dvfac)
         il = int((lamb+1.0)/dlamb)
         il = min(il,nlgrd-1)
         if(ie.gt.(negrd-2))then
            mynobge = mynobge+1
            goto 200
         end if
         we(0) = 1.-fave+(float(ie+1)-vfac/dvfac)*fave
         we(1) = 1.-we(0)
         wp(0) = 1-favl+(float(il+1)-(lamb+1.)/dlamb)*favl
         wp(1) = 1.-wp(0)

         do j1 = 0,1
            do j2 = 0,1
               do j3 = 0,1
                  do j4 = 0,1
                     do j5 = 0,1
                        w(j1,j2,j3,j4,j5) = wx(j1)*wy(j2)*wz(j3)*we(j4)*wp(j5)
                     end do
                  end do
               end do
            end do
         end do

         dum = 0.
         gdum = 0.
!         g = tg(i,ie,il) !exp(-vfac)/(pidum*vte**3)!tg(i,ie,il)

         do j1 = 0,1
            do j2 = 0,1
               do j3 = 0,1
                  do j4 = 0,1
                     do j5 = 0,1
                      dum=dum+h1(i+j1,j+j2,j3,ie+j4,il+j5)*w(j1,j2,j3,j4,j5)
                      gdum=gdum+tg(i+j1,j+j2,j3,ie+j4,il+j5)*w(j1,j2,j3,j4,j5)
                     end do
                  end do
               end do
            end do
         end do
         wght(m) = dum/gdum

 200          continue
         myavewe = myavewe+abs(wght(m))
         mydwcell(i,j) = mydwcell(i,j)+eprs*(wh3(m) -wght(m))
         mydwecell(i,j) = mydwecell(i,j)+eprs*(wh3(m) -wght(m))*vfac
         mydwpcell(i,j) = mydwpcell(i,j)+eprs*(wh3(m) -wght(m))*uh3(m)
         wh3(m) = wh3(m)*(1.-eprs)+wght(m)*eprs
         wh2(m) = wh3(m)
      end do

      call MPI_ALLREDUCE(myavewe, avwh, &
 1        ,MPI_REAL8,       &
                MPI_SUM,MPI_COMM_WORLD,ierr)
      avwh = avwh/tmm(1)

      call MPI_ALLREDUCE(mydwcell, dwcell, &
              imx*jmx,MPI_REAL8,       &
                MPI_SUM,GRID_COMM,ierr)
      call MPI_ALLREDUCE(mydwecell, dwecell, &
              imx*jmx,MPI_REAL8,       &
                MPI_SUM,GRID_COMM,ierr)
      call MPI_ALLREDUCE(mydwpcell, dwpcell, &
              imx*jmx,MPI_REAL8,       &
                MPI_SUM,GRID_COMM,ierr)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine rstpi
      use gem_com
      use equil
      implicit none
      INTEGER :: m,n,ip,i,j,k,j1,j2,j3,j4,j5,nobge
      real(8) :: wx0,wx1,wz0,wz1,bfldp,ter,gdum
      real(8) :: favx=1,favy=1,favz=1,fave=0,favl=1
      REAL(8) :: wx(0:1),wy(0:1),wz(0:1),we(0:1),wp(0:1),vte
      REAL(8) :: r,qr,th,jacp,b,vfac,js,jv,dum,dum1,dum2,pidum,vel,jfnp
      REAL(8) :: xt,yt,zt,emac=10.,dvfac,dlamb,lamb,dely,g,myavewe,avwe
      integer :: isign,il,ie,bfcnt,mynobge
      real(8) :: w(0:1,0:1,0:1,0:1,0:1)
      real(8) :: mytotw(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1), &
                 totw(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1)
      real(8) :: wght(1:mmx),y(1:mmx)
      real(8) :: myg(0:imx,0:jmx,0:1,0:negrd,0:nlgrd), &
                 tg(0:imx,0:jmx,0:1,0:negrd,0:nlgrd)
      integer :: ng(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1), &
                 myng(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1)
      real(8) ::    h(0:imx,0:jmx,0:1,0:negrd,0:nlgrd), &
                    h1(0:imx,0:jmx,0:1,0:negrd,0:nlgrd), &
                    lbfs(0:imx,0:jmx,0:negrd,0:nlgrd), &
                    rbfs(0:imx,0:jmx,0:negrd,0:nlgrd), &
                    lbfr(0:imx,0:jmx,0:negrd,0:nlgrd), &
                    rbfr(0:imx,0:jmx,0:negrd,0:nlgrd)
      real(8) :: myncell(0:imx-1,0:jmx-1),ncell(0:imx-1,0:jmx-1)  
      real(8) :: mydwcell(0:imx-1,0:jmx-1),dwcell(0:imx-1,0:jmx-1)  
      real(8) :: mydwpcell(0:imx-1,0:jmx-1),dwpcell(0:imx-1,0:jmx-1)  
      real(8) :: mydwecell(0:imx-1,0:jmx-1),dwecell(0:imx-1,0:jmx-1)  
      real(8) :: myavwcell(0:imx-1,0:jmx-1),avwcell(0:imx-1,0:jmx-1)  
      real(8) :: dwn(0:imx-1,0:jmx-1),dwe(0:imx-1,0:jmx-1)  
      real(8) :: dwp(0:imx-1,0:jmx-1)

      dvfac = emac/negrd
      dlamb = 2.0/nlgrd
      pidum = sqrt(pi)*pi2
      vte = sqrt(amie)

      if(idg==1)write(*,*)'before 1 loop',myid
      h = 0.
      mytotw = 0.
      mynobge = 0
      myg = 0.
      myng = 0
      myncell = 0.
      myavwcell = 0.

      do m=1,mm(1)
         r=x3(m)-0.5*lx+lr0

         k = int(z3(m)/delz)
         wz0 = ((k+1)*delz-z3(m))/delz
         wz1 = 1-wz0
         th = wz0*thfnz(k)+wz1*thfnz(k+1)

         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         jacp = wx0*wz0*jacob(i,k)+wx0*wz1*jacob(i,k+1) &
                 +wx1*wz0*jacob(i+1,k)+wx1*wz1*jacob(i+1,k+1)
         jacp = jacp*jfnp
         ter = wx0*t0i(i)+wx1*t0i(i+1)
         b=1.-tor+tor*bfldp
         vfac = 0.5*(mims(1)*u3(m)**2 + 2.*mu(m)*b)/ter
         vel = sqrt(2.*vfac*ter/mims(1))
         lamb = u3(m)/(vel+1.e-6)

         xt = x3(m)
         i=int(xt/dx)
         i = min(i,im-1)
         yt = y3(m)
         j = int(y3(m)/dy)
         j = min(j,jmx-1)
         myncell(i,j) = myncell(i,j)+1
         myavwcell(i,j) = myavwcell(i,j)+abs(w3(m))
         k=0
         wx(0)=1.-favx+(float(i+1)-xt/dx)*favx
         wx(1)=1.-wx(0)
         wy(0)=1.-favy+(float(j+1)-yt/dy)*favy
         wy(1)=1.-wy(0)
         wz(0)=1.-favz+(float(gclr*kcnt+k+1)-z3(m)/dz)*favz
         wz(1)=1.-wz(0)
         ie = int(vfac/dvfac)
         il = int((lamb+1.0)/dlamb)
         il = min(il,nlgrd-1)
         if(ie.gt.(negrd-1))then
            mynobge = mynobge+1
            goto 100
         end if
         we(0) = 1.-fave+(float(ie+1)-vfac/dvfac)*fave
         we(1) = 1.-we(0)
         wp(0) = 1-favl+(float(il+1)-(lamb+1.)/dlamb)*favl
         wp(1) = 1.-wp(0)

         do j1 = 0,1
            do j2 = 0,1
               do j3 = 0,1
                  do j4 = 0,1
                     do j5 = 0,1
                        w(j1,j2,j3,j4,j5) = wx(j1)*wy(j2)*wz(j3)*we(j4)*wp(j5)
                     end do
                  end do
               end do
            end do
         end do


         js = jacp
         jv = 1. !sqrt(2.)*vte**3*sqrt(vfac+1.e-3)
         dum = totvol/tmm(1)/(dx*dy*dz*dvfac*dlamb*js*jv)

         myng(i,j,ie,il) = myng(i,j,ie,il)+1
         mytotw(i,j,ie,il) = mytotw(i,j,ie,il)+w3(m)

         do j1 = 0,1
            do j2 = 0,1
               do j3 = 0,1
                  do j4 = 0,1
                     do j5 = 0,1
                        h(i+j1,j+j2,j3,ie+j4,il+j5) =  &
                           h(i+j1,j+j2,j3,ie+j4,il+j5) + &
                           dum*w(j1,j2,j3,j4,j5)*w3(m)

                        myg(i+j1,j+j2,j3,ie+j4,il+j5) =  &
                           myg(i+j1,j+j2,j3,ie+j4,il+j5) + &
                           dum*w(j1,j2,j3,j4,j5)
                     end do
                  end do
               end do
            end do
         end do
 100               continue
      end do

      if(idg==1)write(*,*)'before enforce'
      bfcnt = (imx+1)*(jmx+1)*(negrd+1)*(nlgrd+1)
      rbfs(:,:,:,:)=h(:,:,1,:,:)
      call MPI_SENDRECV(rbfs(0,0,0,0),bfcnt, &
          MPI_REAL8, &
          rngbr,101, &
          lbfr(0,0,0,0),bfcnt, &
          MPI_REAL8, &
          lngbr,101, &
          tube_comm,stat,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      lbfs(:,:,:,:)=h(:,:,0,:,:)
      call MPI_SENDRECV(lbfs(0,0,0,0),bfcnt, &
          MPI_REAL8, &
          lngbr,102, &
          rbfr(0,0,0,0),bfcnt, &
          MPI_REAL8, &
          rngbr,102, &
          tube_comm,stat,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do i = 0,imx
         do j = 0,jmx
            do ie = 0,negrd
               do il = 0,nlgrd
                  h(i,j,0,ie,il) = h(i,j,0,ie,il) &
                +weightp(i)*lbfr(i,jpl(i,j),ie,il)  &
                +weightpn(i)*lbfr(i,jpn(i,j),ie,il)
               end do
            end do
         end do
      enddo

      do i = 0,imx
         do j = 0,jmx
            do ie = 0,negrd
               do il = 0,nlgrd
                  h(i,j,1,ie,il) = h(i,j,1,ie,il) &
                +weightm(i)*rbfr(i,jmi(i,j),ie,il) &
                +weightmn(i)*rbfr(i,jmn(i,j),ie,il)
               end do
            end do
         end do
      enddo

      rbfs(:,:,:,:)=myg(:,:,1,:,:)
      call MPI_SENDRECV(rbfs(0,0,0,0),bfcnt, &
          MPI_REAL8, &
          rngbr,103, &
          lbfr(0,0,0,0),bfcnt, &
          MPI_REAL8, &
          lngbr,103, &
          tube_comm,stat,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      lbfs(:,:,:,:)=myg(:,:,0,:,:)
      call MPI_SENDRECV(lbfs(0,0,0,0),bfcnt, &
          MPI_REAL8, &
          lngbr,104, &
          rbfr(0,0,0,0),bfcnt, &
          MPI_REAL8, &
          rngbr,104, &
          tube_comm,stat,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do i = 0,imx
         do j = 0,jmx
            do ie = 0,negrd
               do il = 0,nlgrd
                  myg(i,j,0,ie,il) = myg(i,j,0,ie,il) &
                +weightp(i)*lbfr(i,jpl(i,j),ie,il)  &
                +weightpn(i)*lbfr(i,jpn(i,j),ie,il)
               end do
            end do
         end do
      enddo

      do i = 0,imx
         do j = 0,jmx
            do ie = 0,negrd
               do il = 0,nlgrd
                  myg(i,j,1,ie,il) = myg(i,j,1,ie,il) &
                +weightm(i)*rbfr(i,jmi(i,j),ie,il) &
                +weightmn(i)*rbfr(i,jmn(i,j),ie,il)
               end do
            end do
         end do
      enddo

      if(idg==1)write(*,*)'before reduce'
      call MPI_ALLREDUCE(h,  &
                h1,             &
                2*bfcnt,MPI_REAL8,       &
                MPI_SUM,GRID_COMM,ierr)

      call MPI_ALLREDUCE(mytotw, totw, &
              imx*jmx*negrd*nlgrd,MPI_REAL8,       &
                MPI_SUM,GRID_COMM,ierr)

      call MPI_ALLREDUCE(myg, tg, &
                2*bfcnt,MPI_REAL8,       &
                MPI_SUM,GRID_COMM,ierr)

      call MPI_ALLREDUCE(myng, ng, &
              imx*jmx*negrd*nlgrd,MPI_INTEGER,       &
                MPI_SUM,GRID_COMM,ierr)

      call MPI_ALLREDUCE(mynobge, nobge, &
              1,MPI_integer,       &
                MPI_SUM,MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(myavwcell, avwcell, &
              imx*jmx,MPI_REAL8,       &
                MPI_SUM,GRID_COMM,ierr)

      call MPI_ALLREDUCE(myncell, ncell, &
              imx*jmx,MPI_REAL8,       &
                MPI_SUM,GRID_COMM,ierr)

      if(idg==1)write(*,*)'before 2 loop',myid
      myavewe = 0.
      mydwcell = 0.
      mydwecell = 0.
      mydwpcell = 0.
      do m=1,mm(1)
         wght(m) = w3(m)
         r=x3(m)-0.5*lx+lr0

         k = int(z3(m)/delz)
         wz0 = ((k+1)*delz-z3(m))/delz
         wz1 = 1-wz0
         th = wz0*thfnz(k)+wz1*thfnz(k+1)

         i = int((r-rin)/dr)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1.-wz0
         bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
                 +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
         jacp = wx0*wz0*jac(i,k)+wx0*wz1*jac(i,k+1) &
                 +wx1*wz0*jac(i+1,k)+wx1*wz1*jac(i+1,k+1)
         ter = wx0*t0i(i)+wx1*t0i(i+1)
         b=1.-tor+tor*bfldp
         vfac = 0.5*(mims(1)*u3(m)**2 + 2.*mu(m)*b)/ter
         vel = sqrt(2.*vfac*ter/mims(1))
         lamb = u3(m)/(vel+1.e-6)

         xt = x3(m)
         i=int(xt/dx)
         i = min(i,im-1)
         yt = y3(m)
         j = int(y3(m)/dy)
         j = min(j,jmx-1)
         k=0
         wx(0)=1.-favx+(float(i+1)-xt/dx)*favx
         wx(1)=1.-wx(0)
         wy(0)=1.-favy+(float(j+1)-yt/dy)*favy
         wy(1)=1.-wy(0)
         wz(0)=1.-favz+(float(gclr*kcnt+k+1)-z3(m)/dz)*favz
         wz(1)=1.-wz(0)
         ie = int(vfac/dvfac)
         il = int((lamb+1.0)/dlamb)
         il = min(il,nlgrd-1)
         if(ie.gt.(negrd-2))then
            mynobge = mynobge+1
            goto 200
         end if
         we(0) = 1.-fave+(float(ie+1)-vfac/dvfac)*fave
         we(1) = 1.-we(0)
         wp(0) = 1-favl+(float(il+1)-(lamb+1.)/dlamb)*favl
         wp(1) = 1.-wp(0)

         do j1 = 0,1
            do j2 = 0,1
               do j3 = 0,1
                  do j4 = 0,1
                     do j5 = 0,1
                        w(j1,j2,j3,j4,j5) = wx(j1)*wy(j2)*wz(j3)*we(j4)*wp(j5)
                     end do
                  end do
               end do
            end do
         end do

         dum = 0.
         gdum = 0.
!         g = tg(i,ie,il) !exp(-vfac)/(pidum*vte**3)!tg(i,ie,il)

         do j1 = 0,1
            do j2 = 0,1
               do j3 = 0,1
                  do j4 = 0,1
                     do j5 = 0,1
                      dum=dum+h1(i+j1,j+j2,j3,ie+j4,il+j5)*w(j1,j2,j3,j4,j5)
                      gdum=gdum+tg(i+j1,j+j2,j3,ie+j4,il+j5)*w(j1,j2,j3,j4,j5)
                     end do
                  end do
               end do
            end do
         end do
         wght(m) = dum/gdum

 200          continue
         myavewe = myavewe+abs(wght(m))
         mydwcell(i,j) = mydwcell(i,j)+eprs*(w3(m) -wght(m))
         mydwecell(i,j) = mydwecell(i,j)+eprs*(w3(m) -wght(m))*vfac
         mydwpcell(i,j) = mydwpcell(i,j)+eprs*(w3(m) -wght(m))*u3(m)
         w3(m) = w3(m)*(1.-eprs)+wght(m)*eprs
         w2(m) = w3(m)
      end do

      call MPI_ALLREDUCE(myavewe, avwe, &
 1        ,MPI_REAL8,       &
                MPI_SUM,MPI_COMM_WORLD,ierr)
      avwe = avwe/tmm(1)

      call MPI_ALLREDUCE(mydwcell, dwcell, &
              imx*jmx,MPI_REAL8,       &
                MPI_SUM,GRID_COMM,ierr)
      call MPI_ALLREDUCE(mydwecell, dwecell, &
              imx*jmx,MPI_REAL8,       &
                MPI_SUM,GRID_COMM,ierr)
      call MPI_ALLREDUCE(mydwpcell, dwpcell, &
              imx*jmx,MPI_REAL8,       &
                MPI_SUM,GRID_COMM,ierr)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sorth
      use gem_com
      use equil
      implicit none
      INTEGER :: np_old,np_new

      np_old=mm(3) 
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
      call pmove(muh2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mui,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xih,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z0h,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(pzh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ekh,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(vih,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(index,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(isrl,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call end_pmove(ierr)
      mm(3)=np_new

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gkps1(nstep,ip)

      use gem_com
      use equil
      use fft_wrapper
      implicit none
      REAL(8) :: b2,gam0,gam1,th,delyz,bf,rkper,r,qr,shat,ter !,b
      REAL(8) :: kx,ky
      complex(8) :: lbfr(0:imx-1,0:jcnt-1)
      complex(8) :: lbfs(0:imx-1,0:jcnt-1)
      complex(8) :: rbfr(0:imx-1,0:jcnt-1)
      complex(8) :: rbfs(0:imx-1,0:jcnt-1)
      REAL(8),dimension(:),allocatable :: dely,aky
      complex(8),dimension(:,:,:,:),allocatable:: nab1,nab2
      complex(8),dimension(:,:,:,:),allocatable :: mx
      integer,dimension(:,:,:,:),allocatable :: ipiv
      REAL(8),dimension(:,:,:),allocatable :: formphi
      REAL(8) :: sgnx,sgny,sz,myfe
      INTEGER :: i,i1,j,j1,k,k1,l,m,n,ifirst,nstep,ip,INFO
      INTEGER :: l1,m1,myk,myj,izonal,ix,ikx,iext
      COMPLEX(8) :: temp3dxy(0:imx-1,0:jmx-1,0:1),&
            v(0:imx-1,0:jcnt-1,0:1)
      COMPLEX(8) :: sbuf(0:imx*jcnt*2-1),rbuf(0:imx*jmx*2-1)
      COMPLEX(8) :: sl(1:imx-1,0:jcnt-1,0:1)
      REAL(8) :: dum,u(0:imx,0:jmx,0:1)
      complex(8) :: ua(0:imx-1,0:jcnt-1,0:1),calph,cbeta
      REAL(8) :: grp,gthp,gxdgyp,dydrp,qhatp,grdgtp,bfldp,g2xp,g2yp, &
                 radiusp,g2zp,gzp,gxdgzp,gydgzp
      REAL(8) :: wx0,wx1,wz0,wz1
      character*70 fname
      character*4 holdmyid

      save formphi,ifirst,dely,aky,mx,ipiv
      calph = 1.
      cbeta = 0.
      write(holdmyid,'(I4.4)') MyId
      fname='./matrix/'//'mx_phi_'//holdmyid
      if (ifirst.ne.-99) then
         allocate(dely(0:imx-1),aky(0:jmx-1),nab1(0:imx-1,0:jcnt-1,0:imx-1,0:1))
         allocate(nab2(0:imx-1,0:jcnt-1,0:imx-1,0:1))
         allocate(mx(imx-1,imx-1,0:jcnt-1,0:1),formphi(0:imx-1,0:jcnt-1,0:1))
         allocate(ipiv(imx-1,imx-1,0:jcnt-1,0:1))

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
                     ter = wx0*t0i(i1)+wx1*t0i(i1+1)
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
                     radiusp = wx0*wz0*radius(i1,k1)+wx0*wz1*radius(i1,k1+1) &
                     +wx1*wz0*radius(i1+1,k1)+wx1*wz1*radius(i1+1,k1+1)

                     m1 = mstart+int((float(m)+1.0)/2)
                     if(m==0)m1=0
                     sgny = isgnft(m)

                     ky=sgny*2.*pi*float(m1)/ly
                     kx=pi*float(l)/lx
                     bf=bfldp
!                     b=mims(1)*(kx*kx*grp**2 + &
!                        ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
!                        +2*dydrp*lr0/q0*qhatp*grdgtp) &
!                        +2*kx*ky*gxdgyp)/bf/bf

                     nab1(l,m,i,n) = kx**2*grp**2+ky**2*  &
                        (dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                        +2*dydrp*lr0/q0*qhatp*grdgtp)

                     nab2(l,m,i,n) = -IU*ky*kx*2*gxdgyp
!   formfactor in ezamp
                     formphi(l,m,n) = 1./jmx
                     if(abs(ky)>kycut)formphi(l,m,n) = 0.
                     if(m1.ne.mlk.and.onemd==1)formphi(l,m,n) = 0.
 53                                 continue
 54                                             continue
 52                                                      continue
 51                                                            continue

         do k = 0,1
            do j = 0,jcnt-1
               do i = 1,imx-1
                  do ix = 1,imx-1
                     mx(i,ix,j,k) = 0.
!                     if(i==ix)mx(i,ix,j,k) = beta*amie*(1+nh*ish)
                     do ikx = 1,imx-1
                        mx(i,ix,j,k) = mx(i,ix,j,k) &
                             +mims(1)/q(1)*nab1(ikx,j,i,k) &
                             *sin(ix*ikx*pi/imx)*sin(i*ikx*pi/imx)*2.0/imx*gn0e(i)/(bmag(i,k)*bmag(i,k)) &
                             +mims(1)/q(1)*nab2(ikx,j,i,k) &
                              *sin(ix*ikx*pi/imx)*cos(i*ikx*pi/imx)*2.0/imx*gn0e(i)/(bmag(i,k)*bmag(i,k))

                     end do
                  end do
               end do
            end do
         end do
         if(gclr==1.and.tclr==0.and.nstep==0)then
            open(20,file="mx",status='unknown')
            j = 0
            k = 0
            do i = 1,imx-1
               do ix = 1,imx-1
                  write(20,10)i,ix,mx(i,ix,j,k),mx(ix,i,j,k)
               end do
            end do
 10                  format(1x,i5,1x,i5,2(2x,e12.5,2x,e12.5))
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

 200          ifirst=-99
      endif

!   now do field solve...
!  find rho(kx,ky)
      u = rho+hden*ishgk+bden*isbgk
      call dcmpy(u(0:imx-1,0:jmx-1,0:1),v)
      if(idg==1)write(*,*)'pass first fft'

      temp3dxy = 0.
      do k = 0,1
         do j = 0,jcnt-1
            myj=jft(j)
            sl(1:imx-1,j,k) = v(1:imx-1,j,k)
            call ZGETRS('N',imx-1,1,mx(:,:,j,k),imx-1,ipiv(:,:,j,k), &
                         sl(:,j,k),imx-1,INFO)
            temp3dxy(1:imx-1,myj,k) = sl(1:imx-1,j,k)
            temp3dxy(0,myj,k) = 0.
         end do
      end do

!  from rho(kx,ky) to phi(kx,ky)
      do k = 0,1
         do j = 0,jmx-1
            do i = 0,imx-1
               temp3dxy(i,j,k)=temp3dxy(i,j,k)/jmx
            end do
         end do
      end do

!  from phi(kx,ky) to phi(x,y)
      do k=0,mykm
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3dxy(i,j,k)
            end do
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3dxy(i,j,k) = tmpy(j)  ! phi(x,y)
            end do
         end do
      end do

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

      return
      end

!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ftcamp
      use gem_com
      use fft_wrapper
      implicit none

      complex(8) :: mycampf(0:mynf-1)
      real(8) :: aomega(0:nfreq-1)
      integer :: i,j,k,thisf,nsize,ir
      real(8) :: domega,om,dum1,dum2,dum3,dum4,dum5,x,gam
      complex(8) :: sbuf(0:mynf-1),rbuf(0:nfreq-1)

      if(idg==1)write(*,*)'enter ftcamp'
      nsize=nm/ifskp
      domega = 2*frmax/nfreq
      do i = 0,nfreq-1
         aomega(i) = -frmax+i*domega
      end do

!      ir = irlk
      do ir = 0,6
!make the data stationary
         dum1 = 0.
         do i = 100,500
            dum1 = dum1+abs(camp(ir,i))
         end do
         dum1 = dum1/400
         dum2 = 0.
         do i = nsize-401,nsize-1
            dum2 = dum2+abs(camp(ir,i))
         end do
         dum2 = dum2/400
         gam = dlog(dum2/dum1)/(nsize-400)
         do i = 0,nsize-1
!         camp(:,i) = camp(:,i)/exp(gam*i)
!         camp(:,i) = exp(IU*1.5e-3*dt*ifskp*i)*exp(gam*i)  !test FT effect
         end do

         do j = 0,mynf-1
            thisf = gclr*mynf+j
            om = aomega(thisf)
            mycampf(j) = 0.
            do i = 0,nsize-1
               mycampf(j) = mycampf(j)+camp(ir,i)*exp(-IU*om*dt*ifskp*i)
            end do
         end do

         do j = 0,mynf-1
            sbuf(j) = mycampf(j)
         end do
         
         call mpi_allgather(sbuf(0),mynf,mpi_double_complex,rbuf,mynf, &
                         mpi_double_complex,tube_comm,ierr)

         do j = 0,nfreq-1
            campf(j) = rbuf(j)
         end do

!find 5 peaks
         dum1 = 0.
         do i=0,nfreq-1
            x = abs(campf(i))
            if(x>dum1)then
               dum1 = x
            end if
         end do
         dum2 = 0.
         do i=0,nfreq-1
            x = abs(campf(i))
            if(x>dum2.and. x .ne. dum1)then
               dum2 = x
            end if
         end do
         dum3 = 0.
         do i=0,nfreq-1
            x = abs(campf(i))
            if(x>dum3 .and. x.ne.dum1 .and. x.ne.dum2)then
               dum3 = x
            end if
         end do
         dum4 = 0.
         do i=0,nfreq-1
            x = abs(campf(i))
            if(x>dum4 .and. x.ne.dum1 .and. x.ne.dum2 .and. x.ne.dum3)then
               dum4 = x
            end if
         end do
         dum5 = 0.
         do i=0,nfreq-1
            x = abs(campf(i))
            if(x>dum5 .and. x.ne.dum1 .and. x.ne.dum2 .and. x.ne.dum3 .and. x.ne.dum4)then
               dum5 = x
            end if
         end do
         

         if(myid==0)then
            do i = 0,nfreq-1
               x = abs(campf(i))
               if(x==dum1) write(15,10)i,aomega(i),abs(campf(i))**2
               if(x==dum2) write(15,10)i,aomega(i),abs(campf(i))**2
               if(x==dum3) write(15,10)i,aomega(i),abs(campf(i))**2
               if(x==dum4) write(15,10)i,aomega(i),abs(campf(i))**2
               if(x==dum5) write(15,10)i,aomega(i),abs(campf(i))**2
            end do
            
            write(15,*)'frequency, ir=', ir
            do i = 0,nfreq-1
               write(15,10)i,aomega(i),abs(campf(i))**2
            end do
            write(15,*)'camp(ir,t)'
            do i = 0,nsize-1
!            write(15,10)i,real(camp(ir,i)),aimag(camp(ir,i)),abs(camp(ir,i))
            end do
         end if
 10      format(1x,i6, 3(2x,e12.5))
      end do

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sweep
      use gem_com
      use fft_wrapper
      implicit none

      complex(8) :: mycampf(0:mynf-1)
      real(8) :: aomega(0:nfreq-1)
      integer :: i,j,k,thisf,nsize,ir,isnap,nbegin,nend
      real(8) :: domega,om,dum1,dum2,dum3,dum4,dum5,x,gam
      complex(8) :: sbuf(0:mynf-1),rbuf(0:nfreq-1)

!      write(*,*)'enter sweep'
      nsize=nm/ifskp
      domega = 2*frmax/nfreq
      do i = 0,nfreq-1
         aomega(i) = -frmax+i*domega
      end do

!      ir = irlk
      do ir = 0,6
         do isnap = 1,nsnap-2
            nbegin = (isnap-1)*nsize/nsnap
            nend = nbegin+3*nsize/nsnap
            if(nend>(nsize-1))nend = nsize
            do j = 0,mynf-1
               thisf = gclr*mynf+j
               om = aomega(thisf)
               mycampf(j) = 0.
               do i = nbegin,nend-1
                  mycampf(j) = mycampf(j)+camp(ir,i)*exp(-IU*om*dt*ifskp*i)
               end do
            end do

            do j = 0,mynf-1
               sbuf(j) = mycampf(j)
            end do
         
            call mpi_allgather(sbuf(0),mynf,mpi_double_complex,rbuf,mynf, &
                         mpi_double_complex,tube_comm,ierr)

            do j = 0,nfreq-1
               campft(ir,isnap,j) = rbuf(j)/(nend-nbegin)
            end do
         end do
         campft(ir,0,:) = campft(ir,1,:)
         campft(ir,nsnap-1,:) = campft(ir,nsnap-2,:)
      end do

      if(myid==0)then
         do k = 0,6
            do i = 0,nsnap-1
               do j = 0,nfreq-1
                  write(21,10)i,j,aomega(j),real(campft(k,i,j)),aimag(campft(k,i,j)),(abs(campft(k,i,j)))**2
               end do
            end do
         end do
      end if
 10      format(1x,i6, 1x,i6,4(2x,e12.5))


      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine jie(n)

!    source quantities are are calculated: n_i
!    right now only ion quantitities are calculated...

      use gem_com
      use equil
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: enerb,vxdum,dum,xdot,ydot,xdt,ydt,zdot,pidum,dum1,dum2,dumv
      INTEGER :: m,n,i,j,k,l,ns,ip,nonfi,nonfe
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,vte
      REAL(8) :: sz,wght,wght0,wght1,wght2,r,th,cost,sint,b,qr,dv,kap,ter
      REAL(8) :: kapnp,kaptp,xnp,kaphp,nhp,pzp,psp,nhip,v,edot,dist,bstar
      REAL(8) :: xt,yt,rhog,vpar,xs,dely,vfac,vp0
      REAL(8) :: lbfs(0:imx,0:jmx)
      REAL(8) :: rbfs(0:imx,0:jmx)
      REAL(8) :: lbfr(0:imx,0:jmx)
      REAL(8) :: rbfr(0:imx,0:jmx)
      real(8) :: myjpar(0:imx,0:jmx,0:1),myjpex(0:imx,0:jmx,0:1)
      real(8) :: myjpey(0:imx,0:jmx,0:1),mydnidt(0:imx,0:jmx,0:1)
      real(8) :: myhpar(0:imx,0:jmx,0:1),myhpex(0:imx,0:jmx,0:1)
      real(8) :: myhpey(0:imx,0:jmx,0:1),mydnhdt(0:imx,0:jmx,0:1)
      real(8) :: myreynx(0:imx,0:jmx,0:1),myreyny(0:imx,0:jmx,0:1)
      real(8) :: mymaxwx(0:imx,0:jmx,0:1),mymaxwy(0:imx,0:jmx,0:1)
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,psip2p,cvbzp
      REAL(8) :: grp,gxdgyp,rhox(4),rhoy(4),vncp,vparspp
      real(8) :: vperp,plam,vpatcl,auglam,dplam,dlam2,y,dydv,dlamdv,ddlamdv

      if(idg==1)write(*,*)'enter jie'
      nonfi = 1 
      nonfe = 1 

      pidum = 1./(pi*2)**1.5*(vwidth)**3
      if(isuni.eq.0)pidum = 1.

      ns = 1
      myjpar = 0.
      myjpex = 0.
      myjpey = 0.
      mydnidt = 0.
      myreynx = 0.
      myreyny = 0.
      mymaxwx = 0.
      mymaxwy = 0.
      do m=1,mm(ns)
         dv=float(lr(ns))*(dx*dy*dz)

         r=x3(m)-0.5*lx+lr0

         k = int(z3(m)/delz)
         wz0 = ((k+1)*delz-z3(m))/delz
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
         ter = wx0*t0i(i)+wx1*t0i(i+1)        
         kaptp = wx0*capti(i)+wx1*capti(i+1)        
         kapnp = wx0*capni(i)+wx1*capni(i+1)        
         xnp = wx0*xn0i(i)+wx1*xn0i(i+1)        

!         b=1.-lr0/br0*cost
         b=1.-tor+tor*bfldp

         rhog=sqrt(2.*b*mu(m)*mims(ns))/(q(ns)*b)*iflr
         vfac = 0.5*(mims(ns)*u3(m)**2 + 2.*mu(m)*b)

         vpar = u3(m)
         kap = kapnp - (1.5-vfac/ter)*kaptp
         wght=w3(m)/dv
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
         do 100 l=1,lr(1)
            xs=x3(m)+rhox(l) !rwx(1,l)*rhog
            yt=y3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            include 'cpushli.h'
 100     continue

         exp1=exp1/4.
         eyp=eyp/4.
         delbxp=delbxp/4.
         delbyp=delbyp/4.
         aparp = aparp/4.

         vpar = u3(m)
         enerb=(mu(m)+mims(ns)*vpar*vpar/b)/q(ns)*tor
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1

         xdot = vxdum*nonlin  &
              -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         xdt = -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlin  &
             +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)
         ydt = iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)
         zdot =  vpar*(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         wght1 = wght0*(vxdum*kap+q(ns)*(xdt*exp1/ter+ydt*eyp/ter))*xnp*(tload/ter)**1.5*exp(vfac*(1/tload-1./ter))

!    now do 1,2,4 point average, where lr is the no. of points...
         do 101 l=1,lr(1)
            xs=x3(m)+rhox(l) !rwx(1,l)*rhog
            yt=y3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            include 'jieli.h'

!            myreynx(i,j,k) = myreynx(i,j,k)+wght*eyp/b*dum1
!            myreyny(i,j,k) = myreyny(i,j,k)-wght*exp1/b*dum1
!            mymaxwx(i,j,k) = mymaxwx(i,j,k)+wght*vpar/b*delbxp*dum1
!            mymaxwy(i,j,k) = mymaxwy(i,j,k)+wght*vpar/b*delbyp*dum1
 101     continue
! subtract the result if there is no FLR         
         if(icncl==1)then
            xs=x3(m)
            yt=y3(m)
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            i=int(xt/dx+0.5)
            j=int(yt/dy+0.5)
            k=int(z3(m)/dz+0.5)-gclr*kcnt
            exp1=ex(i,j,k)
            eyp=ey(i,j,k)
            vxdum = eyp/b*dum1
            xdt = -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
            ydt = iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)
!         xdt = -2./bfldp**3*fp/radiusp*dbdtp*grcgtp
!         ydt = 2./bfldp**3*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)

            wght1 = 4*wght0*(vxdum*kap+q(ns)*(xdt*exp1+ydt*eyp)/ter)*xnp*(tload/ter)**1.5*exp(vfac*(1/tload-1./ter))
            mydnidt(i,j,k) = mydnidt(i,j,k)-wght1
         end if
      enddo

!   enforce periodicity
      call enforce(myjpar)
      call enforce(myjpex)
      call enforce(myjpey)
      call enforce(mydnidt)
      call enforce(myreynx)
      call enforce(myreyny)
      call enforce(mymaxwx)
      call enforce(mymaxwy)
!      call filter(myjpar)
!      call filter(mydnidt)

      do 110 i=0,im
         do 120 j=0,jm
            do 130 k=0,mykm
               jpar(i,j,k) = q(ns)*myjpar(i,j,k)/n0/jac(i,k)*cn0i
               jpex(i,j,k) = q(ns)*myjpex(i,j,k)/n0/jac(i,k)*cn0i
               jpey(i,j,k) = q(ns)*myjpey(i,j,k)/n0/jac(i,k)*cn0i
               dnidt(i,j,k) = q(ns)*mydnidt(i,j,k)/n0/jac(i,k)*cn0i
               reynix(i,j,k) = q(ns)*myreynx(i,j,k)/n0/jac(i,k)*cn0i
               reyniy(i,j,k) = q(ns)*myreyny(i,j,k)/n0/jac(i,k)*cn0i
               maxwix(i,j,k) = q(ns)*mymaxwx(i,j,k)/n0/jac(i,k)*cn0i
               maxwiy(i,j,k) = q(ns)*mymaxwy(i,j,k)/n0/jac(i,k)*cn0i
 130        continue
 120     continue
 110  continue

!for now delete electrons. When kinetic electron closure is used the electron code might be needed

! energetic particles
      if(ishgk==0)goto 299
      myhpar = 0.
      myhpex = 0.
      myhpey = 0.
      mydnhdt = 0.
      myreynx = 0.
      myreyny = 0.
      mymaxwx = 0.
      mymaxwy = 0.
      do m=1,mm(3)
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
         cvbzp = wx0*wz0*curvbz(i,k)+wx0*wz1*curvbz(i,k+1) &
                 +wx1*wz0*curvbz(i+1,k)+wx1*wz1*curvbz(i+1,k+1) 
         fp = wx0*f(i)+wx1*f(i+1)        
         psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         psp = wx0*psi(i)+wx1*psi(i+1)
         pzp = mims(3)*uh3(m)/b*fp/br0-q(3)*psp/br0
         kaphp = wx0*kaphi(i)+wx1*kaphi(i+1)
         nhip = wx0*nhi(i)+wx1*nhi(i+1)        
         xnp = wx0*xn0h(i)+wx1*xn0h(i+1)        

!         b=1.-lr0/br0*cost
         b=1.-tor+tor*bfldp

         rhog=sqrt(2.*b*muh(m)*mims(3))/(q(3)*b)*iflrh
         vfac = 0.5*(mims(3)*uh3(m)**2 + 2.*muh(m)*b)
         v = sqrt(2.*vfac/mims(3))

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
         ezp = 0.
         delbxp = 0.
         delbyp = 0.
         aparp = 0.
         do 200 l=1,lr(1)
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
            ezp=ezp + ez(i,j,k)
            delbxp = delbxp+delbx(i,j,k)
            delbyp = delbyp+delby(i,j,k)
            aparp = aparp+apar(i,j,k)
 200     continue

         exp1=exp1/4.
         eyp=eyp/4.
         ezp=ezp/4.
         delbxp=delbxp/4.
         delbyp=delbyp/4.
         aparp = aparp/4.

         vpar = uh3(m)
         bstar = b*(1+mims(3)*vpar/(b*q(3))*fp/(b**2*radiusp)*(psip2p*grp**2/radiusp+cvbzp))
         enerb=(muh(m)+mims(3)*vpar*vpar/b)/q(3)*tor
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         dumv = dum1*b/bstar
         vxdum = (eyp/b+vpar/b*delbxp)*dum1*b/bstar

         xdot = vxdum*nonlinh -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp             
         xdt = -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp             
         ydot = (-exp1/b+vpar/b*delbyp)*dum1*b/bstar*nonlinh     &
             +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp) &
             -mims(3)*vpar**2/(q(3)*bstar*b)*(psip2p*grp**2/radiusp+cvbzp)*lr0/(radiusp*q0)
         ydt = tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp) &
             -mims(3)*vpar**2/(q(3)*bstar*b)*(psip2p*grp**2/radiusp+cvbzp)*lr0/(radiusp*q0)
         zdot =  (vpar)*b/bstar &
                 *(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         edot = q(3)*(xdt*exp1+ydt*eyp+vpar*ezp )

         dist = nh/(v**3+vi**3)*nhip
         dum = v**3+vi**3
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1*b/bstar

!    now do 1,2,4 point average, where lr is the no. of points...
         wght=wh3(m)/dv
         wght0 = exp(-vfac)/dv
         if(isuni.eq.0)wght0 = 1./dv
         do 201 l=1,lr(1)
            xs=xh3(m)+rhox(l) !rwx(1,l)*rhog
            yt=yh3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            i=int(xt/dx+0.5)
            j=int(yt/dy+0.5)
            k=int(zh3(m)/dz+0.5)-gclr*kcnt

            wght1 = wght0*(vxdum*kaphp*dist+edot*3*v/mims(3)*dist/(v**3+vi**3))*dum
            myhpar(i,j,k) = myhpar(i,j,k)+wght*zdot
            myhpex(i,j,k) = myhpex(i,j,k)+wght*xdot
            myhpey(i,j,k) = myhpey(i,j,k)+wght*ydot
            mydnhdt(i,j,k) = mydnhdt(i,j,k)+wght1
            myreynx(i,j,k) = myreynx(i,j,k)+wght*eyp/b*dumv
            myreyny(i,j,k) = myreyny(i,j,k)-wght*exp1/b*dumv
            mymaxwx(i,j,k) = mymaxwx(i,j,k)+wght*vpar/b*delbxp*dumv
            mymaxwy(i,j,k) = mymaxwy(i,j,k)+wght*vpar/b*delbyp*dumv
 201     continue
! subtract the result if there is no FLR         
         if(icncl==1)then
            xs=xh3(m)
            yt=yh3(m)
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            i=int(xt/dx+0.5)
            j=int(yt/dy+0.5)
            k=int(zh3(m)/dz+0.5)-gclr*kcnt
            exp1=ex(i,j,k)
            eyp=ey(i,j,k)
            vxdum = eyp/b*dum1
            xdt = -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
            ydt = tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)
            edot = q(3)*(xdt*exp1+ydt*eyp)
            wght1 = 4*wght0*(vxdum*kaphp*dist+edot*3*v/mims(3)*dist/(v**3+vi**3))*dum
            mydnhdt(i,j,k) = mydnhdt(i,j,k)-wght1
         end if
      enddo

!   enforce periodicity
      call enforce(myhpar)
      call enforce(myhpex)
      call enforce(myhpey)
      call enforce(mydnhdt)
      call enforce(myreynx)
      call enforce(myreyny)
      call enforce(mymaxwx)
      call enforce(mymaxwy)
!      call filter(myjpar)
!      call filter(mydnidt)

      do 210 i=0,im
         do 220 j=0,jm
            do 230 k=0,mykm
               hpar(i,j,k) = q(3)*myhpar(i,j,k)/n0/jac(i,k)*cv
               hpex(i,j,k) = q(3)*myhpex(i,j,k)/n0/jac(i,k)*cv
               hpey(i,j,k) = q(3)*myhpey(i,j,k)/n0/jac(i,k)*cv
               dnhdt(i,j,k) = q(3)*mydnhdt(i,j,k)/n0/jac(i,k)*pidum*cv
               reynhx(i,j,k) = q(ns)*myreynx(i,j,k)/n0/jac(i,k)*cn0i
               reynhy(i,j,k) = q(ns)*myreyny(i,j,k)/n0/jac(i,k)*cn0i
               maxwhx(i,j,k) = q(ns)*mymaxwx(i,j,k)/n0/jac(i,k)*cn0i
               maxwhy(i,j,k) = q(ns)*mymaxwy(i,j,k)/n0/jac(i,k)*cn0i
 230        continue
 220     continue
 210  continue

 299  continue
!Beam particles
      if(isbgk==0)goto 399
      myhpar = 0.
      myhpex = 0.
      myhpey = 0.
      mydnhdt = 0.
      myreynx = 0.
      myreyny = 0.
      mymaxwx = 0.
      mymaxwy = 0.

      do m=1,mm(4)
         dv=float(lr(1))*(dx*dy*dz)

         r=xb3(m)-0.5*lx+lr0

         k = int(zb3(m)/delz)
         wz0 = ((k+1)*delz-zb3(m))/delz
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
         psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         psp = wx0*psi(i)+wx1*psi(i+1)
         pzp = mims(4)*ub3(m)/b*fp/br0-q(4)*psp/br0
         kaphp = wx0*kapbi(i)+wx1*kapbi(i+1)
         nhip = wx0*nbi(i)+wx1*nbi(i+1)        
         xnp = wx0*xn0b(i)+wx1*xn0b(i+1)        

!         b=1.-lr0/br0*cost
         b=1.-tor+tor*bfldp

         rhog=sqrt(2.*b*mub(m)*mims(4))/(q(4)*b)*iflrh
         vfac = 0.5*(mims(4)*ub3(m)**2 + 2.*mub(m)*b)
         v = sqrt(2.*vfac/mims(4))

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
         ezp = 0.
         delbxp = 0.
         delbyp = 0.
         aparp = 0.
         do 300 l=1,lr(1)
            xs=xb3(m)+rhox(l) !rwx(1,l)*rhog
            yt=yb3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            i=int(xt/dx+0.5)
            j=int(yt/dy+0.5)
            k=int(zb3(m)/dz+0.5)-gclr*kcnt
            exp1=exp1 + ex(i,j,k)
            eyp=eyp + ey(i,j,k)
            ezp=ezp + ez(i,j,k)
            delbxp = delbxp+delbx(i,j,k)
            delbyp = delbyp+delby(i,j,k)
            aparp = aparp+apar(i,j,k)
 300     continue

         exp1=exp1/4.
         eyp=eyp/4.
         ezp=ezp/4.
         delbxp=delbxp/4.
         delbyp=delbyp/4.
         aparp = aparp/4.

         vpar = ub3(m)
         bstar = b*(1+mims(4)*vpar/(b*q(4))*fp/(b**2*radiusp)*(psip2p*grp**2/radiusp+cvbzp))
         enerb=(mub(m)+mims(4)*vpar*vpar/b)/q(4)*tor
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1*b/bstar
         dumv = dum1*b/bstar

         xdot = vxdum*nonlinh -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp             
         xdt = -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp             
         ydot = (-exp1/b+vpar/b*delbyp)*dum1*b/bstar*nonlinh     &
             +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp) &
             -mims(4)*vpar**2/(q(4)*bstar*b)*(psip2p*grp**2/radiusp+cvbzp)*lr0/(radiusp*q0)
         ydt = tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp) &
             -mims(4)*vpar**2/(q(4)*bstar*b)*(psip2p*grp**2/radiusp+cvbzp)*lr0/(radiusp*q0)
         zdot =  (vpar)*b/bstar &
                 *(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         edot = q(4)*(xdt*exp1+ydt*eyp+vpar*ezp )

         vperp = sqrt(2*mub(m)*b/mims(4))
!         v = sqrt(vperp**2+vpar**2)
         plam = vperp**2/v**2/b
         vpatcl = v/vbeam
         auglam=(vpatcl**3+vcrit**3)/(vpatcl**3*(1.0+vcrit**3)) 
         dlam2=dplam0**2+coedlam*(1.0-plam0)*dlog(auglam)
         dplam=sqrt(dlam2)
         y = auglam
         dydv = -3*vbeam**3/(vbeam**3+vibeam**3)*vibeam**3/v**4
         dlamdv = -mub(m)/vfac**2*mims(4)*v
         ddlamdv = coedlam*(1-plam0)/y*dydv/(2*dplam)
         dist = nbeam*nhip/(v**3+vibeam**3)*exp(-(plam-plam0)**2/dlam2)
         dum = v**3+vibeam**3
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1*b/bstar

!    now do 1,2,4 point average, where lr is the no. of points...
         wght=wb3(m)/dv
         wght0 = exp(-vfac)/dv
         if(isuni.eq.0)wght0 = 1./dv
         if(ipass(m)==1 .and. vpar>0)wght0 = 0.
         do 301 l=1,lr(1)
            xs=xb3(m)+rhox(l) !rwx(1,l)*rhog
            yt=yb3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            i=int(xt/dx+0.5)
            j=int(yt/dy+0.5)
            k=int(zb3(m)/dz+0.5)-gclr*kcnt

            wght1 = wght0*(kaphp*vxdum*dist &
                                 -edot/(v*mims(4))*dist*(-3*v**2/(v**3+vibeam**3)-2*(plam-plam0)/dlam2*dlamdv+2*(plam-plam0)**2/dplam**3*ddlamdv))*dum 
            myhpar(i,j,k) = myhpar(i,j,k)+wght*zdot
            myhpex(i,j,k) = myhpex(i,j,k)+wght*xdot
            myhpey(i,j,k) = myhpey(i,j,k)+wght*ydot
            mydnhdt(i,j,k) = mydnhdt(i,j,k)+wght1
            myreynx(i,j,k) = myreynx(i,j,k)+wght*eyp/b*dumv
            myreyny(i,j,k) = myreyny(i,j,k)-wght*exp1/b*dumv
            mymaxwx(i,j,k) = mymaxwx(i,j,k)+wght*vpar/b*delbxp*dumv
            mymaxwy(i,j,k) = mymaxwy(i,j,k)+wght*vpar/b*delbyp*dumv
 301     continue
! subtract the result if there is no FLR         
         if(icncl==1)then
            xs=xb3(m)
            yt=yb3(m)
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            i=int(xt/dx+0.5)
            j=int(yt/dy+0.5)
            k=int(zb3(m)/dz+0.5)-gclr*kcnt
            exp1=ex(i,j,k)
            eyp=ey(i,j,k)
            vxdum = eyp/b*dum1
            xdt = -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
            ydt = tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)
            edot = q(4)*(xdt*exp1+ydt*eyp)
            wght1 = 4*wght0*(kaphp*vxdum*dist &
                                 -edot/(v*mims(4))*dist*(-3*v**2/(v**3+vibeam**3)-2*(plam-plam0)/dlam2*dlamdv+2*(plam-plam0)**2/dplam**3*ddlamdv))*dum
            mydnhdt(i,j,k) = mydnhdt(i,j,k)-wght1
         end if
      enddo

!   enforce periodicity
      call enforce(myhpar)
      call enforce(myhpex)
      call enforce(myhpey)
      call enforce(mydnhdt)
      call enforce(myreynx)
      call enforce(myreyny)
      call enforce(mymaxwx)
      call enforce(mymaxwy)
!      call filter(myjpar)
!      call filter(mydnidt)

      do 310 i=0,im
         do 320 j=0,jm
            do 330 k=0,mykm
               bpar(i,j,k) = q(4)*myhpar(i,j,k)/n0/jac(i,k)*cvbeam
               bpex(i,j,k) = q(4)*myhpex(i,j,k)/n0/jac(i,k)*cvbeam
               bpey(i,j,k) = q(4)*myhpey(i,j,k)/n0/jac(i,k)*cvbeam
               dnbdt(i,j,k) = q(4)*mydnhdt(i,j,k)/n0/jac(i,k)*pidum*cvbeam
               reynbx(i,j,k) = q(ns)*myreynx(i,j,k)/n0/jac(i,k)*cn0i
               reynby(i,j,k) = q(ns)*myreyny(i,j,k)/n0/jac(i,k)*cn0i
               maxwbx(i,j,k) = q(ns)*mymaxwx(i,j,k)/n0/jac(i,k)*cn0i
               maxwby(i,j,k) = q(ns)*mymaxwy(i,j,k)/n0/jac(i,k)*cn0i
 330        continue
 320     continue
 310  continue
 399  continue

      if(iscgk==0)goto 499
      myjpar = 0.
      myjpex = 0.
      myjpey = 0.
      mydnidt = 0.
      myreynx = 0.
      myreyny = 0.
      mymaxwx = 0.
      mymaxwy = 0.
      do m=1,mm(2)
         dv=float(lr(1))*(dx*dy*dz)

         r=xc3(m)-0.5*lx+lr0

         k = int(zc3(m)/delz)
         wz0 = ((k+1)*delz-zc3(m))/delz
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
         ter = wx0*t0c(i)+wx1*t0c(i+1)        
         kaptp = wx0*captc(i)+wx1*captc(i+1)        
         kapnp = wx0*capnc(i)+wx1*capnc(i+1)        
         xnp = wx0*xn0c(i)+wx1*xn0c(i+1)        

!         b=1.-lr0/br0*cost
         b=1.-tor+tor*bfldp

         rhog=sqrt(2.*b*muc(m)*mims(2))/(q(2)*b)*iflr
         vfac = 0.5*(mims(2)*uc3(m)**2 + 2.*muc(m)*b)

         vpar = uc3(m)
         kap = kapnp - (1.5-vfac/ter)*kaptp
         wght=wc3(m)/dv
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
         do 400 l=1,lr(1)
            xs=xc3(m)+rhox(l) !rwx(1,l)*rhog
            yt=yc3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            include 'mcushli.h'
 400     continue

         exp1=exp1/4.
         eyp=eyp/4.
         delbxp=delbxp/4.
         delbyp=delbyp/4.
         aparp = aparp/4.

         vpar = uc3(m)
         enerb=(muc(m)+mims(2)*vpar*vpar/b)/q(2)*tor
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
         dumv = dum1*b/bstar

         xdot = vxdum*nonlinh  &
              -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         xdt = -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlinh  &
             +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)
         ydt = iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)
         zdot =  vpar*(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         wght1 = wght0*(vxdum*kap+q(2)*(xdt*exp1/ter+ydt*eyp/ter))*xnp*(tloadc/ter)**1.5*exp(vfac*(1/tloadc-1./ter))

!    now do 1,2,4 point average, where lr is the no. of points...
         do 401 l=1,lr(1)
            xs=xc3(m)+rhox(l) !rwx(1,l)*rhog
            yt=yc3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            include 'jielic.h'
!            myreynx(i,j,k) = myreynx(i,j,k)+wght*eyp/b*dum1
!            myreyny(i,j,k) = myreyny(i,j,k)-wght*exp1/b*dum1
!            mymaxwx(i,j,k) = mymaxwx(i,j,k)+wght*vpar/b*delbxp*dum1
!            mymaxwy(i,j,k) = mymaxwy(i,j,k)+wght*vpar/b*delbyp*dum1
 401     continue
! subtract the result if there is no FLR         
         if(icncl==1)then
            xs=xc3(m)
            yt=yc3(m)
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            i=int(xt/dx+0.5)
            j=int(yt/dy+0.5)
            k=int(zc3(m)/dz+0.5)-gclr*kcnt
            exp1=ex(i,j,k)
            eyp=ey(i,j,k)
            vxdum = eyp/b*dum1
            xdt = -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
            ydt = iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)
!         xdt = -2./bfldp**3*fp/radiusp*dbdtp*grcgtp
!         ydt = 2./bfldp**3*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)

            wght1 = 4*wght0*(vxdum*kap+q(2)*(xdt*exp1+ydt*eyp)/ter)*xnp*(tloadc/ter)**1.5*exp(vfac*(1/tloadc-1./ter))
            mydnidt(i,j,k) = mydnidt(i,j,k)-wght1
         end if
      enddo

!   enforce periodicity
      call enforce(myjpar)
      call enforce(myjpex)
      call enforce(myjpey)
      call enforce(mydnidt)
      call enforce(myreynx)
      call enforce(myreyny)
      call enforce(mymaxwx)
      call enforce(mymaxwy)
!      call filter(myjpar)
!      call filter(mydnidt)

      do 410 i=0,im
         do 420 j=0,jm
            do 430 k=0,mykm
               cpar(i,j,k) = q(2)*myjpar(i,j,k)/n0/jac(i,k)*cn0c
               cpex(i,j,k) = q(2)*myjpex(i,j,k)/n0/jac(i,k)*cn0c
               cpey(i,j,k) = q(2)*myjpey(i,j,k)/n0/jac(i,k)*cn0c
               dncdt(i,j,k) = q(2)*mydnidt(i,j,k)/n0/jac(i,k)*cn0c
               reyncx(i,j,k) = q(2)*myreynx(i,j,k)/n0/jac(i,k)*cn0c
               reyncy(i,j,k) = q(2)*myreyny(i,j,k)/n0/jac(i,k)*cn0c
               maxwcx(i,j,k) = q(2)*mymaxwx(i,j,k)/n0/jac(i,k)*cn0c
               maxwcy(i,j,k) = q(2)*mymaxwy(i,j,k)/n0/jac(i,k)*cn0c
 430        continue
 420     continue
 410  continue
 499  continue
 999  continue
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine drdt(n)
      use gem_com
      use equil
      implicit none
      integer :: i,j,k,n,isdndt
      real(8) :: lbfr(0:imx,0:jmx)
      real(8) :: lbfs(0:imx,0:jmx)
      real(8) :: rbfr(0:imx,0:jmx)
      real(8) :: rbfs(0:imx,0:jmx)
      real(8) :: dum
      real(8) :: djdx(0:imx,0:jmx,0:1),djdy(0:imx,0:jmx,0:1)
      real*8 :: dmnl1(0:imx,0:jmx,0:1),dmnl2(0:imx,0:jmx,0:1),dmnl3(0:imx,0:jmx,0:1),dmnl4(0:imx,0:jmx,0:1)
      real*8 :: ux(0:imx,0:jmx,0:1),uy(0:imx,0:jmx,0:1)
      real(8) :: myjaca(0:imx-1),jaca(0:imx-1),myava(0:imx-1),ava(0:imx-1)

!      call gradx(reynix+reynhx*ishgk+reynbx*isbgk+reyncx*iscgk,ux)
      call gradx(reynix+reyncx*iscgk,ux)
!      call grady(reyniy+reynhy*ishgk+reynby*isbgk+reyncy*iscgk,uy)
      call grady(reyniy+reyncy*iscgk,uy)
      reyn = -(ux+uy)
      call gradx(reyncx,ux)
      call grady(reyncy,uy)
      reynh = -(ux+uy)

!      call gradx(maxwix+maxwhx*ishgk+maxwbx*isbgk+maxwcx*iscgk,ux)
      call gradx(maxwix+maxwcx*iscgk,ux)
!      call grady(maxwiy+maxwhy*ishgk+maxwby*isbgk+maxwcy*iscgk,uy)
      call grady(maxwiy+maxwcy*iscgk,uy)
      maxw = -(ux+uy)
      call gradx(maxwcx,ux)
      call grady(maxwcy,uy)
      maxwh = -(ux+uy)

      if(idg==1)write(*,*)'enter drdt'
!electron ddedt
      do i = 0,im
         do j = 0,jm
            rbfs(i,j)=upar(i,j,0)/bmag(i,0)
         end do   
      end do
      call MPI_SENDRECV(rbfs(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8,rngbr,404, &
          lbfr(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8,lngbr,404, &
          tube_comm,stat,ierr)
      do i = 0,im
         do j = 0,jm
            lbfs(i,j)=upar(i,j,1)/bmag(i,1)
         end do   
      end do
      call MPI_SENDRECV(lbfs(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8,lngbr,405, &
          rbfr(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8,rngbr,405, &
          tube_comm,stat,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      if(iske==0)then
         call gradu(dene*ex,ux,uy)
         dmnl1 = uy
         call gradu(dene*ey,ux,uy)
         dmnl2 = ux
         call gradu(delby*upar,ux,uy)
         dmnl3 = uy
         call gradu(delbx*upar,ux,uy)
         dmnl4 = ux
      end if

      if(iske==1)then
         call gradu(denek*ex,ux,uy)
         dmnl1 = uy
         call gradu(denek*ey,ux,uy)
         dmnl2 = ux
         call gradu(delby*upark,ux,uy)
         dmnl3 = uy
         call gradu(delbx*upark,ux,uy)
         dmnl4 = ux
      end if

!dene_p for PC, denes for RK
      do i = 0,im-1
         do j = 0,jm-1
            dum =  weightp(i)*lbfr(i,jpl(i,j)) &
       	+weightpn(i)*lbfr(i,jpn(i,j))
            ddedt(i,j,0) = -(upar(i,j,1)/bmag(i,1)-dum) &
                /(2.*dz)*bmag(i,0)*bdgrzn(i,0)*gn0e(i) &
                +gn0e(i)*gcpne(i)*ey(i,j,0)/bmag(i,0)*bdgxcgy(i,0)*icmprs*(1-icncl) &
                +2.*(cfx(i,0)*(-dnedx(i,j,0)*gt0e(i)/gn0e(i)*ispre-ex(i,j,0)*(1-icncl)) &
                +cfy(i,0)*(-dnedy(i,j,0)*gt0e(i)/gn0e(i)*ispre-ey(i,j,0)*(1-icncl)))/br0*icmprs*gn0e(i) &
                -(delby(i,j,0)*gnuobx(i,0)+delbx(i,j,0)*gnuoby(i,0))    &
                -bdgxcgy(i,0)*(-dmnl1(i,j,0)+dmnl2(i,j,0))/bmag(i,0)*nonline &
                -bdgxcgy(i,0)*(dmnl3(i,j,0)+dmnl4(i,j,0))*gn0e(i)/bmag(i,0)*nonline*iflut
            reyn(i,j,0) = reyn(i,j,0)+bdgxcgy(i,0)*(-dmnl1(i,j,0)+dmnl2(i,j,0))/bmag(i,0)/ntube
            maxw(i,j,0) = maxw(i,j,0)+bdgxcgy(i,0)*(dmnl3(i,j,0)+dmnl4(i,j,0))*gn0e(i)/bmag(i,0)/ntube
            maxwe(i,j,0) = bdgxcgy(i,0)*(dmnl3(i,j,0)+dmnl4(i,j,0))*gn0e(i)/bmag(i,0)/ntube

            dum = weightm(i)*rbfr(i,jmi(i,j)) &
                +weightmn(i)*rbfr(i,jmn(i,j))
            ddedt(i,j,1) = -(dum-upar(i,j,0)/bmag(i,0)) &
                /(2.*dz)*bmag(i,1)*bdgrzn(i,1)*gn0e(i) &
                +gn0e(i)*gcpne(i)*ey(i,j,1)/bmag(i,1)*bdgxcgy(i,1)*icmprs*(1-icncl) &
                +2.*(cfx(i,1)*(-dnedx(i,j,1)*gt0e(i)/gn0e(i)*ispre-ex(i,j,1)*(1-icncl)) &
                +cfy(i,1)*(-dnedy(i,j,1)*gt0e(i)/gn0e(i)*ispre-ey(i,j,1)*(1-icncl)))/br0*icmprs*gn0e(i) &
                -(delby(i,j,1)*gnuobx(i,1)+delbx(i,j,1)*gnuoby(i,1))    &
                -bdgxcgy(i,1)*(-dmnl1(i,j,1)+dmnl2(i,j,1))/bmag(i,1)*nonline &
                -bdgxcgy(i,1)*(dmnl3(i,j,1)+dmnl4(i,j,1))*gn0e(i)/bmag(i,1)*nonline*iflut
            reyn(i,j,1) = reyn(i,j,1)+bdgxcgy(i,1)*(-dmnl1(i,j,1)+dmnl2(i,j,1))/bmag(i,1)/ntube
            maxw(i,j,1) = maxw(i,j,1)+bdgxcgy(i,1)*(dmnl3(i,j,1)+dmnl4(i,j,1))*gn0e(i)/bmag(i,1)/ntube
            maxwe(i,j,1) = bdgxcgy(i,1)*(dmnl3(i,j,1)+dmnl4(i,j,1))*gn0e(i)/bmag(i,1)/ntube
         end do
      end do

      call gradx(jpex*ision+hpex*ishgk+bpex*isbgk+cpex*iscgk,djdx)
      call grady(jpey*ision+hpey*ishgk+bpey*isbgk+cpey*iscgk,djdy)
!      djdx = 0.
!      djdy = 0.
      isdndt = 1
!      upa0(:,:,:) = apar(ip,:,:,:)
      do i = 0,im
         do j = 0,jm
            rbfs(i,j)=jpar(i,j,0)*ision+hpar(i,j,0)*ishgk+bpar(i,j,0)*isbgk+cpar(i,j,0)*iscgk
         end do   
      end do
      call MPI_SENDRECV(rbfs,(imx+1)*(jmx+1),  &
          MPI_REAL8,rngbr,404, &
          lbfr,(imx+1)*(jmx+1), &
          MPI_REAL8,lngbr,404, &
          tube_comm,stat,ierr)
      do i = 0,im
         do j = 0,jm
            lbfs(i,j)=jpar(i,j,1)*ision+hpar(i,j,1)*ishgk+bpar(i,j,1)*isbgk+cpar(i,j,1)*iscgk
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
            drhodt(i,j,0) = -(jpar(i,j,1)*ision+hpar(i,j,1)*ishgk+bpar(i,j,1)*isbgk+cpar(i,j,1)*iscgk-dum)/(2.*dz) &
                -djdx(i,j,0)-djdy(i,j,0) &
                +(dnidt(i,j,0)*ision+dnhdt(i,j,0)*ishgk+dnbdt(i,j,0)*isbgk+dncdt(i,j,0)*iscgk-ddedt(i,j,0)/ntube)

            dum = weightm(i)*rbfr(i,jmi(i,j)) &
                +weightmn(i)*rbfr(i,jmn(i,j))
            drhodt(i,j,1) = -(dum-(jpar(i,j,0)*ision+hpar(i,j,0)*ishgk+bpar(i,j,0)*isbgk+cpar(i,j,0)*iscgk))/(2.*dz)  &
                -djdx(i,j,1)-djdy(i,j,1) &
                +(dnidt(i,j,1)*ision+dnhdt(i,j,1)*ishgk+dnbdt(i,j,1)*isbgk+dncdt(i,j,1)*iscgk-ddedt(i,j,1)/ntube)
         end do
      end do
      call enfxy(drhodt(:,:,:))
      call enfz(drhodt(:,:,:))
      call enfxy(reyn)
      call enfz(reyn)
      call enfxy(maxw)
      call enfz(maxw)
!      call filter(drhodt(:,:,:))

      do i=0,im-1
         myava(i) = 0.
         myjaca(i) = 0.
         do j=0,jm-1
            myava(i)=myava(i)+reyn(i,j,0)*jac(i,0)
            myjaca(i)=myjaca(i)+jac(i,0)
         enddo
      enddo
      call MPI_ALLREDUCE(myava,ava,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,TUBE_COMM,ierr)
      call MPI_ALLREDUCE(myjaca,jaca,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,TUBE_COMM,ierr)

      ava(0:imx-1)=ava(0:imx-1)/jaca(0:imx-1)
      call MPI_ALLREDUCE(ava,reyn00,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,GRID_COMM,ierr)

      do i=0,im-1
         myava(i) = 0.
         myjaca(i) = 0.
         do j=0,jm-1
            myava(i)=myava(i)+reynh(i,j,0)*jac(i,0)
            myjaca(i)=myjaca(i)+jac(i,0)
         enddo
      enddo
      call MPI_ALLREDUCE(myava,ava,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,TUBE_COMM,ierr)
      call MPI_ALLREDUCE(myjaca,jaca,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,TUBE_COMM,ierr)

      ava(0:imx-1)=ava(0:imx-1)/jaca(0:imx-1)
      call MPI_ALLREDUCE(ava,reynh00,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,GRID_COMM,ierr)

      do i=0,im-1
         myava(i) = 0.
         myjaca(i) = 0.
         do j=0,jm-1
            myava(i)=myava(i)+maxw(i,j,0)*jac(i,0)
            myjaca(i)=myjaca(i)+jac(i,0)
         enddo
      enddo
      call MPI_ALLREDUCE(myava,ava,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,TUBE_COMM,ierr)
      call MPI_ALLREDUCE(myjaca,jaca,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,TUBE_COMM,ierr)

      ava(0:imx-1)=ava(0:imx-1)/jaca(0:imx-1)
      call MPI_ALLREDUCE(ava,maxw00,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,GRID_COMM,ierr)

      do i=0,im-1
         myava(i) = 0.
         myjaca(i) = 0.
         do j=0,jm-1
            myava(i)=myava(i)+maxwh(i,j,0)*jac(i,0)
            myjaca(i)=myjaca(i)+jac(i,0)
         enddo
      enddo
      call MPI_ALLREDUCE(myava,ava,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,TUBE_COMM,ierr)
      call MPI_ALLREDUCE(myjaca,jaca,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,TUBE_COMM,ierr)

      ava(0:imx-1)=ava(0:imx-1)/jaca(0:imx-1)
      call MPI_ALLREDUCE(ava,maxwh00,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,GRID_COMM,ierr)

      do i=0,im-1
         myava(i) = 0.
         myjaca(i) = 0.
         do j=0,jm-1
            myava(i)=myava(i)+maxwe(i,j,0)*jac(i,0)
            myjaca(i)=myjaca(i)+jac(i,0)
         enddo
      enddo
      call MPI_ALLREDUCE(myava,ava,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,TUBE_COMM,ierr)
      call MPI_ALLREDUCE(myjaca,jaca,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,TUBE_COMM,ierr)

      ava(0:imx-1)=ava(0:imx-1)/jaca(0:imx-1)
      call MPI_ALLREDUCE(ava,maxwe00,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,GRID_COMM,ierr)

      do i=0,im-1
         myava(i) = 0.
         myjaca(i) = 0.
         do j=0,jm-1
            myava(i)=myava(i)+drhodt(i,j,0)*jac(i,0)
            myjaca(i)=myjaca(i)+jac(i,0)
         enddo
      enddo
      call MPI_ALLREDUCE(myava,ava,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,TUBE_COMM,ierr)
      call MPI_ALLREDUCE(myjaca,jaca,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,TUBE_COMM,ierr)

      ava(0:imx-1)=ava(0:imx-1)/jaca(0:imx-1)
      call MPI_ALLREDUCE(ava,drdt00,imx, &
                 MPI_REAL8,                               &
                 MPI_SUM,GRID_COMM,ierr)

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dpdt(n)   
      use gem_com
      use equil
      use fft_wrapper
      implicit none
      REAL(8) :: b,b1,b2,b1c,b2c,gam0,gam1,delyz,th,bf,dum,r,qr,shat,ter,terc
      REAL(8) :: kx1,kx2,ky
      REAL(8),dimension(:),allocatable :: akx,aky
      real(8),dimension(:,:,:,:),allocatable:: gamb1,gamb2,gamb1c,gamb2c
      complex(8),dimension(:,:,:,:),allocatable :: mx
      integer,dimension(:,:,:,:),allocatable :: ipiv
      REAL(8),dimension(:,:,:),allocatable :: formdpt
      REAL(8) :: sgnx,sgny,sz,myfe,u(0:imx,0:jmx,0:1)
      INTEGER :: i,i1,j,j1,k,k1,l,m,n,ifirst,nstep,INFO
      INTEGER :: l1,m1,myk,myj,ix,ikx
      COMPLEX(8) :: temp3dxy(0:imx-1,0:jmx-1,0:1),v(0:imx-1,0:jcnt-1,0:1)
      COMPLEX(8) :: sbuf(0:imx*jcnt*2-1),rbuf(0:imx*jmx*2-1)
      COMPLEX(8) :: sl(1:imx-1,0:jcnt-1,0:1)
      complex(8) :: cdum
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
      REAL(8) :: grp,gxdgyp,grdgtp,gthp
      REAL(8) :: wx0,wx1,wz0,wz1
      character*70 fname
      character*5 holdmyid

      save formdpt,ifirst,akx,aky,mx,gamb1,gamb2,ipiv
      if(idg==1)write(*,*)'enter dpdt'

!     form factors....
      write(holdmyid,'(I5.5)') MyId
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
            close(30000+myid)
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

                     m1 = mstart+int((float(m)+1.0)/2)
                     if(m==0)m1=0
                     sgny = isgnft(m)

                     ky=sgny*2.*pi*float(m1)/ly
                     kx1=pi*float(l)/lx
                     kx2=-pi*float(l)/lx
                     bf=bfldp
                     b1=mims(1)*(kx1*kx1*grp**2 + &
                       ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                       +2*dydrp*lr0/q0*qhatp*grdgtp) &
                       +2*kx1*ky*gxdgyp)/(bf*bf)*ter/(q(1)*q(1))

                     b2=mims(1)*(kx2*kx2*grp**2 + &
                       ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                       +2*dydrp*lr0/q0*qhatp*grdgtp) &
                       +2*kx2*ky*gxdgyp)/(bf*bf)*ter/(q(1)*q(1))

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
                        mx(i,ix,j,k) = mx(i,ix,j,k)+q(1)*sin(ix*ikx*pi/imx)* & 
                             ((1-gamb1(ikx,j,i,k))*exp(IU*ikx*i*pi/imx)- &
                              (1-gamb2(ikx,j,i,k))*exp(-IU*ikx*i*pi/imx)) &
                             /ter*gn0e(i)/(IU*imx)
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
            close(30000+myid)
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
      if(idg==1)write(*,*)'pass first fft in dpdt'

      temp3dxy = 0.
      do k = 0,1
         do j = 0,jcnt-1
            myj=jft(j)
            sl(1:imx-1,j,k) = v(1:imx-1,j,k)
            call ZGETRS('N',imx-1,1,mx(:,:,j,k),imx-1,ipiv(:,:,j,k), &
                         sl(:,j,k),imx-1,INFO) 
            temp3dxy(1:imx-1,myj,k) = sl(1:imx-1,j,k)
            temp3dxy(0,myj,k) = 0.
         end do
      end do
      if(idg==1)write(*,*)'pass solving in dpdt'

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
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
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
      call fltx(dphidt(:,:,:),0,0,0,0)
      call filter(dphidt(:,:,:))
      if(idg==1)write(*,*)'pass enfz', myid

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine constrphi
      use gem_com
      use equil
      implicit none
      integer :: n,i,j,k,m,ip
      real(8) :: ky1,ky2

      ky1 = pi2/ly
      ky2 = -pi2/ly
      phiti=0.
      do i = 0,im
         do j = 0,jm
            do k = 0,1
               do m = 0,nom-1
                  phiti(i,j,k) = phiti(i,j,k)+(phiom(i,0,k,m)*exp(-IU*ky1*yg(j))+phiom(i,1,k,m)*exp(-IU*ky2*yg(j))) * &
                                               exp(IU*omlk(m)*ipg*nplot*dt)/nm*iomlk(m)
               end do
            end do
         end do
      end do
      phiti = phiti*nm*dt/pi*2*omu/nom
      call enfxy(phiti(:,:,:))
      call enfz(phiti(:,:,:))
      call fltx(phiti(:,:,:),0,0,1,0)
      call filter(phiti(:,:,:))
      call dump3d(phiti(:,:,:),'phiti',52,n)
      call pol2d

end subroutine constrphi      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine constrdte
      use gem_com
      use equil
      implicit none
      integer :: n,i,j,k,m,ip,nskip
      real(8) :: ky1,ky2,dum,dum1
      real(8) :: mydbr(0:imx-1),myjac(0:imx-1),v(0:imx-1),dtav(1:8)

      ky1 = pi2/ly
      ky2 = -pi2/ly
      dtav = 0.

      nskip = 1000
      do ip = 1,nm/nskip
         dteti=0.
         do i = 0,im
            do j = 0,jm
               do k = 0,1
                  do m = 0,nom-1
                     dteti(i,j,k) = dteti(i,j,k)+(dteom(i,0,k,m)*exp(-IU*ky1*yg(j))+dteom(i,1,k,m)*exp(-IU*ky2*yg(j))) * &
                                               exp(IU*omlk(m)*ip*nskip*dt)/nm*iomlk(m)
                  end do
               end do
            end do
         end do
         dteti = dteti*nm*dt/pi*omu/nom/jmx
         call enfxy(dteti(:,:,:))
         call enfz(dteti(:,:,:))
         call fltx(dteti(:,:,:),0,0,1,0)
         call filter(dteti(:,:,:))

!compute delter
         do i = 0,nxpp-1
            dum = 0.
            dum1 = 0.
            do j = 0,jm-1
               dum = dum+(dteti(i,j,0))**2*jac(i,0)
               dum1 = dum1+jac(i,0)
            end do
            mydbr(i) = dum
            myjac(i) = dum1
         end do
         call MPI_ALLREDUCE(mydbr,dtr,nxpp,  &
             MPI_REAL8,MPI_SUM,           &
             tube_comm,ierr)
         call MPI_ALLREDUCE(myjac,v,nxpp,  &
             MPI_REAL8,MPI_SUM,           &
             tube_comm,ierr)
         dtr(0:imx-1) = sqrt(dtr(0:imx-1)/v(0:imx-1))

 12      format(1x,i7,7(2x,e10.3))
         m = imx/8
         do i = 1,7
            dtav(i) = dtav(i)+dtr(i*m)
         end do
         if(myid==master)then
            write(*,12)ip,dtr(m),dtr(2*m),dtr(3*m),dtr(4*m),dtr(5*m),dtr(6*m),dtr(7*m)
            write(22,12)ip,dtr(m),dtr(2*m),dtr(3*m),dtr(4*m),dtr(5*m),dtr(6*m),dtr(7*m)
         end if
      end do
      dtav = dtav/(nm/nskip)
 11   format(8(2x,e10.3))
      if(myid==0)then
         write(*,11)delom,dtav(1),dtav(2),dtav(3),dtav(4),dtav(5),dtav(6),dtav(7)
         write(22,11)delom,dtav(1),dtav(2),dtav(3),dtav(4),dtav(5),dtav(6),dtav(7)
      end if
      call dump3d(dteti(:,:,:),'dteti',52,n)
      call pol2d

end subroutine constrdte      
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dpdt1(n)   
      use gem_com
      use equil
      use fft_wrapper
      implicit none
      REAL(8) :: b,b1,b2,b1c,b2c,gam0,gam1,delyz,th,bf,dum,r,qr,shat,ter,terc
      REAL(8) :: kx1,kx2,ky,kx
      REAL(8),dimension(:),allocatable :: akx,aky
      real(8),dimension(:,:,:,:),allocatable:: nab1,nab2
      complex(8),dimension(:,:,:,:),allocatable :: mx
      integer,dimension(:,:,:,:),allocatable :: ipiv
      REAL(8),dimension(:,:,:),allocatable :: formdpt
      REAL(8) :: sgnx,sgny,sz,myfe,u(0:imx,0:jmx,0:1)
      INTEGER :: i,i1,j,j1,k,k1,l,m,n,ifirst,nstep,INFO
      INTEGER :: l1,m1,myk,myj,ix,ikx
      COMPLEX(8) :: temp3dxy(0:imx-1,0:jmx-1,0:1),v(0:imx-1,0:jcnt-1,0:1)
      COMPLEX(8) :: sbuf(0:imx*jcnt*2-1),rbuf(0:imx*jmx*2-1)
      COMPLEX(8) :: sl(1:imx-1,0:jcnt-1,0:1)
      complex(8) :: cdum
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
      REAL(8) :: grp,gxdgyp,grdgtp,gthp
      REAL(8) :: wx0,wx1,wz0,wz1
      character*70 fname
      character*5 holdmyid

      save formdpt,ifirst,akx,aky,mx,ipiv
      if(idg==1)write(*,*)'enter dpdt'

!     form factors....
      write(holdmyid,'(I5.5)') MyId
      fname='./matrix/'//'mx_dpt_'//holdmyid
      if (ifirst.ne.-99) then
         allocate(akx(0:imx-1),aky(0:jcnt-1))
         allocate(nab1(0:imx-1,0:jcnt-1,0:imx-1,0:1),nab2(0:imx-1,0:jcnt-1,0:imx-1,0:1))
         allocate(mx(imx-1,imx-1,0:jcnt-1,0:1),formdpt(0:imx-1,0:jcnt-1,0:1))
         allocate(ipiv(imx-1,imx-1,0:jcnt-1,0:1))

         if(iget.eq.1) then
            open(30000+MyId,file=fname,form='unformatted',status='old')
            read(30000+MyId)mx,ipiv
            close(30000+myid)
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

                     m1 = mstart+int((float(m)+1.0)/2)
                     if(m==0)m1=0
                     sgny = isgnft(m)

                     ky=sgny*2.*pi*float(m1)/ly
                     kx=pi*float(l)/lx
                     bf=bfldp

                     nab1(l,m,i,n) = kx**2*grp**2+ky**2*  &
                        (dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                        +2*dydrp*lr0/q0*qhatp*grdgtp)

                     nab2(l,m,i,n) = -IU*ky*kx*2*gxdgyp

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
                        mx(i,ix,j,k) = mx(i,ix,j,k) &
                             +mims(1)/q(1)*nab1(ikx,j,i,k) &
                             *sin(ix*ikx*pi/imx)*sin(i*ikx*pi/imx)*2.0/imx*gn0e(i)/(bmag(i,k)*bmag(i,k)) &
                             +mims(1)/q(1)*nab2(ikx,j,i,k) &
                              *sin(ix*ikx*pi/imx)*cos(i*ikx*pi/imx)*2.0/imx*gn0e(i)/(bmag(i,k)*bmag(i,k))
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
            close(30000+myid)
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
      if(idg==1)write(*,*)'pass first fft in dpdt'

      temp3dxy = 0.
      do k = 0,1
         do j = 0,jcnt-1
            myj=jft(j)
            sl(1:imx-1,j,k) = v(1:imx-1,j,k)
            call ZGETRS('N',imx-1,1,mx(:,:,j,k),imx-1,ipiv(:,:,j,k), &
                         sl(:,j,k),imx-1,INFO) 
            temp3dxy(1:imx-1,myj,k) = sl(1:imx-1,j,k)
            temp3dxy(0,myj,k) = 0.
         end do
      end do
      if(idg==1)write(*,*)'pass solving in dpdt'

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
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
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
      call fltx(dphidt(:,:,:),0,0,0,0)
      call filter(dphidt(:,:,:))
      if(idg==1)write(*,*)'pass enfz', myid

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine bpush(n)

      use gem_com
      use equil
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum,dum1
      INTEGER :: m,i,j,k,l,n,mynopz
      INTEGER :: np_old,np_new
      REAL(8) :: rhog,vfac,kap,vpar,pidum,kaphp
      REAL(8) :: b,th,r,enerb,cost,sint,qr,laps,sz,ter,bstar
      REAL(8) :: xt,xs,yt,xdot,ydot,xdt,ydt,zdot,pzdot,edot,pzd0,vp0
      REAL(8) :: pzd1,dpzdt,dthdt,dpsidt,grdgtp,psip2p
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,kaphip,nhip
      REAL(8) :: grp,gxdgyp,rhox(4),rhoy(4),psp,pzp,cvbzp
      REAL(8) :: v,dist,s,s0,rpsi,drdpsi,dfdpsi
      real(8) :: vperp,plam,vpatcl,auglam,dplam,dlam2,y,dydv,dlamdv,ddlamdv

      pidum = 1./(pi*2)**1.5*vwidth**3
      do m=1,mm(4)
         r=xb2(m)-0.5*lx+lr0
         k = int(zb2(m)/delz)
         wz0 = ((k+1)*delz-zb2(m))/delz
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
         grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
                 +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1) 
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
         psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         psp = wx0*psi(i)+wx1*psi(i+1)        
         kaphp = wx0*kapbi(i)+wx1*kapbi(i+1)
         nhip = wx0*nbi(i)+wx1*nbi(i+1)        

         b=1.-tor+tor*bfldp
         pzp = mims(4)*ub2(m)/b*fp/br0-q(4)*psp/br0

         rhog=sqrt(2.*b*mub(m)*mims(4))/(q(4)*b)*iflrh

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
            xs=xb2(m)+rhox(l) !rwx(1,l)*rhog
            yt=yb2(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
!
!   particle can go out of bounds during gyroavg...
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            include "bpushli.h"
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
         vfac = 0.5*(mims(4)*ub2(m)**2 + 2.*mub(m)*b)
         v = sqrt(2.*vfac/mims(4))

         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))

         vpar = ub2(m)
         bstar = b*(1+mims(4)*vpar/(b*q(4))*fp/(b**2*radiusp)*(psip2p*grp**2/radiusp+cvbzp))
         enerb=(mub(m)+mims(4)*vpar*vpar/b)/q(4)*b/bstar*tor
         
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1*b/bstar
         dum = 1/(b*b)*fp/radiusp*cvbzp
         xdot = vxdum*nonlinh -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         xdt = -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1/b+vpar/b*delbyp)*dum1*b/bstar*nonlinh &
             +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0  &
             -mims(4)*vpar**2/(q(4)*bstar*b)*(psip2p*grp**2/radiusp+cvbzp)*lr0/(radiusp*q0) 
         ydt = tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0  &
             -mims(4)*vpar**2/(q(4)*bstar*b)*(psip2p*grp**2/radiusp+cvbzp)*lr0/(radiusp*q0) 
         zdot =  (vpar)*b/bstar &
                 *(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         pzd0 = tor*(-mub(m)/mims(4)/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar
         pzd1 = q(4)/mims(4)*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp 
         pzdot = pzd0 + pzd1*iparah

         edot = q(4)*(xdt*exp1+ydt*eyp+vpar*ezp )

         xb3(m) = xb2(m) + 0.5*dt*xdot
         yb3(m) = yb2(m) + 0.5*dt*ydot
         zb3(m) = zb2(m) + 0.5*dt*zdot
         ub3(m) = ub2(m) + 0.5*dt*pzdot

         vperp = sqrt(2*mub(m)*b/mims(4))
!         v = sqrt(vperp**2+vpar**2)
         plam = vperp**2/v**2/b
         vpatcl = v/vbeam
         auglam=(vpatcl**3+vcrit**3)/(vpatcl**3*(1.0+vcrit**3)) 
         dlam2=dplam0**2+coedlam*(1.0-plam0)*dlog(auglam)
         dplam=sqrt(dlam2)
         y = auglam
         dydv = -3*vbeam**3/(vbeam**3+vibeam**3)*vibeam**3/v**4
         dlamdv = -mub(m)/vfac**2*mims(4)*v
         ddlamdv = coedlam*(1-plam0)/y*dydv/(2*dplam)
         dist = nbeam*nhip/(v**3+vibeam**3)*exp(-(plam-plam0)**2/dlam2)
         dum = v**3+vibeam**3
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1*b/bstar

         wb3(m)=wb2(m) + 0.5*dt*(kaphp*vxdum*dist &
                                 -edot/(v*mims(4))*dist*(-3*v**2/(v**3+vibeam**3)-2*(plam-plam0)/dlam2*dlamdv+2*(plam-plam0)**2/dplam**3*ddlamdv))*dum &
                        -0.5*dt*gambm*wb2(m)

!         if(xb3(m)>lx .or. xb3(m)<0.)wb3(m) = 0.

         go to 333

 333     continue
         laps=anint((zb3(m)/lz)-.5)*(1-peritr)
         r=xb3(m)-0.5*lx+lr0
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         i = max(i,0)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         qr = wx0*sf(i)+wx1*sf(i+1)
         yb3(m)=dmod(yb3(m)-laps*2*pi*qr*lr0/q0*dsign(1.d0,q0)+800.*ly,ly)
         if(xb3(m)>lx)then
            xb3(m) = lx-1.e-8
            zb3(m)=lz-zb3(m)
            zb2(m)=zb3(m)
            xb2(m) = xb3(m)
            yb3(m) = ly*ran2(iseed)
            yb2(m) = yb3(m)
!            wb2(m) = 0.
!            wb3(m) = 0.
         end if
         if(xb3(m)<0.)then
            xb3(m) = 1.e-8
            zb3(m)=lz-zb3(m)
            zb2(m)=zb3(m)
            xb2(m) = xb3(m)
            yb3(m) = ly*ran2(iseed)
            yb2(m) = yb3(m)
!            wb2(m) = 0.
!            wb3(m) = 0.
         end if
         zb3(m)=dmod(zb3(m)+8.*lz,lz)
         xb3(m)=dmod(xb3(m)+8.*lx,lx)         
         xb3(m) = min(xb3(m),lx-1.0e-8)
         yb3(m) = min(yb3(m),ly-1.0e-8)
         zb3(m) = min(zb3(m),lz-1.0e-2)
         
      enddo

!      call MPI_ALLREDUCE(mynopz,nopz,1,MPI_integer, &
!          MPI_SUM, MPI_COMM_WORLD,ierr)

      np_old=mm(4)
      call init_pmove(zb3,np_old,lz,ierr)

      call pmove(xb2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xb3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yb2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yb3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zb2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zb3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ub2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ub3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wb2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wb3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mub,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mub2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ipass,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call end_pmove(ierr)
      mm(4)=np_new

      return
      end
!-----------------------------------------------------------------------

      subroutine bcush(n)

      use gem_com
      use equil
      implicit none
      INTEGER :: n
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,w3old
      INTEGER :: m,i,j,k,ns,l
      INTEGER :: np_old,np_new
      REAL(8) :: vfac,vpar,vxdum,dum,xdot,ydot,xdt,ydt,zdot,pzdot,edot,pzd0,vp0
      REAL(8) :: pzd1,dpzdt,dthdt,dpsidt,grdgtp,psip2p
      REAL(8) :: rhog,xt,yt,zt,kap,xs,pidum,dum1,kaphp
      REAL(8) :: b,th,r,enerb,cost,sint,qr,laps,sz,ter,bstar
      REAL(8) :: myke,mypfl,myavewh
      REAL(8) :: myefl,mynos
      REAL(8) :: ketemp,pfltemp
      REAL(8) :: efltemp,nostemp
      REAL(8) :: sbuf(10),rbuf(10)
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,kaphip,nhip
      REAL(8) :: grp,gxdgyp,rhox(4),rhoy(4),psp,pzp,cvbzp
      REAL(8) :: v,dist,s,s0,rpsi,drdpsi,dfdpsi
      real(8) :: vperp,plam,vpatcl,auglam,dplam,dlam2,y,dydv,dlamdv,ddlamdv

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

      do m=1,mm(4)
         r=xb3(m)-0.5*lx+lr0

         k = int(zb3(m)/delz)
         wz0 = ((k+1)*delz-zb3(m))/delz
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
         grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
                 +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1) 
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
         psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         psp = wx0*psi(i)+wx1*psi(i+1)        
         b=1.-tor+tor*bfldp
         pzp = mims(4)*ub3(m)/b*fp/br0-q(4)*psp/br0
         kaphp = wx0*kapbi(i)+wx1*kapbi(i+1)
         nhip = wx0*nbi(i)+wx1*nbi(i+1)        

         rhog=sqrt(2.*b*mub(m)*mims(4))/(q(4)*b)*iflrh

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
            xs=xb3(m)+rhox(l) !rwx(1,l)*rhog
            yt=yb3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
!   BOUNDARY
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            include "bcushli.h"
 200     continue

         exp1=exp1/4.
         eyp=eyp/4.
         ezp=ezp/4.
         delbxp=delbxp/4.
         delbyp=delbyp/4.
         dpdzp = dpdzp/4.
         dadzp = dadzp/4.
         aparp = aparp/4.


         vfac = 0.5*(mims(4)*ub3(m)**2 + 2.*mub(m)*b)
         v = sqrt(2.*vfac/mims(4))

         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))

         vpar = ub3(m)
         bstar = b*(1+mims(4)*vpar/(b*q(4))*fp/(b**2*radiusp)*(psip2p*grp**2/radiusp+cvbzp))
         enerb=(mub(m)+mims(4)*vpar*vpar/b)/q(4)*b/bstar*tor

         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1*b/bstar
         dum = 1/(b*b)*fp/radiusp*cvbzp
         xdot = vxdum*nonlinh -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         xdt = -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1/b+vpar/b*delbyp)*dum1*b/bstar*nonlinh     &
             +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
             -mims(4)*vpar**2/(q(4)*bstar*b)*(psip2p*grp**2/radiusp+cvbzp)*lr0/(radiusp*q0) 
         ydt = tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
             -mims(4)*vpar**2/(q(4)*bstar*b)*(psip2p*grp**2/radiusp+cvbzp)*lr0/(radiusp*q0) 
         zdot =  (vpar)*b/bstar &
                 *(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         pzd0 = tor*(-mub(m)/mims(4)/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar
         pzd1 = q(4)/mims(4)*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp 
         pzdot = pzd0 + pzd1*iparah

         edot = q(4)*(xdot*exp1+ydot*eyp+vpar*ezp )

         xb3(m) = xb2(m) + dt*xdot
         yb3(m) = yb2(m) + dt*ydot
         zb3(m) = zb2(m) + dt*zdot
         ub3(m) = ub2(m) + dt*pzdot

         vperp = sqrt(2*mub(m)*b/mims(4))
!         v = sqrt(vperp**2+vpar**2)
         plam = vperp**2/v**2/b
         vpatcl = v/vbeam
         auglam=(vpatcl**3+vcrit**3)/(vpatcl**3*(1.0+vcrit**3)) 
         dlam2=dplam0**2+coedlam*(1.0-plam0)*dlog(auglam)
         dplam=sqrt(dlam2)
         y = auglam
         dydv = -3*vbeam**3/(vbeam**3+vibeam**3)*vibeam**3/v**4
         dlamdv = -mub(m)/vfac**2*mims(4)*v
         ddlamdv = coedlam*(1-plam0)/y*dydv/(2*dplam)
         dist = nbeam*nhip/(v**3+vibeam**3)*exp(-(plam-plam0)**2/dlam2)
         dum = v**3+vibeam**3
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1*b/bstar

         w3old = wb3(m)
         wb3(m)=wb2(m) + dt*(kaphp*vxdum*dist &
                                 -edot/(v*mims(4))*dist*(-3*v**2/(v**3+vibeam**3)-2*(plam-plam0)/dlam2*dlamdv+2*(plam-plam0)**2/dplam**3*ddlamdv))*dum &
                        -dt*gambm*wb2(m)

         if(abs(wb3(m)).gt.nbeam.and.nonlinh==1)then
!            wb3(m) = 0.
!            wb2(m) = 0.
         end if


         laps=anint((zb3(m)/lz)-.5)*(1-peritr)
         r=xb3(m)-0.5*lx+lr0
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         i = max(i,0)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         qr = wx0*sf(i)+wx1*sf(i+1)
         yb3(m)=dmod(yb3(m)-laps*2*pi*qr*lr0/q0*dsign(1.d0,q0)+800.*ly,ly)
         if(xb3(m)>lx)then
            xb3(m) = lx-1.e-8
            zb3(m)=lz-zb3(m)
            zb2(m)=zb3(m)
            xb2(m) = xb3(m)
            yb3(m) = ly*ran2(iseed)
            yb2(m) = yb3(m)
!            wb2(m) = 0.
!            wb3(m) = 0.
         end if
         if(xb3(m)<0.)then
            xb3(m) = 1.e-8
            zb3(m)=lz-zb3(m)
            zb2(m)=zb3(m)
            xb2(m) = xb3(m)
            yb3(m) = ly*ran2(iseed)
            yb2(m) = yb3(m)
!            wb2(m) = 0.
!            wb3(m) = 0.
         end if
         zb3(m)=dmod(zb3(m)+8.*lz,lz)
         xb3(m)=dmod(xb3(m)+8.*lx,lx)         
         xb3(m) = min(xb3(m),lx-1.0e-8)
         yb3(m) = min(yb3(m),ly-1.0e-8)
         zb3(m) = min(zb3(m),lz-1.0e-2)
         if(abs(zb3(m)-zb2(m)).gt.lz/2)ipass(m)=1
         if(ipass(m)==1 .and. vpar>0)then
            wb2(m) = 0.
            wb3(m) = 0.
         end if

!     particle diagnostics done here because info is available...
         mypfl=mypfl + w3old*(eyp+vpar*delbxp/b) 
         myefl=myefl + vfac*w3old*(eyp+vpar*delbxp/b)
         myke=myke + vfac !*wb3(m)
         mynos=mynos + wb3(m)
         myavewh = myavewh+abs(wb3(m))

!     xn+1 becomes xn...
         ub2(m)=ub3(m)
         xb2(m)=xb3(m)
         yb2(m)=yb3(m)
         zb2(m)=zb3(m)
         wb2(m)=wb3(m)

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
      avewb(n) = rbuf(5)/( float(tmm(1)) )
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      nos(4,n)=nostemp/( float(tmm(1)) )
      pfl(4,n)=pfltemp/( float(tmm(1)) )
      efl(4,n)=efltemp/( float(tmm(1)) )
      ke(4,n)=ketemp/(float(tmm(1)))
      np_old=mm(4) 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call init_pmove(zb3,np_old,lz,ierr)

      call pmove(xb2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xb3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yb2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yb3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zb2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zb3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ub2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ub3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wb2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wb3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mub,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mub2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ipass,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call end_pmove(ierr)
      mm(4)=np_new
!     write(*,*)MyId,mm(1)

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine loadb

      use gem_com
      use equil
      implicit none
      INTEGER :: i,j,k,m,idum,ns,m1,nvgrd=100
      INTEGER :: np_old,np_new
      REAL(8) :: vpar,vperp2,r,qr,th,b,cost,ter,x
      REAL(8) :: avgv,myavgv,avgw,myavgw
      real(8) :: dumx,dumy,dumz,jacp
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,psp
      REAL(8) :: grp,gxdgyp,zoldp
      REAL(8) :: wx0,wx1,wz0,wz1
      REAL(8) :: vx,vy,v,vg(0:10000),dv,f1,f2

      cnt=int(tmm(1)/numprocs)
!   write(*,*)'in loader cnt, mm(1)=',cnt,mm(1)

      myavgv=0.
      avgv=0.
      avgw = 0.
      myavgw = 0.

!find a set of v grids for slowing-down distribution
      dv = (vcbeam-vmin)/nvgrd
      vg(0) = vmin
      vg(1) = vmin+dv
      do i=1,1000000
         f1 = 1./(vg(i-1)**3+vibeam**3)
         f2 = 1./(vg(i)**3+vibeam**3)
         vg(i+1) = vg(i)+(vg(i)-vg(i-1))*vg(i-1)**2*f1/(vg(i)**2*f2)
         if(vg(i+1)>vcbeam)then
            nvgrd = i
            goto 150
         end if
      end do
 150  continue
!compute x=int v^4/(v^3+vi^3)
      x = 0.
      dv = (vcbeam-vmin)/nvgrd
      do i = 0, nvgrd-1
         v = vmin+(i+0.5)*dv
         x = x+v**4/(v**3+vibeam**3)*dv
      end do
!      if(myid==0)write(*,*)'nvgrd, x= ', nvgrd,x*pi*mims(4)*nbeam

      mm(4)=int(tmm(1)/numprocs)
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
         v = vg(int(ran2(iseed)*nvgrd))+(ran2(iseed)-0.5)*dv
         vpar = 2.*v*(ran2(iseed)-0.5)

         if((ran2(iseed)<(0.5*jacp/jacmax)).and.(v<vcbeam).and.(v>vmin))then
         m = m+1
         if(m>mm(4))goto 170
         xb2(m)=min(dumx,lx-1.d-8)
         yb2(m)=min(dumy,ly-1.d-8)

         k = int((th+pi)/dth)
         wz0 = (-pi+(k+1)*dth-th)/dth
         wz1 = 1-wz0
         zb2(m) = wz0*zfnth(k)+wz1*zfnth(k+1)
         zb2(m)=min(zb2(m),lz-1.d-8)

!   normalizations will be done in following loop...

         r=xb2(m)-0.5*lx+lr0
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
         fp = wx0*f(i)+wx1*f(i+1)        
!         ter = tets(4) !wx0*t0i(i)+wx1*t0i(i+1)
         b=1.-tor+tor*bfldp

         ub2(m)=vpar
         mub(m)=0.5*mims(4)*(v**2-vpar**2)/b
         mub2(m) = mub(m)
         myavgv=myavgv+ub2(m)

!    LINEAR: perturb w(m) to get linear growth...
         wb2(m)=2.*amp*nbeam/cvbeam*(revers(MyId*cnt+m,13)-0.5) !(ran2(iseed) - 0.5 )
         myavgw=myavgw+wb2(m)
         ipass(m) = 0
         end if
 160  continue
 170  continue
      myavgw = myavgw/mm(4)
!      write(*,*)'myid ', myid,x2(10),y2(20),u2(20),w2(20)
!    subtract off avg. u...
      call MPI_ALLREDUCE(myavgv,avgv,1, &
          MPI_REAL8, &
          MPI_SUM,MPI_COMM_WORLD,ierr)
      if(idg.eq.1)write(*,*)'all reduce'
      avgv=avgv/float(tmm(1))
      do 180 m=1,mm(4)
         ub2(m)=ub2(m)-avgv
         xb3(m)=xb2(m)
         yb3(m)=yb2(m)
         zb3(m)=zb2(m)
         ub3(m)=ub2(m)
         wb3(m)=wb2(m)
 180  continue

      np_old=mm(4)
      call init_pmove(zb3,np_old,lz,ierr)
      if(idg.eq.1)write(*,*)'pass init_pmove'
!     
      call pmove(xb2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xb3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yb2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yb3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zb2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zb3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ub2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ub3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wb2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wb3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mub,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mub2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ipass,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
!     
      call end_pmove(ierr)
      mm(4)=np_new
!     write(*,*)MyId,j,mm(4)

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine collb(n)
      use gem_com
      use equil
      implicit none
      integer :: i,ip,k,m,n,ncol,icol,isdum,np_old,np_new
      real(8) :: edum,vdum,dum,dum1,ptch,vte,r,qr,th,cost,b,psp
      real(8) :: h_x,h_coll,x,eps,dtcol,uold,hee,nue,ter
      real(8) :: wx0,wx1,wz0,wz1

      ncol = 1
      dtcol = dt/ncol
      if(rneu==0.0)return
      do k = 1,mm(4)
         r=xb3(k)-0.5*lx+lr0

         m = int(zb3(k)/delz)
         wz0 = ((m+1)*delz-zb3(k))/delz
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
         psp = wx0*psi(i)+wx1*psi(i+1)
         uold = ub3(k)
         edum = b*mub(k)+0.5*mims(4)*ub3(k)*ub3(k)
         vdum = sqrt(2.*edum/mims(4))
         ptch = ub3(k)/vdum

! collision frequency for experimental profiles

         nue=rneu*vi**3/(vdum**3+1.0)
         dum = mims(1)/mims(4)*dtcol*nue

         do icol = 1,ncol
            ptch =ptch-dum*ptch+sqrt(dum*(1.-ptch*ptch)) &
                *dsign(1.d0,ran2(iseed)-0.5)
            ptch = dmin1(ptch,0.999d0)
            ptch = dmax1(ptch,-0.999d0)

!            vdum = vdum-rneu*(vdum+vi**3/vdum**2)*dtcol*islow
         end do
         ub3(k) = vdum*ptch
         mub(k) = 0.5*mims(4)*vdum*vdum*(1.-ptch*ptch)/b
         mub2(k) = mub(k)
         ub2(k) = ub3(k)
      end do

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pint
      use gem_com
      use equil
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dgdtp,dpdzp,dadzp,aparp,vncp
      REAL(8) :: dum,vxdum,dum1,dum2,eps,x,h_x,h_coll,hee,nue
      INTEGER :: m,i,j,k,l,n,ipover,ieover
      INTEGER :: np_old,np_new
      REAL(8) :: vfac,kap,kapnp,kaptp,sz,vpar,ppar,vpdum,pidum,xnp
      REAL(8) :: b,th,r,enerb,cost,sint,qr,laps,ter
      REAL(8) :: xt,xs,yt,zt,xdot,ydot,zdot,pzdot,edot,pzd0,pzd1,vp0
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp
      REAL(8) :: grp,gxdgyp,psp,pzp,psip2p,bdcrvbp,curvbzp,dipdrp,bstar
      REAL(8) :: myavptch,myaven
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

         curvbzp = wx0*wz0*curvbz(i,k)+wx0*wz1*curvbz(i,k+1) &
                 +wx1*wz0*curvbz(i+1,k)+wx1*wz1*curvbz(i+1,k+1) 
         bdcrvbp = wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
                 +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1) 
         grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
                 +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1) 

         fp = wx0*f(i)+wx1*f(i+1)       
         jfnp = wz0*jfn(k)+wz1*jfn(k+1) 
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         psp = wx0*psi(i)+wx1*psi(i+1)        
         ter = wx0*t0e(i)+wx1*t0e(i+1)        
         kaptp = wx0*capte(i)+wx1*capte(i+1)        
         kapnp = wx0*capne(i)+wx1*capne(i+1)        
         xnp = wx0*xn0e(i)+wx1*xn0e(i+1)        
!         vncp = wx0*phincp(i)+wx1*phincp(i+1)        
         psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
         dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)        
         b=1.-tor+tor*bfldp
         pzp = emass*u2e(m)/b+psp/br0

         xt = x2e(m)
         yt = y2e(m)

         include 'ppushlie.h'

         vfac = 0.5*(emass*u2e(m)**2 + 2.*mue2(m)*b)
!         vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
!         vp0 = vp0*vncp*vexbsw
         kap = kapnp - (1.5-vfac/ter)*kaptp

         ppar = u2e(m)
         vpar = u2e(m)-qel/emass*aparp*ipara
         bstar = b*(1+emass*vpar/(qel*b)*bdcrvbp)
         enerb=(mue2(m)+emass*vpar*vpar/b)/qel*b/bstar*tor
         dum2 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp*iflut)*dum2
         xdot = vxdum*nonline  &
             -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1/b+vpar/b*delbyp*iflut)*dum2*nonline &
             +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0*0. &
             +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
             -emass*vpar**2/(qel*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
             -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*grcgtp*lr0/q0*qhatp  

         zdot =  vpar*b/bstar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
                 -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*q0*br0*grcgtp/jfnp

         pzd0 = tor*(-mue2(m)/emass/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar &
                +mue2(m)*vpar/(qel*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp
         pzd1 = qel/emass*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
             +qel/emass*(-xdot*delbyp+ydot*delbxp+zdot*dadzp) 
         pzdot = pzd0+pzd1*ipara
                      
         vpdum = vpar
         edot = qel*(xdot*exp1+(ydot-vp0)*eyp+vpar*ezp)

         x3e(m) = x2e(m) + 0.5*dte*xdot
         y3e(m) = y2e(m) + 0.5*dte*ydot
         z3e(m) = z2e(m) + 0.5*dte*zdot
         u3e(m) = u2e(m) + 0.5*dte*pzdot
         mue3(m) = mue2(m)

         eps = (b*mue2(m)+0.5*emass*u2e(m)*u2e(m))/ter
         x = sqrt(eps)
         h_x    = 4.0*x*x/(3.0*sqrt(pi))
         h_coll = h_x/sqrt(1.0+h_x**2)
!         h_coll = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*DERF(x)
! collision frequency for experimental profiles
         hee = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*erf(x)
         nue=rneu*(wx0*nue0(i)+wx1*nue0(i+1))*(wx0*zeff(i)+wx1*zeff(i+1)+hee)
!         dum1 = 0.5*rneu/eps**1.5*(1+h_coll)
         dum1 = 0.5*nue/(eps+0.1)**1.5*(1+h_coll)
!         if(x<0.3)dum1=0.0

         dum = 1-w2e(m)*nonline*0.
         if(isuni.eq.1)dum = exp(-vfac)*pidum
         vxdum = (eyp/b+vpdum/b*delbxp)*dum2
         w3e(m)=w2e(m) + 0.5*dte*(vxdum*kap + edot/ter)*xnp*(tload/ter)**1.5*exp(vfac*(1/tload-1./ter))-0.5*dt*gamele*w2e(m)

!         if(x3e(m)>lx .or. x3e(m)<0.)w3e(m) = 0. 

         laps=anint((z3e(m)/lz)-.5)*(1-peritr)
         r=x3e(m)-0.5*lx+lr0
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         i = max(i,0)
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
         x3e(m)=dmod(x3e(m)+80.*lx,lx)         
         x3e(m) = min(x3e(m),lx-1.0e-8)
         y3e(m) = min(y3e(m),ly-1.0e-8)
         z3e(m) = min(z3e(m),lz-1.0e-8)

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
      call pmove(u0e,np_old,np_new,ierr)
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
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dgdtp,dpdzp,dadzp,aparp,vncp
      REAL(8) :: dum,vxdum,dum1,dum2,eps,x,h_x,h_coll,hee,nue
      INTEGER :: m,i,j,k,l,n,mynowe
      INTEGER :: np_old,np_new
      REAL(8) :: vfac,kap,kapnp,kaptp,sz,vpar,ppar,vpdum,pidum,xnp
      REAL(8) :: b,th,r,enerb,cost,sint,qr,laps,ter
      REAL(8) :: xt,xs,yt,zt,xdot,ydot,zdot,pzdot,edot,pzd0,pzd1,vp0
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,w3old
      REAL(8) :: myke,mypfl(1:nsubd),myptrp
      REAL(8) :: myefl(1:nsubd),mynos,myavewe
      REAL(8) :: ketemp,pfltemp,ptrptmp
      REAL(8) :: efltemp,nostemp
      REAL(8) :: sbuf(10),rbuf(10)
      REAL(8) :: mytotn,mytrap,totn,ttrap
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp
      REAL(8) :: grp,gxdgyp,psp,pzp,psip2p,bdcrvbp,curvbzp,dipdrp,bstar
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

         curvbzp = wx0*wz0*curvbz(i,k)+wx0*wz1*curvbz(i,k+1) &
                 +wx1*wz0*curvbz(i+1,k)+wx1*wz1*curvbz(i+1,k+1) 
         bdcrvbp = wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
                 +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1) 
         grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
                 +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1) 

         fp = wx0*f(i)+wx1*f(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         psp = wx0*psi(i)+wx1*psi(i+1)        
         ter = wx0*t0e(i)+wx1*t0e(i+1)        
         kaptp = wx0*capte(i)+wx1*capte(i+1)        
         kapnp = wx0*capne(i)+wx1*capne(i+1)        
         xnp = wx0*xn0e(i)+wx1*xn0e(i+1)        
!         vncp = wx0*phincp(i)+wx1*phincp(i+1)        
         psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
         dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)        
         b=1.-tor+tor*bfldp
         pzp = emass*u3e(m)/b+psp/br0

         xt = x3e(m)
         yt = y3e(m)

         include 'cpushlie.h'

         vfac = 0.5*(emass*u3e(m)**2 + 2.*mue3(m)*b)
!         vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
!         vp0 = vp0*vncp*vexbsw
         kap = kapnp - (1.5-vfac/ter)*kaptp

         ppar = u3e(m)
         vpar = u3e(m)-qel/emass*aparp*ipara
         bstar = b*(1+emass*vpar/(qel*b)*bdcrvbp)
         enerb=(mue3(m)+emass*vpar*vpar/b)/qel*b/bstar*tor
         dum2 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp*iflut)*dum2
         xdot = vxdum*nonline &
             -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1/b+vpar/b*delbyp*iflut)*dum2*nonline  &
             +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0*0. &
             +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
             -emass*vpar**2/(qel*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
             -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*grcgtp*lr0/q0*qhatp  

         zdot =  vpar*b/bstar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
                 -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*q0*br0*grcgtp/jfnp

         pzd0 = tor*(-mue3(m)/emass/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar &
                +mue3(m)*vpar/(qel*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp
         pzd1 = qel/emass*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
             +qel/emass*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)
         pzdot = pzd0+pzd1*ipara

         vpdum = vpar
         edot = qel*(xdot*exp1+(ydot-vp0)*eyp+vpar*ezp)

         x3e(m) = x2e(m) + dte*xdot
         y3e(m) = y2e(m) + dte*ydot
         z3e(m) = z2e(m) + dte*zdot
         u3e(m) = u2e(m) + dte*pzdot
         mue3(m) = mue2(m)

         eps = (b*mue3(m)+0.5*emass*u3e(m)*u3e(m))/ter
         x = sqrt(eps)
         h_x    = 4.0*x*x/(3.0*sqrt(pi))
         h_coll = h_x/sqrt(1.0+h_x**2)
!         h_coll = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*DERF(x)
! collision frequency for experimental profiles
         hee = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*erf(x)
         nue=rneu*(wx0*nue0(i)+wx1*nue0(i+1))*(wx0*zeff(i)+wx1*zeff(i+1)+hee)
!         dum1 = 0.5*rneu/eps**1.5*(1+h_coll)
         dum1 = 0.5*nue/(eps+0.1)**1.5*(1+h_coll)
!         if(x<0.3)dum1=0.0

         dum = 1-w3e(m)*nonline*0.
         if(isuni.eq.1)dum = exp(-vfac)*pidum
         vxdum = (eyp/b+vpdum/b*delbxp)*dum2
!             -(2*u3e(m)*aparp+qel/emass*aparp*aparp)/(b*b)/br0*sint
         w3old = w3e(m)
         w3e(m)=w2e(m) + dte*(vxdum*kap + edot/ter)*xnp*(tload/ter)**1.5*exp(vfac*(1/tload-1./ter))-dt*gamele*w2e(m)

         if(abs(w3e(m)).gt.1.0.and.nonline==1)then
            w3e(m) = 0.
            w2e(m) = 0.
            mynowe = mynowe+1
         end if


         mytotn = mytotn+1

         laps=anint((z3e(m)/lz)-.5)*(1-peritr)
         r=x3e(m)-0.5*lx+lr0
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         i = max(i,0)
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
         x3e(m)=dmod(x3e(m)+80.*lx,lx)         
         x3e(m) = min(x3e(m),lx-1.0e-8)
         y3e(m) = min(y3e(m),ly-1.0e-8)
         z3e(m) = min(z3e(m),lz-1.0e-8)

!     particle diagnostics done here because info is available...
         k = int(x3e(m)/(lx/nsubd))
         k = min(k,nsubd-1)
         k = k+1
         mypfl(k)=mypfl(k) + w3old*(eyp+vpar*delbxp/b) 
         myefl(k)=myefl(k) + vfac*w3old*(eyp+vpar*delbxp/b)
         myke=myke + vfac*w3e(m)
         mynos=mynos + w3e(m)
         myavewe = myavewe+abs(w3e(m))

         u2e(m)=u3e(m)
         x2e(m)=x3e(m)
         y2e(m)=y3e(m)
         z2e(m)=z3e(m)
         w2e(m)=w3e(m)

      enddo
      call MPI_ALLREDUCE(mynowe,nowe,1,MPI_integer, &
          MPI_SUM, MPI_COMM_WORLD,ierr)

      sbuf(1)=myke
      sbuf(2)=myefl(nsubd/2)
      sbuf(3)=mypfl(nsubd/2)
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
      avewe(n) = rbuf(8)/( float(tmme) )
!      ke(2,n)=ketemp/( 2.*float(tmm(1))*mims(ns) )
      ftrap = ttrap/totn
!      nos(2,n)=nostemp/( float(tmm(1)) )

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      sbuf(1:nsubd) = myefl(1:nsubd)
      call MPI_ALLREDUCE(sbuf,rbuf,10,  &
          MPI_REAL8,MPI_SUM,  &
          MPI_COMM_WORLD,ierr)
      do k = 1,nsubd   
         efle(k,n)=rbuf(k)/( float(tmme) )*totvol/vol(k)*cn0e
      end do
      sbuf(1:nsubd) = mypfl(1:nsubd)
      call MPI_ALLREDUCE(sbuf,rbuf,10,  &
          MPI_REAL8,MPI_SUM,  &
          MPI_COMM_WORLD,ierr)
      do k = 1,nsubd   
         pfle(k,n)=rbuf(k)/( float(tmme) )*totvol/vol(k)*cn0e
      end do

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
      call pmove(u0e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call end_pmove(ierr)
      mme=np_new

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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

      cnt=int(tmme/numprocs)
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

         u2e(m)=vpar*sqrt(amie*tload)
         mue(m)=0.5*vperp2/b*tload
         eke(m) = mue(m)*b+0.5*emass*u2e(m)**2
         pze(m) = emass*u2e(m)/b+psp/br0
         z0e(m) = z2e(m)
         xie(m) = x2e(m)
         u0e(m) = u2e(m)
         myavgv=myavgv+u2e(m)

!    LINEAR: perturb w(m) to get linear growth...
         w2e(m)=2.*amp*(revers(MyId*cnt+m,13)-0.5) !(ran2(iseed) - 0.5 )
!         w2e(m) = amp*u2e(m)/(abs(u2e(m))+0.1)*sin(x2e(m)*pi2/lx)
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
      call pmove(u0e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

!     
      call end_pmove(ierr)
      mme=np_new

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine setw(ip,n)
      use gem_com
      use equil
      implicit none
      INTEGER :: m,n,ip,i,j,k
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,vte
      REAL(8) :: wght
      REAL(8) :: xt,yt,zt

      vte = sqrt(amie*t0e(nr/2))
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
!         if(abs(u3e(m)/vte).gt.vcut)wght = 0.

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
      subroutine laplace(u)   

      use gem_com
      use equil
      use fft_wrapper
      implicit none
      REAL(8) :: b2,gam0,gam1,th,delyz,bf,rkper,r,qr,shat,ter !,b
      REAL(8) :: kx,ky
      complex(8) :: lbfr(0:imx-1,0:jcnt-1)
      complex(8) :: lbfs(0:imx-1,0:jcnt-1)
      complex(8) :: rbfr(0:imx-1,0:jcnt-1)
      complex(8) :: rbfs(0:imx-1,0:jcnt-1)
      REAL(8),dimension(:),allocatable :: dely,aky
      complex(8),dimension(:,:,:,:),allocatable:: nab1,nab2
      complex(8),dimension(:,:,:,:),allocatable :: mx
      REAL(8),dimension(:,:,:),allocatable :: formlap
      REAL(8) :: sgnx,sgny,sz,myfe
      INTEGER :: i,i1,j,j1,k,k1,l,m,n,ifirst,nstep,ip,INFO
      INTEGER :: l1,m1,myk,myj,izonal,ix,ikx,iext
      COMPLEX(8) :: temp3dxy(0:imx-1,0:jmx-1,0:1),va(0:imx-1,0:jmx-1,0:1),&
            v(0:imx-1,0:jcnt-1,0:1)
      COMPLEX(8) :: sbuf(0:imx*jcnt*2-1),rbuf(0:imx*jmx*2-1)
      COMPLEX(8) :: sl(1:imx-1,0:jcnt-1,0:1)
      REAL(8) :: dum,u(0:imx,0:jmx,0:1)
      complex(8) :: ua(0:imx-1,0:jcnt-1,0:1),calph,cbeta
      REAL(8) :: grp,gthp,gxdgyp,dydrp,qhatp,grdgtp,bfldp,g2xp,g2yp, &
                 radiusp,g2zp,gzp,gxdgzp,gydgzp
      REAL(8) :: wx0,wx1,wz0,wz1
      character*70 fname
      character*4 holdmyid

      save formlap,ifirst,dely,aky,mx
      calph = 1.
      cbeta = 0.
      write(holdmyid,'(I4.4)') MyId
      fname='./matrix/'//'mx_lap_'//holdmyid
      if (ifirst.ne.-99) then
         allocate(dely(0:imx-1),aky(0:jmx-1),nab1(0:imx-1,0:jcnt-1,0:imx-1,0:1))
         allocate(nab2(0:imx-1,0:jcnt-1,0:imx-1,0:1))
         allocate(mx(imx-1,imx-1,0:jcnt-1,0:1),formlap(0:imx-1,0:jcnt-1,0:1))

         if(iget.eq.1) then
            open(20000+MyId,file=fname,form='unformatted',status='old')
            read(20000+MyId)mx
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
                     ter = wx0*t0i(i1)+wx1*t0i(i1+1)
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
                     radiusp = wx0*wz0*radius(i1,k1)+wx0*wz1*radius(i1,k1+1) &
                     +wx1*wz0*radius(i1+1,k1)+wx1*wz1*radius(i1+1,k1+1)

                     m1 = mstart+int((float(m)+1.0)/2)
                     if(m==0)m1=0
                     sgny = isgnft(m)

                     ky=sgny*2.*pi*float(m1)/ly
                     kx=pi*float(l)/lx
                     bf=bfldp
!                     b=mims(1)*(kx*kx*grp**2 + &
!                        ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
!                        +2*dydrp*lr0/q0*qhatp*grdgtp) &
!                        +2*kx*ky*gxdgyp)/bf/bf
                  
                     nab1(l,m,i,n) = kx**2*grp**2+ky**2*  &
                        (dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                        +2*dydrp*lr0/q0*qhatp*grdgtp)

                     nab2(l,m,i,n) = -IU*ky*kx*2*gxdgyp
!   formfactor in ezamp
                     formlap(l,m,n) = -1./jmx 
                     if(abs(ky)>kycut)formlap(l,m,n) = 0.
                     if(m1.ne.mlk.and.onemd==1)formlap(l,m,n) = 0.
 53               continue
 54            continue
 52         continue
 51      continue

         do k = 0,1
            do j = 0,jcnt-1
               do i = 1,imx-1
                  do ix = 1,imx-1
                     mx(i,ix,j,k) = 0.
!                     if(i==ix)mx(i,ix,j,k) = beta*amie*(1+nh*ish)
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

         if(iget.eq.0) then
            open(20000+MyId,file=fname,form='unformatted',status='unknown')
            write(20000+MyId)mx
            goto 200
         end if

 200     ifirst=-99
      endif

!   now do field solve...

      call dcmpy(u(0:imx-1,0:jmx-1,0:1),v)

      temp3dxy = 0.
      do k = 0,1
         do j = 0,jcnt-1  !no n=0 !!!
            myj=jft(j)
            call ZGEMV('N',imx-1,imx-1,calph,mx(:,:,j,k),imx-1, &
                       v(1:imx-1,j,k), 1,cbeta,sl(1:imx-1,j,k), 1)
            temp3dxy(1:imx-1,myj,k) = sl(1:imx-1,j,k)
            temp3dxy(0,myj,k) = 0.
         end do
      end do

!  from rho(kx,ky) to phi(kx,ky)
      do k = 0,1
         do j = 0,jmx-1
            do i = 0,imx-1
               temp3dxy(i,j,k)=-temp3dxy(i,j,k)/jmx
            end do
         end do
      end do

!  from phi(kx,ky) to phi(x,y)
      do k=0,mykm
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3dxy(i,j,k)
            end do
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3dxy(i,j,k) = tmpy(j)  ! phi(x,y)
            end do
         end do
      end do

      do i = 0,nxpp-1
         do j = 0,jm-1
            do k = 0,mykm
                  u(i,j,k) = temp3dxy(i,j,k)
            end do
         end do
      end do

!    x-y boundary points 
      call enfxy(u(:,:,:))
      call enfz(u(:,:,:))
!      call filter(u(:,:,:))

      return
      end

!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine jie814(ip,n)

!    source quantities are are calculated: n_i
!    right now only ion quantitities are calculated...

      use gem_com
      use equil
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: enerb,vxdum,dum,xdot,ydot,zdot,pzd0,pidum,dum1,dum2
      INTEGER :: m,n,i,j,k,l,ns,ip,nonfi,nonfe
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,vte
      REAL(8) :: sz,wght,wght0,wght1,wght2,r,th,cost,sint,b,qr,dv,kap,ter
      REAL(8) :: kapnp,kaptp,xnp,psip2p,bdcrvbp,curvbzp,dipdrp,bstar
      REAL(8) :: xt,yt,rhog,vpar,xs,dely,vfac,vp0
      REAL(8) :: lbfs(0:imx,0:jmx)
      REAL(8) :: rbfs(0:imx,0:jmx)
      REAL(8) :: lbfr(0:imx,0:jmx)
      REAL(8) :: rbfr(0:imx,0:jmx)
      real(8) :: myupar(0:imx,0:jmx,0:1)
      real(8) :: myupex(0:imx,0:jmx,0:1),myupey(0:imx,0:jmx,0:1)
      real(8) :: myupazd(0:imx,0:jmx,0:1),mydene(0:imx,0:jmx,0:1)
      real(8) :: mydnedt(0:imx,0:jmx,0:1)
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp
      REAL(8) :: grp,gxdgyp,rhox(4),rhoy(4),vncp,vparspp

      nonfe = 1 

! electrons current
      pidum = 1./(pi*2)**1.5*(vwidthe)**3
      if(isuni.eq.0)pidum = 1.

      vte = sqrt(amie*t0e(nr/2))
      ppex = 0.
      ppey = 0.
      djedt = 0.
      myupex = 0.
      myupey = 0.
      myupazd = 0.
      mydnedt = 0.
      mydene = 0.
      myupar = 0.

      if(iske==0 .and. ieqmo814==1)then
         do  i=0,im
            do  j=0,jm
               do  k=0,mykm
                  ppazd(i,j,k) = gt0e(i)*dene(i,j,k)*bdgrzn(i,k)/emass/ntube 
               end do
            end do
         end do
         return
      end if

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

         curvbzp = wx0*wz0*curvbz(i,k)+wx0*wz1*curvbz(i,k+1) &
                 +wx1*wz0*curvbz(i+1,k)+wx1*wz1*curvbz(i+1,k+1) 
         bdcrvbp = wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
                 +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1) 
         grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
                 +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1) 

         fp = wx0*f(i)+wx1*f(i+1)        
         jfnp = wz0*jfn(k)+wz1*jfn(k+1)
         psipp = wx0*psip(i)+wx1*psip(i+1)        
         ter = wx0*t0e(i)+wx1*t0e(i+1)        
         kaptp = wx0*capte(i)+wx1*capte(i+1)        
         kapnp = wx0*capne(i)+wx1*capne(i+1)        
         xnp = wx0*xn0e(i)+wx1*xn0e(i+1)        
         b=1.-tor+tor*bfldp
         psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
         dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)        

         vfac = 0.5*(emass*u3e(m)**2 + 2.*mue3(m)*b)
         kap = kapnp - (1.5-vfac/ter)*kaptp
!!!! attention, ifrzt=1 assumed
!         kap = kapnp
         wght=w3e(m)/dv
         wght0=exp(-vfac)/dv
         if(isuni.eq.0)wght0=1./dv
         wght0 = wght0
         vpar = u3e(m)
         if(abs(vfac/ter).gt.vcut)then
            wght = 0.
            wght0 = 0.
         endif

         bstar = b*(1+emass*vpar/(qel*b)*bdcrvbp)
         enerb=(mue3(m)+emass*vpar*vpar/b)/qel*b/bstar*tor

         xdot = -tor*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = +tor*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
             +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
             -emass*vpar**2/(qel*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
             -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*grcgtp*lr0/q0*qhatp  

         zdot =  u3e(m)*b/bstar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
                 -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*q0*br0*grcgtp/jfnp

         pzd0 = tor*(-mue2(m)/emass/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar &
                +mue2(m)*vpar/(qel*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp

         xt = x3e(m)
         yt = y3e(m)
         i=int(xt/dx)
         j=int(yt/dy)
         k=0 !int(z3e(m)/dz)-gclr*kcnt

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
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
         xdot = xdot+vxdum*nonline
         ydot = ydot+(-exp1/b+vpar/b*delbyp)*dum1*nonline

         wght1 = vpar*wght
         wght2 = wght*vpar*zdot

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
         vxdum = (vpar/b*delbxp)*dum2 
         wght0 = wght0*vpar*(vxdum*kap)*xnp*(tload/ter)**1.5*exp(vfac*(1/tload-1./ter))+(wght)*pzd0 !*(1-ifrzt)

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

      do  i=0,im
         do  j=0,jm
            do  k=0,mykm
               ppex(i,j,k) = myupex(i,j,k)/n0e/jac(i,k)*cn0e
               ppey(i,j,k) = myupey(i,j,k)/n0e/jac(i,k)*cn0e
               ppazd(i,j,k) = myupazd(i,j,k)/n0e/jac(i,k)*cn0e
               djedt(i,j,k) = mydnedt(i,j,k)/n0e/jac(i,k)*cn0e
            end do
         end do
      end do

      if(ifrzt==0)goto 999
!ifrzt==1 assumed
      if(idenek==0)then
         do  i=0,im
            do  j=0,jm
               do  k=0,mykm
                  ppazd(i,j,k) = gt0e(i)*dene(i,j,k)*bdgrzn(i,k)/emass/ntube 
               end do
            end do
         end do
      end if

      if(idenek==1)then
         do  i=0,im
            do  j=0,jm
               do  k=0,mykm
                  ppazd(i,j,k) = gt0e(i)*denek(i,j,k)*bdgrzn(i,k)/emass/ntube 
               end do
            end do
         end do
      end if

 999  continue
      call filter(ppazd(:,:,:))

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine drdt814(ip)
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
      real(8) :: dumpp(0:imx,0:jmx,0:1),dumx(0:imx,0:jmx,0:1),dumy(0:imx,0:jmx,0:1)
      real(8) :: u(0:imx,0:jmx,0:1)

      do i = 0,im
         do j = 0,jm
            do k = 0,1
               u(i,j,k) = ezs(i,j,k) !*bdgrzn(i,k)
               dumpp(i,j,k) = gn0e(i)*phi(i,j,k)-gt0e(i)*dene(i,j,k)-gn0e(i)*delte(i,j,k)
            end do
         end do
      end do
!      call laplace(u)
!      call neq0(dumpp)

      call gradx(-ppex,djdx)
      call grady(-ppey,djdy)
      call gradx(delbx*dumpp,dumx)
      call grady(delby*dumpp,dumy)
      call neq0(dumx)
      call neq0(dumy)

      isdndt = 1
      do i = 0,im
         do j = 0,jm
            rbfs(i,j)=(-ppazd(i,j,0))*jac(i,0) !+1./(amie*beta)*u(i,j,0)
         end do   
      end do
      call MPI_SENDRECV(rbfs,(imx+1)*(jmx+1),  &
          MPI_REAL8,rngbr,404, &
          lbfr,(imx+1)*(jmx+1), &
          MPI_REAL8,lngbr,404, &
          tube_comm,stat,ierr)
      do i = 0,im
         do j = 0,jm
            lbfs(i,j)=(-ppazd(i,j,1))*jac(i,1) !+1./(amie*beta)*u(i,j,1)
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
            djdt(i,j,0) = ((-ppazd(i,j,1))*jac(i,1)-dum)/(2.*dz)/jac(i,0) &
                +djdx(i,j,0)+djdy(i,j,0) &
                +djedt(i,j,0) &
                +gn0e(i)*cn0e*gt0e(i)*delbx(i,j,0)*(gcpne(i)+gcpte(i)*0.)*bdgxcgy(i,0)/emass/ntube*ieqmo814 &
                +(dumx(i,j,0)+dumy(i,j,0))/bmag(i,0)*bdgxcgy(i,0)/emass/ntube*nonline*ieqmo814

            dum = weightm(i)*rbfr(i,jmi(i,j)) &
                +weightmn(i)*rbfr(i,jmn(i,j))
            djdt(i,j,1) = (dum-(-ppazd(i,j,0))*jac(i,0))/(2.*dz)/jac(i,1) &
                +djdx(i,j,1)+djdy(i,j,1) &
                +(djedt(i,j,1)) &
                +gn0e(i)*cn0e*gt0e(i)*delbx(i,j,1)*(gcpne(i)+gcpte(i)*0.)*bdgxcgy(i,1)/emass/ntube*ieqmo814 &
                +(dumx(i,j,1)+dumy(i,j,1))/bmag(i,1)*bdgxcgy(i,1)/emass/ntube*nonline*ieqmo814
         end do
      end do
      do i = 0,im-1
         do j = 0,jm-1
            do k = 0,1
               djdt(i,j,k) = djdt(i,j,k)*emass-gn0e(i)*cn0e*ezs(i,j,k)/ntube
            end do
         end do
      end do
      call enfxy(djdt(:,:,:))
      call enfz(djdt(:,:,:))
!      call filter(drhodt(:,:,:))
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine eqmo814(nstep,ip)   

      use gem_com
      use equil
      use fft_wrapper
      implicit none
      REAL(8) :: b2,gam0,gam1,th,delyz,bf,rkper,r,qr,shat,ter !,b
      REAL(8) :: kx,ky
      complex(8) :: lbfr(0:imx-1,0:jcnt-1)
      complex(8) :: lbfs(0:imx-1,0:jcnt-1)
      complex(8) :: rbfr(0:imx-1,0:jcnt-1)
      complex(8) :: rbfs(0:imx-1,0:jcnt-1)
      REAL(8),dimension(:),allocatable :: dely,aky
      complex(8),dimension(:,:,:,:),allocatable:: nab1,nab2
      complex(8),dimension(:,:,:,:),allocatable :: mx
      integer,dimension(:,:,:,:),allocatable :: ipiv
      REAL(8),dimension(:,:,:),allocatable :: form814
      REAL(8) :: sgnx,sgny,sz,myfe
      INTEGER :: i,i1,j,j1,k,k1,l,m,n,ifirst,nstep,ip,INFO
      INTEGER :: l1,m1,myk,myj,izonal,ix,ikx,iext
      COMPLEX(8) :: temp3dxy(0:imx-1,0:jmx-1,0:1),va(0:imx-1,0:jmx-1,0:1),&
            v(0:imx-1,0:jcnt-1,0:1)
      COMPLEX(8) :: sbuf(0:imx*jcnt*2-1),rbuf(0:imx*jmx*2-1)
      COMPLEX(8) :: sl(1:imx-1,0:jcnt-1,0:1)
      REAL(8) :: dum,u(0:imx,0:jmx,0:1)
      complex(8) :: ua(0:imx-1,0:jcnt-1,0:1),calph,cbeta
      REAL(8) :: grp,gthp,gxdgyp,dydrp,qhatp,grdgtp,bfldp,g2xp,g2yp, &
                 radiusp,g2zp,gzp,gxdgzp,gydgzp
      REAL(8) :: wx0,wx1,wz0,wz1
      character*70 fname
      character*4 holdmyid

      save form814,ifirst,dely,aky,mx,ipiv
      calph = 1.
      cbeta = 0.
      write(holdmyid,'(I4.4)') MyId
      fname='./matrix/'//'mx_814_'//holdmyid
      if (ifirst.ne.-99) then
         allocate(dely(0:imx-1),aky(0:jmx-1),nab1(0:imx-1,0:jcnt-1,0:imx-1,0:1))
         allocate(nab2(0:imx-1,0:jcnt-1,0:imx-1,0:1))
         allocate(mx(imx-1,imx-1,0:jcnt-1,0:1),form814(0:imx-1,0:jcnt-1,0:1))
         allocate(ipiv(imx-1,imx-1,0:jcnt-1,0:1))

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
                     ter = wx0*t0i(i1)+wx1*t0i(i1+1)
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
                     radiusp = wx0*wz0*radius(i1,k1)+wx0*wz1*radius(i1,k1+1) &
                     +wx1*wz0*radius(i1+1,k1)+wx1*wz1*radius(i1+1,k1+1)

                     m1 = mstart+int((float(m)+1.0)/2)
                     if(m==0)m1=0
                     sgny = isgnft(m)

                     ky=sgny*2.*pi*float(m1)/ly
                     kx=pi*float(l)/lx
                     bf=bfldp
!                     b=mims(1)*(kx*kx*grp**2 + &
!                        ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
!                        +2*dydrp*lr0/q0*qhatp*grdgtp) &
!                        +2*kx*ky*gxdgyp)/bf/bf
                  
                     nab1(l,m,i,n) = kx**2*grp**2+ky**2*  &
                        (dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                        +2*dydrp*lr0/q0*qhatp*grdgtp)

                     nab2(l,m,i,n) = -IU*ky*kx*2*gxdgyp
!   formfactor in ezamp
                     form814(l,m,n) = beta/jmx 
                     if(abs(ky)>kycut)form814(l,m,n) = 0.
                     if(m1.ne.mlk.and.onemd==1)form814(l,m,n) = 0.
 53               continue
 54            continue
 52         continue
 51      continue

         do k = 0,1
            do j = 0,jcnt-1
               do i = 1,imx-1
                  do ix = 1,imx-1
                     mx(i,ix,j,k) = 0.
                     if(i==ix)mx(i,ix,j,k) = gn0e(i)*cn0e
                     do ikx = 1,imx-1
                        mx(i,ix,j,k) = mx(i,ix,j,k) &
                             +nab1(ikx,j,i,k) &
                             *sin(ix*ikx*pi/imx)*sin(i*ikx*pi/imx)*2.0/(imx*beta*amie) &
                             +nab2(ikx,j,i,k) &
                              *sin(ix*ikx*pi/imx)*cos(i*ikx*pi/imx)*2.0/(imx*beta*amie)

                     end do
                  end do
               end do
            end do
         end do

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
      do i = 0,imx
         do j = 0,jmx
            do k = 0,1
               u(i,j,k) = djdt(i,j,k)+etaohm*djpa(i,j,k)
            end do
         end do
      end do
      call dcmpy(u(0:imx-1,0:jmx-1,0:1),v)

      temp3dxy = 0.
      do k = 0,1
         do j = 0,jcnt-1
            myj=jft(j)
            sl(1:imx-1,j,k) = v(1:imx-1,j,k)
            call ZGETRS('N',imx-1,1,mx(:,:,j,k),imx-1,ipiv(:,:,j,k), &
                         sl(:,j,k),imx-1,INFO) 
            temp3dxy(1:imx-1,myj,k) = sl(1:imx-1,j,k)
            temp3dxy(0,myj,k) = 0.
         end do
      end do

!  from rho(kx,ky) to phi(kx,ky)
      do k = 0,1
         do j = 0,jmx-1
            do i = 0,imx-1
               temp3dxy(i,j,k)=temp3dxy(i,j,k)/jmx
            end do
         end do
      end do

!  from phi(kx,ky) to phi(x,y)
      do k=0,mykm
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3dxy(i,j,k)
            end do
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3dxy(i,j,k) = tmpy(j)  ! phi(x,y)
            end do
         end do
      end do

      do i = 0,nxpp-1
         do j = 0,jm-1
            do k = 0,mykm
                  ez(i,j,k) = temp3dxy(i,j,k)+ezs(i,j,k)
            end do
         end do
      end do

!    x-y boundary points 
      call enfxy(ez(:,:,:))
      call enfz(ez(:,:,:))
!      call filter(upar(:,:,:))

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine gam(u,v,isdamp)
      use gem_com
      use equil
      use fft_wrapper
      implicit none
!                                                                                                                                                                                                                  
      complex(8) :: u(0:imx-1),v(0:imx-1)
      INTEGER :: n,i,j,k,k1,i1,isdamp

      call ccfft('z',0,kmx,1.d0,tmpz,coefz,workz,0)

      do i1 = 0,imx-1
         tmpz(0) = u(i1)

         if(GCLR.ne.master)then
            call MPI_SEND(tmpz(0),mykm,MPI_DOUBLE_COMPLEX,master, &
               gclr,tube_comm,ierr)
         end if

         if(gclr.eq.master)then
            do i = 1,GLST
               call MPI_RECV(tmpz(i*mykm),mykm,MPI_DOUBLE_COMPLEX,i, &
                   i,tube_comm,stat,ierr)
            end do
         end if

         if(GCLR.eq.master) then
            call ccfft('z',1,kmx,1.d0,tmpz,coefz,workz,0)
         end if

         call MPI_BCAST(tmpz,kmx,MPI_DOUBLE_COMPLEX,master, &
            tube_comm,ierr)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         do k = 0,kmx-1
            k1 = k
            if(k>kmx/2)k1=kmx-k
!            if(k1>1)tmpz(k) = 0.
            tmpz(k) = (1-nugam*dt*isdamp)*tmpz(k)
         end do
         call ccfft('z',-1,kmx,1.d0,tmpz,coefz,workz,0)
         tmpz = tmpz/kmx

         v(i1) = tmpz(gclr)
      end do

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                                                               
      subroutine gam1(u,v,isdamp)
      use gem_com
      use equil
      use fft_wrapper
      implicit none
!                                                                                                                                                                                                                  
      complex(8) :: u(0:imx-1),v(0:imx-1)
      INTEGER :: n,i,j,k,k1,i1,isdamp

      call ccfft('z',0,kmx,1.d0,tmpz,coefz,workz,0)

      do i1 = 0,imx-1
         tmpz(0) = u(i1)

         if(GCLR.ne.master)then
            call MPI_SEND(tmpz(0),mykm,MPI_DOUBLE_COMPLEX,master, &
               gclr,tube_comm,ierr)
         end if

         if(gclr.eq.master)then
            do i = 1,GLST
               call MPI_RECV(tmpz(i*mykm),mykm,MPI_DOUBLE_COMPLEX,i, &
                   i,tube_comm,stat,ierr)
            end do
         end if

         if(GCLR.eq.master) then
            call ccfft('z',1,kmx,1.d0,tmpz,coefz,workz,0)
         end if

         call MPI_BCAST(tmpz,kmx,MPI_DOUBLE_COMPLEX,master, &
            tube_comm,ierr)
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         do k = 0,kmx-1
            k1 = k
            if(k>kmx/2)k1=kmx-k
            if(k1>1)tmpz(k) = 0.
         end do
         call ccfft('z',-1,kmx,1.d0,tmpz,coefz,workz,0)
         tmpz = tmpz/kmx

         v(i1) = tmpz(gclr)
      end do

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                                                                               
      subroutine mdamp(u,n)

!     calculate mode histories. calculate modes for u at timestep n
!     and store in modehis(mode,n).

      use gem_com
      use equil
      use fft_wrapper
      implicit none
!     
      REAL(8) :: u(0:imx,0:jmx,0:1)
      INTEGER :: n,i,j,k,m
      COMPLEX(8) :: tmp3d(0:imx,0:jmx,0:1),cdum

      call ccfft('z',0,kmx,1.d0,tmpz,coefz,workz,0)

!      if(n.eq.0) return

      tmp3d = 0.
      do j=0,jm-1
         do i=0,im-1
            do k = 0,mykm
               tmp3d(i,j,k)=u(i,j,k)
            enddo
         end do
      end do

      do k = 0,mykm
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = tmp3d(i,j,k)
            end do
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               tmp3d(i,j,k) = tmpy(j)
            end do
         end do
      end do   
      tmp3d = tmp3d/jmx

      do m = 0,jcnt-1
         cdum = 0.
         do i = 1,imx-1
            cdum = cdum+abs(tmp3d(i,jft(m),0)**2)
         end do
         cdum = sqrt(cdum/(imx-1))
         mdhis(m) = cdum

         do i = 0,6
            phihis(i,m) = tmp3d(imx/8*(i+1),jft(m),0)
         end do

      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mdampa(u,n)

!     calculate mode histories. calculate modes for u at timestep n
!     and store in modehis(mode,n).

      use gem_com
      use equil
      use fft_wrapper
      implicit none
!     
      REAL(8) :: u(0:imx,0:jmx,0:1)
      INTEGER :: n,i,j,k,m
      COMPLEX(8) :: tmp3d(0:imx,0:jmx,0:1),cdum

      call ccfft('z',0,kmx,1.d0,tmpz,coefz,workz,0)

!      if(n.eq.0) return

      tmp3d = 0.
      do j=0,jm-1
         do i=0,im-1
            do k = 0,mykm
               tmp3d(i,j,k)=u(i,j,k)
            enddo
         end do
      end do

      do k = 0,mykm
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = tmp3d(i,j,k)
            end do
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               tmp3d(i,j,k) = tmpy(j)
            end do
         end do
      end do   
      tmp3d = tmp3d/jmx

      do m = 0,jcnt-1
         cdum = 0.
         do i = 1,imx-1
            cdum = cdum+abs(tmp3d(i,jft(m),0)**2)
         end do
         cdum = sqrt(cdum/(imx-1))
         mdhisa(m) = cdum

         do i = 0,6
            aparhis(i,m) = tmp3d(imx/8*(i+1),jft(m),0)
         end do

      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mdampd(u,v)

!     calculate mode histories. calculate modes for u at timestep n
!     and store in modehis(mode,n).

      use gem_com
      use equil
      use fft_wrapper
      implicit none
!     
      REAL(8) :: u(0:imx,0:jmx,0:1),v(0:100)
      INTEGER :: n,i,j,k,m
      COMPLEX(8) :: tmp3d(0:imx,0:jmx,0:1),cdum

      call ccfft('z',0,kmx,1.d0,tmpz,coefz,workz,0)

!      if(n.eq.0) return

      tmp3d = 0.
      do j=0,jm-1
         do i=0,im-1
            do k = 0,mykm
               tmp3d(i,j,k)=u(i,j,k)
            enddo
         end do
      end do

      do k = 0,mykm
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = tmp3d(i,j,k)
            end do
            call ccfft('y',1,jmx,1.d0,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               tmp3d(i,j,k) = tmpy(j)
            end do
         end do
      end do   
      tmp3d = tmp3d/jmx

      do m = 0,4
         cdum = 0.
         do i = 1,imx-1
            cdum = cdum+abs(tmp3d(i,m,0)**2)
         end do
         cdum = sqrt(cdum/(imx-1))
         v(m) = cdum
      end do

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine neq0(u)   
      use gem_com
      use fft_wrapper
      implicit none
      INTEGER :: i,j,k,k1,l,m,n,ifirst,nstep,ip,INFO,isbl,idigit,NNI1
      INTEGER :: l1,m1,myk,myj,izonal,ix,ikx,ism,isdamp
      real(8) :: u(0:imx,0:jmx,0:1)
      real :: kx,ky,kx0,th,shat,sgny,dum,dum1

      do k = 0,1
         do i = 0,imx
            dum = 0.
            do j = 0,jmx-1
               dum = dum+u(i,j,k)
            end do
            dum = dum/jmx
            u(i,:,k) = dum
         end do
      end do

      call enfxy(u(:,:,:))
      call enfz(u)
      return

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mpush(n)

      use gem_com
      use equil
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum,dum1
      INTEGER :: m,i,j,k,l,n
      INTEGER :: np_old,np_new
      REAL(8) :: rhog,vfac,kap,vpar,pidum,kaptp,kapnp,xnp
      REAL(8) :: b,th,r,enerb,cost,sint,qr,laps,sz,ter
      REAL(8) :: xt,xs,yt,xdot,ydot,zdot,xdt,ydt,pzdot,edot,pzd0,vp0
      REAL(8) :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
      REAL(8) :: grp,gxdgyp,rhox(4),rhoy(4),psp,pzp

      pidum = 1./(pi*2)**1.5*vwidth**3
      do m=1,mm(2)
         r=xc2(m)-0.5*lx+lr0
         k = int(zc2(m)/delz)
         wz0 = ((k+1)*delz-zc2(m))/delz
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
         ter = wx0*t0c(i)+wx1*t0c(i+1)        
         kaptp = wx0*captc(i)+wx1*captc(i+1)
         kapnp = wx0*capnc(i)+wx1*capnc(i+1)
         xnp = wx0*xn0c(i)+wx1*xn0c(i+1)

         b=1.-tor+tor*bfldp
         pzp = mims(2)*uc2(m)/b*fp/br0-q(2)*psp/br0

         rhog=sqrt(2.*b*muc(m)*mims(2))/(q(2)*b)*iflr

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
            xs=xc2(m)+rhox(l) !rwx(1,l)*rhog
            yt=yc2(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
!
!   particle can go out of bounds during gyroavg...
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            include "mpushli.h"
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
         vfac = 0.5*(mims(2)*uc2(m)**2 + 2.*muc(m)*b)
         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))
         kap = kapnp - (1.5-vfac/ter)*kaptp

         vpar = uc2(m)-q(2)/mims(2)*aparp*nonlinh*0.
         enerb=(muc(m)+mims(2)*vpar*vpar/b)/q(2)*tor
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
         xdot = vxdum*nonlinh -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         xdt = -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlinh &
             +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0
         ydt = iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0
         zdot =  vpar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         pzd0 = tor*(-muc(m)/mims(2)/radiusp/bfldp*psipp*dbdtp*grcgtp)
         pzdot = pzd0 + (q(2)/mims(2)*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
             +q(2)/mims(2)*(-xdot*delbyp+ydot*delbxp+zdot*dadzp))*ipara

         edot = q(2)*(xdt*exp1+ydt*eyp+vpar*ezp)

         xc3(m) = xc2(m) + 0.5*dt*xdot
         yc3(m) = yc2(m) + 0.5*dt*ydot
         zc3(m) = zc2(m) + 0.5*dt*zdot
         uc3(m) = uc2(m) + 0.5*dt*pzdot

         dum = 1-wc2(m)*nonlinh*0.
         if(isuni.eq.1)dum = exp(-vfac)*pidum
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
!         vxdum = eyp+vpar/b*delbxp
         wc3(m)=wc2(m) + 0.5*dt*(vxdum*kap + edot/ter)*dum*xnp*(tloadc/ter)**1.5*exp(vfac*(1/tloadc-1./ter))-0.5*dt*gamion*wc2(m)
         

 333     continue
         laps=anint((zc3(m)/lz)-.5)*(1-peritr)
         r=xc3(m)-0.5*lx+lr0
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         i = max(i,0)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         qr = wx0*sf(i)+wx1*sf(i+1)
         yc3(m)=dmod(yc3(m)-laps*2*pi*qr*lr0/q0*dsign(1.d0,q0)+800.*ly,ly)
         if(xc3(m)>lx)then
            xc3(m) = lx-1.e-8
            zc3(m)=lz-zc3(m)
            xc2(m) = xc3(m)
            zc2(m) = zc3(m)
            wc2(m) = 0.
            wc3(m) = 0.
         end if
         if(xc3(m)<0.)then
            xc3(m) = 1.e-8
            zc3(m)=lz-zc3(m)
            xc2(m) = xc3(m)
            zc2(m) = zc3(m)
            wc2(m) = 0.
            wc3(m) = 0.
         end if
         zc3(m)=dmod(zc3(m)+8.*lz,lz)
         xc3(m)=dmod(xc3(m)+8.*lx,lx)         
         xc3(m) = min(xc3(m),lx-1.0e-8)
         yc3(m) = min(yc3(m),ly-1.0e-8)
         zc3(m) = min(zc3(m),lz-1.0e-2)
         
      enddo

      np_old=mm(2)
      call init_pmove(zc3,np_old,lz,ierr)

      call pmove(xc2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xc3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yc2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yc3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zc2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zc3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(uc2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(uc3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wc2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wc3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(muc,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit


      call end_pmove(ierr)
      mm(2)=np_new

      return
      end

!-----------------------------------------------------------------------

      subroutine mcush(n)

      use gem_com
      use equil
      implicit none
      INTEGER :: n
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,w3old
      INTEGER :: m,i,j,k,ns,l
      INTEGER :: np_old,np_new
      REAL(8) :: vfac,vpar,vxdum,dum,xdot,ydot,xdt,ydt,zdot,pzdot,edot,pzd0,vp0
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

      do m=1,mm(2)
         r=xc3(m)-0.5*lx+lr0

         k = int(zc3(m)/delz)
         wz0 = ((k+1)*delz-zc3(m))/delz
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
         ter = wx0*t0c(i)+wx1*t0c(i+1)        
         kaptp = wx0*captc(i)+wx1*captc(i+1)
         kapnp = wx0*capnc(i)+wx1*capnc(i+1)
         xnp = wx0*xn0c(i)+wx1*xn0c(i+1)
         b=1.-tor+tor*bfldp
         pzp = mims(2)*uc3(m)/b*fp/br0-q(2)*psp/br0

         rhog=sqrt(2.*b*muc(m)*mims(2))/(q(2)*b)*iflr

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
            xs=xc3(m)+rhox(l) !rwx(1,l)*rhog
            yt=yc3(m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
!   BOUNDARY
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            xt = min(xt,lx-1.0e-8)
            yt = min(yt,ly-1.0e-8)

            include "mcushli.h"
 200     continue

         exp1=exp1/4.
         eyp=eyp/4.
         ezp=ezp/4.
         delbxp=delbxp/4.
         delbyp=delbyp/4.
         dpdzp = dpdzp/4.
         dadzp = dadzp/4.
         aparp = aparp/4.


         vfac = 0.5*(mims(2)*uc3(m)**2 + 2.*muc(m)*b)
         vp0 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vp0 = vp0*(vt0+vpp*(r-lr0)+yd0*sin(kxz*(r+0.5*lx-lr0)*2*pi/lx))
         kap = kapnp -(1.5-vfac/ter)*kaptp

         vpar = uc3(m)-q(2)/mims(2)*aparp*nonlinh*0.
         enerb=(muc(m)+mims(2)*vpar*vpar/b)/q(2)*tor
         dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
         xdot = vxdum*nonlinh -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         xdt = -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
         ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlinh     &
             +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
                 (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0
         ydt = iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp*(-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0
         zdot =  vpar*(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
                 +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp

         pzd0 = tor*(-muc(m)/mims(2)/radiusp/bfldp*psipp*dbdtp*grcgtp)
         pzdot = pzd0 + (q(2)/mims(2)*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
             +q(2)/mims(2)*(-xdot*delbyp+ydot*delbxp+zdot*dadzp))*ipara

         edot = q(2)*(xdt*exp1+ydt*eyp+vpar*ezp)

         xc3(m) = xc2(m) + dt*xdot
         yc3(m) = yc2(m) + dt*ydot
         zc3(m) = zc2(m) + dt*zdot
         uc3(m) = uc2(m) + dt*pzdot

         dum = 1-wc3(m)*nonlinh*0.
         if(isuni.eq.1)dum = exp(-vfac)*pidum
!         vxdum = eyp+vpar/b*delbxp
         vxdum = (eyp/b+vpar/b*delbxp)*dum1
         w3old = wc3(m)
         wc3(m) = wc2(m) + dt*(vxdum*kap+edot/ter)*dum*xnp*(tloadc/ter)**1.5*exp(vfac*(1/tloadc-1./ter))-dt*gamion*wc2(m)

         if(abs(wc3(m)).gt.1..and.nonlinh==1)then
            wc3(m) = 0.
            wc2(m) = 0.
         end if


         laps=anint((zc3(m)/lz)-.5)*(1-peritr)
         r=xc3(m)-0.5*lx+lr0
         i = int((r-rin)/dr)
         i = min(i,nr-1)
         wx0 = (rin+(i+1)*dr-r)/dr
         wx1 = 1.-wx0
         qr = wx0*sf(i)+wx1*sf(i+1)
         yc3(m)=dmod(yc3(m)-laps*2*pi*qr*lr0/q0*dsign(1.d0,q0)+800.*ly,ly)
         if(xc3(m)>lx)then
            xc3(m) = lx-1.e-8
            zc3(m)=lz-zc3(m)
            xc2(m) = xc3(m)
            zc2(m) = zc3(m)
            wc2(m) = 0.
            wc3(m) = 0.
         end if
         if(xc3(m)<0.)then
            xc3(m) = 1.e-8
            zc3(m)=lz-zc3(m)
            xc2(m) = xc3(m)
            zc2(m) = zc3(m)
            wc2(m) = 0.
            wc3(m) = 0.
         end if
         zc3(m)=dmod(zc3(m)+8.*lz,lz)
         xc3(m)=dmod(xc3(m)+8.*lx,lx)         
         xc3(m) = min(xc3(m),lx-1.0e-8)
         yc3(m) = min(yc3(m),ly-1.0e-8)
         zc3(m) = min(zc3(m),lz-1.0e-2)

!     particle diagnostics done here because info is available...
         mypfl=mypfl + w3old*(eyp+vpar*delbxp/b) 
         myefl=myefl + vfac*w3old*(eyp+vpar*delbxp/b)
         myke=myke + vfac*wc3(m)
         mynos=mynos + wc3(m)
         myavewi = myavewi+abs(wc3(m))

!     xn+1 becomes xn...
         uc2(m)=uc3(m)
         xc2(m)=xc3(m)
         yc2(m)=yc3(m)
         zc2(m)=zc3(m)
         wc2(m)=wc3(m)

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
      avewc(n) = rbuf(5)/( float(tmm(1)) )
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      nos(2,n)=nostemp/( float(tmm(1)) )
      pfl(2,n)=pfltemp/( float(tmm(1)) )
      efl(2,n)=efltemp/( float(tmm(1)) )
      ke(2,n)=ketemp/(float(tmm(1)))
      np_old=mm(2) 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call init_pmove(zc3,np_old,lz,ierr)
!     
      call pmove(xc2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(xc3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yc2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(yc3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zc2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(zc3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(uc2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(uc3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wc2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(wc3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(muc,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call end_pmove(ierr)
      mm(2)=np_new
!     write(*,*)MyId,mm(1)

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
