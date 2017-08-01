	 subroutine outd(n)
	 use gem_com
	 implicit none
	 INTEGER :: n,np
	 integer :: i,j,k

!02/17/2009 debugging when removing "contu,pskip,wmax,imovie' from gem.in
	 imovie = 1000000000
!	 if ( (mod(n,imovie).eq.0) ) then
!	   call mphxy(phi(:,:,:),'mphxy',25)
!	   call mphxz(phi(:,:,:),'mphxz',26)
!	 endif

	 if (mod(n,nplot).eq.0) then
!	   call phixy(rho,'rhoxy',27,n)
!	   call phixz(rho,'rhoxz',28,n)

	   call aphir(denh(:,:,:),'denhr',30,n)
	   call aphir(phi(:,:,:),'aphir',52,n)
	   call aphir(apar(:,:,:),'aparr',53,n)

	   call phixy(phi(:,:,:),'phixy',31,n)
	   call phixz(phi(:,:,:),'phixz',32,n)

	   call aphir(hden0(:,:,:),'hden0r',35,n)
	   call phixy(apar(:,:,:),'apaxy',36,n)
	   call phixz(apar(:,:,:),'apaxz',37,n)

!	   call aphir(upar(:,:,:),'aupar',38,n)
!	   call phixy(upar(:,:,:),'upaxy',39,n)
!	   call phixz(upar(:,:,:),'upaxz',40,n)

!	   call phixy(dti(:,:,:),'tpixy',29,n)
	   call dump3d(phi(:,:,:),'phi3d',50,n)
	   call dump3d(apar(:,:,:),'apa3d',51,n)

	 endif

	 if (n.eq.nm-1) call histout(34)

	 return
	 end

!--------------------------------------------------------------

	 subroutine phixy(grd,fl,unt,n)

	 use gem_com
	 implicit none
	 REAL(8) :: grd(0:imx,0:jmx,0:1)
	 integer :: i,j,k,n
	 INTEGER :: unt,tind,flag,oproc,outk
	 character*5 fl
	 character*70 flnm
!	 save flag

	 oproc=int(km/2/kcnt*ntube)
	 if (MyId.eq.oproc) then

	 flnm=outdir//fl//'.out'
	 open(unt,file=flnm,form='formatted', &
     	   status='unknown',position='append')

	 if (n.eq.0) then
!      if (flag.ne.-1) then
	   tind=int(nm/nplot)+1
	   if(iget.eq.0)write(unt,110)im,jm,tind,nplot,int(dt)
	 endif
!       flag=-1
	 outk=0  !km/2-MyId*kcnt
	 write(unt,99)n
 99	 format('time step= ',I6)
	 do j=0,jm
	   do i=0,im
	     write(unt,100)grd(i,j,outk)
	   enddo
	 enddo


 100	 format (e10.3)
 110	 format (5I5)
	 endfile unt
	 close(unt)
	 endif

	 return
	 end

!--------------------------------------------------------------

	 subroutine phixz(grd,fl,unt,n)

	 use gem_com
	 implicit none
	 integer :: i,j,k,n
	 REAL(8) :: grd(0:imx,0:jmx,0:1)
	 REAL(8) :: sbuf(0:(imx+1)*kmx)
	 REAL(8) :: rbuf(0:(imx+1)*kmx)
	 INTEGER :: unt,tind,flag
	 INTEGER :: lst,m,pkm
	 character*5 fl
	 character*70 flnm
!      save flag

	 flnm=outdir//fl//'.out'

!      if ( (MyId.eq.Last).and.(flag.ne.-1) ) then
	 if ( (MyId.eq.Last).and.(n.eq.0) ) then
	   open(unt,file=flnm,form='formatted', &
     	     status='unknown',position='append')
	   tind=int(nm/nplot)+1
	   if(iget.eq.0)write(unt,110)im,km,tind,nplot,int(dt)
	   endfile unt
	   close(unt)
	 endif
!      flag=-1

	 pkm=int(ntube*km/numprocs)
	 do k=0,pkm-1
	   do i=0,im
	     m=k*(im+1)+i
	     sbuf(m)=grd(i,jm/2,k)
	   enddo
	 enddo
	 cnt=pkm*(im+1)
	 lst=cnt-1
	 call MPI_GATHER(sbuf(0:lst),cnt, &
               MPI_REAL8, &
               rbuf,cnt,MPI_REAL8, &
               glst,tube_comm,ierr)
! 

	 if (MyId.eq.Last) then
	    open(unt,file=flnm,form='formatted', &
                status='unknown',position='append')
	    
	    write(unt,99)n
 99	    format('time step= ',I6)

	    do i=0,((im+1)*pkm)*numprocs/ntube-1
	       write(unt,100)rbuf(i)
	    enddo
	    do k=pkm,mykm
	       do i=0,im
		  write(unt,100)grd(i,jm/2,k)
	       enddo
	    enddo
	    endfile unt
	    close(unt)
	 endif

 100	 format (e10.3)
 110	 format (5I5)

	 return
	 end

!--------------------------------------------------------------

	 subroutine mphxy(grd,fl,unt)

	 use gem_com
	 implicit none
	 integer :: i,j,k
	 REAL(8) :: grd(0:imx,0:jmx,0:1)
	 INTEGER :: unt,tind,flag,oproc,outk
	 character*5 fl
	 character*70 flnm
	 save flag

	 oproc=int(km/2/kcnt)
	 if (MyId.eq.oproc) then

	 flnm=outdir//fl//'.out'
	 open(unt,file=flnm,form='formatted', &
     	   status='unknown',position='append')
	 if (flag.ne.-1) then
	   tind=int(nm/imovie)+1
	   if(iget.eq.0)write(unt,110)im,jm,tind,imovie,int(dt)
	 endif
	 flag=-1
	 outk=km/2-MyId*kcnt
	 do j=0,jm
	   do i=0,im
	     write(unt,100)grd(i,j,outk)
	   enddo
	 enddo


 100	 format (e9.2)
 110	 format (5I5)
	 endfile unt
	 close(unt)
	 endif

	 return
	 end

!--------------------------------------------------------------

	 subroutine mphxz(grd,fl,unt)

	 use gem_com
	 implicit none
	 integer :: i,j,k,n
	 REAL(8) :: grd(0:imx,0:jmx,0:1)
	 REAL(8) :: sbuf(0:(imx+1)*kmx)
	 REAL(8) :: rbuf(0:(imx+1)*kmx)
	 INTEGER :: unt,tind,flag
	 INTEGER :: lst,m,pkm
	 character*5 fl
	 character*70 flnm
	 save flag

	 flnm=outdir//fl//'.out'

	 if ( (MyId.eq.Last).and.(flag.ne.-1) ) then
	   open(unt,file=flnm,form='formatted', &
     	     status='unknown',position='append')
	   tind=int(nm/imovie)+1
	   if(iget.eq.0)write(unt,110)im,km,tind,imovie,int(dt)
	   endfile unt
	   close(unt)
	 endif
	 flag=-1

	 pkm=int(km/numprocs)
	 do k=0,pkm-1
	   do i=0,im
	     m=k*(im+1)+i
	     sbuf(m)=grd(i,jm/2,k)
	   enddo
	 enddo
	 cnt=pkm*(im+1)
	 lst=cnt-1
	 call MPI_GATHER(sbuf(0:lst),cnt, &
               MPI_REAL8, &
               rbuf,cnt,MPI_REAL8, &
               Last,MPI_COMM_WORLD,ierr)
 

	 if (MyId.eq.Last) then
         open(unt,file=flnm,form='formatted', &
          status='unknown',position='append')
	 do i=0,((im+1)*pkm)*numprocs-1
	   write(unt,100)rbuf(i)
         enddo
	 do k=pkm,mykm
	   do i=0,im
	     write(unt,100)grd(i,jm/2,k)
	   enddo
	 enddo
	 endfile unt
	 close(unt)
	 endif

 100	 format (e9.2)
 110	 format (5I5)

	 return
	 end

!--------------------------------------------------------------

	 subroutine histout(unt)

	 use gem_com
	 implicit none
	 integer :: i,j,k
	 INTEGER :: unt,n,m
	 character*70 flnm
      
      
	 if (MyId.eq.Master) then

	 flnm=outdir//'hist.out'
	 open(unt,file=flnm,form='formatted', &
     	   status='unknown',position='rewind')
	 write(unt,80)nsm,nm,modem,int(dt)
	 do m=1,modem
	   write(unt,90)lmode(m),mmode(m),nmode(m)
	 enddo
	 do m=1,modem
	   do n=1,nm
	     write(unt,100)pmodehis(m,n)
	   enddo
	 enddo
	 do m=1,nsm
	   do n=1,nm
	     write(unt,110)ke(m,n)
	   enddo
	 enddo
	 do m=1,nsm
	   do n=1,nm
	     write(unt,110)pfl(m,n)
	   enddo
	 enddo
	 do m=1,nsm
	   do n=1,nm
	     write(unt,110)efl(m,n)
	   enddo
	 enddo
	 do n=1,nm
	   write(unt,110)fe(n)
	 enddo
	 do n=1,nm
	   write(unt,110)rmsphi(n)
	 enddo

 80	 format (4I6)
 90	 format (3I6)
 100	 format (2e10.3)
 110	 format (1e10.3)

	 endfile unt
	 close(unt)
	 endif

	 return
	 end

!--------------------------------------------------------------

	 subroutine aphir(grd,fl,unt,n)

	 use gem_com
	 use fft_wrapper
	 implicit none
	 integer :: i,j,k,n
	 REAL(8) :: grd(0:imx,0:jmx,0:1)
	 REAL(8) :: aph(0:imx)
	 REAL(8) :: myaph(0:imx)
	 REAL(8) :: myjac(0:imx),v(0:imx),dum,dum1
	 INTEGER :: unt,tind,flag
	 character*5 fl
	 character*70 flnm
!      save flag

	 flnm=outdir//fl//'.out'

!      if ( (MyId.eq.Master).and.(flag.ne.-1) ) then
	 if ( (MyId.eq.Master).and.(n.eq.0) ) then
	   open(unt,file=flnm,form='formatted', &
     	     status='unknown',position='append')
	   tind=int(nm/nplot)+1
	   if(iget.eq.0)write(unt,110)im,tind,nplot,int(dt)
	 endif
!      flag=-1

	 do i=0,im
	   myaph(i)=0.
	   aph(i)=0.
	   myjac(i) = 0.
	   v(i) = 0.
	 enddo
	 do k=0,mykm-1
	   do j=0,jm-1
	     do i=0,im
	       myaph(i)=myaph(i)+grd(i,j,k)
	       myjac(i) = myjac(i)+jac(i,k)
	     enddo
	   enddo
	 enddo

	 cnt=imx+1
	 call MPI_REDUCE(myaph(0:imx),aph(0:imx),cnt, &
     			MPI_REAL8, &
     			MPI_SUM,Master, &
     			MPI_COMM_WORLD,ierr)

	 call MPI_REDUCE(myjac(0:imx),v(0:imx),cnt, &
     			MPI_REAL8, &
     			MPI_SUM,Master, &
     			MPI_COMM_WORLD,ierr)
!      
	 if (MyId.eq.Master) then
	 open(unt,file=flnm,form='formatted', &
     	   status='unknown',position='append')

	 write(unt,99)n
 99	 format('time step= ',I6)

	 do i=0,im
	   aph(i)=aph(i)/v(i)
	   write(unt,100)i,aph(i)
	 enddo

	 if(izon.eq.5)then
	    do i=0,im-1
	       tmpx(i) = aph(i) 
	    enddo
	    call ccfft('x',1,imx,1.d0,tmpx,tmpx,coefx,workx,0)
	    write(*,98)n,aimag(tmpx(0)),aimag(tmpx(1)),aimag(tmpx(2)), &
	    aimag(tmpx(3)),aimag(tmpx(4))

!        write(*,98)n,real(tmpx(0)),real(tmpx(1)),real(tmpx(2)),
!     1              real(tmpx(3)),real(tmpx(4))

	    write(15,98)n,aimag(tmpx(0)),aimag(tmpx(1)),aimag(tmpx(2)),&
	    aimag(tmpx(3)),aimag(tmpx(4))
 98	    format(1x,i5, 5(2x, e12.3))
	 end if
	 endfile unt
	 close(unt)
	 endif

 100	 format (i5, 2x, e10.3)
 110	 format (4I5)

	 return
	 end

!ccccccccccccccccccccccccc

	 subroutine dump3d(grd,fl,unt,n)

	 use gem_com
	 implicit none
	 integer :: i,j,k,n
	 REAL(8) :: grd(0:imx,0:jmx,0:1)
	 REAL(8) :: sbuf(0:(imx+1)*(jmx+1))
	 REAL(8) :: rbuf(0:(imx+1)*(jmx+1)*kmx)
	 INTEGER :: unt,tind,flag
	 INTEGER :: lst,m,pkm
	 character*5 fl
	 character*70 flnm
!      save flag

	 flnm=outdir//fl//'.out'

!      if ( (MyId.eq.Last).and.(flag.ne.-1) ) then
	 if ( (MyId.eq.Last).and.(n.eq.0) ) then
	   open(unt,file=flnm,form='formatted', &
     	     status='unknown',position='append')
	   tind=int(nm/nplot)+1
!	   if(iget.eq.0)write(unt,110)im,jm,km
	   endfile unt
	   close(unt)
	 endif
!      flag=-1

	 pkm=int(ntube*km/numprocs)
	 do k=0,pkm-1
	    do j = 0,jm
	       do i=0,im
		  m=k*(im+1)*(jm+1)+j*(im+1)+i
		  sbuf(m)=grd(i,j,k)
	       end do
	   enddo
	 enddo
	 cnt=pkm*(im+1)*(jm+1)
	 lst=cnt-1
	 call MPI_GATHER(sbuf(0:lst),cnt, &
               MPI_REAL8, &
               rbuf,cnt,MPI_REAL8, &
               glst,tube_comm,ierr)
! 

	 if (MyId.eq.Last) then
	    open(unt,file=flnm,form='formatted', &
                status='unknown',position='append')
	    
!         write(unt,99)n
! 99	    format('time step= ',I6)

	    do i=0,((im+1)*(jm+1)*pkm)*numprocs/ntube-1
	       write(unt,100)rbuf(i)
	    enddo
	    do k=pkm,mykm
	       do j = 0,jm
		  do i=0,im
		     write(unt,100)grd(i,j,k)
		  end do
	       enddo
	    enddo
	    endfile unt
	    close(unt)
	 endif

 100	 format (e10.3)
 110	 format (3I5)

	 return
	 end

!--------------------------------------------------------------
	 subroutine pol2d

	 use gem_com
	 use equil
	 implicit none
	 integer :: i,j,i1,j1,k,n,ir,ith,n1,n2, numx,numz,&
                    unt=71,unt1=72,unt2=73,upol=74,upolxz=75
	 parameter (n1=5000,n2=10000,numx=400,numz=600)
         integer :: m,mpol(1:25),npts(0:numx,0:numz)
	 REAL(8) :: grd(0:imx,0:jmx,0:kmx),yxz(0:numx,0:numz)
	 real(8) :: r,th,thf,qhatp,qr,dum,xt,yt,zt,dum1,dum2,wx0,wx1,wy0,wy1,wz0,wz1
	 real(8) :: xmin,xmax,zmin,zmax,dumx,dumz
	 real(8) :: hghtp,radiusp,wr0,wr1,rplt(0:n1),qplt(0:n1)
	 complex(8) :: grdm(1:25,0:n1)
	 character*5 fl
	 character*70 flnm,polout,bdout,mpolout,polxzout

	 do i = 1,25
	    mpol(i) = i
	 end do

	 if(isphi==0)flnm=outdir//'apa3d'//'.out'
	 if(isphi==1)flnm=outdir//'phi3d'//'.out'
	 if(isphi==2)flnm=outdir//'dteti'//'.out'
	 polout=outdir//'pol1d'//'.out'
	 polxzout=outdir//'polxz'//'.out'
	 bdout=outdir//'bd'//'.out'
	 mpolout=outdir//'mpol'//'.out'

	 if (MyId.eq.Last) then
	   open(unt,file=flnm,form='formatted', &
     	     status='old',action='read')
!	   read(unt,110)im,jm,km
! if iphi=2 then read from phiti, one page only
	   i1 = 1
	   if(isphi==2)i1=ipg
	   do n=i1,ipg
	      do k = 0,kmx
		 do j = 0,jmx
		    do i = 0,imx
		       read(unt,100)grd(i,j,k)
		    end do
		 end do
	      end do
	   end do
!	   endfile unt
	   close(unt)

	   open(unt2,file=bdout,form='formatted', &
     	     status='unknown',position='rewind')

	   xmin = 100000.
	   xmax = -100000.
	   zmin = 100000.
	   zmax = -100000.
	   do ir = 0,1
	      r = rin+(rout-rin)*ir
	      i1 = int((r-rin)/dr)
	      i1 = min(i1,nr-1)
	      wx0 = (rin+(i1+1)*dr-r)/dr
	      wx1 = 1.-wx0
	      do ith = 0,n2
		 th = -pi+pi2/float(n2)*ith
		 j1 = int((th+pi)/dth)
		 j1 = min(j1,ntheta-1)
		 wz0 = (-pi+(j1+1)*dth-th)/dth
		 wz1 = 1.-wz0

		 radiusp = wx0*wz0*radius(i1,j1)+wx0*wz1*radius(i1,j1+1) &
                     +wx1*wz0*radius(i1+1,j1)+wx1*wz1*radius(i1+1,j1+1) 
		 hghtp = wx0*wz0*hght(i1,j1)+wx0*wz1*hght(i1,j1+1) &
                     +wx1*wz0*hght(i1+1,j1)+wx1*wz1*hght(i1+1,j1+1) 
		     
		 write(unt2,130)radiusp,hghtp
		 if(radiusp>=xmax)xmax = radiusp
		 if(radiusp<=xmin)xmin = radiusp
		 if(hghtp>=zmax)zmax = hghtp
		 if(hghtp<=zmin)zmin = hghtp
	      end do
	   end do
	   endfile unt2
	   close(unt2)

	   dumx = (xmax-xmin)/numx
	   dumz = (zmax-zmin)/numz
	   npts = 0
	   yxz = 0.
	   open(unt1,file=polout,form='formatted', &
     	     status='unknown',position='rewind')

	   do ir = 0,n1
	      r = rin+(rout-rin)/float(n1)*ir
	      i1 = int((r-rin)/dr)
	      i1 = min(i1,nr-1)
	      wr0 = (rin+(i1+1)*dr-r)/dr
	      wr1 = 1.-wr0
	      rplt(ir) = rhogem(i1)*wr0+rhogem(i1+1)*wr1  !r/a
	      qplt(ir) = sf(i1)*wr0+sf(i1+1)*wr1
	      qr = wr0*sf(i1)+wr1*sf(i1+1)
	      grdm(:,ir) = 0.
	      do ith = 0,n2
		 th = -pi+pi2/float(n2)*ith
		 j1 = int((th+pi)/dth)
		 j1 = min(j1,ntheta-1)
		 wz0 = (-pi+(j1+1)*dth-th)/dth
		 wz1 = 1.-wz0

		 radiusp = wr0*wz0*radius(i1,j1)+wr0*wz1*radius(i1,j1+1) &
                     +wr1*wz0*radius(i1+1,j1)+wr1*wz1*radius(i1+1,j1+1) 
		 hghtp = wr0*wz0*hght(i1,j1)+wr0*wz1*hght(i1,j1+1) &
                     +wr1*wz0*hght(i1+1,j1)+wr1*wz1*hght(i1+1,j1+1) 
		 thf = wr0*wz0*thflx(i1,j1)+wr0*wz1*thflx(i1,j1+1) &
                     +wr1*wz0*thflx(i1+1,j1)+wr1*wz1*thflx(i1+1,j1+1) 
		 qhatp = wr0*wz0*qhat(i1,j1)+wr0*wz1*qhat(i1,j1+1) &
                     +wr1*wz0*qhat(i1+1,j1)+wr1*wz1*qhat(i1+1,j1+1) 
		     

		 yt = wr0*wz0*yfn(i1,j1)+wr0*wz1*yfn(i1,j1+1) &
                     +wr1*wz0*yfn(i1+1,j1)+wr1*wz1*yfn(i1+1,j1+1) 
		 yt = dmod(yt+800.*ly,ly)    
		 zt = wz0*zfnth(j1)+wz1*zfnth(j1+1)
		 zt = min(zt,lz-1.e-2)
		 xt = r-rin

		 i=int(xt/dx)
		 i = min(i,imx-1)
		 j=int(yt/dy)
		 j = min(j,jmx-1)
		 k=int(zt/dz)
		 k = min(k,kmx-1)

		 wx0=float(i+1)-xt/dx 
		 wx1=1.-wx0
		 wy0=float(j+1)-yt/dy
		 wy1=1.-wy0
		 wz0=float(k+1)-zt/dz
		 wz1=1.-wz0
		 dum = wx0*wy0*wz0*grd(i,j,k)  &
		 + wx1*wy0*wz0*grd(i+1,j,k) &
		 + wx0*wy1*wz0*grd(i,j+1,k) &
		 + wx1*wy1*wz0*grd(i+1,j+1,k) &
		 + wx0*wy0*wz1*grd(i,j,k+1) &
		 + wx1*wy0*wz1*grd(i+1,j,k+1) &
		 + wx0*wy1*wz1*grd(i,j+1,k+1) &
		 + wx1*wy1*wz1*grd(i+1,j+1,k+1)
		 do m = 1,25
		    grdm(m,ir) = grdm(m,ir)+dum*exp(IU*mpol(m)*thf)*qhatp/qr*pi2/n2
		 end do
		 if(ith==n2/2)write(unt1,120)r/a,dum

		 i=int((radiusp-xmin)/dumx)
		 i = min(i,numx-1)
		 k = int((hghtp-zmin)/dumz)
		 k = min(k,numz-1)
		 npts(i,k) = npts(i,k)+1
		 yxz(i,k) = yxz(i,k)+dum
                 npts(i+1,k) = npts(i+1,k)+1
                 yxz(i+1,k) = yxz(i+1,k)+dum
                 npts(i,k+1) = npts(i,k+1)+1
                 yxz(i,k+1) = yxz(i,k+1)+dum
                 npts(i+1,k+1) = npts(i+1,k+1)+1
                 yxz(i+1,k+1) = yxz(i+1,k+1)+dum

	      end do
	   end do
	   endfile unt1
	   close(unt1)


           open(upolxz,file=polxzout,form='formatted', &
             status='unknown',position='rewind')
           yxz = yxz/(npts+0.01)
	   do i = 0,numx
	      do k = 0,numz
		 write(upolxz,140)npts(i,k),xmin+i*dumx,zmin+k*dumz,yxz(i,k)
	      end do
	   end do
	   close(upolxz)

	   open(upol,file=mpolout,form='formatted', &
     	     status='unknown',position='rewind')
	   do m = 1,25
	      do ir = 0,n1
		 write(upol,140)mpol(m),rplt(ir),abs(grdm(m,ir)),real(grdm(m,ir)),aimag(grdm(m,ir))
	      end do
	   end do
	   endfile upol
	   close(upol)

	 endif

 100	 format (e12.5)
 110	 format (3I5)
 120	 format (3(2x,e12.5))
 130	 format (2(2x,e12.5))
 140	 format (1x,I5,4(2x,e12.5))
	 return
	 end

!--------------------------------------------------------------
	 subroutine tor3d

	 use gem_com
	 use equil
	 implicit none
	 integer :: i,j,i1,j1,k,n,ir,ith,izeta,nzeta,n1,n2, numx,numz,&
                    unt=81,unt1=82,unt2=83,upol=84,upolxz=85,uprof=86
	 parameter (n1=128,n2=512,numx=64,numz=64,nzeta=120)
         integer :: m,mpol(1:25),npts(0:numx,0:numz)
	 REAL(8) :: grd(0:imx,0:jmx,0:kmx),yxz(0:numx,0:numz)
	 real(8) :: r,th,thf,qhatp,qr,dum,xt,yt,zt,dum1,dum2,wx0,wx1,wy0,wy1,wz0,wz1,rmajr,elonr,triar
	 real(8) :: xmin,xmax,zmin,zmax,dumx,dumz,rtor,thtor
	 real(8) :: hghtp,radiusp,wr0,wr1,rplt(0:n1),qplt(0:n1)
	 complex(8) :: grdm(1:25,0:n1)
	 character*5 fl
	 character*70 flnm,polout,bdout,mpolout,polxzout,profiles

	 do i = 1,25
	    mpol(i) = 5+i
	 end do

	 if(isphi==0)flnm=outdir//'apa3d'//'.out'
	 if(isphi==1)flnm=outdir//'phi3d'//'.out'
	 if(isphi==2)flnm=outdir//'dteti'//'.out'
	 polout=outdir//'tor3d'//'.out'
	 bdout=outdir//'bd'//'.out'
	 profiles=outdir//'profiles'//'.out'

	 if (MyId.eq.Last) then
	   open(unt,file=flnm,form='formatted', &
     	     status='old',action='read')
!	   read(unt,110)im,jm,km
! if iphi=2 then read from phiti, one page only
	   i1 = 1
	   if(isphi==2)i1=ipg
	   do n=i1,ipg
	      do k = 0,kmx
		 do j = 0,jmx
		    do i = 0,imx
		       read(unt,100)grd(i,j,k)
		    end do
		 end do
	      end do
	   end do
 100	   format (e12.5)
!	   endfile unt
	   close(unt)

	   open(unt2,file=bdout,form='formatted', &
     	     status='unknown',position='rewind')

	   xmin = 100000.
	   xmax = -100000.
	   zmin = 100000.
	   zmax = -100000.
	   do ir = 0,1
	      r = rin+(rout-rin)*ir
	      i1 = int((r-rin)/dr)
	      i1 = min(i1,nr-1)
	      wx0 = (rin+(i1+1)*dr-r)/dr
	      wx1 = 1.-wx0
	      do ith = 0,n2
		 th = -pi+pi2/float(n2)*ith
		 j1 = int((th+pi)/dth)
		 j1 = min(j1,ntheta-1)
		 wz0 = (-pi+(j1+1)*dth-th)/dth
		 wz1 = 1.-wz0

		 radiusp = wx0*wz0*radius(i1,j1)+wx0*wz1*radius(i1,j1+1) &
                     +wx1*wz0*radius(i1+1,j1)+wx1*wz1*radius(i1+1,j1+1) 
		 hghtp = wx0*wz0*hght(i1,j1)+wx0*wz1*hght(i1,j1+1) &
                     +wx1*wz0*hght(i1+1,j1)+wx1*wz1*hght(i1+1,j1+1) 
		     
		 write(unt2,130)radiusp,hghtp
 130		 format (2(2x,e12.5))
		 if(radiusp>=xmax)xmax = radiusp
		 if(radiusp<=xmin)xmin = radiusp
		 if(hghtp>=zmax)zmax = hghtp
		 if(hghtp<=zmin)zmin = hghtp
	      end do
	   end do
	   endfile unt2
	   close(unt2)

	   dumx = (xmax-xmin)/numx
	   dumz = (zmax-zmin)/numz
	   npts = 0
	   yxz = 0.
	   open(unt1,file=polout,form='formatted', &
     	     status='unknown',position='rewind')
	   open(uprof,file=profiles,form='formatted', &
     	     status='unknown',position='rewind')

	   do ir = 0,n1
	      r = rin+(rout-rin)/float(n1)*ir
	      i1 = int((r-rin)/dr)
	      i1 = min(i1,nr-1)
	      wr0 = (rin+(i1+1)*dr-r)/dr
	      wr1 = 1.-wr0
	      rmajr = rmaj(i1)*wr0+rmaj(i1+1)*wr1
	      elonr = elon(i1)*wr0+elon(i1+1)*wr1
	      triar = tria(i1)*wr0+tria(i1+1)*wr1
	      rplt(ir) = r/(2*r0)
	      qplt(ir) = sf(i1)*wr0+sf(i1+1)*wr1
	      qr = wr0*sf(i1)+wr1*sf(i1+1)
	      write(uprof,120)qr,rmajr,elonr,triar
	      grdm(:,ir) = 0.
	      do ith = 0,n2
		 th = -pi+pi2/float(n2)*ith
		 j1 = int((th+pi)/dth)
		 j1 = min(j1,ntheta-1)
		 wz0 = (-pi+(j1+1)*dth-th)/dth
		 wz1 = 1.-wz0

		 radiusp = wr0*wz0*radius(i1,j1)+wr0*wz1*radius(i1,j1+1) &
                     +wr1*wz0*radius(i1+1,j1)+wr1*wz1*radius(i1+1,j1+1) 
		 hghtp = wr0*wz0*hght(i1,j1)+wr0*wz1*hght(i1,j1+1) &
                     +wr1*wz0*hght(i1+1,j1)+wr1*wz1*hght(i1+1,j1+1) 
		 thf = wr0*wz0*thflx(i1,j1)+wr0*wz1*thflx(i1,j1+1) &
                     +wr1*wz0*thflx(i1+1,j1)+wr1*wz1*thflx(i1+1,j1+1) 
		 qhatp = wr0*wz0*qhat(i1,j1)+wr0*wz1*qhat(i1,j1+1) &
                     +wr1*wz0*qhat(i1+1,j1)+wr1*wz1*qhat(i1+1,j1+1) 
		     
		 rtor = sqrt((radiusp-rmaj0)**2+hghtp**2)
		 thtor = asin(hghtp/rtor)
		 if(radiusp < rmaj0)then
		    if(thtor .ge. 0.0)thtor=pi-thtor
		    if(thtor .lt. 0.0)thtor=-pi-thtor
		 end if

		 do izeta = 0,nzeta
		    
		    yt = wr0*wz0*yfn(i1,j1)+wr0*wz1*yfn(i1,j1+1) &
                          +wr1*wz0*yfn(i1+1,j1)+wr1*wz1*yfn(i1+1,j1+1) - r0/q0*pi2/nzeta*izeta
		    yt = dmod(yt+800.*ly,ly)    
		    zt = wz0*zfnth(j1)+wz1*zfnth(j1+1)
		    zt = min(zt,lz-1.e-2)
		    xt = r-rin

		    i=int(xt/dx)
		    i = min(i,imx-1)
		    j=int(yt/dy)
		    j = min(j,jmx-1)
		    k=int(zt/dz)
		    k = min(k,kmx-1)

		    wx0=float(i+1)-xt/dx 
		    wx1=1.-wx0
		    wy0=float(j+1)-yt/dy
		    wy1=1.-wy0
		    wz0=float(k+1)-zt/dz
		    wz1=1.-wz0
		    dum = wx0*wy0*wz0*grd(i,j,k)  &
         		 + wx1*wy0*wz0*grd(i+1,j,k) &
		         + wx0*wy1*wz0*grd(i,j+1,k) &
		         + wx1*wy1*wz0*grd(i+1,j+1,k) &
		         + wx0*wy0*wz1*grd(i,j,k+1) &
		         + wx1*wy0*wz1*grd(i+1,j,k+1) &
		         + wx0*wy1*wz1*grd(i,j+1,k+1) &
		         + wx1*wy1*wz1*grd(i+1,j+1,k+1)

!		    write(unt1,120)r,th,pi2/nzeta*izeta,dum
		    write(unt1,120)radiusp,hghtp,pi2/nzeta*izeta,dum
 120		    format (4(2x,e12.5))
		 enddo
	      end do
	   end do
	   endfile unt1
	   close(unt1)
	   endfile uprof
	   close(uprof)

	 endif

	 return
	 end

!--------------------------------------------------------------
