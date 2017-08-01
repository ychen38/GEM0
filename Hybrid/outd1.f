	 subroutine outd(n)
	 use gem_com
	 implicit none
	 INTEGER :: n,np
	 integer :: i,j,k
	 REAL(8) :: ur(0:imx,0:jmx,0:1),ui(0:imx,0:jmx,0:1)

	 if ( (mod(n,imovie).eq.0) ) then
	   call mphxy(phi(:,:,:),'mphxy',25)
	   call mphxz(phi(:,:,:),'mphxz',26)
	 endif

	 if (mod(n,nplot).eq.0) then
	   call phixy(rho,'rhoxy',27,n)
	   call phixz(rho,'rhoxz',28,n)

	   call aphir(phi(:,:,:),'aphir',30,n)
!	   call ffty(phi(:,:,:),ur,ui)
	   call phixy(phi(:,:,:),'phixy',31,n)
	   call phixz(phi(:,:,:),'phixz',32,n)
!	   call phixz(ur(:,:,:),'prexz',33,n)
!	   call phixz(ui(:,:,:),'pimxz',34,n)

	   call aphir(apar(:,:,:),'aapar',35,n)
!	   call ffty(apar(:,:,:),ur,ui)
	   call phixy(apar(:,:,:),'apaxy',36,n)
	   call phixz(apar(:,:,:),'apaxz',37,n)
!	   call phixz(ur(:,:,:),'arexz',53,n)
!	   if (n==nm-nplot) call xinxz(ur(:,:,:),n)
!	   call phixz(ui(:,:,:),'aimxz',54,n)

!	   call aphir(upar(:,:,:),'aupar',38,n)
!	   call phixy(upar(:,:,:),'upaxy',39,n)
!	   call phixz(upar(:,:,:),'upaxz',40,n)

	   call phixy(dti(:,:,:),'tpixy',29,n)
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
 99	 format('time step= ',I5)
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
	     sbuf(m)=grd(i,1,k)
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
 99	    format('time step= ',I5)

	    do i=0,((im+1)*pkm)*numprocs/ntube-1
	       write(unt,100)rbuf(i)
	    enddo
	    do k=pkm,mykm
	       do i=0,im
		  write(unt,100)grd(i,1,k)
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
	 enddo
	 do k=0,mykm-1
	   do j=0,jm-1
	     do i=0,im
	       myaph(i)=myaph(i)+grd(i,j,k) 
	     enddo
	   enddo
	 enddo

	 cnt=imx+1
	 call MPI_REDUCE(myaph(0:imx),aph(0:imx),cnt, &
     			MPI_REAL8, &
     			MPI_SUM,Master, &
     			MPI_COMM_WORLD,ierr)
!      
	 if (MyId.eq.Master) then
	 open(unt,file=flnm,form='formatted', &
     	   status='unknown',position='append')

	 write(unt,99)n
 99	 format('time step= ',I5)

	 do i=0,im
	   aph(i)=aph(i)/(ntube*km*jm)
	   write(unt,100)aph(i)
	 enddo

	 if(izon.eq.1)then
	    do i=0,im-1
	       tmpx(i) = aph(i) 
	    enddo
	    call ccfft('x',1,imx,1.,tmpx,tmpx,coefx,workx,0)
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

 100	 format (e10.3)
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
	   if(iget.eq.0)write(unt,110)im,jm,km
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
! 99	    format('time step= ',I5)

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
	 integer :: i,j,i1,j1,k,n,ir,ith,n1,n2, &
                    unt=71,unt1=72,unt2=73,upol=74
         integer :: m,mpol(1:4) 
	 parameter (n1=300,n2=2000)
	 REAL(8) :: grd(0:imx,0:jmx,0:kmx)
	 real(8) :: r,th,dum,xt,yt,zt,dum1,dum2,wx0,wx1,wy0,wy1,wz0,wz1
	 real(8) :: hghtp,radiusp,wr0,wr1,rplt(0:n1),qplt(0:n1)
	 complex(8) :: grdm(1:4,0:n1)
	 character*5 fl
	 character*70 flnm,polout,bdout,mpolout1,mpolout2,mpolout3,mpolout4

	 mpol(1) = 1
	 mpol(2) = 2
	 mpol(3) = 3
	 mpol(4) = 4

	 if(isphi==0)flnm=outdir//'apa3d'//'.out'
	 if(isphi==1)flnm=outdir//'phi3d'//'.out'
	 polout=outdir//'pol'//'.out'
	 bdout=outdir//'bd'//'.out'
	 mpolout1=outdir//'mpol1'//'.out'
	 mpolout2=outdir//'mpol2'//'.out'
	 mpolout3=outdir//'mpol3'//'.out'
	 mpolout4=outdir//'mpol4'//'.out'

	 if (MyId.eq.Last) then
	   open(unt,file=flnm,form='formatted', &
     	     status='old')
	   read(unt,110)im,jm,km
	   do n=1,ipg
	      do k = 0,kmx
		 do j = 0,jmx
		    do i = 0,imx
		       read(unt,100)grd(i,j,k)
		    end do
		 end do
	      end do
	   end do
	   endfile unt
	   close(unt)

	   open(unt1,file=polout,form='formatted', &
     	     status='unknown')

	   do ir = 0,n1
	      r = rin+(rout-rin)/float(n1)*ir
	      i1 = int((r-rin)/dr)
	      i1 = min(i1,nr-1)
	      wr0 = (rin+(i1+1)*dr-r)/dr
	      wr1 = 1.-wr0
	      rplt(ir) = r/(2*r0)
	      qplt(ir) = sf(i1)*wr0+sf(i1+1)*wr1
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
		 do m = 1,4
		    grdm(m,ir) = grdm(m,ir)+dum*exp(IU*mpol(m)*th)*pi2/n2
		 end do
		 write(unt1,120)radiusp,hghtp,dum
	      end do
	   end do
	   endfile unt1
	   close(unt1)

	   open(upol,file=mpolout1,form='formatted', &
     	     status='unknown')
	   do ir = 0,n1
	      write(upol,140)rplt(ir),abs(grdm(1,ir)),aimag(grdm(1,ir)),qplt(ir)
	   end do
	   endfile upol
	   close(upol)

	   open(upol,file=mpolout2,form='formatted', &
     	     status='unknown')
	   do ir = 0,n1
	      write(upol,140)rplt(ir),abs(grdm(2,ir)),aimag(grdm(2,ir))
	   end do
	   endfile upol
	   close(upol)

	   open(upol,file=mpolout3,form='formatted', &
     	     status='unknown')
	   do ir = 0,n1
	      write(upol,140)rplt(ir),abs(grdm(3,ir)),aimag(grdm(3,ir))
	   end do
	   endfile upol
	   close(upol)

	   open(upol,file=mpolout4,form='formatted', &
     	     status='unknown')
	   do ir = 0,n1
	      write(upol,140)rplt(ir),abs(grdm(4,ir)),aimag(grdm(4,ir))
	   end do
	   endfile upol
	   close(upol)

	   open(unt2,file=bdout,form='formatted', &
     	     status='unknown')

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
	      end do
	   end do
	   endfile unt2
	   close(unt2)
	 endif

 100	 format (e10.3)
 110	 format (3I5)
 120	 format (3(2x,e10.3))
 130	 format (2(2x,e10.3))
 140	 format (4(2x,e10.3))
	 return
	 end

!--------------------------------------------------------------
      subroutine ffty(u,ur,ui)

      use gem_com
      use equil
      use fft_wrapper
      implicit none
!     
      REAL(8) :: u(0:imx,0:jmx,0:1),ur(0:imx,0:jmx,0:1),ui(0:imx,0:jmx,0:1)
      REAL(8) :: ur1(0:imx-1,0:jmx-1,0:1),ui1(0:imx-1,0:jmx-1,0:1)	
      INTEGER :: n,i,j,k
      COMPLEX(8) :: tmp3d(0:imx-1,0:jmx-1,0:1)

     do j=0,jmx
        do i=0,imx
            do k = 0,1
               ur(i,j,k)=0.
	       ui(i,j,k)=0.
            enddo
         end do
      end do

      do j=0,jm-1
         do i=0,im-1
            do k = 0,1
               tmp3d(i,j,k)=u(i,j,k)
            enddo
         end do
      end do

      do k = 0,1
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
      ur1 = real(tmp3d)/jmx
      ui1 = aimag(tmp3d)/jmx

      do j=0,jmx-1
         do i=0,imx-1
            do k = 0,1
               ur(i,j,k)=ur1(i,j,k)
	       ui(i,j,k)=ui1(i,j,k)
            enddo
         end do
      end do

      return
      end
!--------------------------------------------------------------

	 subroutine xinxz(grd,n)

	 use gem_com
	 implicit none
	 integer :: i,j,k,n
	 REAL(8) :: grd(0:imx,0:jmx,0:1)
	 REAL(8) :: sbuf(0:(imx+1)*kmx)
	 REAL(8) :: rbuf(0:(imx+1)*kmx)
	 INTEGER :: unt,tind,flag
	 INTEGER :: lst,m,pkm

	 pkm=int(ntube*km/numprocs)
	 do k=0,pkm-1
	   do i=0,im
	     m=k*(im+1)+i
	     sbuf(m)=grd(i,1,k)
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
	    open(UNIT=129,file=outdir//'azl.out',form='formatted', &
	        status='unknown',position='append')
            do i=0,((im+1)*pkm)*numprocs/ntube-1
	       if(mod(i,(im+1))==44) write(129,*)rbuf(i)
	    enddo
            do k=pkm,mykm
	       write(129,*)grd(44,1,k)	
	    enddo	
            endfile 129
	    close(129)
	 endif

 100	 format (e10.3)
 110	 format (5I5)

	 return
	 end
