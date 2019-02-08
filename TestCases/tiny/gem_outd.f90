subroutine outd(n)

  use gem_com

  implicit none
  INTEGER :: n,np
  integer :: i,j,k

  if ( (mod(n,imovie).eq.0) ) then
     call mphxy(phi(:,:,:),'mphxy',25)
     call mphxz(phi(:,:,:),'mphxz',26)
  endif

  if (mod(n,nplot).eq.0) then
     !         call phixy(rho,'rhoxy',27,n)
     !         call phixz(rho,'rhoxz',28,n)

     call aphir(phi(:,:,:),'aphir',30,n)
     call phixy(phi(:,:,:),'phixy',31,n)
     call phixz(phi(:,:,:),'phixz',32,n)

     call aphir(apar(:,:,:),'aapar',35,n)
     call phixy(apar(:,:,:),'apaxy',36,n)
     call phixz(apar(:,:,:),'apaxz',37,n)

     !         call aphir(upar(:,:,:),'aupar',38,n)
     !         call phixy(upar(:,:,:),'upaxy',39,n)
     !         call phixz(upar(:,:,:),'upaxz',40,n)

     !         call phixy(dti(1,:,:,:),'tpixy',29,n)
     !         call dump3d(phi(:,:,:),'phi3d',50,n)
     !         call dump3d(apar(:,:,:),'apa3d',51,n)

  endif

  !       if (n.eq.nm-1) call histout(34)

  !         call timephi(phi(:,:,:),57,'tmphi',43,n)

  !       return
end subroutine outd

!--------------------------------------------------------------

subroutine phixy(grd,fl,unt,n)

  use gem_com

  implicit none

  REAL :: grd(0:imx,0:jmx,0:1)
  integer :: i,j,k,n
  INTEGER :: unt,tind,flag,oproc,outk
  character(len=5) fl
  character(len=70) flnm
  !       save flag

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
99   format('time step= ',I6)
     do j=0,jm
        do i=0,im
           write(unt,100)grd(i,j,outk)
        enddo
     enddo


100  format (e10.3)
110  format (5I5)
     endfile unt
     close(unt)
  endif

  !       return
end subroutine phixy

!--------------------------------------------------------------

subroutine phixz(grd,fl,unt,n)

  use gem_com

  implicit none

  integer :: i,j,k,n
  REAL :: grd(0:imx,0:jmx,0:1)
  REAL :: sbuf(0:(imx+1)*kmx)
  REAL :: rbuf(0:(imx+1)*kmx)
  INTEGER :: unt,tind,flag
  INTEGER :: lst,m,pkm
  character(len=5) fl
  character(len=70) flnm
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
99   format('time step= ',I6)

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

100 format (e10.3)
110 format (5I5)

  !       return
end subroutine phixz

!--------------------------------------------------------------

subroutine mphxy(grd,fl,unt)

  use gem_com

  implicit none

  integer :: i,j,k
  REAL :: grd(0:imx,0:jmx,0:1)
  INTEGER :: unt,tind,flag,oproc,outk
  character(len=5) fl
  character(len=70) flnm
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


100  format (e9.2)
110  format (5I5)
     endfile unt
     close(unt)
  endif

  !       return
end subroutine mphxy

!--------------------------------------------------------------

subroutine mphxz(grd,fl,unt)

  use gem_com

  implicit none

  integer :: i,j,k,n
  REAL :: grd(0:imx,0:jmx,0:1)
  REAL :: sbuf(0:(imx+1)*kmx)
  REAL :: rbuf(0:(imx+1)*kmx)
  INTEGER :: unt,tind,flag
  INTEGER :: lst,m,pkm
  character(len=5) fl
  character(len=70) flnm
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

100 format (e9.2)
110 format (5I5)

  !       return
end subroutine mphxz

!--------------------------------------------------------------

subroutine histout(unt)

  use gem_com
  use gem_equil

  implicit none

  integer :: i,j,k
  INTEGER :: unt,n,m
  character(len=70) flnm


  if (MyId.eq.Master) then

     flnm=outdir//'hist.out'
     open(unt,file=flnm,form='formatted', &
          status='unknown',position='rewind')
     write(unt,80)nsm,nm,modem,dt
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
           write(unt,110)pfl_es(m,1,n)
        enddo
     enddo
     do m=1,nsm
        do n=1,nm
           write(unt,110)efl_es(m,1,n)
        enddo
     enddo
     do n=1,nm
        write(unt,110)fe(n)
     enddo
     do n=1,nm
        write(unt,110)rmsphi(n)
     enddo

80   format (3I6,e10.3)
90   format (3I6)
100  format (2e10.3)
110  format (1e10.3)

     endfile unt
     close(unt)

     ! output particle and energy fluxes to file
     flnm=outdir//'pefluxes.out'
     open(18,file=flnm,form='formatted', &
          status='unknown',position='rewind')
     write(18,2) nsm,nsubd,nm,tir0,xnir0,xu,frequ,vu,rina,routa
2    format(' no. ion species, no. radial slices, no. timesteps, '/&
          'energy unit (eV), ref. density (cm**-3), length unit (cm), '/&
          'frequency unit (1/sec), veloc. unit (cm/sec),'/&
          ' r-in/a, r-out/a'/3i5,5x,1p3e20.6/4e20.6)
     write(18,1)
1    format(/' species 0 (electrons)')
     write(18,3)
3    format(/' particle flux(1<=iradius<=nsubd,1<=timestep<=nm)'/)
     write(18,4) (n,(pfle_es(i,n),i=1,nsubd),n=1,nm)
4    format((i5,1p8e13.5))
     write(18,5)
5    format(/' energy flux(1<=iradius<=nsubd,1<=timestep<=nm)'/)
     write(18,4) (n,(efle_es(i,n),i=1,nsubd),n=1,nm)
     do m=1,nsm
        write(18,6) m
6       format(/' ion species',i2)
        write(18,3)
        write(18,4) (n,(pfl_es(m,i,n),i=1,nsubd),n=1,nm)
        write(18,5)
        write(18,4) (n,(efl_es(m,i,n),i=1,nsubd),n=1,nm)
     end do

  endif

  !       return
end subroutine histout

!--------------------------------------------------------------

subroutine aphir(grd,fl,unt,n)

  use gem_com
  use gem_fft_wrapper

  implicit none

  integer :: i,j,k,n
  REAL :: grd(0:imx,0:jmx,0:1)
  REAL :: aph(0:imx)
  REAL :: myaph(0:imx)
  INTEGER :: unt,tind,flag
  character(len=5) fl
  character(len=70) flnm
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
99   format('time step= ',I6)

     do i=0,im
        aph(i)=aph(i)/(ntube*km*jm)
        write(unt,100)aph(i)
     enddo

     if(izonal.eq.2)then
        do i=0,im-1
           tmpx(i) = aph(i) 
        enddo
        call ccfft('x',1,imx,1.0,tmpx,coefx,workx,0)
        write(*,98)n,aimag(tmpx(0)),aimag(tmpx(1)),aimag(tmpx(2)), &
             aimag(tmpx(3)),aimag(tmpx(4))

        !        write(*,98)n,real(tmpx(0)),real(tmpx(1)),real(tmpx(2)),
        !     1              real(tmpx(3)),real(tmpx(4))

        write(15,98)n,aimag(tmpx(0)),aimag(tmpx(1)),aimag(tmpx(2)),&
             aimag(tmpx(3)),aimag(tmpx(4))
98      format(1x,i5, 5(2x, e12.3))
     end if
     endfile unt
     close(unt)
  endif

100 format (e10.3)
110 format (4I5)

  !       return
end subroutine aphir

!ccccccccccccccccccccccccc

subroutine dump3d(grd,fl,unt,n)

  use gem_com
  implicit none
  integer :: i,j,k,n
  REAL :: grd(0:imx,0:jmx,0:1)
  REAL :: sbuf(0:(imx+1)*(jmx+1))
  REAL :: rbuf(0:(imx+1)*(jmx+1)*kmx)
  INTEGER :: unt,tind,flag
  INTEGER :: lst,m,pkm
  character(len=5) fl
  character(len=70) flnm
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
     ! 99          format('time step= ',I5)

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

100 format (e10.3)
110 format (3I5)

  !       return
end subroutine dump3d

!--------------------------------------------------------------
subroutine pol2d

  use gem_com
  use gem_equil

  implicit none

  integer :: i,j,i1,j1,k,n,ir,ith,n1,n2, numx,numz,&
       unt=71,unt1=72,unt2=73,upol=74,upolxz=75
  parameter (n1=5000,n2=10000,numx=400,numz=600)
  integer :: m,mpol(1:25),npts(0:numx,0:numz)
  REAL :: grd(0:imx,0:jmx,0:kmx),yxz(0:numx,0:numz)
  real :: r,th,thf,qhatp,qr,dum,xt,yt,zt,dum1,dum2,wx0,wx1,wy0,wy1,wz0,wz1
  real :: xmin,xmax,zmin,zmax,dumx,dumz
  real :: hghtp,radiusp,wr0,wr1,rplt(0:n1),qplt(0:n1)
  complex :: grdm(1:25,0:n1)
  character(len=5) fl
  character(len=70) flnm,polout,bdout,mpolout,polxzout

  do i = 1,25
     mpol(i) = -1+i
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
     !         read(unt,110)im,jm,km
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
     !         endfile unt
     close(unt)

     !remove n=0
     do k = 0,kmx
        do i = 0,imx
           dum = 0.
           do j = 0,jmx
              dum = dum+grd(i,j,k)
           end do
           dum = dum/jmx
           grd(i,:,k) = grd(i,:,k)-dum
        end do
     end do


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
        rplt(ir) = r/a
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
           yt = mod(yt+800.0*ly,ly)    
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

100 format (e12.5)
110 format (3I5)
120 format (3(2x,e12.5))
130 format (2(2x,e12.5))
140 format (1x,I5,4(2x,e12.5))
  !       return
end subroutine pol2d

!--------------------------------------------------------------

subroutine timephi(grd,rim,fl,unt,n)

  use gem_com

  implicit none

  integer :: i,j,k,n
  REAL :: grd(0:imx,0:jmx,0:1)
  REAL :: aph(0:imx),yph(0:imx)
  REAL :: myaph(0:imx)
  INTEGER :: rim,unt,oproc,tnow
  character(len=5) fl
  character(len=70) flnm

  flnm=outdir//fl//'.out'

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

  oproc=int(km/2/kcnt*ntube)
  ! get grd(r,Ly/2,Lz/2)...         
  if(myid==oproc) then
     do i=0,im
        yph(i) = grd(i,jm/2,0)
     enddo
  endif
  call MPI_BCAST(yph,im+1,MPI_REAL8,oproc, &
       MPI_COMM_WORLD,ierr)

  if (MyId.eq.Master) then
     open(unt,file=flnm,form='formatted', &
          status='unknown',position='append')

     tnow=tcurr-dt

     do i=0,im
        aph(i)=aph(i)/(ntube*km*jm)
     enddo
     write(unt,99)tnow,aph(rim),yph(rim)
99   format(1x,I5,2(2x,e10.3))

     endfile unt
     close(unt)
  endif

  !       return
end subroutine timephi

!--------------------------------------------------------------
