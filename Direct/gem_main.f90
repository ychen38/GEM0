!       3D Flux Tube Toroidal Electromagnetic GK Delta-f Code
!   global variables...
program gem_main

  use gem_com
  use gem_equil
  use gem_fft_wrapper

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

  if(isft==1)then
     call ftcamp
     goto 100
  end if

  do  timestep=ncurr,nm
     tcurr = tcurr+dt

     call accumulate(timestep-1,0)
     call ampere(timestep-1,0)
     call poisson(timestep-1,0)
     call field(timestep-1,0)
     call split_weight(timestep-1,0)
     call diagnose(timestep-1)
     call reporter(timestep-1)

     call push_wrapper(timestep,1)

     call accumulate(timestep,1)
     call ampere(timestep,1)
     call poisson(timestep,1)
     call field(timestep,1)
     call split_weight(timestep,1)

     call push_wrapper(timestep,0)
     if(mod(timestep,1000)==0)then
        do i=0,last 
           !                 if(myid==i)write(*,*)myid,mm(1),mme
           call MPI_BARRIER(MPI_COMM_WORLD,ierr)
        end do
     end if

  end do
  call ftcamp
  lasttm=MPI_WTIME()
  tottm=lasttm-starttm
  !  write(*,*)'ps time=',pstm,'tot time=',tottm
  do i=0,last 
     !            if(myid==i)write(*,*)myid,mm(1),mme
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  end do

100 call MPI_FINALIZE(ierr)

end program gem_main

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine hybinit

  use gem_com
  use gem_equil
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
  !      return
end subroutine hybinit

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine init

  use gem_com
  use gem_equil
  implicit none
  character(len=62) dumchar
  INTEGER :: i,j,k,m,n,ns,idum,i1,j1,k1
  INTEGER :: mm1,mm2,lr1
  real :: mims1,tets1,q1,kappan,kappat,r,qr,th,cost,dum,zdum
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp
  real :: grp,gxdgyp,jacp,jfnp,gn0ep,gt0ep,gt0ip,grdgtp,gthp
  real,DIMENSION (1:5) :: gn0sp
  real :: wx0,wx1,wz0,wz1,b

  !jycheng
  namelist /primary_parameters/ itube,mimp,mcmp,chgi,chgc,imx,jmx,kmx,mmx,mmxe,nmx,nsmx,ntube,lxa, &
                                lymult,jcnt,dt,nm,nsm,amp,fradi,r0a,vpp,vt0,yd0,ifluid,isg,amie,rneui, &
                                betai,nonlin1,nonlin2,nonline,ipara,vwidth,vwidthe,vcut,micell,mecell, &
                                mims3,q3
  namelist /control_parameters/ iperi,iperidf,ibunit,modemx,delra,delri,delre,delrn,nlow,xshape, &
                                yshape,zshape,iput,iget,ision,isiap,peritr,llk,mlk,onemd,izonal, &
                                adiabatic_electron,ineq0,iflut,nzcrt,npze,npzi,npzc,npzb,iphbf,iapbf, &
                                idpbf,cut,kxcut,kycut,bcut,c4,vexbsw,vparsw,mach,gamma_E,isuni, &
                                lr1,iflr,iorb
  namelist /diagnosis_parameters/ icrs_sec,ipg,isphi,nplot,xnplt,modem,isft,mynf,frmax,ifskp,idg
  namelist /fluxtube/ Rovera,elon0,selon0,tria0,stria0,rmaj0p,q0,shat0,teti,tcti,rhoia,Rovlni,Rovlti, &
                       Rovlne,Rovlte,Rovlnc,Rovltc,ncne,nuacs
  namelist /others/ nrst,eprs,tor,ishift,width

  IU=cmplx(0.,1.)
  pi=4.0*atan(1.0)
  pi2 = pi*2.

  call ppinit_mpi(myid,numprocs)
  last=numprocs-1
  !the initial timestep index
  timestep=0
  !the initial timestep 
  tcurr=0.
  open(unit=115,file='gem.in',status='old',action='read')
  read(115,nml=primary_parameters)
  read(115,nml=control_parameters)
  read(115,nml=diagnosis_parameters)
  read(115,nml=fluxtube)
  read(115,nml=others)
  close(115)
  nonlin(1)=nonlin1
  nonlin(2)=nonlin2
  ntube=int(numprocs/kmx)
  if((mod(numprocs,kmx) .ne. 0) .and. myid==0)then
    write(*,*)'WARNING: MPI number is a multiple of kmx, ntube is not an integer'
    stop
  endif

  !mm1=int(imx,8)*int(jmx,8)*int(kmx,8)*int(micell,8)
  mm1=imx*jmx*kmx*micell
  !mm2=int(imx,8)*int(jmx,8)*int(kmx,8)*int(mecell,8)
  mm2=imx*jmx*kmx*mecell
  mmx=int(real(mm1/int(kmx*ntube))*1.5)
  mmxe=int(real(mm2/int(kmx*ntube))*1.5)
  call ppinit_decomp(myid,numprocs,ntube,tube_comm,grid_comm)
  call hybinit
  call mpi_barrier(mpi_comm_world,ierr)
  
  im=imx;jm=jmx;km=kmx

  if(isft==1) iget = 1
  nfreq = kmx*mynf
  call new_gem_com()
  !      rin = r0-lx/2
  !      rout = r0+lx/2
  rina = r0a-lxa/2.
  routa = r0a+lxa/2.

  call new_equil()
  lx=lxa*a
  ly=2.*pi*r0/q0abs/lymult
  !beta has been used all the time, so the definition is fixed, declared in gem_com.f. betai is declared in equil.f
  !in gem.in betai is intended to be GYRO's beta_e (within a factor of 2). It is defined to be (mu0 ne Te)/Bunit**2
  beta=betai 
  br0 = rmaj0
  lr0 = r0
  qp = q0p
  lz = pi2*q0abs*rmaj0
  delz = lz/ntheta
  rneu = nuacs
  if(onemd==1)then
     kxcut=pi2*shat0*pi2/ly*2.1
     kycut=pi2/ly*1.1
  end if
  achii=0.
  achie=0.
  addi=0.

  if(myid.eq.master)then
     open(19, file='eqdat', status='unknown')
     write(19,*) 'r0a,lxa,beta,a,rmaj0=', r0a,lxa,beta,a,rmaj0
     write(19,*)'elon0,tria0,rmaj0p,stria0,selon0,q0,shat=', elon0,tria0,rmaj0p,stria0,selon0,q0,r0*q0p/q0
     write(19,*)'xn0s(1,nr/2),xn0s(2,nr/2),t0s(1,nr/2),t0s(2,nr/2)=',xn0s(1,nr/2),xn0s(2,nr/2),t0s(1,nr/2),t0s(2,nr/2)
     write(19,*)'capns(1,nr/2),capns(2,nr/2),capts(1,nr/2),capts(2,nr/2)=', capns(1,nr/2),capns(2,nr/2),capts(1,nr/2),capts(2,nr/2)
     write(19,*)'capne(nr/2),capte(nr/2)=', capne(nr/2),capte(nr/2)
     write(19,*)'t0i0,t0e0,ni0,ne0=', t0i(nr/2),t0e(nr/2),xn0i(nr/2),xn0e(nr/2)
     write(19,*)'a/cs=', a/sqrt(t0e(nr/2)/mimp)
     write(19,*)'i,   sf,   ni,   ne,   nc,   ti,   tc,   capni,   capnc, captc'
     do i = 0,nr
        write(19,99)i,sf(i),xn0i(i),xn0e(i),t0i(i),t0e(i),capni(i),capne(i),capti(i),capte(i)
     end do
     write(19,99)i,cn0e,cn0s(1),cn0s(2),cn0s(3)
     close(19)
  end if
98 format(7(1x,e14.7))
99 format(1x,i3,2x,15(1x,e10.3))
  ! equilibrium data for transport coefficients calculation by mablab with flux data
  if(myid.eq.master)then
     open(912, file='eqflux', status='unknown')
     do j = 1, nsubd
        i = int(nr*(2*j-1)/(2*nsubd))
        write(912,9991) xn0e(i)*cn0e,capne(i),xn0c(i)*cn0c,capnc(i),t0i(i),capti(i)
     enddo
     close(912)
  endif
9991 format(6(2x,e12.4)) 
  !
  !      write(*,*)'br0,lr0,q0,qp = ', br0,lr0,q0,qp

  if(isuni.eq.0)vwidthe=vwidth
  dte = dt
  iadi = 0
  if(isg.gt.0.)fradi = isg
  if(ifluid.eq.0)then
     iadi = 1
     fradi = cn0e*xn0e(nr/2)/t0e(nr/2)
  end if

  !     begin reading species info, ns=1,nsm...
  ns = 1
  mims(3)=mims3
  q(3)=q3
  mims(1)=mimp;mims(2)=mcmp;q(1)=chgi;q(2)=chgc
  tmm(1)=mm1
  tmm(2)=mm2
  mm(:)=int(mm1/numprocs)
  mme = int(mm2/numprocs)
  if (MyId.eq.Last) mm(ns)=mm1-Last*mm(ns)
  !     write(*,*)'in init  ',Myid,mm(ns)
  tets(1)=1
  lr(1)=lr1
  lr(2)=lr1
  lr(3)=lr1

  pzcrite = abs(psi(nr)-psi(0))/br0/npze
  encrit = 1.
  pzcrit(1) = q(1)*abs(psi(nr)-psi(0))/br0/npzi
  pzcrit(2) = q(2)*abs(psi(nr)-psi(0))/br0/npzc
  pzcrit(3) = q(3)*abs(psi(nr)-psi(0))/br0/npzb

  kapn(ns)=kappan
  kapt(ns)=kappat

  emass = 1./amie
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
  do i=0,nxpp
     xg(i)=i*dx !dx*(tclr*nxpp+i)
  enddo
  do j=0,jm
     yg(j)=dy*float(j)
  enddo
  kcnt=1

  do k=0,mykm
     n=GCLR*kcnt+k 
     zg(k)=dz*float(n)
  enddo

  !      jcnt = 3  !jmx/ntube                                                                                                                                                                         
  mstart = 0
  ntor0 = mstart+1
  do m = 0,jcnt-1
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
127     continue
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
128     continue
     end do
  end if

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

        grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
             +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1) 
        gthp = wx0*wz0*gth(i,k)+wx0*wz1*gth(i,k+1) &
             +wx1*wz0*gth(i+1,k)+wx1*wz1*gth(i+1,k+1) 

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
        gt0ep = wx0*t0e(i)+wx1*t0e(i+1)        
        do ns = 1, nsm
           gn0sp(ns) = wx0*xn0s(ns,i)+wx1*xn0s(ns,i+1)
           gn0s(ns,i1) = gn0sp(ns)
        enddo
        gt0ip = wx0*t0s(1,i)+wx1*t0s(1,i+1)        
        b=1.-tor+tor*bfldp
        cfx(i1,k1) = br0/b**3*fp/radiusp*dbdtp*grcgtp
        cfy(i1,k1) = br0/b**3*fp/radiusp* &
             (dydrp*dbdtp-lr0/q0*qhatp*dbdrp)*grcgtp
        bdgxcgy(i1,k1) = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp             
        bmag(i1,k1) = b
        jac(i1,k1) = jacp*jfnp
        bdgrzn(i1,k1) = q0*br0/radiusp/b*psipp*grcgtp/jfnp
        gn0e(i1) = gn0ep
        gt0e(i1) = gt0ep
        gt0i(i1) = gt0ip

        ggxdgy(i1,k1) = gxdgyp
        ggy2(i1,k1) = dydrp**2*grp**2 + (r0/q0*qhatp*gthp)**2 + 2*dydrp*r0/q0*qhatp*grdgtp
        ggx(i1,k1) = grp
     end do
  end do

  iseed = -(1777+myid*13)
  idum = ran2(iseed)
  phi = 0.
  apar = 0.
  dene = 0.
  upar = 0.
  upa0 = 0.
  camp = 0.

  if(myid.eq.master)then
     write(*,*)zfnth(ntheta),thfnz(ntheta/2),thfnz(ntheta/2+1)
     if(myid.eq.master)open(9, file='plot', &
          status='unknown',position='append')
     write(9,*)'dt,beta= ',dt, beta
     write(9,*)'amp,vpp,yd0 = ',amp,vpp,yd0
     write(9,*)'peritr,ifluid= ',peritr,ifluid
     write(9,*)'tor,nonlin= ',tor,nonlin(1),nonlin(2)
     write(9,*)'isuni= ',isuni, 'amie= ',amie
     write(9,*)'kxcut,kycut,bcut= ',kxcut,kycut,bcut
     write(9,*)'fradi,isg= ',fradi,isg
     write(9,*)'llk,mlk,onemd =',llk,mlk,onemd
     write(9,*)'vwidth (i,e), vcut= ', vwidth,vwidthe,vcut
     write(9,*)'rneu= ', rneu
     write(9,*)'V-ExB switch= ', vexbsw
     write(9,*)'V-parallel switch= ', vparsw
     write(9,*)'mm1= ',mm1
     write(9,*)'pzcrite,encrit = ',pzcrite,encrit
     write(9,*) 'lxa,lymult,delra,r0a,rina,routa=',lxa,lymult,delra,r0a,rina,routa
     write(9,*) 'a,r0,rmaj0,q0,lx,ly,lz=',a,r0,rmaj0,q0,lx,ly,lz
     write(9,*) 't0,kyrhoi_local=',t0i(nr/2),2*pi*sqrt(mims(1))*sqrt(t0i(nr/2))/ly
     close(9)
  end if

  if(myid.eq.master)then
     write(*,*)'dt,beta= ',dt, beta
     write(*,*)'amp,vpp,yd0 = ',amp,vpp,yd0
     write(*,*)'peritr,ifluid= ',peritr,ifluid
     write(*,*)'tor,nonlin= ',tor,nonlin(1),nonlin(2)
     write(*,*)'isuni= ',isuni, 'amie= ',amie
     write(*,*)'kxcut,kycut,bcut= ',kxcut,kycut,bcut
     write(*,*)'fradi,isg= ',fradi,isg
     write(*,*)'llk,mlk,onemd =',llk,mlk,onemd
     write(*,*)'vwidth (i,e), vcut= ', vwidth,vwidthe,vcut
     write(*,*)'rneu= ', rneu
     write(*,*)'V-ExB switch= ', vexbsw
     write(*,*)'V-parallel switch= ', vparsw
     write(*,*)'mm1= ',mm1
     write(*,*)'pzcrite,encrit = ',pzcrite,encrit
     write(*,*)'nue0 = ',nue0(1),nue0(nr/2),nue0(nr-1)
     write(*,*)'xn0e(1),xnir0 = ',xn0e(1),xnir0
     write(*,*)'frequ, eru = ', frequ, eru
     write(*,*) 'lxa,lymult,delra,r0a,rina,routa=',lxa,lymult,delra,r0a,rina,routa
     write(*,*) 'a,r0,rmaj0,q0,lx,ly,lz=',a,r0,rmaj0,q0,lx,ly,lz
     write(*,*) 't0,kyrhoi_local=',t0i(nr/2),2*pi*sqrt(mims(1))*sqrt(t0i(nr/2))/ly
     write(*,*) 'coefu = ', xu**2*frequ
     write(*,*) 'ktheta*rhos = ',2*pi*sqrt(mims(1))*sqrt(t0e(nr/2))/ly
     write(*,*) 'cs/a, q0, q0p, s^hat = ',sqrt(t0e(nr/2)/2.)/a, q0, q0p, q0p/q0*r0
     write(*,*) 'rho* = rhos/a = ', sqrt(mims(1))*sqrt(t0e(nr/2))/a
     write(*,*) 'f0p,psip(nr/2),Bunit,candyf0p = ',f0p,psip(nr/2),bunit,candyf0p
     write(*,*) 'lxa min = ', ly*q0/(2*pi*r0*q0p)/a
     write(*,*) 't0i(nr/2)= ', t0i(nr/2)
     write(*,*) 'Gyrokrs = ', 2*pi*sqrt(mims(1))*sqrt(t0e(nr/2))/ly/bunit
  end if
  close(115)
  !      return
end subroutine init

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine ppush(n,ns)

  use gem_com
  use gem_equil
  implicit none
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
  real :: wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum,dum1,bstar
  INTEGER :: m,i,j,k,l,n,ns
  INTEGER :: np_old,np_new
  real :: rhog,vfac,kap,vpar,pidum,kaptp,kapnp,xnp
  real :: b,th,r,enerb,cost,sint,qr,laps,sz,ter
  real :: xt,xs,yt,xdot,ydot,zdot,pzdot,edot,pzd0,vp0
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp
  real :: grp,gxdgyp,rhox(4),rhoy(4),psp,pzp,vncp,vparspp,psip2p,bdcrvbp,curvbzp,dipdrp
  integer :: mynopi

  pidum = 1./(pi*2)**1.5*vwidth**3
  mynopi = 0
  nopi(ns) = 0
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
     ter = wx0*t0s(ns,i)+wx1*t0s(ns,i+1)        
     kaptp = wx0*capts(ns,i)+wx1*capts(ns,i+1)        
     kapnp = wx0*capns(ns,i)+wx1*capns(ns,i+1)        
     xnp = wx0*xn0s(ns,i)+wx1*xn0s(ns,i+1)        
     vncp = wx0*phincp(i)+wx1*phincp(i+1)        
     vparspp = wx0*vparsp(ns,i)+wx1*vparsp(ns,i+1)        
     psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
     dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)        
     b=1.-tor+tor*bfldp
     pzp = mims(ns)*u2(ns,m)/b-q(ns)*psp/br0

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
     do l=1,lr(1)
        !
        xs=x2(ns,m)+rhox(l) !rwx(1,l)*rhog
        yt=y2(ns,m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
        !
        !   particle can go out of bounds during gyroavg...
        xt=mod(xs+800.*lx,lx)
        yt=mod(yt+800.*ly,ly)
        xt = min(xt,lx-1.0e-8)
        yt = min(yt,ly-1.0e-8)

        include "ppushngp.h"
     enddo
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
     vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
     vp0 = vp0*vncp*vexbsw

     vpar = u2(ns,m)-q(ns)/mims(ns)*aparp*nonlin(ns)*0.
     bstar = b*(1+mims(ns)*vpar/(q(ns)*b)*bdcrvbp)
     enerb=(mu(ns,m)+mims(ns)*vpar*vpar/b)/q(ns)*b/bstar*tor

     kap = kapnp - (1.5-vfac/ter)*kaptp-vpar*mims(ns)/ter*vparspp*vparsw
     dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
     vxdum = (eyp/b+vpar/b*delbxp)*dum1
     xdot = vxdum*nonlin(ns) -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
     ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlin(ns) &
          +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
          (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0   &
          +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
          -mims(ns)*vpar**2/(q(ns)*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
          -dipdrp/radiusp*mims(ns)*vpar**2/(q(ns)*bstar*b)*grcgtp*lr0/q0*qhatp 
     zdot =  vpar*b/bstar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
          +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
          -1./b**2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp &
          -dipdrp/radiusp*mims(ns)*vpar**2/(q(ns)*bstar*b)*q0*br0*grcgtp/jfnp

     pzd0 = tor*(-mu(ns,m)/mims(ns)/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar &
          +mu(ns,m)*vpar/(q(ns)*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp
     pzdot = pzd0 + (q(ns)/mims(ns)*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
          +q(ns)/mims(ns)*(-xdot*delbyp+ydot*delbxp+zdot*dadzp))*ipara

     edot = q(ns)*(xdot*exp1+(ydot-vp0)*eyp+zdot*ezp)                      &
          +q(ns)*pzdot*aparp*tor     &
          +q(ns)*vpar*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)    &
          -q(ns)*vpar*delbxp*vp0

     x3(ns,m) = x2(ns,m) + 0.5*dt*xdot
     y3(ns,m) = y2(ns,m) + 0.5*dt*ydot
     z3(ns,m) = z2(ns,m) + 0.5*dt*zdot
     u3(ns,m) = u2(ns,m) + 0.5*dt*pzdot

     dum = 1-w2(ns,m)*nonlin(ns)*0.
     if(ildu.eq.1)dum = (tgis(ns)/ter)**1.5*exp(vfac*(1/tgis(ns)-1./ter))
     vxdum = (eyp/b+vpar/b*delbxp)*dum1
     !         vxdum = eyp+vpar/b*delbxp
     w3(ns,m)=w2(ns,m) + 0.5*dt*(vxdum*kap + edot/ter)*dum*xnp

     !         if(x3(ns,m)>lx .or. x3(ns,m)<0.)w3(ns,m) = 0.


     if(itube==1)go to 333
     if(abs(pzp-pzi(ns,m))>pzcrit(ns).or.abs(vfac-eki(ns,m))>0.2*eki(ns,m))then
        mynopi = mynopi+1
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
        u3(ns,m) = u0i(ns,m)     
        u2(ns,m) = u3(ns,m)
        w3(ns,m) = 0.
        w2(ns,m) = 0.
        x2(ns,m) = x3(ns,m)
        z2(ns,m) = z3(ns,m)
     end if

333  continue
     laps=anint((z3(ns,m)/lz)-.5)*(1-peritr)
     r=x3(ns,m)-0.5*lx+lr0
     i = int((r-rin)/dr)
     i = min(i,nr-1)
     i = max(i,0)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     qr = wx0*sf(i)+wx1*sf(i+1)
     y3(ns,m)=mod(y3(ns,m)-laps*2*pi*qr*lr0/q0*sign(1.0,q0)+8000.*ly,ly)
     if(x3(ns,m)>lx.and.iperidf==0)then
        x3(ns,m) = lx-1.e-8
        z3(ns,m)=lz-z3(ns,m)
        x2(ns,m) = x3(ns,m)
        z2(ns,m) = z3(ns,m)
        w2(ns,m) = 0.
        w3(ns,m) = 0.
     end if
     if(x3(ns,m)<0..and.iperidf==0)then
        x3(ns,m) = 1.e-8
        z3(ns,m)=lz-z3(ns,m)
        x2(ns,m) = x3(ns,m)
        z2(ns,m) = z3(ns,m)
        w2(ns,m) = 0.
        w3(ns,m) = 0.
     end if
     z3(ns,m)=mod(z3(ns,m)+8.*lz,lz)
     x3(ns,m)=mod(x3(ns,m)+800.*lx,lx)         
     x3(ns,m) = min(x3(ns,m),lx-1.0e-8)
     y3(ns,m) = min(y3(ns,m),ly-1.0e-8)
     z3(ns,m) = min(z3(ns,m),lz-1.0e-8)

  enddo
  call MPI_ALLREDUCE(mynopi,nopi(ns),1,MPI_integer, &
       MPI_SUM, MPI_COMM_WORLD,ierr)

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
  call pmove(u0i(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call end_pmove(ierr)
  mm(ns)=np_new

  !      return
end subroutine ppush

!-----------------------------------------------------------------------

subroutine cpush(n,ns)

  use gem_com
  use gem_equil
  implicit none
  INTEGER :: n
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
  real :: wx0,wx1,wy0,wy1,wz0,wz1,w3old
  INTEGER :: m,i,j,k,ns,l
  INTEGER :: np_old,np_new
  real :: vfac,vpar,vxdum,dum,xdot,ydot,zdot,pzdot,edot,pzd0,vp0
  real :: rhog,xt,yt,zt,kap,xs,pidum,dum1,kaptp,kapnp,xnp
  real :: b,th,r,enerb,cost,sint,qr,laps,sz,ter,bstar
  real :: myke,mypfl_es(1:nsubd),mypfl_em(1:nsubd),myavewi
  real :: myefl_es(1:nsubd),myefl_em(1:nsubd),mynos
  real :: ketemp,pfltemp
  real :: efltemp,nostemp
  real :: sbuf(10),rbuf(10)
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp
  real :: grp,gxdgyp,rhox(4),rhoy(4),psp,pzp,vncp,vparspp,psip2p,bdcrvbp,curvbzp,dipdrp

  sbuf(1:10) = 0.
  rbuf(1:10) = 0.
  myavewi = 0.
  myke=0.    
  mypfl_es=0.    
  mypfl_em=0.    
  myefl_es=0. 
  myefl_em=0. 
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
     ter = wx0*t0s(ns,i)+wx1*t0s(ns,i+1)        
     kaptp = wx0*capts(ns,i)+wx1*capts(ns,i+1)        
     kapnp = wx0*capns(ns,i)+wx1*capns(ns,i+1)        
     xnp = wx0*xn0s(ns,i)+wx1*xn0s(ns,i+1)        
     b=1.-tor+tor*bfldp
     psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
     dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)        
     pzp = mims(ns)*u3(ns,m)/b-q(ns)*psp/br0
     vncp = wx0*phincp(i)+wx1*phincp(i+1)        
     vparspp = wx0*vparsp(ns,i)+wx1*vparsp(ns,i+1)        

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
     do l=1,lr(1)
        xs=x3(ns,m)+rhox(l) !rwx(1,l)*rhog
        yt=y3(ns,m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
        !   BOUNDARY
        xt=mod(xs+800.*lx,lx)
        yt=mod(yt+800.*ly,ly)
        xt = min(xt,lx-1.0e-8)
        yt = min(yt,ly-1.0e-8)

        include "cpushngp.h"
     enddo

     exp1=exp1/4.
     eyp=eyp/4.
     ezp=ezp/4.
     delbxp=delbxp/4.
     delbyp=delbyp/4.
     dpdzp = dpdzp/4.
     dadzp = dadzp/4.
     aparp = aparp/4.


     vfac = 0.5*(mims(ns)*u3(ns,m)**2 + 2.*mu(ns,m)*b)
     vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
     vp0 = vp0*vncp*vexbsw

     vpar = u3(ns,m)-q(ns)/mims(ns)*aparp*nonlin(ns)*0.
     bstar = b*(1+mims(ns)*vpar/(q(ns)*b)*bdcrvbp)
     enerb=(mu(ns,m)+mims(ns)*vpar*vpar/b)/q(ns)*b/bstar*tor

     kap = kapnp - (1.5-vfac/ter)*kaptp-vpar*mims(ns)/ter*vparspp*vparsw
     dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
     vxdum = (eyp/b+vpar/b*delbxp)*dum1
     xdot = vxdum*nonlin(ns) -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
     ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlin(ns)     &
          +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
          (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
          +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
          -mims(ns)*vpar**2/(q(ns)*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
          -dipdrp/radiusp*mims(ns)*vpar**2/(q(ns)*bstar*b)*grcgtp*lr0/q0*qhatp  
     zdot =  vpar*b/bstar*(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
          +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
          -1./b**2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp &
          -dipdrp/radiusp*mims(ns)*vpar**2/(q(ns)*bstar*b)*q0*br0*grcgtp/jfnp

     pzd0 = tor*(-mu(ns,m)/mims(ns)/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar &
          +mu(ns,m)*vpar/(q(ns)*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp
     pzdot = pzd0 + (q(ns)/mims(ns)*ezp*q0*br0/radiusp/b*psipp*grcgtp/jfnp  &
          +q(ns)/mims(ns)*(-xdot*delbyp+ydot*delbxp+zdot*dadzp))*ipara

     edot = q(ns)*(xdot*exp1+(ydot-vp0)*eyp+zdot*ezp)                      &
          +q(ns)*pzdot*aparp*tor     &
          +q(ns)*vpar*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)   &
          -q(ns)*vpar*delbxp*vp0

     x3(ns,m) = x2(ns,m) + dt*xdot
     y3(ns,m) = y2(ns,m) + dt*ydot
     z3(ns,m) = z2(ns,m) + dt*zdot
     u3(ns,m) = u2(ns,m) + dt*pzdot

     dum = 1-w3(ns,m)*nonlin(ns)*0.
     if(ildu.eq.1)dum = (tgis(ns)/ter)**1.5*exp(vfac*(1/tgis(ns)-1./ter))
     !         vxdum = eyp+vpar/b*delbxp
     vxdum = (eyp/b+vpar/b*delbxp)*dum1
     w3old = w3(ns,m)
     w3(ns,m) = w2(ns,m) + dt*(vxdum*kap+edot/ter)*dum*xnp

     if(abs(w3(ns,m)).gt.1.0.and.nonlin(ns)==1)then
        w3(ns,m) = 0.
        w2(ns,m) = 0.
     end if


     laps=anint((z3(ns,m)/lz)-.5)*(1-peritr)
     r=x3(ns,m)-0.5*lx+lr0
     i = int((r-rin)/dr)
     i = min(i,nr-1)
     i = max(i,0)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     qr = wx0*sf(i)+wx1*sf(i+1)
     y3(ns,m)=mod(y3(ns,m)-laps*2*pi*qr*lr0/q0*sign(1.0,q0)+8000.*ly,ly)
     if(x3(ns,m)>lx.and.iperidf==0)then
        x3(ns,m) = lx-1.e-8
        z3(ns,m)=lz-z3(ns,m)
        x2(ns,m) = x3(ns,m)
        z2(ns,m) = z3(ns,m)
        w2(ns,m) = 0.
        w3(ns,m) = 0.
     end if
     if(x3(ns,m)<0..and.iperidf==0)then
        x3(ns,m) = 1.e-8
        z3(ns,m)=lz-z3(ns,m)
        x2(ns,m) = x3(ns,m)
        z2(ns,m) = z3(ns,m)
        w2(ns,m) = 0.
        w3(ns,m) = 0.
     end if
     z3(ns,m)=mod(z3(ns,m)+8.*lz,lz)
     x3(ns,m)=mod(x3(ns,m)+800.*lx,lx)         
     x3(ns,m) = min(x3(ns,m),lx-1.0e-8)
     y3(ns,m) = min(y3(ns,m),ly-1.0e-8)
     z3(ns,m) = min(z3(ns,m),lz-1.0e-8)

     !     particle diagnostics done here because info is available...
     k = int(x3(ns,m)/(lx/nsubd))
     k = min(k,nsubd-1)
     k = k+1
     mypfl_es(k)=mypfl_es(k) + w3old*(eyp) 
     mypfl_em(k)=mypfl_em(k) + w3old*(vpar*delbxp/b) 
     myefl_es(k)=myefl_es(k) + vfac*w3old*(eyp)
     myefl_em(k)=myefl_em(k) + vfac*w3old*(vpar*delbxp/b)
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
  sbuf(2)=myefl_es(nsubd/2)
  sbuf(3)=mypfl_es(nsubd/2)
  sbuf(4)=mynos
  sbuf(5)=myavewi
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,           &
       MPI_COMM_WORLD,ierr)

  ketemp=rbuf(1)
  efltemp=rbuf(2)
  pfltemp=rbuf(3)
  nostemp=rbuf(4)
  avewi(ns,n) = rbuf(5)/( float(tmm(1)) )
  nos(1,n)=nostemp/( float(tmm(1)) )
  ke(1,n)=ketemp/( 2.*float(tmm(1))*mims(ns) )

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  sbuf(1:nsubd) = myefl_es(1:nsubd)
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,  &
       MPI_COMM_WORLD,ierr)
  do k = 1,nsubd
     efl_es(ns,k,n)=rbuf(k)/( float(tmm(1)) )*totvol/vol(k)*cn0s(ns)
  end do

  sbuf(1:nsubd) = myefl_em(1:nsubd)
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,  &
       MPI_COMM_WORLD,ierr)
  do k = 1,nsubd
     efl_em(ns,k,n)=rbuf(k)/( float(tmm(1)) )*totvol/vol(k)*cn0s(ns)
  end do

  sbuf(1:nsubd) = mypfl_es(1:nsubd)
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,  &
       MPI_COMM_WORLD,ierr)
  do k = 1,nsubd
     pfl_es(ns,k,n)=rbuf(k)/( float(tmm(1)) )*totvol/vol(k)*cn0s(ns)
  end do

  sbuf(1:nsubd) = mypfl_em(1:nsubd)
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,  &
       MPI_COMM_WORLD,ierr)
  do k = 1,nsubd
     pfl_em(ns,k,n)=rbuf(k)/( float(tmm(1)) )*totvol/vol(k)*cn0s(ns)
  end do

  !      pfl(1,n)=pfltemp/( float(tmm(1)) )
  !      efl(1,n)=mims(ns)/tets(1)*efltemp/( float(tmm(1)) )

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
  call pmove(u0i(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call end_pmove(ierr)
  mm(ns)=np_new
  !     write(*,*)MyId,mm(ns)

  !      return
end subroutine cpush
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine grad(ip)

  !  currently set up for periodic in x,y,z

  use gem_com
  use gem_equil
  implicit none
  INTEGER :: i,j,k,ip
  real :: ux(0:imx,0:jmx,0:1),uy(0:imx,0:jmx,0:1)
  real :: tmp(0:imx,0:jmx,0:1)

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

  !      return
end subroutine grad

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine grid1(ip,n)

  !    source quantities are are calculated: n_i
  !    right now only ion quantitities are calculated...

  use gem_com
  use gem_equil
  implicit none
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
  real :: enerb,vxdum,dum,xdot,ydot
  INTEGER :: m,n,i,j,k,l,ns,ip
  real :: wx0,wx1,wy0,wy1,wz0,wz1,vte
  real :: sz,wght,r,th,cost,sint,b,qr,dv
  real :: xt,yt,rhog,pidum,vpar,xs,dely,vfac
  real :: lbfs(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
  real :: lbfr(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: myden(0:imx,0:jmx,0:1),myjpar(0:imx,0:jmx,0:1)
  real :: mydene(0:imx,0:jmx,0:1),myupar(0:imx,0:jmx,0:1)
  real :: mydti(0:imx,0:jmx,0:1),mydte(0:imx,0:jmx,0:1)
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
  real :: grp,gxdgyp,rhox(4),rhoy(4)

  rho=0.
  jion = 0.
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
        do l=1,lr(1)
           xs=x3(ns,m)+rhox(l) !rwx(1,l)*rhog
           yt=y3(ns,m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
           xt=mod(xs+800.*lx,lx)
           yt=mod(yt+800.*ly,ly)
           xt = min(xt,lx-1.e-8)
           yt = min(yt,ly-1.e-8)

           i=int(xt/dx+0.5)
           j=int(yt/dy+0.5)
           k=int(z3(ns,m)/dz+0.5)-gclr*kcnt

           myden(i,j,k) = myden(i,j,k) + wght
           myjpar(i,j,k) = myjpar(i,j,k)+wght*vpar
           mydti(i,j,k) = mydti(i,j,k)+wght*vfac

        enddo
     enddo
     if(idg.eq.1)write(*,*)myid,'pass ion grid1'
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     !   enforce periodicity
     call enforce(myden(:,:,:))
     call enforce(myjpar)
     call enforce(mydti)
     !      call filter(myden(:,:,:))
     !      call filter(myjpar(:,:,:))

     do i=0,im
        do j=0,jm
           do k=0,mykm
              den(ns,i,j,k)=q(ns)*myden(i,j,k)/n0/jac(i,k)*cn0s(ns)
              jpar(ns,i,j,k) = q(ns)*myjpar(i,j,k)/n0/jac(i,k)*cn0s(ns)
              mydti(i,j,k) = mydti(i,j,k)/n0/jac(i,k)*cn0s(ns)
           enddo
        enddo
     enddo

     do i = 0,im
        do j = 0,jm
           do k = 0,1
              mydti(i,j,k) = (mydti(i,j,k)-gt0i(i)*den(ns,i,j,k))/gn0s(1,i)
           end do
        end do
     end do
     call MPI_ALLREDUCE(den(ns,0:im,0:jm,0:1),  &
          dti(ns,0:im,0:jm,0:1),             &
          (imx+1)*(jmx+1)*2,MPI_REAL8,       &
          MPI_SUM,GRID_COMM,ierr)


     do i=0,im
        do j=0,jm
           do k=0,mykm
              rho(i,j,k)=rho(i,j,k)+den(ns,i,j,k)
              jion(i,j,k) = jion(i,j,k)+jpar(ns,i,j,k)
           enddo
        enddo
     enddo
  end do

  if(ifluid==0)return
  ! electrons density and current
  vte = sqrt(amie*t0e(nr/2))
  dene(:,:,:) = 0.
  upar(:,:,:) = 0.
  mydene = 0.
  myupar = 0.
  if(idg.eq.1)write(*,*)'enter electron grid1'
  do m=1,mme
     dv=(dx*dy*dz)
     wght=w3e(m)/dv
     vpar = u3e(m) !linearly correct
     !         if(abs(vpar/vte).gt.vcut)wght = 0.

     xt=x3e(m)
     yt=y3e(m)

     i=int(xt/dx)
     j=int(yt/dy)
     k=0 !int(z3e(m)/dz)-gclr*kcnt

     wx0=float(i+1)-xt/dx 
     wx1=1.-wx0
     wy0=float(j+1)-yt/dy
     wy1=1.-wy0
     wz0=float(gclr*kcnt+k+1)-z3e(m)/dz
     wz1=1.-wz0

     mydene(i,j,k)      =mydene(i,j,k)+wght*w000(m)
     mydene(i+1,j,k)    =mydene(i+1,j,k)+wght*w100(m)
     mydene(i,j+1,k)    =mydene(i,j+1,k)+wght*w010(m)
     mydene(i+1,j+1,k)  =mydene(i+1,j+1,k)+wght*w110(m)
     mydene(i,j,k+1)    =mydene(i,j,k+1)+wght*w001(m)
     mydene(i+1,j,k+1)  =mydene(i+1,j,k+1)+wght*w101(m)
     mydene(i,j+1,k+1)  =mydene(i,j+1,k+1)+wght*w011(m)
     mydene(i+1,j+1,k+1)=mydene(i+1,j+1,k+1)+wght*w111(m)

     myupar(i,j,k)      =myupar(i,j,k)+wght*vpar*w000(m)
     myupar(i+1,j,k)    =myupar(i+1,j,k)+wght*vpar*w100(m)
     myupar(i,j+1,k)    =myupar(i,j+1,k)+wght*vpar*w010(m)
     myupar(i+1,j+1,k)  =myupar(i+1,j+1,k)+wght*vpar*w110(m)
     myupar(i,j,k+1)    =myupar(i,j,k+1)+wght*vpar*w001(m)
     myupar(i+1,j,k+1)  =myupar(i+1,j,k+1)+wght*vpar*w101(m)
     myupar(i,j+1,k+1)  =myupar(i,j+1,k+1)+wght*vpar*w011(m)
     myupar(i+1,j+1,k+1)=myupar(i+1,j+1,k+1)+wght*vpar*w111(m)

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
           dene(i,j,k)= mydene(i,j,k)/n0e/jac(i,k)*cn0e*ifluid
           upar(i,j,k) = myupar(i,j,k)/n0e/jac(i,k)*cn0e*ifluid
        end do
     end do
  end do
999 continue

  call MPI_ALLREDUCE(dene(0:im,0:jm,0:1),  &
       delte(0:im,0:jm,0:1),             &
       (imx+1)*(jmx+1)*2,MPI_REAL8,       &
       MPI_SUM,GRID_COMM,ierr)
  call MPI_ALLREDUCE(upar(0:im,0:jm,0:1),  &
       upart(0:im,0:jm,0:1),             &
       (imx+1)*(jmx+1)*2,MPI_REAL8,       &
       MPI_SUM,GRID_COMM,ierr)

  do i = 0,im
     do j = 0,jm
        do k = 0,mykm
           delte(i,j,k) = delte(i,j,k) + phi(i,j,k)*gn0e(i)/gt0e(i)
        enddo
     enddo
  enddo

  do i = 0,im
     do j = 0,jm
        do k = 0,mykm
           rho(i,j,k) = ision*rho(i,j,k) + dene(i,j,k)*qel
        enddo
     enddo
  enddo

  !      return
end subroutine grid1

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine modes2(u,modehis,n)

  !     calculate mode histories. calculate modes for u at timestep n
  !     and store in modehis(mode,n).

  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none
  !     
  real :: u(0:imx,0:jmx,0:1)
  COMPLEX :: modehis(modemx,0:nmx)
  COMPLEX :: modebuf
  INTEGER :: n,i,j,k,l,m,ii

  INTEGER :: mode,jj,thek,oproc,ifirst

  !     

  if(n.eq.0) return

  do mode=1,modem
     oproc=int(nmode(mode)/kcnt*ntube)

     if (MyId.eq.oproc) then
        thek=0
        do j=0,jm-1
           do i=0,im-1
              tmpx(i)=u(i,j,thek)
           enddo

           !     FT in x....
           call ccfft('x',1,imx,1.0,tmpx,coefx,workx,0)
           ii=lmode(mode) !+1
           if(lmode(mode).lt.0) write(*,*) 'lmode < 0, error'
           tmpy(j)=tmpx(ii)/float(im)
        enddo

        !     FT in y....
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
        jj=mmode(mode)  !+1
        if(mmode(mode).lt.0) write(*,*) 'mmode < 0, error'
        modebuf=tmpy(jj)/float(jm)

     endif

     call MPI_BCAST(modebuf,1,MPI_DOUBLE_COMPLEX,oproc, &
          MPI_COMM_WORLD,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     !         write(*,*)myid,modebuf
     modehis(mode,n)=modebuf
  enddo
  !     
  !      return
end subroutine modes2

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine restart(iflag,n)

  use gem_com
  use gem_equil
  implicit none
  INTEGER :: m,ns,iflag,n,i,j,k,ip
  character(len=70) fname
  character(len=5) holdmyid

  write(holdmyid,'(I5.5)') MyId
  fname=directory//'dump_'//holdmyid//'.b'
  if(iflag.eq.1) then
     open(139+MyId,file=fname,form='unformatted',status='old')
     read(139+MyId)ncurr
     read(139+MyId)tcurr,rmpp,rmaa
     do ns = 1,nsm
        read(139+MyId)mm(ns)
        do m=1,mm(ns)
           read(139+MyId) mu(ns,m)
           read(139+MyId) x2(ns,m),y2(ns,m),z2(ns,m),u2(ns,m),w2(ns,m)
           read(139+MyId) xii(ns,m),z0i(ns,m),pzi(ns,m),eki(ns,m),u0i(ns,m)
           w2(ns,m)=w2(ns,m)/cut
           x3(ns,m)=x2(ns,m)
           y3(ns,m)=y2(ns,m)
           z3(ns,m)=z2(ns,m)
           u3(ns,m)=u2(ns,m)
           w3(ns,m)=w2(ns,m)
        enddo
     end do
120  continue

     read(139+MyId)mme
     do  m=1,mme
        read(139+MyId) x2e(m),y2e(m),z2e(m),u2e(m),w2e(m),mue2(m),ipass(m)
        read(139+MyId) xie(m),z0e(m),pze(m),eke(m),mue(m),u0e(m)
        w2e(m)=w2e(m)/cut
        x3e(m)=x2e(m)
        y3e(m)=y2e(m)
        z3e(m)=z2e(m)
        u3e(m)=u2e(m)
        w3e(m)=w2e(m)
        mue3(m)=mue2(m)
     end do

     do i = 0,6
        do j = 0,50000
           read(139+myid)camp(i,j)
        end do
     end do
     close(139+MyId)
  endif

  if(iflag.eq.2) then
     open(139+MyId,file=fname,form='unformatted',status='unknown')
     write(139+MyId)n+1
     write(139+MyId)tcurr-dt,rmpp,rmaa
     do ns = 1,nsm
        write(139+MyId)mm(ns)
        do m=1,mm(ns)
           write(139+MyId) mu(ns,m)
           write(139+MyId) x2(ns,m),y2(ns,m),z2(ns,m),u2(ns,m),w2(ns,m)
           write(139+MyId) xii(ns,m),z0i(ns,m),pzi(ns,m),eki(ns,m),u0i(ns,m)
        enddo
     end do

     write(139+MyId)mme
     do  m=1,mme
        write(139+MyId) x2e(m),y2e(m),z2e(m),u2e(m),w2e(m),mue2(m),ipass(m)
        write(139+MyId) xie(m),z0e(m),pze(m),eke(m),mue(m),u0e(m)
     end do

     do i = 0,6
        do j = 0,50000
           write(139+myid)camp(i,j)
        end do
     end do
     close(139+MyId)
  endif

end subroutine restart
!-----------------------------------------------------------------------

!        Normal distribution random no. generator, stand. dev. = 1.
!        Version 2 does it Willy's way...

subroutine parperp(vpar,vperp2,m,pi,cnt,MyId)

  INTERFACE
     real function revers(num,n)
       integer :: num,n
     end function revers
  END INTERFACE

  real :: vpar,vperp2,r1,r2,t,pi
  INTEGER :: m,iflag,cnt,MyId
  real :: c0,c1,c2
  real :: d1,d2,d3
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
110 continue
  if(r1.ge.1.e-6) then
     t=sqrt(log(1.0/(r1*r1)))
  else
     t=5.0
     write(*,*)'parperp2 warning  m= ',m
  endif

  vpar=t-(c0+c1*t+c2*t**2)/(1.+d1*t+d2*t**2+d3*t**3)
  vpar=vpar*iflag

  vperp2=-2.0*log(r2)

  !        return
end subroutine parperp

!---------------------------------------------------------------------- 
!    calculate weights and delta j for proper periodicity 
!    in the toroidal direction, weights and delta j are
!    dependant only on minor radius, hence x, hence
!    weight is a vector in 0:imx

subroutine weight

  use gem_com
  use gem_equil
  implicit none
  INTEGER :: i,j
  real :: dely,aweight,r,qr,wx0,wx1

  !         peritr = 0
  if (GCLR.eq.Master.and.peritr.eq.0) then
     do i=0,im
        r=float(i)*dx-0.5*lx+lr0
        j = int((r-rin)/dr)
        j = min(j,nr-1)
        wx0 = (rin+(j+1)*dr-r)/dr
        wx1 = 1.-wx0
        qr = wx0*sf(j)+wx1*sf(j+1)
        dely=mod(2.*pi*lr0/q0*qr*sign(1.0,q0)+800.0*ly,ly)
        deljp(i)=int(dely/dy)
        deljm(i)=0
        aweight=mod(dely,dy)
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
        dely=mod(2.*pi*lr0/q0*qr*sign(1.0,q0)+800.*ly,ly)
        deljm(i)=int(dely/dy)
        deljp(i)=0
        aweight=mod(dely,dy)
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

  !      return
end subroutine weight

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine gkps(nstep,ip)   

  use gem_com
  use gem_equil
  use gem_fft_wrapper

  implicit none

  INTEGER :: ns
  real :: b,gam0,gam1,delyz,th,bf,dum,r,qr,shat
  real,DIMENSION(1:5) :: b1,b2,ter
  real :: kx1,kx2,ky
  real,dimension(:),allocatable :: akx,aky
  real,dimension(:,:,:,:,:),allocatable:: gamb1,gamb2
  complex,dimension(:,:,:,:),allocatable :: mx
  integer,dimension(:,:,:,:),allocatable :: ipiv
  real,dimension(:,:,:),allocatable :: formphi,formfe
  real :: sgnx,sgny,sz,myfe
  INTEGER :: i,i1,j,j1,k,k1,l,m,n,ifirst,nstep,ip,INFO
  INTEGER :: l1,m1,myk,myj,ix,ikx
  COMPLEX :: temp3dxy(0:imx-1,0:jmx-1,0:1),v(0:imx-1,0:jcnt-1,0:1)
  COMPLEX :: sbuf(0:imx*jcnt*2-1),rbuf(0:imx*jmx*2-1)
  COMPLEX :: sl(1:imx-1,0:jcnt-1,0:1)
  COMPLEX :: aphik(0:nxpp-1),myaphik(0:nxpp-1)
  real :: myaph(0:nxpp),aph(0:nxpp),u(0:imx,0:jmx,0:1)
  complex :: cdum
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
  real :: grp,gxdgyp,grdgtp,gthp
  real :: wx0,wx1,wz0,wz1
  character(len=70) fname
  character(len=5) holdmyid

  save formphi,formfe,ifirst,akx,aky,mx,gamb1,gamb2,ipiv
  if(idg==1)write(*,*)'enter gkps'

  write(holdmyid,'(I5.5)') MyId
  fname='./matrix/'//'mx_phi_'//holdmyid
  !     form factors....
  if (ifirst.ne.-99) then
     allocate(akx(0:imx-1),aky(0:jcnt-1), &
          gamb1(1:5,0:imx-1,0:jcnt-1,0:imx-1,0:1), &
          gamb2(1:5,0:imx-1,0:jcnt-1,0:imx-1,0:1))
     allocate(mx(imx-1,imx-1,0:jcnt-1,0:1),formphi(0:imx-1,0:jcnt-1,0:1))
     allocate(formfe(0:imx-1,0:jcnt-1,0:1),ipiv(imx-1,imx-1,0:jcnt-1,0:1))

     if(iget.eq.1) then
        open(10000+MyId,file=fname,form='unformatted',status='old')
        read(10000+MyId)mx,ipiv
        close(10000+myid)
        goto 200
     end if

     do l=0,im-1
        do m=0,jcnt-1
           j = tclr*jcnt+m
           do i=0,im-1 
              r = lr0-lx/2+i*dx
              i1 = int((r-rin)/dr)
              i1 = min(i1,nr-1)
              wx0 = (rin+(i1+1)*dr-r)/dr
              wx1 = 1.-wx0
              do ns = 1, nsm
                 ter(ns) = wx0*t0s(ns,i1)+wx1*t0s(ns,i1+1)
              enddo
              do n = 0,1
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

                 do ns = 1, nsm
                    b1(ns)=mims(ns)*(kx1*kx1*grp**2 + &
                         ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                         +2*dydrp*lr0/q0*qhatp*grdgtp) &
                         +2*kx1*ky*gxdgyp)/(bf*bf)*ter(ns)/(q(ns)*q(ns))

                    b2(ns)=mims(ns)*(kx2*kx2*grp**2 + &
                         ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                         +2*dydrp*lr0/q0*qhatp*grdgtp) &
                         +2*kx2*ky*gxdgyp)/(bf*bf)*ter(ns)/(q(ns)*q(ns))

                    call srcbes(b1(ns),gam0,gam1)
                    gamb1(ns,l,m,i,n)=gam0
                    call srcbes(b2(ns),gam0,gam1)
                    gamb2(ns,l,m,i,n)=gam0
                 enddo

                 !   formfactor in gkps
                 formfe(l,m,n) = 1.-gam0
                 formphi(l,m,n) = 1./jmx 
                 if(abs(ky)>kycut)formphi(l,m,n) = 0.
                 !                    if(abs(ky)==0.)formphi(l,m,n) = 0.
                 if(m1.ne.mlk.and.onemd==1)formphi(l,m,n) = 0.
              enddo
           enddo
        enddo
     enddo

     do k = 0,1
        do j = 0,jcnt-1
           do i = 1,imx-1
              r = lr0-lx/2+i*dx
              i1 = int((r-rin)/dr)
              i1 = min(i1,nr-1)
              wx0 = (rin+(i1+1)*dr-r)/dr
              wx1 = 1.-wx0
              do ns = 1, nsm
                 ter(ns) = wx0*t0s(ns,i1)+wx1*t0s(ns,i1+1)
              enddo
              do ix = 1,imx-1
                 mx(i,ix,j,k) = 0.0
                 if(i==ix)mx(i,ix,j,k) = fradi*gn0e(i)*cn0e/gt0e(i)
                 do ikx = 0,imx-1
                    do ns = 1, nsm
                       mx(i,ix,j,k) = mx(i,ix,j,k)+q(ns)*sin(ix*ikx*pi/imx)* &
                            ((1-gamb1(ns,ikx,j,i,k))*exp(IU*ikx*i*pi/imx)- &
                            (1-gamb2(ns,ikx,j,i,k))*exp(-IU*ikx*i*pi/imx)) &
                            /ter(ns)*cn0s(ns)*gn0s(ns,i)/(IU*imx)
                    end do
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
        close(10000+myid)
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
200  ifirst=-99
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
           u(i,j,k) = rho(i,j,k)+den0apa(i,j,k)+fradi*(phi(i,j,k)/ntube*gn0e(i)*cn0e/gt0e(i)-den0(i,j,k))
        end do
     end do
  end do
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
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3dxy(i,j,k) = tmpy(j)  ! phi(x,y)
        end do
     end do
  end do

100 continue
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

  !      return
end subroutine gkps

!      End of gkps....

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine eqmo(ip)

  use gem_com
  use gem_equil

  implicit none
  real :: lbfr(0:imx,0:jmx)
  real :: lbfs(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
  integer :: i,j,k,ip
  real :: eta

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

  !      return
end subroutine eqmo
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine spec(n)
  use gem_com
  use gem_equil
  implicit none
  integer :: i,j,k,l,m,n
  real :: pf,efe,efi,pfi,efc,pfc,efb,pfb,pflxgb,eflxgb,x
  real :: pf_em,efe_em,efi_em,pfi_em,efc_em,pfc_em,efb_em,pfb_em
  real :: tdum,v(0:imx)

!  call zon(upart+amie*upa0t,v)
!  call joule

  eflxgb = xn0e(nr2)*cn0e*t0e(nr2)*sqrt(t0e(nr2)/mimp)*rhoia**2
  pflxgb = eflxgb/t0e(nr2)
  i = tcurr-dt
  pf = 0.
  efe = 0.
  efi = 0.
  pfi = 0.
  efc = 0.
  pfc = 0.
  efb = 0.
  pfb = 0.

  pf_em = 0.
  efe_em = 0.
  efi_em = 0.
  pfi_em = 0.
  efc_em = 0.
  pfc_em = 0.
  efb_em = 0.
  pfb_em = 0.
  k = 2
  x = float(nsubd)/float(nsubd-2*k)
  do j = 1+k,nsubd-k
     pf = pf+pfle_es(j,n)*vol(j)/totvol*x
     efe = efe+efle_es(j,n)*vol(j)/totvol*x
     efi = efi+efl_es(1,j,n)*vol(j)/totvol*x
     pfi = pfi+pfl_es(1,j,n)*vol(j)/totvol*x
     efc = efc+efl_es(2,j,n)*vol(j)/totvol*x
     pfc = pfc+pfl_es(2,j,n)*vol(j)/totvol*x
     efb = efb+efl_es(3,j,n)*vol(j)/totvol*x
     pfb = pfb+pfl_es(3,j,n)*vol(j)/totvol*x
  end do

  do j = 1+k,nsubd-k
     pf_em = pf_em+pfle_em(j,n)*vol(j)/totvol*x
     efe_em = efe_em+efle_em(j,n)*vol(j)/totvol*x
     efi_em = efi_em+efl_em(1,j,n)*vol(j)/totvol*x
     pfi_em = pfi_em+pfl_em(1,j,n)*vol(j)/totvol*x
     efc_em = efc_em+efl_em(2,j,n)*vol(j)/totvol*x
     pfc_em = pfc_em+pfl_em(2,j,n)*vol(j)/totvol*x
     efb_em = efb_em+efl_em(3,j,n)*vol(j)/totvol*x
     pfb_em = pfb_em+pfl_em(3,j,n)*vol(j)/totvol*x
  end do

  tdum = tcurr-dt
  if(myid.eq.master)then
     open(9, file='plot', status='unknown',position='append')
     open(11, file='flux', status='unknown',position='append')
     open(17, file='yyre', status='unknown',position='append')

     write(*,10)i,rmsphi(n),rmsapa(n),pf,efe,pfi,efi,avewi(1,n),&
          avewe(n),yyre(1,0),yyim(1,0),yyamp(1,0)

10   format(1x,i6,12(2x,e10.3))
11   format(6x,5(2x,e12.5))
12   format(1x,i6,12(2x,e12.5))
13   format(1x,i6,4(2x,e10.3),2x,i7,2x,i7)
15   format(1x,i6,8(2x,e10.3))

     write(9,10)i,rmsphi(n),rmsapa(n),pf,efe,pfi,efi,avewi(1,n),&
          avewe(n),yyre(1,0),yyim(1,0),yyamp(1,0)

     write(11,12)i,pf/pflxgb,pfi/pflxgb,pfc/pflxgb,efe/eflxgb,efi/eflxgb,&
          efc/eflxgb,pf_em/pflxgb,pfi_em/pflxgb,pfc_em/pflxgb,efe_em/eflxgb,&
          efi_em/eflxgb,efc_em/eflxgb

     write(17,12)i,yyre(1,0),yyre(1,1),yyre(1,2),yyre(1,3),yyre(1,4)
     close(9)
     close(11)
     close(17)
  end if

  return
  if(gclr==kmx/2 .and. tclr==0)then
     open(22, file='yyre2', status='unknown',position='append')
     open(23, file='mdhis', status='unknown',position='append')
     open(24, file='mdhisa', status='unknown',position='append')
     open(25, file='stress', status='unknown',position='append')

     write(23,16)tdum,mdhis(0),mdhis(1),mdhis(2),mdhis(3),mdhis(4),&
          mdhisa(0),mdhisa(1),mdhisa(2),mdhisa(3),mdhisa(4)
     write(24,16)tdum,mdhisb(0),mdhisb(1),mdhisc(0),mdhisc(1),&
          mdhisd(0),mdhisd(1)
     write(25,17)tdum,(v(i),i = 0,imx-1)

     do  i = 0,6
        write(22,14)tdum,i,real(phihis(i,0)),(real(phihis(i,j)), &
             aimag(phihis(i,j)), j = 1,jcnt-2,2)
        write(22,14)tdum,i,real(aparhis(i,0)),(real(aparhis(i,j)), &
             aimag(aparhis(i,j)), j = 1,jcnt-2,2)
     end do
     close(22)
     close(23)
     close(24)
     close(25)
14   format(1x,f10.1,1x,i2,10(2x,e12.5))
16   format(1x,f10.1,1x,10(2x,e12.5))
17   format(1x,f10.1,256(2x,e12.5))
  end if

end subroutine spec

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine ezamp(nstep,ip)   

  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none
  real :: b,b2,gam0,gam1,th,delyz,bf,rkper,r,qr,shat
  real :: kx,ky
  real,dimension(:),allocatable :: akx,aky
  complex,dimension(:,:,:,:),allocatable:: nab1,nab2
  complex,dimension(:,:,:,:),allocatable :: mx
  real,dimension(:,:,:),allocatable :: formapa
  real :: sgnx,sgny,sz,myfe
  integer,dimension(:,:,:,:),allocatable :: ipiv
  INTEGER :: i,i1,j,j1,k,k1,l,m,n,ifirst,nstep,ip,INFO
  INTEGER :: l1,m1,myk,myj,ix,ikx,iext
  COMPLEX :: temp3dxy(0:imx-1,0:jmx-1,0:1),v(0:imx-1,0:jcnt-1,0:1)
  COMPLEX :: sbuf(0:imx*jcnt*2-1),rbuf(0:imx*jmx*2-1)
  COMPLEX :: sl(1:imx-1,0:jcnt-1,0:1)
  real :: dum,u(0:imx,0:jmx,0:1)
  real :: grp,gthp,gxdgyp,dydrp,qhatp,grdgtp,bfldp
  real :: wx0,wx1,wz0,wz1
  character(len=70) fname
  character(len=5) holdmyid

  save formapa,ifirst,akx,aky,mx,IPIV

  write(holdmyid,'(I5.5)') MyId
  fname='./matrix/'//'mx_apa_'//holdmyid
  if (ifirst.ne.-99) then
     allocate(akx(0:imx-1),aky(0:jcnt-1),nab1(0:imx-1,0:jcnt-1,0:imx-1,0:1))
     allocate(nab2(0:imx-1,0:jcnt-1,0:imx-1,0:1),ipiv(imx-1,imx-1,0:jcnt-1,0:1))
     allocate(mx(imx-1,imx-1,0:jcnt-1,0:1),formapa(0:imx-1,0:jcnt-1,0:1))

     if(iget.eq.1) then
        open(20000+MyId,file=fname,form='unformatted',status='old')
        read(20000+MyId)mx,ipiv
        close(20000+myid)
        goto 200
     end if

     do l=0,im-1
        do m=0,jcnt-1
           j = tclr*jcnt+m
           do i=0,im-1
              r = lr0-lx/2+i*dx
              do n = 0,1
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

                 m1 = mstart+int((float(m)+1.0)/2)   
                 if(m==0)m1=0       
                 sgny = isgnft(m)
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
              enddo
           enddo
        enddo
     enddo

     do k = 0,1
        do j = 0,jcnt-1
           do i = 1,imx-1
              do ix = 1,imx-1
                 mx(i,ix,j,k) = 0.
                 if(i==ix)mx(i,ix,j,k) = beta*(amie*gn0e(i)*cn0e+q(1)*q(1)/mims(1)*gn0s(1,i)*cn0s(1)*isiap)
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
        close(20000+myid)
        goto 200
     end if

200  ifirst=-99
  endif

  !   now do field solve...

  !  find jtot(kx,ky)
  iext = 1
  if(nstep>200)iext = 0
  do i = 0,imx
     do j = 0,jmx
        do k = 0,1
           u(i,j,k) = jion(i,j,k)*ision-upar(i,j,k) &
                +amie*(apar(i,j,k)/ntube*gn0e(i)*cn0e-upa00(i,j,k)) &
                +0.0*sin(mlk*yg(j)-0.01*tcurr)*exp(-(xg(i)-lx/2)**2/(lx/3)**2) &
                *exp(-(zg(k)-lz/2)**2/(lz/2)**2)*iext
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
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
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

  !      return
end subroutine ezamp

!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine filter(u)

  use gem_com
  use gem_equil
  implicit none
  integer :: i,j,k,l,m,n,ip,NMDZ,NNI,ikz,ndum
  parameter(NMDZ=2)
  real :: u(0:imx,0:jmx,0:1)
  real :: temp(0:imx,0:jmx,0:1)
  real :: cc0(0:imx,0:jmx,0:NMDZ-1),ss0(0:imx,0:jmx,0:NMDZ-1)
  real :: cc1(0:imx,0:jmx,0:NMDZ-1),ss1(0:imx,0:jmx,0:NMDZ-1)
  real ::  rkz,DW1,DW2
  parameter(NNI = 3)
  parameter(DW1 = 0.5)
  parameter(DW2 = -1./6.)
  real :: ttemp3d(0:1,0:imx-1,0:jmx-1)
  real :: htemp3d(0:1,0:imx-1,0:jmx-1)
  real :: lbfs(0:imx-1,0:jmx-1) 
  real :: lbfr(0:imx-1,0:jmx-1)
  real :: rbfs(0:imx-1,0:jmx-1)
  real :: rbfr(0:imx-1,0:jmx-1)
  real :: tmpbc(0:imx,0:jmx)
  real :: uz0(0:imx,0:jmx),uz1(0:imx,0:jmx)
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
        enddo
     enddo
  enddo

  do n=1,NNI

     htemp3d=ttemp3d

     !  build buffers and send grid info to neighbors

     do j=0,jm-1
        do i=0,im-1
           rbfs(i,j)=htemp3d(mykm-1,i,j)
           lbfs(i,j)=htemp3d(0,i,j)
        enddo
     enddo

     call MPI_SENDRECV(rbfs(0:imx-1,0:jmx-1),imx*jmx, &
          MPI_REAL8,rngbr,9, &
          lbfr(0:imx-1,0:jmx-1),imx*jmx, &
          MPI_REAL8,lngbr,9,  &
          TUBE_COMM,stat,ierr)

     call MPI_SENDRECV(lbfs(0:imx-1,0:jmx-1),imx*jmx, &
          MPI_REAL8,lngbr,8,  &
          rbfr(0:imx-1,0:jmx-1),imx*jmx,  &
          MPI_REAL8,rngbr,8,  &
          TUBE_COMM,stat,ierr)
     !	 write(*,*)'past buff'

     do j=0,jm-1
        do i=0,im-1
           ttemp3d(0,i,j)=( htemp3d(0,i,j)   &
                + DW1*(  &
                + weightp(i)*lbfr(i,jpl(i,j)) &
                + weightpn(i)*lbfr(i,jpn(i,j)) &
                + weightm(i)*rbfr(i,jmi(i,j)) &
                + weightmn(i)*rbfr(i,jmn(i,j)) ) &
                ) &
                / (1.0+2*DW1)
        enddo
     enddo

     htemp3d=ttemp3d

     !  build buffers and send grid info to neighbors

     do j=0,jm-1
        do i=0,im-1
           rbfs(i,j)=htemp3d(mykm-1,i,j)
           lbfs(i,j)=htemp3d(0,i,j)
        enddo
     enddo

     call MPI_SENDRECV(rbfs(0:imx-1,0:jmx-1),imx*jmx, &
          MPI_REAL8,rngbr,9, &
          lbfr(0:imx-1,0:jmx-1),imx*jmx, &
          MPI_REAL8,lngbr,9, &
          TUBE_COMM,stat,ierr) 

     call MPI_SENDRECV(lbfs(0:imx-1,0:jmx-1),imx*jmx, &
          MPI_REAL8,lngbr,8,  &
          rbfr(0:imx-1,0:jmx-1),imx*jmx, &
          MPI_REAL8,rngbr,8, &
          TUBE_COMM,stat,ierr)
     !	 write(*,*)'past buff'

     do j=0,jm-1
        do i=0,im-1
           ttemp3d(0,i,j)=(  htemp3d(0,i,j) &
                + DW2*( nflag*htemp3d(1,i,j) &
                + weightp(i)*lbfr(i,jpl(i,j)) &
                + weightpn(i)*lbfr(i,jpn(i,j)) &
                + weightm(i)*rbfr(i,jmi(i,j)) &
                + weightmn(i)*rbfr(i,jmn(i,j)) ) &
                ) &
                /(1.0+2*DW2)
        enddo
     enddo

  enddo

  !  assign ttemp3d(mykm,:,:)

  do i = 0,im-1
     do j = 0,jm-1
        lbfs(i,j) = ttemp3d(0,i,j)
     end do
  end do

  call MPI_SENDRECV(lbfs(0:imx-1,0:jmx-1),imx*jmx, &
       MPI_REAL8,lngbr,10, &
       rbfr(0:imx-1,0:jmx-1),imx*jmx, &
       MPI_REAL8,rngbr,10, &
       TUBE_COMM,stat,ierr)

  do i = 0,im-1
     do j = 0,jm-1
        ttemp3d(mykm,i,j) = weightm(i)*rbfr(i,jmi(i,j)) &
             + weightmn(i)*rbfr(i,jmn(i,j)) 
     end do
  end do

  do k=0,mykm
     do j=0,jm-1
        do i=0,im-1
           temp(i,j,k)=ttemp3d(k,i,j)
        enddo
     enddo
  enddo

  call enfxy(temp)

100 continue
  do i = 0,im
     do j = 0,jm
        do k = 0,mykm
           u(i,j,k) = temp(i,j,k)
        end do
     end do
  end do
  return

200 continue

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

  !      return
end subroutine filter

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine yveck(u,n)

  !     calculate mode histories. calculate modes for u at timestep n
  !     and store in modehis(mode,n).

  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none
  !     
  real :: u(0:imx,0:jmx,0:1)
  INTEGER :: n,i,j,k
  complex :: tmp3d(0:imx-1,0:jmx-1,0:1-1)

  call ccfft('z',0,kmx,1.0,tmpz,coefz,workz,0)

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
        call ccfft('x',1,imx,1.0,tmpx,coefx,workx,0)
        do i = 0,imx-1
           tmp3d(i,j,k) = tmpx(i)
        end do
     end do
     do i = 0,imx-1
        do j = 0,jmx-1
           tmpy(j) = tmp3d(i,j,k)
        end do
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
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
        call ccfft('z',1,kmx,1.0,tmpz,coefz,workz,0)
     end if

     call MPI_BCAST(tmpz,kmx,MPI_DOUBLE_COMPLEX,master, &
          tube_comm,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)

     do i = 0,4
        yyamp(j,i) = abs(tmpz(i)) !cabs
        yyre(j,i) = real(tmpz(i))
        yyim(j,i) = aimag(tmpz(i))
     end do
  end do

  !      return
end subroutine yveck
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine yveck1(u,n)

  !     calculate mode histories. calculate modes for u at timestep n
  !     and store in modehis(mode,n).

  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none
  !     
  real :: u(0:imx,0:jmx,0:1)
  INTEGER :: n,i,j,k
  COMPLEX :: tmp3d(0:imx-1,0:jmx-1,0:1-1)

  call ccfft('z',0,kmx,1.0,tmpz,coefz,workz,0)

  !      if(n.eq.0) return

  do j=0,jm-1
     do i=0,im-1
        do k = 0,mykm-1
           tmp3d(i,j,k)=u(i,j,k)
        enddo
     end do
  end do

  do k = 0,mykm-1
     do i = 0,imx-1
        do j = 0,jmx-1
           tmpy(j) = tmp3d(i,j,k)
        end do
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           tmp3d(i,j,k) = tmpy(j)
        end do
     end do
  end do

  j = int(n/ifskp)
  if(mod(n,ifskp)==0)then
     camp(:,j)=0.
     if(j>0)camp(:,j-1)=camp(:,j-1)/ifskp
  end if
  do i = 0,6
     camp(i,j) = camp(i,j)+tmp3d(imx/8*(i+1),mlk,0)
  end do
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

end subroutine yveck1

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

real function ran2(idum)

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
  real :: temp

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
end function ran2

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine loadi(ns)

  use gem_com
  use gem_equil
  implicit none
  INTEGER :: i,j,k,m,idum,ns,m1
  INTEGER :: np_old,np_new
  real :: vpar,vperp2,r,qr,th,b,cost,ter
  real :: avgv,myavgv,avgw,myavgw
  real :: dumx,dumy,dumz,jacp
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,psp
  real :: grp,gxdgyp,zoldp
  real :: wx0,wx1,wz0,wz1

  cnt=int(tmm(1)/numprocs)
  !   write(*,*)'in loader cnt, mm(1)=',cnt,mm(1)

  myavgv=0.
  avgv=0.
  avgw = 0.
  myavgw = 0.

  !      ns = 1
  m = 0
  do j = 1,100000000

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
        x2(ns,m)=min(dumx,lx-1e-8)
        y2(ns,m)=min(dumy,ly-1e-8)

        k = int((th+pi)/dth)
        wz0 = (-pi+(k+1)*dth-th)/dth
        wz1 = 1-wz0
        z2(ns,m) = wz0*zfnth(k)+wz1*zfnth(k+1)
        z2(ns,m)=min(z2(ns,m),lz-1e-8)

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
        if(ildu==1)ter = tgis(ns)
        b=1.-tor+tor*bfldp

        u2(ns,m)=vpar/sqrt(mims(ns)/ter)
        mu(ns,m)=0.5*vperp2/b*ter
        eki(ns,m) = mu(ns,m)*b+0.5*mims(ns)*u2(ns,m)**2
        pzi(ns,m) = mims(ns)*u2(ns,m)/b-q(ns)*psp/br0
        z0i(ns,m) = z2(ns,m)
        xii(ns,m) = x2(ns,m)
        u0i(ns,m) = u2(ns,m)
        myavgv=myavgv+u2(ns,m)

        !    LINEAR: perturb w(ns,m) to get linear growth...
        !         w2(ns,m)=2.*amp*(revers(MyId*cnt+m,13)-0.5) !(ran2(iseed) - 0.5 )
        w2(ns,m)=2.*amp*sin(pi2/ly*y2(ns,m))*exp(-(z2(ns,m)-lz/2)**2/(lz/8)**2)*exp(-(x2(ns,m)-0.4*lx)**2/(lx/8)**2)
        if(ns==2)w2(ns,m) = 0.
        myavgw=myavgw+w2(ns,m)
     end if
  enddo
170 continue
  myavgw = myavgw/mm(ns)
  !      write(*,*)'myid ', myid,x2(10),y2(20),u2(20),w2(20)
  !    subtract off avg. u...
  call MPI_ALLREDUCE(myavgv,avgv,1, &
       MPI_REAL8, &
       MPI_SUM,MPI_COMM_WORLD,ierr)
  if(idg.eq.1)write(*,*)'all reduce'
  avgv=avgv/float(tmm(1))
  do m=1,mm(ns)
     u2(ns,m)=u2(ns,m)-avgv
     x3(ns,m)=x2(ns,m)
     y3(ns,m)=y2(ns,m)
     z3(ns,m)=z2(ns,m)
     u3(ns,m)=u2(ns,m)
     !         w2(ns,m) = w2(ns,m)-myavgw
     w3(ns,m)=w2(ns,m)
  enddo

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
  call pmove(u0i(ns,:),np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  !     
  call end_pmove(ierr)
  mm(ns)=np_new

  !      return
end subroutine loadi
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine enforce(u)
  use gem_com
  use gem_equil
  implicit none
  INTEGER :: m,i,j,k,ns,jj
  real :: u(0:imx,0:jmx,0:1)      
  real :: lbfs(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
  real :: lbfr(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: dum,dum1,dely,th,wy1,ydum

  do j=0,jm-1
     do k=0,mykm
        u(0,j,k) = u(0,j,k)+u(im,j,k)
     enddo
  enddo

  do i=0,im-1 
     do k=0,mykm
        u(i,0,k) = u(i,0,k)+u(i,jm,k)
        u(i,jm,k) = u(i,0,k)
     enddo
  enddo

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
  !      return
end subroutine enforce
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine enfxy(u)
  use gem_com
  use gem_equil
  implicit none
  real :: u(0:imx,0:jmx,0:1)
  integer :: i,j,k,l,m,n,jj
  real :: ydum,th,dely,wy1

  !    periodic bc in y...
  do k=0,mykm
     do i=0,im-1
        u(i,jm,k)=u(i,0,k)
     enddo
  enddo

  !   bc for x
  do k=0,mykm
     do j=0,jm
        u(im,j,k)=u(0,j,k)
     enddo
  enddo

  !      return
end subroutine enfxy
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine gradu(u,ux,uy)
  use gem_com
  use gem_equil
  implicit none
  real :: u(0:imx,0:jmx,0:1)
  real :: ux(0:imx,0:jmx,0:1),uy(0:imx,0:jmx,0:1)
  integer :: i,j,k,l,m,n,jj,ju,jl
  real :: ydum,th,dely,wy1,ul

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
  !      return
end subroutine gradu

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine gradx(u,ux)
  use gem_com
  use gem_equil
  implicit none
  real :: u(0:imx,0:jmx,0:1)
  real :: ux(0:imx,0:jmx,0:1)
  integer :: i,j,k,l,m,n,jj,ju,jl
  real :: ydum,th,dely,wy1,ul

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

end subroutine gradx

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine grady(u,uy)
  use gem_com
  use gem_equil
  implicit none
  real :: u(0:imx,0:jmx,0:1)
  real :: uy(0:imx,0:jmx,0:1)
  integer :: i,j,k,l,m,n,jj,ju,jl
  real :: ydum,th,dely,wy1,ul

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

end subroutine grady

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine enfz(u)
  use gem_com
  use gem_equil
  implicit none
  real :: u(0:imx,0:jmx,0:1)
  real :: lbfr(0:imx,0:jmx)
  real :: lbfs(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
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

  !      return
end subroutine enfz
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine pint
  use gem_com
  use gem_equil
  implicit none
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dgdtp,dpdzp,dadzp,aparp,vncp
  real :: dum,vxdum,dum1,dum2,eps,x,h_x,h_coll,hee,nue
  INTEGER :: m,i,j,k,l,n,ipover,ieover
  INTEGER :: np_old,np_new
  real :: vfac,kap,kapnp,kaptp,sz,vpar,ppar,vpdum,pidum,xnp
  real :: b,th,r,enerb,cost,sint,qr,laps,ter
  real :: xt,xs,yt,zt,xdot,ydot,zdot,pzdot,edot,pzd0,pzd1,vp0
  real :: wx0,wx1,wy0,wy1,wz0,wz1
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp
  real :: grp,gxdgyp,psp,pzp,psip2p,bdcrvbp,curvbzp,dipdrp,bstar
  real :: myavptch,myaven
  integer :: myopz,myoen
  real :: x000,x001,x010,x011,x100,x101,x110,x111

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
     vncp = wx0*phincp(i)+wx1*phincp(i+1)        
     psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
     dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)        
     b=1.-tor+tor*bfldp
     pzp = emass*u2e(m)/b+psp/br0

     xt = x2e(m)
     yt = y2e(m)

     include 'ppushlie.h'

     vfac = 0.5*(emass*u2e(m)**2 + 2.*mue2(m)*b)
     vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
     vp0 = vp0*vncp*vexbsw
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
          (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
          +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
          -emass*vpar**2/(qel*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
          -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*grcgtp*lr0/q0*qhatp  

     zdot =  vpar*b/bstar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
          +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
          -1./b**2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp &
          -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*q0*br0*grcgtp/jfnp

     pzd0 = tor*(-mue2(m)/emass/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar &
          +mue2(m)*vpar/(qel*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp
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

     eps = (b*mue2(m)+0.5*emass*u2e(m)*u2e(m))/ter
     x = sqrt(eps)
     h_x    = 4.0*x*x/(3.0*sqrt(pi))
     h_coll = h_x/sqrt(1.0+h_x**2)
     !         h_coll = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*DERF(x)
     ! collision frequency for experimental profiles
     hee = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*erf(x)
     nue=rneu*(wx0*nue0(i)+wx1*nue0(i+1))*(wx0*zeff(i)+wx1*zeff(i+1)+hee)
     !         dum1 = 0.5*rneu/eps**1.5*(1+h_coll)
     dum1 = nue/(2*(eps+0.1))**1.5  !*(1+h_coll)
     !         if(x<0.3)dum1=0.0

     dum = 1-w2e(m)*nonline*0.
     if(eldu.eq.1)dum = (tge/ter)**1.5*exp(vfac*(1/tge-1./ter))
     vxdum = (eyp/b+vpdum/b*delbxp)*dum2
     w3e(m)=w2e(m) + 0.5*dte*(  &
          (vxdum*kap + edot/ter-dum1*ppar*aparp/ter)*xnp     &
          +isg*(-dgdtp-zdot*dpdzp+xdot*exp1+ydot*eyp)/ter*xnp)*dum

     !         if(x3e(m)>lx .or. x3e(m)<0.)w3e(m) = 0. 

     !         go to 333

     ieover = 0
     ipover = 0
     if(abs(pzp-pze(m))>pzcrite)then
        myopz = myopz+1
        ipover = 1
     end if
     if(abs(vfac-eke(m))>encrit.and.abs(x2e(m)-lx/2)<(lx/2-1))then
        myoen = myoen+1
        ieover = 1
        myaven = myaven+eke(m)
        myavptch = myavptch+abs(vpar)/sqrt(2/emass*vfac)
     end if
     if(itube==1)goto 333
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
        u3e(m) = u0e(m)
        u2e(m) = u3e(m)
        w3e(m) = 0.
        w2e(m) = 0.
        x2e(m) = x3e(m)
        z2e(m) = z3e(m)
        mue3(m) = mue(m)
        mue2(m) = mue(m)
     end if
333  continue
     laps=anint((z3e(m)/lz)-.5)*(1-peritr)
     r=x3e(m)-0.5*lx+lr0
     i = int((r-rin)/dr)
     i = min(i,nr-1)
     i = max(i,0)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     qr = wx0*sf(i)+wx1*sf(i+1)
     y3e(m)=mod(y3e(m)-laps*2*pi*qr*lr0/q0*sign(1.0,q0)+8000.*ly,ly)
     if(x3e(m)>lx.and.iperidf==0)then
        x3e(m) = lx-1.e-8
        z3e(m)=lz-z3e(m)
        x2e(m) = x3e(m)
        z2e(m) = z3e(m)
        w2e(m) = 0.
        w3e(m) = 0.
     end if
     if(x3e(m)<0..and.iperidf==0)then
        x3e(m) = 1.e-8
        z3e(m)=lz-z3e(m)
        x2e(m) = x3e(m)
        z2e(m) = z3e(m)
        w2e(m) = 0.
        w3e(m) = 0.
     end if
     z3e(m)=mod(z3e(m)+8.*lz,lz)
     x3e(m)=mod(x3e(m)+800.*lx,lx)         
     x3e(m) = min(x3e(m),lx-1.0e-8)
     y3e(m) = min(y3e(m),ly-1.0e-8)
     z3e(m) = min(z3e(m),lz-1.0e-8)

  end do
  call MPI_ALLREDUCE(myopz,nopz,1,MPI_integer, &
       MPI_SUM, MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(myoen,noen,1,MPI_integer, &
       MPI_SUM, MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(myaven,aven,1,MPI_REAL8, &
       MPI_SUM, MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(myavptch,avptch,1,MPI_REAL8, &
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
  call pmove(u0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call end_pmove(ierr)
  mme=np_new

  !      return
end subroutine pint
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine cint(n)

  use gem_com
  use gem_equil

  implicit none

  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dgdtp,dpdzp,dadzp,aparp,vncp
  real :: dum,vxdum,dum1,dum2,eps,x,h_x,h_coll,hee,nue
  INTEGER :: m,i,j,k,l,n,mynowe
  INTEGER :: np_old,np_new
  real :: vfac,kap,kapnp,kaptp,sz,vpar,ppar,vpdum,pidum,xnp
  real :: b,th,r,enerb,cost,sint,qr,laps,ter
  real :: xt,xs,yt,zt,xdot,ydot,zdot,pzdot,edot,pzd0,pzd1,vp0
  real :: wx0,wx1,wy0,wy1,wz0,wz1,w3old
  real :: myke,mypfl_es(1:nsubd),mypfl_em(1:nsubd),myptrp
  real :: myefl_es(1:nsubd),myefl_em(1:nsubd),mynos,myavewe
  real :: ketemp,pfltemp,ptrptmp
  real :: efltemp,nostemp
  real :: sbuf(10),rbuf(10)
  real :: mytotn,mytrap,totn,ttrap
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp
  real :: grp,gxdgyp,psp,pzp,psip2p,bdcrvbp,curvbzp,dipdrp,bstar
  real :: x000,x001,x010,x011,x100,x101,x110,x111

  sbuf(1:10) = 0.
  rbuf(1:10) = 0.
  myavewe = 0.
  myke=0.    
  mypfl_es=0.    
  mypfl_em=0.    
  myptrp=0.
  myefl_es=0. 
  myefl_em=0. 
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
     vncp = wx0*phincp(i)+wx1*phincp(i+1)        
     psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
     dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)        
     b=1.-tor+tor*bfldp
     pzp = emass*u3e(m)/b+psp/br0

     xt = x3e(m)
     yt = y3e(m)

     include 'cpushlie.h'

     vfac = 0.5*(emass*u3e(m)**2 + 2.*mue3(m)*b)
     vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
     vp0 = vp0*vncp*vexbsw
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
          (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
          +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
          -emass*vpar**2/(qel*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
          -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*grcgtp*lr0/q0*qhatp  

     zdot =  vpar*b/bstar*(1-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
          +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
          -1./b**2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp &
          -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*q0*br0*grcgtp/jfnp

     pzd0 = tor*(-mue3(m)/emass/radiusp/bfldp*psipp*dbdtp*grcgtp)*b/bstar &
          +mue3(m)*vpar/(qel*bstar*b)*dipdrp/radiusp*dbdtp*grcgtp
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

     eps = (b*mue3(m)+0.5*emass*u3e(m)*u3e(m))/ter
     x = sqrt(eps)
     h_x    = 4.0*x*x/(3.0*sqrt(pi))
     h_coll = h_x/sqrt(1.0+h_x**2)
     !         h_coll = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*DERF(x)
     ! collision frequency for experimental profiles
     hee = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*erf(x)
     nue=rneu*(wx0*nue0(i)+wx1*nue0(i+1))*(wx0*zeff(i)+wx1*zeff(i+1)+hee)
     !         dum1 = 0.5*rneu/eps**1.5*(1+h_coll)
     dum1 = nue/(2*(eps+0.1))**1.5  !*(1+h_coll)
     !         if(x<0.3)dum1=0.0

     dum = 1-w3e(m)*nonline*0.
     if(eldu.eq.1)dum = (tge/ter)**1.5*exp(vfac*(1/tge-1./ter))
     vxdum = (eyp/b+vpdum/b*delbxp)*dum2
     !             -(2*u3e(m)*aparp+qel/emass*aparp*aparp)/(b*b)/br0*sint
     w3old = w3e(m)
     w3e(m)=w2e(m) + dte*(  &
          (vxdum*kap + edot/ter  -dum1*ppar*aparp/ter)*xnp    & 
          +isg*(-dgdtp-zdot*dpdzp+xdot*exp1+ydot*eyp)/ter*xnp)*dum

     if(abs(w3e(m)).gt.4.0.and.nonline==1)then
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
     i = max(i,0)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     qr = wx0*sf(i)+wx1*sf(i+1)
     y3e(m)=mod(y3e(m)-laps*2*pi*qr*lr0/q0*sign(1.0,q0)+8000.*ly,ly)
     if(x3e(m)>lx.and.iperidf==0)then
        x3e(m) = lx-1.e-8
        z3e(m)=lz-z3e(m)
        x2e(m) = x3e(m)
        z2e(m) = z3e(m)
        w2e(m) = 0.
        w3e(m) = 0.
     end if
     if(x3e(m)<0..and.iperidf==0)then
        x3e(m) = 1.e-8
        z3e(m)=lz-z3e(m)
        x2e(m) = x3e(m)
        z2e(m) = z3e(m)
        w2e(m) = 0.
        w3e(m) = 0.
     end if
     z3e(m)=mod(z3e(m)+8.*lz,lz)
     x3e(m)=mod(x3e(m)+800.*lx,lx)         
     x3e(m) = min(x3e(m),lx-1.0e-8)
     y3e(m) = min(y3e(m),ly-1.0e-8)
     z3e(m) = min(z3e(m),lz-1.0e-8)

     !     particle diagnostics done here because info is available...
     k = int(x3e(m)/(lx/nsubd))
     k = min(k,nsubd-1)
     k = k+1
     mypfl_es(k)=mypfl_es(k) + w3old*(eyp) 
     mypfl_em(k)=mypfl_em(k) + w3old*(vpar*delbxp/b) 
     myptrp=myptrp + w3old*eyp*(1-ipass(m)) 
     myefl_es(k)=myefl_es(k) + vfac*w3old*(eyp)
     myefl_em(k)=myefl_em(k) + vfac*w3old*(vpar*delbxp/b)
     myke=myke + vfac*w3e(m)
     mynos=mynos + w3e(m)
     myavewe = myavewe+abs(w3e(m))

     if(abs(z3e(m)-z2e(m)).gt.lz/2)ipass(m)=1

     ! Turning off trapped electrons...
     !         if(tcurr>2000.and.ipass(m)==0) then
     !            w3e(m)=0
     !         endif
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
  sbuf(2)=myefl_es(nsubd/2)
  sbuf(3)=mypfl_es(nsubd/2)
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
  avewe(n) = rbuf(8)/( float(tmm(2)) )
  !      ke(2,n)=ketemp/( 2.*float(tmm(2))*mims(ns) )
  ftrap = ttrap/totn
  !      nos(2,n)=nostemp/( float(tmm(2)) )

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  sbuf(1:nsubd) = myefl_es(1:nsubd)
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,  &
       MPI_COMM_WORLD,ierr)
  do k = 1,nsubd   
     efle_es(k,n)=rbuf(k)/( float(tmm(2)) )*totvol/vol(k)*cn0e
  end do

  sbuf(1:nsubd) = myefl_em(1:nsubd)
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,  &
       MPI_COMM_WORLD,ierr)
  do k = 1,nsubd   
     efle_em(k,n)=rbuf(k)/( float(tmm(2)) )*totvol/vol(k)*cn0e
  end do

  sbuf(1:nsubd) = mypfl_es(1:nsubd)
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,  &
       MPI_COMM_WORLD,ierr)
  do k = 1,nsubd   
     pfle_es(k,n)=rbuf(k)/( float(tmm(2)) )*totvol/vol(k)*cn0e
  end do

  sbuf(1:nsubd) = mypfl_em(1:nsubd)
  call MPI_ALLREDUCE(sbuf,rbuf,10,  &
       MPI_REAL8,MPI_SUM,  &
       MPI_COMM_WORLD,ierr)
  do k = 1,nsubd   
     pfle_em(k,n)=rbuf(k)/( float(tmm(2)) )*totvol/vol(k)*cn0e
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
  call pmove(u0e,np_old,np_new,ierr)
  if (ierr.ne.0) call ppexit

  call end_pmove(ierr)
  mme=np_new

  !      return
end subroutine cint
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine drdt(ip)
  use gem_com
  use gem_equil
  implicit none
  integer :: i,j,k,ip,isdndt
  real :: lbfr(0:imx,0:jmx)
  real :: lbfs(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
  real :: dum
  real :: djdx(0:imx,0:jmx,0:1),djdy(0:imx,0:jmx,0:1)

  call gradx(jionx*ision-upex,djdx)
  call grady(jiony*ision-upey,djdy)
  !      djdx = 0.
  !      djdy = 0.
  isdndt = 1
  !      upa0(:,:,:) = apar(ip,:,:,:)
  do i = 0,im
     do j = 0,jm
        rbfs(i,j)=(jion(i,j,0)*ision-q(1)*q(1)/mims(1)*gn0s(1,i)*cn0s(1)*&
             apar(i,j,0)*bdgrzn(i,0)/ntube*isiap-(upazd(i,j,0)+amie*upa0(i,j,0)))*jac(i,0)
     end do
  end do
  call MPI_SENDRECV(rbfs,(imx+1)*(jmx+1),  &
       MPI_REAL8,rngbr,404, &
       lbfr,(imx+1)*(jmx+1), &
       MPI_REAL8,lngbr,404, &
       tube_comm,stat,ierr)
  do i = 0,im
     do j = 0,jm
        lbfs(i,j)=(jion(i,j,1)*ision-q(1)*q(1)/mims(1)*gn0s(1,i)*cn0s(1)*apar(i,j,1)*&
             bdgrzn(i,1)/ntube*isiap-(upazd(i,j,1)+amie*upa0(i,j,1)))*jac(i,1)
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
        drhodt(i,j,0) = -((jion(i,j,1)*ision-q(1)*q(1)/mims(1)*gn0s(1,i)*cn0s(1)*apar(i,j,1)*&
             bdgrzn(i,1)/ntube*isiap-upazd(i,j,1) &
             -amie*upa0(i,j,1))*jac(i,1)-dum)/(2.*dz)/jac(i,0) &
             -djdx(i,j,0)-djdy(i,j,0) &
             +(drhoidt(i,j,0)*ision-dnedt(i,j,0))

        dum = weightm(i)*rbfr(i,jmi(i,j)) &
             +weightmn(i)*rbfr(i,jmn(i,j))
        drhodt(i,j,1) = -(dum-(jion(i,j,0)*ision-q(1)*q(1)/mims(1)*gn0s(1,i)*cn0s(1)*apar(i,j,0)*&
             bdgrzn(i,0)/ntube*isiap-upazd(i,j,0)  &
             -amie*upa0(i,j,0))*jac(i,0))/(2.*dz)/jac(i,1)  &
             -djdx(i,j,1)-djdy(i,j,1) &
             +(drhoidt(i,j,1)*ision-dnedt(i,j,1))
     end do
  end do
  call enfxy(drhodt(:,:,:))
  call enfz(drhodt(:,:,:))
  !      call filter(drhodt(:,:,:))
  !      return
end subroutine drdt
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine dpdt(ip)   
  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none
  INTEGER :: ns
  real :: b,gam0,gam1,delyz,th,bf,dum,r,qr,shat
  real,DIMENSION(1:5) :: b1,b2,ter
  real :: kx1,kx2,ky
  real,dimension(:),allocatable :: akx,aky
  real,dimension(:,:,:,:,:),allocatable:: gamb1,gamb2
  complex,dimension(:,:,:,:),allocatable :: mx
  integer,dimension(:,:,:,:),allocatable :: ipiv
  real,dimension(:,:,:),allocatable :: formdpt
  real :: sgnx,sgny,sz,myfe,u(0:imx,0:jmx,0:1)
  INTEGER :: i,i1,j,j1,k,k1,l,m,n,ifirst,nstep,ip,INFO
  INTEGER :: l1,m1,myk,myj,ix,ikx
  complex :: temp3dxy(0:imx-1,0:jmx-1,0:1),v(0:imx-1,0:jcnt-1,0:1)
  complex :: sbuf(0:imx*jcnt*2-1),rbuf(0:imx*jmx*2-1)
  complex :: sl(1:imx-1,0:jcnt-1,0:1)
  complex :: cdum
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
  real :: grp,gxdgyp,grdgtp,gthp
  real :: wx0,wx1,wz0,wz1
  character(len=70) fname
  character(len=5) holdmyid

  save formdpt,ifirst,akx,aky,mx,gamb1,gamb2,ipiv
  if(idg==1)write(*,*)'enter gkps'

  !     form factors....
  write(holdmyid,'(I5.5)') MyId
  fname='./matrix/'//'mx_dpt_'//holdmyid
  if (ifirst.ne.-99) then
     allocate(akx(0:imx-1),aky(0:jcnt-1),&
          gamb1(1:5,0:imx-1,0:jcnt-1,0:imx-1,0:1),&
          gamb2(1:5,0:imx-1,0:jcnt-1,0:imx-1,0:1))
     allocate(mx(imx-1,imx-1,0:jcnt-1,0:1),formdpt(0:imx-1,0:jcnt-1,0:1))
     allocate(ipiv(imx-1,imx-1,0:jcnt-1,0:1))

     if(iget.eq.1) then
        open(30000+MyId,file=fname,form='unformatted',status='old')
        read(30000+MyId)mx,ipiv
        close(30000+myid)
        goto 200
     end if

     do l=0,im-1
        do m=0,jcnt-1
           j = tclr*jcnt+m
           do i=0,im-1 
              r = lr0-lx/2+i*dx
              i1 = int((r-rin)/dr)
              i1 = min(i1,nr-1)
              wx0 = (rin+(i1+1)*dr-r)/dr
              wx1 = 1.-wx0
              do ns = 1, nsm
                 ter(ns) = wx0*t0s(ns,i1)+wx1*t0s(ns,i1+1)
              enddo
              do n = 0,1
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
                 do ns = 1, nsm
                    b1(ns)=mims(ns)*(kx1*kx1*grp**2 + &
                         ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                         +2*dydrp*lr0/q0*qhatp*grdgtp) &
                         +2*kx1*ky*gxdgyp)/(bf*bf)*ter(ns)/(q(ns)*q(ns))

                    b2(ns)=mims(ns)*(kx2*kx2*grp**2 + &
                         ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                         +2*dydrp*lr0/q0*qhatp*grdgtp) &
                         +2*kx2*ky*gxdgyp)/(bf*bf)*ter(ns)/(q(ns)*q(ns))

                    call srcbes(b1(ns),gam0,gam1)
                    gamb1(ns,l,m,i,n)=gam0
                    call srcbes(b2(ns),gam0,gam1)
                    gamb2(ns,l,m,i,n)=gam0
                 end do

                 !   formfactor in gkps
                 formdpt(l,m,n) = 1./jmx 
                 if(abs(ky)>kycut)formdpt(l,m,n) = 0.
                 !                    if(abs(ky)==0.)formdpt(l,m,n) = 0.
                 if(m1.ne.mlk.and.onemd==1)formdpt(l,m,n) = 0.
              enddo
           enddo
        enddo
     enddo

     do k = 0,1
        do j = 0,jcnt-1
           do i = 1,imx-1
              r = lr0-lx/2+i*dx
              i1 = int((r-rin)/dr)
              i1 = min(i1,nr-1)
              wx0 = (rin+(i1+1)*dr-r)/dr
              wx1 = 1.-wx0
              do ns = 1, nsm
                 ter(ns) = wx0*t0s(ns,i1)+wx1*t0s(ns,i1+1)
              enddo
              do ix = 1,imx-1
                 mx(i,ix,j,k) = 0.
                 do ikx = 0,imx-1
                    do ns = 1, nsm
                       mx(i,ix,j,k) = mx(i,ix,j,k)+q(ns)*sin(ix*ikx*pi/imx)* &
                            ((1-gamb1(ns,ikx,j,i,k))*exp(IU*ikx*i*pi/imx)- &
                            (1-gamb2(ns,ikx,j,i,k))*exp(-IU*ikx*i*pi/imx)) &
                            /ter(ns)*cn0s(ns)*gn0s(ns,i)/(IU*imx)
                    end do
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
200  ifirst=-99
  endif
  if(idg==1)write(*,*)'pass form factors'

  call dcmpy(drhodt(0:imx-1,0:jmx-1,0:1),v)

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
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3dxy(i,j,k) = tmpy(j)  ! phi(x,y)
        end do
     end do
  end do

100 continue
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

  if(idpbf==0)call fltx(dphidt(:,:,:),0,1)
  if(idpbf==1)call fltx(dphidt(:,:,:),1,1)
  if(idpbf==0)call filter(dphidt(:,:,:))

  if(idg==1)write(*,*)'pass enfz', myid

  !      return
end subroutine dpdt

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine jie(ip,n)

  !    source quantities are are calculated: n_i
  !    right now only ion quantitities are calculated...

  use gem_com
  use gem_equil
  implicit none
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
  real :: enerb,vxdum,dum,xdot,ydot,zdot,pidum,dum1,dum2
  INTEGER :: m,n,i,j,k,l,ns,ip,nonfi,nonfe
  real :: wx0,wx1,wy0,wy1,wz0,wz1,vte
  real :: sz,wght,wght0,wght1,wght2,r,th,cost,sint,b,qr,dv,kap,ter
  real :: kapnp,kaptp,xnp,psip2p,bdcrvbp,curvbzp,dipdrp,bstar
  real :: xt,yt,rhog,vpar,xs,dely,vfac,vp0
  real :: lbfs(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
  real :: lbfr(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: myjpar(0:imx,0:jmx,0:1),myjpex(0:imx,0:jmx,0:1)
  real :: myjpey(0:imx,0:jmx,0:1),myupar(0:imx,0:jmx,0:1)
  real :: myupex(0:imx,0:jmx,0:1),myupey(0:imx,0:jmx,0:1)
  real :: myupazd(0:imx,0:jmx,0:1)
  real :: mydnidt(0:imx,0:jmx,0:1),mydnedt(0:imx,0:jmx,0:1)
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp
  real :: grp,gxdgyp,rhox(4),rhoy(4),vncp,vparspp

  nonfi = 1 
  nonfe = 1 
  jion = 0.
  jionx = 0.
  jiony = 0.
  drhoidt = 0.
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

        curvbzp = wx0*wz0*curvbz(i,k)+wx0*wz1*curvbz(i,k+1) &
             +wx1*wz0*curvbz(i+1,k)+wx1*wz1*curvbz(i+1,k+1) 
        bdcrvbp = wx0*wz0*bdcrvb(i,k)+wx0*wz1*bdcrvb(i,k+1) &
             +wx1*wz0*bdcrvb(i+1,k)+wx1*wz1*bdcrvb(i+1,k+1) 
        grdgtp = wx0*wz0*grdgt(i,k)+wx0*wz1*grdgt(i,k+1) &
             +wx1*wz0*grdgt(i+1,k)+wx1*wz1*grdgt(i+1,k+1) 

        fp = wx0*f(i)+wx1*f(i+1)        
        jfnp = wz0*jfn(k)+wz1*jfn(k+1)
        psipp = wx0*psip(i)+wx1*psip(i+1)        
        ter = wx0*t0s(ns,i)+wx1*t0s(ns,i+1)        
        kaptp = wx0*capts(ns,i)+wx1*capts(ns,i+1)        
        kapnp = wx0*capns(ns,i)+wx1*capns(ns,i+1)        
        xnp = wx0*xn0s(ns,i)+wx1*xn0s(ns,i+1)        

        vncp = wx0*phincp(i)+wx1*phincp(i+1)        
        vparspp = wx0*vparsp(ns,i)+wx1*vparsp(ns,i+1)        
        psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
        dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)        

        !         b=1.-lr0/br0*cost
        b=1.-tor+tor*bfldp

        rhog=sqrt(2.*b*mu(ns,m)*mims(ns))/(q(ns)*b)*iflr
        vfac = 0.5*(mims(ns)*u3(ns,m)**2 + 2.*mu(ns,m)*b)
        if(ildu.eq.1)xnp = xnp*(tgis(ns)/ter)**1.5*exp(vfac*(1/tgis(ns)-1./ter))
        vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
        vp0 = vp0*vncp*vexbsw

        vpar = u3(ns,m)
        kap = kapnp - (1.5-vfac/ter)*kaptp-vpar*mims(ns)/ter*vparspp*vparsw
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
        do l=1,lr(1)
           xs=x3(ns,m)+rhox(l) !rwx(1,l)*rhog
           yt=y3(ns,m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
           xt=mod(xs+800.*lx,lx)
           yt=mod(yt+800.*ly,ly)
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
        enddo

        exp1=exp1/4.
        eyp=eyp/4.
        delbxp=delbxp/4.
        delbyp=delbyp/4.
        aparp = aparp/4.

        vpar = u3(ns,m)-q(ns)/mims(ns)*aparp*0.
        bstar = b*(1+mims(ns)*vpar/(q(ns)*b)*bdcrvbp)
        enerb=(mu(ns,m)+mims(ns)*vpar*vpar/b)/q(ns)*b/bstar*tor
        dum1 = 1./b*lr0/q0*qhatp*fp/radiusp*grcgtp
        vxdum = (eyp/b+vpar/b*delbxp)*dum1

        xdot = vxdum*nonlin(ns)  &
             -iorb*enerb/bfldp/bfldp*fp/radiusp*dbdtp*grcgtp
        ydot = (-exp1/b+vpar/b*delbyp)*dum1*nonlin(ns)  &
             +iorb*enerb/bfldp/bfldp*fp/radiusp*grcgtp* &
             (-dydrp*dbdtp+r0/q0*qhatp*dbdrp)+vp0 &
             +enerb/(bfldp**2)*psipp*lr0/q0/radiusp**2*(dbdrp*grp**2+dbdtp*grdgtp) &
             -mims(ns)*vpar**2/(q(ns)*bstar*b)*(psip2p*grp**2/radiusp+curvbzp)*lr0/(radiusp*q0) &
             -dipdrp/radiusp*mims(ns)*vpar**2/(q(ns)*bstar*b)*grcgtp*lr0/q0*qhatp  
        zdot =  vpar*b/bstar*(1.-tor+tor*q0*br0/radiusp/b*psipp*grcgtp)/jfnp &
             +q0*br0*enerb/(b*b)*fp/radiusp*dbdrp*grcgtp/jfnp &
             -1./b**2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp &
             -dipdrp/radiusp*mims(ns)*vpar**2/(q(ns)*bstar*b)*q0*br0*grcgtp/jfnp

        !    now do 1,2,4 point average, where lr is the no. of points...
        do l=1,lr(1)
           xs=x3(ns,m)+rhox(l) !rwx(1,l)*rhog
           yt=y3(ns,m)+rhoy(l) !(rwy(1,l)+sz*rwx(1,l))*rhog
           xt=mod(xs+800.*lx,lx)
           yt=mod(yt+800.*ly,ly)
           xt = min(xt,lx-1.0e-8)
           yt = min(yt,ly-1.0e-8)

           i=int(xt/dx+0.5)
           j=int(yt/dy+0.5)
           k=int(z3(ns,m)/dz+0.5)-gclr*kcnt

           wght1 = wght0*(vxdum*kap+q(ns)*(xdot*exp1/ter+(ydot-vp0)*eyp/ter))*xnp
           myjpar(i,j,k) = myjpar(i,j,k)+wght*zdot
           myjpex(i,j,k) = myjpex(i,j,k)+wght*xdot
           myjpey(i,j,k) = myjpey(i,j,k)+wght*ydot
           mydnidt(i,j,k) = mydnidt(i,j,k)+wght1
        enddo
     enddo

     !   enforce periodicity
     call enforce(myjpar)
     call enforce(myjpex)
     call enforce(myjpey)
     call enforce(mydnidt)
     !      call filter(myjpar)
     !      call filter(mydnidt)

     do i=0,im
        do j=0,jm
           do k=0,mykm
              jpar(ns,i,j,k) = q(ns)*myjpar(i,j,k)/n0/jac(i,k)*cn0s(ns)
              jpex(ns,i,j,k) = q(ns)*myjpex(i,j,k)/n0/jac(i,k)*cn0s(ns)
              jpey(ns,i,j,k) = q(ns)*myjpey(i,j,k)/n0/jac(i,k)*cn0s(ns)
              dnidt(ns,i,j,k) = q(ns)*mydnidt(i,j,k)/n0/jac(i,k)*cn0s(ns)
           enddo
        enddo
     enddo

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

  vte = sqrt(amie*t0e(nr/2))
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

     vncp = wx0*phincp(i)+wx1*phincp(i+1)        
     !         b=1.-lr0/br0*cost
     b=1.-tor+tor*bfldp
     psip2p = wx0*psip2(i)+wx1*psip2(i+1)        
     dipdrp = wx0*dipdr(i)+wx1*dipdr(i+1)        

     vfac = 0.5*(emass*u3e(m)**2 + 2.*mue3(m)*b)
     if(eldu.eq.1)xnp = xnp*(tge/ter)**1.5*exp(vfac*(1/tge-1./ter))
     vp0 = 1./b**2*lr0/q0*qhatp*fp/radiusp*grcgtp
     vp0 = vp0*vncp*vexbsw
     kap = kapnp - (1.5-vfac/ter)*kaptp

     wght=w3e(m)/dv
     wght0=exp(-vfac)/dv
     if(isuni.eq.0)wght0=1./dv
     wght0 = wght0
     vpar = u3e(m)
     !         if(abs(vpar/vte).gt.vcut)wght = 0.
     !         if(abs(vpar/vte).gt.vcut)wght0 = 0.

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
          -1./b**2*q0*br0*fp/radiusp*grcgtp*vncp*vexbsw/jfnp &
          -dipdrp/radiusp*emass*vpar**2/(qel*bstar*b)*q0*br0*grcgtp/jfnp

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
     vxdum = (eyp/b+vpar/b*delbxp*iflut)*dum1
     xdot = xdot+vxdum*nonline
     ydot = ydot+(-exp1/b+vpar/b*delbyp*iflut)*dum1*nonline

     phip = w000(m)*phi(i,j,k)  &
          + w100(m)*phi(i+1,j,k) &
          + w010(m)*phi(i,j+1,k) &
          + w110(m)*phi(i+1,j+1,k) &
          + w001(m)*phi(i,j,k+1) &
          + w101(m)*phi(i+1,j,k+1) &
          + w011(m)*phi(i,j+1,k+1) &
          + w111(m)*phi(i+1,j+1,k+1)

     wght1 = (wght+1./dv*phip*isg/ter*xnp)
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
     vxdum = (eyp/b+vpar/b*delbxp)*dum2 
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
           upex(i,j,k) = myupex(i,j,k)/n0e/jac(i,k)*cn0e
           upey(i,j,k) = myupey(i,j,k)/n0e/jac(i,k)*cn0e
           upazd(i,j,k) = myupazd(i,j,k)/n0e/jac(i,k)*cn0e
           dnedt(i,j,k) = mydnedt(i,j,k)/n0e/jac(i,k)*cn0e
        end do
     end do
  end do


999 continue
  !      return
end subroutine jie
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine lorentz(ip,n)

  use gem_com
  use gem_equil

  implicit none

  integer :: i,ip,k,m,n,ncol,icol
  real :: edum,vdum,dum,dum1,ptch,vte,r,qr,th,cost,b
  real :: h_x,h_coll,x,eps,dtcol,uold,hee,nue,ter
  real :: wx0,wx1,wz0,wz1

  ncol = 1
  if(ip.eq.1)dtcol = dt/ncol*2
  if(ip.eq.0)dtcol = dt/ncol
  vte = sqrt(amie*t0e(nr/2))
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
     ter = wx0*t0e(i)+wx1*t0e(i+1)        
     b = wx0*wz0*bfld(i,m)+wx0*wz1*bfld(i,m+1) &
          +wx1*wz0*bfld(i+1,m)+wx1*wz1*bfld(i+1,m+1) 
     uold = u3e(k)
     edum = b*mue3(k)+0.5*emass*u3e(k)*u3e(k)
     vdum = sqrt(2.*edum/emass)
     ptch = u3e(k)/vdum
     eps = edum/ter
     x = sqrt(eps)
     h_x    = 4.0*x*x/(3.0*sqrt(pi))
     h_coll = h_x/sqrt(1.0+h_x**2)
     !         h_coll = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*DERF(x)
     ! collision frequency for experimental profiles
     hee = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*erf(x)
     nue=rneu*(wx0*nue0(i)+wx1*nue0(i+1))*(wx0*zeff(i)+wx1*zeff(i+1)+hee)
     dum1 = 1/(2*(eps+0.1))**1.5  !*(1+h_coll)
     !         dum = dtcol*rneu*dum1
     dum = dtcol*nue*dum1
     !         if(x<0.3)dum=0.0
     do icol = 1,ncol
        ptch =ptch-dum*ptch+sqrt(dum*(1.-ptch*ptch)) &
             *sign(1.0,ran2(iseed)-0.5)
        ptch = min(ptch,0.999)
        ptch = max(ptch,-0.999)
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

end subroutine lorentz
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine initialize
  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none
  real :: dum,dum1,dum2,jacp,xndum,r,wx0,wx1
  !        complex(8),dimension(0:1) :: x,y
  real,dimension(0:1) :: x,y
  integer :: n,i,j,k,ip

  call init
  call mpi_barrier(mpi_comm_world,ierr)

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
  n0e=mme*numprocs/totvol
  if(myid==0)then
     write(*,*)'totvol,jacp,dum2=',totvol,jacp,dum2
  end if

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

  call ccfft('x',0,imx,0.0,tmpx,coefx,workx,0)
  call ccfft('y',0,jmx,0.0,tmpy,coefy,worky,0)
  call ccfft('z',0,kmx,0.0,tmpz,coefz,workz,0)
  call dsinf(1,x,1,0,1,0,imx*2,1,1.0,aux1,50000,aux2,20000)

  ncurr = 1
  call blendf
end subroutine initialize
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine loader_wrapper
  use gem_com
  use gem_equil
  implicit none

  integer :: n,i,j,k,ip,ns

  do ns = 1,nsm
     if(isuni.eq.0)call loadi(ns)
  end do
  if(ifluid.eq.1)call ldel
  if(idg.eq.1)write(*,*)'past loader'
end subroutine loader_wrapper
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine accumulate(n,ip)
  use gem_com
  use gem_equil
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
  use gem_equil
  implicit none
  integer :: n,i,i1,j,k,ip,it,iter=1
  real :: myrmsphi,rmp(20),myavap(0:imx-1)
  real :: myjaca(0:imx-1),jaca(0:imx-1)

  do it = 1,iter
     call den0_phi(ip,n,it)
     if(iperi==1)call gkpsL(n,ip)
     if(iperi==0)then
       if(adiabatic_electron==0)then ! jycheng
         call gkps(n,ip)
       elseif(adiabatic_electron==1)then
         call gkps_adiabatic_electron(n,ip)
       else
         write(*,*)'wrong adiabatic_electron option'
         exit
       endif
     endif

     if(idg.eq.1)write(*,*)'pass gkps in poisson'
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
  if(idg.eq.1)write(*,*)'pass iter loop in poisson'
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
100 continue

  if(iperi==0 .and. iphbf==0)call fltx(phi(:,:,:),0,1)
  if(iperi==0 .and. iphbf==1)call fltx(phi(:,:,:),1,1)
  if(iperi==1 .and. iphbf==1)call flty(phi(:,:,:),1)
  if(iphbf==0)call filter(phi(:,:,:))

  !remove zonal field
  if(izonal==0)then
     do i=0,im-1
        myavap(i) = 0.
        myjaca(i) = 0.
        do j=0,jm-1
           myavap(i)=myavap(i)+phi(i,j,0)*jac(i,0)
           myjaca(i)=myjaca(i)+jac(i,0)
        enddo
     enddo
     call MPI_ALLREDUCE(myavap,avap,imx, &
          MPI_REAL8,                               &
          MPI_SUM,TUBE_COMM,ierr)
     call MPI_ALLREDUCE(myjaca,jaca,imx, &
          MPI_REAL8,                               &
          MPI_SUM,TUBE_COMM,ierr)

     avap(0:imx-1)=avap(0:imx-1)/jaca(0:imx-1)

     do i = 0,imx-1
        do j = 0,jmx
           do k = 0,1
              phi(i,j,k) = phi(i,j,k)-avap(i)
           end do
        end do
     end do
  end if

  myrmsphi=0.
  do k=0,mykm-1
     do j=0,jm-1
        do i1=0,im-1
           myrmsphi=myrmsphi+phi(i1,j,k)*phi(i1,j,k)
        enddo
     enddo
  enddo
  call MPI_ALLREDUCE(myrmsphi,rmsphi(n),1, &
       MPI_REAL8,                               &
       MPI_SUM,TUBE_COMM,ierr)
  rmsphi(n)=sqrt(rmsphi(n)/(im*jm*km))

  if(idg.eq.1)write(*,*)'pass poisson'
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine poisson
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine ampere(n,ip)
  use gem_com
  use gem_equil
  implicit none
  integer :: iter=10
  integer :: n,i,i1,j,k,ip
  real :: myrmsapa,rma(20),myavap(0:imx-1)
  real :: myjaca(0:imx-1),jaca(0:imx-1)

  if(ifluid==1.and.beta.gt.1.e-8)then
     do i = 1,iter
        call jpar0(ip,n,i,0)
        if(idg.eq.1)write(*,*)'pass jpar0'
        if(iperi==1)call ezampL(n,ip)
        if(iperi==0)call ezamp(n,ip)
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

  if(iperi==0 .and. iapbf==0)call fltx(apar(:,:,:),0,1)
  if(iperi==0 .and. iapbf==1)call fltx(apar(:,:,:),1,1)
  if(iperi==1 .and. iapbf==1)call flty(apar(:,:,:),1)
  if(iapbf==0)call filter(apar(:,:,:))

  !remove zonal field
  if(izonal==0)then
     do i=0,im-1
        myavap(i) = 0.
        myjaca(i) = 0.
        do j=0,jm-1
           myavap(i)=myavap(i)+apar(i,j,0)*jac(i,0)
           myjaca(i)=myjaca(i)+jac(i,0)
        enddo
     enddo
     call MPI_ALLREDUCE(myavap,avap,imx, &
          MPI_REAL8,                               &
          MPI_SUM,TUBE_COMM,ierr)
     call MPI_ALLREDUCE(myjaca,jaca,imx, &
          MPI_REAL8,                               &
          MPI_SUM,TUBE_COMM,ierr)

     avap(0:imx-1)=avap(0:imx-1)/jaca(0:imx-1)

     do i = 0,imx-1
        do j = 0,jmx
           do k = 0,1
              apar(i,j,k) = apar(i,j,k)-avap(i)
           end do
        end do
     end do
  end if

  myrmsapa=0.
  do k=0,mykm-1
     do j=0,jm-1
        do i1=0,im-1
           myrmsapa=myrmsapa+apar(i1,j,k)*apar(i1,j,k)
        enddo
     enddo
  enddo
  call MPI_ALLREDUCE(myrmsapa,rmsapa(n),1, &
       MPI_REAL8,                               &
       MPI_SUM,TUBE_COMM,ierr)
  rmsapa(n)=sqrt(rmsapa(n)/(im*jm*km))

  if(idg.eq.1)write(*,*)'pass filter(apar)'
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine ampere
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine split_weight(n,ip)
  use gem_com
  use gem_equil
  implicit none

  integer :: n,i,j,k,ip
  if(isg.gt.0..and.ifluid.eq.1)then
     call jie(ip,n)
     if(idg.eq.1)write(*,*)'pass jie'
     call drdt(ip)
     if(idg.eq.1)write(*,*)'pass drdt'
     if(iperi==1)call dpdtL(ip)
     if(iperi==0)call dpdt(ip)
     if(idg.eq.1)write(*,*)'pass dpdt'
  end if
  if(idg.eq.1)write(*,*)'pass split_weight'
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine split_weight
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine field(n,ip)
  use gem_com
  use gem_equil
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
  use gem_equil
  implicit none
  integer :: n,i,j,k,ip,ns

  do ns = 1,nsm
     if(ip.eq.1.and.ision==1)call ppush(n,ns)
     if(ip.eq.0.and.ision==1)call cpush(n,ns)
  end do
  if(ip==0)call colli(ip,n)
  if(idg.eq.1)write(*,*)'pass ppush'

  if(ip.eq.1.and.ifluid==1)call pint
  if(ip.eq.0.and.ifluid==1)call cint(n)
  if(idg.eq.1)write(*,*)'pass pint'
  if(ifluid==1.and.ip==0)call lorentz(ip,n)
  !         if(ifluid==1.and.ip==0)call col(dt)
  !         if(ip.eq.0.and.mod(n+1,10).eq.0)call avgi(1)
  if(eprs>0.and.ip.eq.0.and.mod(n+1,nrst).eq.0)call rstpe

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)          
end subroutine push_wrapper
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine diagnose(n)
  use gem_com
  use gem_equil
  implicit none
  integer :: n,i,j,k,ip

  call modes2(phi,pmodehis,n)
  if(idg.eq.1)write(*,*)'pass modes2'  
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)          
  call yveck(phi(:,:,:),n)
  call yveck1(phi(:,:,:),n)
!  call mdampd(phi(:,:,:),mdhis)
!  call mdampd(apar(:,:,:),mdhisa)
!  call fltd(dti(1,:,:,:))
!  call fltd(dti(2,:,:,:))
!  call fltd(delte(:,:,:))
!  call mdampd(dti(1,:,:,:),mdhisb)
!  call mdampd(dti(2,:,:,:),mdhisc)
!  call mdampd(delte(:,:,:),mdhisd)
  if(idg.eq.1)write(*,*)'pass yvec'    
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)        

end subroutine diagnose
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine reporter(n)
  use gem_com
  use gem_equil
  implicit none
  integer :: n,i,j,k,ip

  if(mod(n,xnplt).eq.0) then
     call spec(n)
  endif
13 format(1x,i6,7(2x,i7))
  if(myid.eq.master)then
     open(16, file='indicator', status='unknown',position='append')
     write(16,13)n,ipred,icorr,jpred,jcorr,nopz,noen,nowe
     close(16)
  end if
  !        if(myid==0)write(*,13)n,ipred,icorr,jpred,jcorr,nopz,noen
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  !     save particle arrays for restart if iput=1...
  !     do this before the code crashes due to graphics problems
  if((iput.eq.1).and.mod(n+1,500).eq.0)call restart(2,n)

  !     periodically make output for plots
  call outd(n)
  if(idg.eq.1)write(*,*)'pass outd'

end subroutine reporter
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine jpar0(ip,n,it,itp)

  !    source quantities are are calculated: n_i
  !    right now only ion quantitities are calculated...

  use gem_com
  use gem_equil
  implicit none
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
  real :: enerb,vxdum,dum,xdot,ydot,avede0
  INTEGER :: m,n,i,j,k,l,ns,ip,it,itp
  real :: wx0,wx1,wy0,wy1,wz0,wz1,vte
  real :: sz,wght,wght0,wght1,r,th,cost,sint,b,qr,dv,xnp
  real :: xt,yt,zt,rhog,pidum,vpar,xs,dely,vfac,ter
  real :: lbfs(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
  real :: lbfr(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: myupa(0:imx,0:jmx,0:1),myupa0(0:imx,0:jmx,0:1),myden0(0:imx,0:jmx,0:1)
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
  real :: grp,gxdgyp
  real :: zdot

  pidum = 1./(pi*2)**1.5*(vwidthe)**3
  if(isuni.eq.0)pidum = 1.
  ns=1

  !    electrons current due to f_M (p_para)
  vte = sqrt(amie*t0e(nr/2))
  myupa = 0.
  myupa0 = 0.
  myden0 = 0.
  !      upa0(:,:,:) = apar(ip,:,:,:)
  !      return
  if(it.eq.1)then
     apar(:,:,:) = 0.
     upa0(:,:,:) = 0.
     upa00(:,:,:) = 0.
     den0apa(:,:,:) = 0.
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
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
     b=1.-tor+tor*bfldp

     vfac = 0.5*(emass*u3e(m)**2 + 2.*mue3(m)*b)
     if(eldu.eq.1)xnp = xnp*(tge/ter)**1.5*exp(vfac*(1/tge-1./ter))
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
     !         if(abs(vpar/vte).gt.vcut)then
     !            wght0 = 0.
     !            wght1 = 0.
     !         end if

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

        wght0 = 1./dv*aparp*vpar/ter*xnp 
        myden0(i,j,k)      =myden0(i,j,k)+wght0*w000(m)
        myden0(i+1,j,k)    =myden0(i+1,j,k)+wght0*w100(m)
        myden0(i,j+1,k)    =myden0(i,j+1,k)+wght0*w010(m)
        myden0(i+1,j+1,k)  =myden0(i+1,j+1,k)+wght0*w110(m)
        myden0(i,j,k+1)    =myden0(i,j,k+1)+wght0*w001(m)
        myden0(i+1,j,k+1)  =myden0(i+1,j,k+1)+wght0*w101(m)
        myden0(i,j+1,k+1)  =myden0(i,j+1,k+1)+wght0*w011(m)
        myden0(i+1,j+1,k+1)=myden0(i+1,j+1,k+1)+wght0*w111(m)
     end if
  enddo

  !   enforce periodicity
  if(itp==1)call enforce(myupa(:,:,:))
  if(itp==1)call enforce(myden0(:,:,:))
  if(itp==0)call enforce(myupa0(:,:,:))
  !      call filter(myupa(:,:,:))

  if(itp==1)then
     do  i=0,im
        do  j=0,jm
           do  k=0,mykm
              upa0(i,j,k)= myupa(i,j,k)/n0e/jac(i,k)*pidum*cn0e
              den0apa(i,j,k)= myden0(i,j,k)/n0e/jac(i,k)*pidum*cn0e
           end do
        end do
     end do
  end if
  if(itp==0)then
     do  i=0,im
        do  j=0,jm
           do  k=0,mykm
              upa00(i,j,k)= myupa0(i,j,k)/n0e/jac(i,k)*pidum*cn0e
           end do
        end do
     end do
  end if

  call MPI_ALLREDUCE(upa00(0:im,0:jm,0:1),  &
       upa0t(0:im,0:jm,0:1),             &
       (imx+1)*(jmx+1)*2,MPI_REAL8,       &
       MPI_SUM,GRID_COMM,ierr)


  upa0(:,:,:) = upa0(:,:,:)+ &
       (isg*phi(:,:,:)+dene(:,:,:))*apar(:,:,:)*nonline*0.

  upa00(:,:,:) = upa00(:,:,:)+ &
       (isg*phi(:,:,:)+dene(:,:,:))*apar(:,:,:)*nonline*0.
999 continue

  !      return
end subroutine jpar0

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine den0_phi(ip,n,it)

  !    source quantities are are calculated: n_i
  !    right now only ion quantitities are calculated...

  use gem_com
  use gem_equil
  real :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
  real :: enerb,vxdum,dum,xdot,ydot,avede0
  INTEGER :: m,n,i,j,k,l,ns,ip,it
  real :: wx0,wx1,wy0,wy1,wz0,wz1,vte
  real :: sz,wght,wght0,r,th,cost,sint,b,qr,dv
  real :: xt,yt,zt,rhog,pidum,vpar,xs,dely,vfac
  real :: lbfs(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
  real :: lbfr(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: myden0(0:imx,0:jmx,0:1)
  real :: xnp,ter,dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp,grdgtp
  real :: grp,gxdgyp,psp,pzp,psip2p,bdcrvbp,curvbzp,dipdrp,bstar

  pidum = 1./(pi*2)**1.5*(vwidthe)**3
  if(isuni.eq.0)pidum = 1.
  ns=1

  ! electrons density due to phi*f_M (p_para)
  vte = sqrt(amie*t0e(nr/2))
  myden0 = 0.
  den0 = 0.
  if(it.eq.1)then
     phi(:,:,:) = 0.
     den0(:,:,:) = 0.
     return
  end if
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

     ter = wx0*t0e(i)+wx1*t0e(i+1)        
     xnp = wx0*xn0e(i)+wx1*xn0e(i+1)        
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
     b=1.-tor+tor*bfldp
     vfac = 0.5*(emass*u3e(m)**2 + 2.*mue3(m)*b)
     if(eldu.eq.1)xnp = xnp*(tge/ter)**1.5*exp(vfac*(1/tge-1./ter))
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

     wght0 = wght0*phip*xnp/ter
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
           den0(i,j,k)= myden0(i,j,k)/n0e/jac(i,k)*pidum*cn0e
        end do
     end do
  end do

  !      return
end subroutine den0_phi
!!!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine countw(n)
  use gem_com
  use gem_equil
  implicit none

  INTEGER :: n,nbin,m,i,j,k
  parameter(nbin=200)
  real :: myavwe,avwe,mymaxw,myminw,maxw,minw
  real :: sbuf(10),rbuf(10)
  real,dimension(:),allocatable :: wghtmax,wghtmin
  real :: dw
  integer :: npbin(0:nbin),mynpbin(0:nbin)
  real :: wpbin(0:nbin),mywpbin(0:nbin)

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
       MPI_REAL8, &
       MPI_SUM, &
       MPI_COMM_WORLD,ierr)

  call MPI_ALLREDUCE(myavwe,avwe,1, &
       MPI_REAL8, &
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

  !      return
end subroutine countw
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine setw(ip,n)
  use gem_com
  use gem_equil
  implicit none
  INTEGER :: m,n,ip,i,j,k
  real :: wx0,wx1,wy0,wy1,wz0,wz1,vte
  real :: wght,r,bfldp,ter,b,vfac,th
  real :: xt,yt,zt

  vte = sqrt(amie*t0e(nr/2))
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
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1) 
     ter = wx0*t0e(i)+wx1*t0e(i+1)        
     b=1.-tor+tor*bfldp

     vfac = 0.5*(emass*u3e(m)**2 + 2.*mue3(m)*b)

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
     if((vfac/ter).gt.vcut)wght = 0.

     w000(m)=wx0*wy0*wz0*wght
     w001(m)=wx0*wy0*wz1*wght
     w010(m)=wx0*wy1*wz0*wght
     w011(m)=wx0*wy1*wz1*wght
     w100(m)=wx1*wy0*wz0*wght
     w101(m)=wx1*wy0*wz1*wght
     w110(m)=wx1*wy1*wz0*wght
     w111(m)=wx1*wy1*wz1*wght
  end do
  !      return
end subroutine setw

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine fltx(u,isbl,ism)   
  use gem_com
  use gem_fft_wrapper
  implicit none
  INTEGER :: i,j,k,k1,l,m,n,ifirst,nstep,ip,INFO,isbl,ism
  INTEGER :: l1,m1,myk,myj,ix,ikx,ierror
  real :: u(0:imx,0:jmx,0:1)
  complex :: temp3dxy(0:imx-1,0:jmx-1,0:1),v(0:imx-1),u1(0:imx-1)
  real :: filterx(0:imx-1,0:jmx-1),gr(0:imx-1),gi(0:imx-1)
  real :: kx,ky,sgny,dum,dum1,b2,myuav(0:imx-1)
  real :: myjaca(0:imx-1),jaca(0:imx-1),uav(0:imx-1)

  xshape=1.0
  yshape=1.0

  do j = 0,jmx-1
     do i = 0,imx-1
        if(j.ge.(jm/2+1)) then
           m1=jm-j
           sgny=-1.0
        else
           m1=j
           sgny=1.0
        endif
        ky=sgny*2.0*pi*float(m1)/ly
        kx = i*pi/lx
        b2 = 0. !xshape**2*kx**2+yshape**2*ky**2

        filterx(i,j) = 2.0/imx*exp(-b2**2)   !be careful with always filtering with hyper-Gaussian

        if(abs(ky)>kycut)filterx(i,j) = 0.0
        if(m1==0.and.ism==0)filterx(i,j) = 0.
        if(kx>kxcut)filterx(i,j) = 0.0
        if(onemd==1.and.m1.ne.mlk)filterx(i,j) = 0.0
        !            if(m1>1)filterx(i,j) = 0.0
        if(j==0.and.i<=nzcrt)filterx(i,j) = 0.0
        if(m1>0.and.m1<nlow)filterx(i,j) = 0.0!2.0/imx*float(m1*m1)/(nlow**2)
        if(j==0 .and. ineq0==0)filterx(i,j) = 0.0
     end do
  end do

  if(ism==2)then
     do i=0,im-1
        myuav(i) = 0.
        myjaca(i) = 0.
        do j=0,jm-1
           myuav(i)=myuav(i)+u(i,j,0)*jac(i,0)
           myjaca(i)=myjaca(i)+jac(i,0)
        enddo
     enddo
     call MPI_ALLREDUCE(myuav,uav,imx, &
          MPI_REAL8,                               &
          MPI_SUM,TUBE_COMM,ierr)
     call MPI_ALLREDUCE(myjaca,jaca,imx, &
          MPI_REAL8,                               &
          MPI_SUM,TUBE_COMM,ierr)

     uav=uav/jaca

     do i = 0,imx-1
        do j = 0,jmx
           do k = 0,1
              u(i,j,k) = u(i,j,k)-uav(i)
           end do
        end do
     end do
  end if

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
        call ccfft('y',-1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3dxy(i,j,k) = tmpy(j)   !rho(ky,x)
        end do
     end do
  enddo

  do k = 0,1
     do j = 0,jmx-1
        gr(:) = real(temp3dxy(:,j,k))
        gi(:) = aimag(temp3dxy(:,j,k))
        !       call dsinf(0,gr(0),1,0,gr(0),1,0,imx*2,1, 1.0,aux1,50000,aux2,20000)
        call dsinf(0,gr,1,0,1,0,imx*2,1,1.0,aux1,50000,aux2,20000)
        !       call dsinf(0,gi(0),1,0,gi(0),1,0,imx*2,1, 1.0,aux1,50000,aux2,20000)
        call dsinf(0,gi,1,0,1,0,imx*2,1,1.0,aux1,50000,aux2,20000)
        if(j==0.and.ism==0)then
           u1(:) = cmplx(gr(:),gi(:))
           call gam(u1(:),v(:))
           gr(:) = real(v(:))
           gi(:) = aimag(v(:))
        end if
        gr(:) = gr(:)*filterx(:,j)
        gi(:) = gi(:)*filterx(:,j)
        !       call dsinf(0,gr(0),1,0,gr(0),1,0,imx*2,1, 1.0,aux1,50000,aux2,20000)
        call dsinf(0,gr,1,0,1,0,imx*2,1,1.0,aux1,50000,aux2,20000)
        !       call dsinf(0,gi(0),1,0,gi(0),1,0,imx*2,1, 1.0,aux1,50000,aux2,20000)
        call dsinf(0,gi,1,0,1,0,imx*2,1,1.0,aux1,50000,aux2,20000)
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
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3dxy(i,j,k) = tmpy(j)  ! phi(x,y)
        end do
     end do
  end do

100 continue
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

  if(ism==2)then
     do i = 0,imx-1
        do j = 0,jmx
           do k = 0,1
              u(i,j,k) = u(i,j,k)+uav(i)
           end do
        end do
     end do
  end if

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

  !      return
end subroutine fltx
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine fltd(u)   
  use gem_com
  use gem_fft_wrapper
  implicit none
  INTEGER :: i,j,k,k1,l,m,n,ifirst,nstep,ip,INFO,isbl,ism
  INTEGER :: l1,m1,myk,myj,ix,ikx,ierror
  real :: u(0:imx,0:jmx,0:1)
  complex :: temp3dxy(0:imx-1,0:jmx-1,0:1),v(0:imx-1),u1(0:imx-1)
  real :: filterx(0:imx-1,0:jmx-1),gr(0:imx-1),gi(0:imx-1)
  real :: kx,ky,sgny,dum,dum1,b2,myuav(0:imx-1)
  real :: myjaca(0:imx-1),jaca(0:imx-1),uav(0:imx-1)

  xshape=1.0
  yshape=1.0

  do j = 0,jmx-1
     do i = 0,imx-1
        if(j.ge.(jm/2+1)) then
           m1=jm-j
           sgny=-1.0
        else
           m1=j
           sgny=1.0
        endif
        ky=sgny*2.0*pi*float(m1)/ly
        kx = i*pi/lx
        b2 = xshape**2*kx**2+yshape**2*ky**2

        filterx(i,j) = 2.0/imx*exp(-b2**2)   !be careful with always filtering with hyper-Gaussian

        if(abs(ky)>kycut)filterx(i,j) = 0.0
        if(kx>kxcut)filterx(i,j) = 0.0
        if(j==0.and.i<=nzcrt)filterx(i,j) = 0.0
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
        call ccfft('y',-1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3dxy(i,j,k) = tmpy(j)   !rho(ky,x)
        end do
     end do
  enddo

  do k = 0,1
     do j = 0,jmx-1
        gr(:) = real(temp3dxy(:,j,k))
        gi(:) = aimag(temp3dxy(:,j,k))
        !       call dsinf(0,gr(0),1,0,gr(0),1,0,imx*2,1, 1.0,aux1,50000,aux2,20000)
        call dsinf(0,gr,1,0,1,0,imx*2,1,1.0,aux1,50000,aux2,20000)
        !       call dsinf(0,gi(0),1,0,gi(0),1,0,imx*2,1, 1.0,aux1,50000,aux2,20000)
        call dsinf(0,gi,1,0,1,0,imx*2,1,1.0,aux1,50000,aux2,20000)
        if(j==0)then
           u1(:) = cmplx(gr(:),gi(:))
           call gam(u1(:),v(:))
           gr(:) = real(v(:))
           gi(:) = aimag(v(:))
        end if
        gr(:) = gr(:)*filterx(:,j)
        gi(:) = gi(:)*filterx(:,j)
        !       call dsinf(0,gr(0),1,0,gr(0),1,0,imx*2,1, 1.0,aux1,50000,aux2,20000)
        call dsinf(0,gr,1,0,1,0,imx*2,1,1.0,aux1,50000,aux2,20000)
        !       call dsinf(0,gi(0),1,0,gi(0),1,0,imx*2,1, 1.0,aux1,50000,aux2,20000)
        call dsinf(0,gi,1,0,1,0,imx*2,1,1.0,aux1,50000,aux2,20000)
        temp3dxy(:,j,k) = cmplx(gr(:),gi(:))
     end do
  end do

  temp3dxy(:,:,:)=temp3dxy(:,:,:)/jmx

  do k=0,mykm
     do i = 0,imx-1
        do j = 0,jmx-1
           tmpy(j) = temp3dxy(i,j,k)
        end do
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3dxy(i,j,k) = tmpy(j)  ! phi(x,y)
        end do
     end do
  end do

100 continue
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

  !      return
end subroutine fltd
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine dcmpy(u,v)
  use gem_com
  use gem_fft_wrapper
  implicit none
  INTEGER :: i,j,k,k1,l,m,n,ifirst,nstep,ip,INFO
  INTEGER :: l1,m1,myk,myj,ix,ikx,id
  INTEGER :: recvcnt(0:ntube-1)
  real :: u(0:imx-1,0:jmx-1,0:1)
  complex :: v(0:imx-1,0:jcnt-1,0:1),myv(0:imx-1,0:jcnt-1,0:1)
  complex :: sbuf(0:imx*jmx*2-1),rbuf(0:imx*jcnt*2-1)
  complex :: temp3d(0:imx-1,0:jmx-1,0:1)
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
        call ccfft('y',-1,jmx,1.0,tmpy,coefy,worky,0)
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

  !      return
end subroutine dcmpy
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
real function en3(s)
  real :: s

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
end function en3
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine blendf
  use gem_com
  use gem_equil
  implicit none

  INTEGER :: m,i,j,k,m1,m2
  complex :: dum,dum1
  complex :: tmp(nb,nb),work(100)
  integer :: IPIV(10),INFO
  real :: r,qr,s1,s2,s3,dth1,wx0,wx1
  real :: aky(0:jmx-1),dely(0:imx-1)

  dth1 = pi2/nb
  do j = 0,im-1
     r = rin+xg(j)
     i = int((r-rin)/dr)
     i = min(i,nr-1)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     qr = wx0*sf(i)+wx1*sf(i+1)
     dely(j) = mod(-pi2*lr0/q0*qr*sign(1.0,q0)+8000.*ly,ly)*tor
  enddo
  do j = 0,jm-1
     if(j.ge.(jm/2+1)) then
        aky(j) = -2.*pi*float(jm-j)/ly
     else
        aky(j) = 2.*pi*float(j)/ly
     end if
  enddo
  do j = 0,jm-1
     do i = 0,im-1
        pfac(i,j) = exp(IU*aky(j)*dely(i))
     end do
  enddo

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
                 dum = dum+conjg(pol(m1,i,j,k))*pol(m2,i,j,k)
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

11 format(12(1x,e10.3))

  !      return
end subroutine blendf
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine filtbl(u)   
  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none
  complex :: lbfs(0:imx-1,0:jmx-1)
  complex :: rbfr(0:imx-1,0:jmx-1)
  complex :: u(0:imx-1,0:jmx-1,0:1)
  complex :: bm(1:nb),cm(1:nb),dum
  real :: qr
  real :: aky(0:jmx-1),dely(0:imx-1),dklz(0:imx-1,0:jmx-1)
  INTEGER :: i,j,k,l,m,n,ind,m1,m2
  INTEGER :: myk,mynum,id
  complex,dimension(:),allocatable :: holdu,rbuf,sbuf

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
             mpi_double_complex,id,10, &
             rbuf(0),mynum, &
             mpi_double_complex,id,10, &
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
           dum = dum+conjg(pol(m1,i,j,k))*tmpz(k)
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
             mpi_double_complex,id,20, &
             rbuf(0),mynum, &
             mpi_double_complex,id,20, &
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
       MPI_double_complex,lngbr,10, &
       rbfr(0:imx-1,0:jmx-1),imx*jmx, &
       MPI_double_complex,rngbr,10, &
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
  !       return
end subroutine filtbl
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine ldel

  use gem_com
  use gem_equil
  implicit none
  INTEGER :: i,j,k,m,idum,ns,m1
  INTEGER :: np_old,np_new
  real :: vpar,vperp2,r,qr,th,b,cost,ter
  real :: avgv,myavgv,avgw,myavgw
  real :: dumx,dumy,dumz,jacp
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,psp
  real :: grp,gxdgyp,zoldp
  real :: wx0,wx1,wz0,wz1

  cnt=mme
  !   write(*,*)'in loader cnt, mm(1)=',cnt,mm(1)

  myavgv=0.
  avgv=0.
  avgw = 0.
  myavgw = 0.

  m = 0
  do j = 1,100000000

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
        x2e(m)=min(dumx,lx-1e-8)
        y2e(m)=min(dumy,ly-1e-8)

        k = int((th+pi)/dth)
        wz0 = (-pi+(k+1)*dth-th)/dth
        wz1 = 1-wz0
        z2e(m) = wz0*zfnth(k)+wz1*zfnth(k+1)
        z2e(m)=min(z2e(m),lz-1e-8)

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
        if(eldu==1)ter = tge
        b=1.-tor+tor*bfldp

        u2e(m)=vpar*sqrt(amie*ter)
        mue(m)=0.5*vperp2/b*ter
        eke(m) = mue(m)*b+0.5*emass*u2e(m)**2
        pze(m) = emass*u2e(m)/b+psp/br0
        z0e(m) = z2e(m)
        xie(m) = x2e(m)
        u0e(m) = u2e(m)
        myavgv=myavgv+u2e(m)

        !    LINEAR: perturb w(m) to get linear growth...
        w2e(m)= 1e-10 !2.*amp*(revers(MyId*cnt+m,13)-0.5) !(ran2(iseed) - 0.5 ) !don't use 0, to have finite apar at n=0 for contour plot
        !         w2e(m) = amp*u2e(m)/(abs(u2e(m))+0.1)*sin(x2e(m)*pi2/lx)
        myavgw=myavgw+w2e(m)
     end if
  enddo
170 continue
  myavgw = myavgw/mme
  !      write(*,*)'myid ', myid,x2(10),y2(20),u2(20),w2(20)
  !    subtract off avg. u...
  call MPI_ALLREDUCE(myavgv,avgv,1, &
       MPI_REAL8, &
       MPI_SUM,MPI_COMM_WORLD,ierr)
  if(idg.eq.1)write(*,*)'all reduce'
  avgv=avgv/float(tmm(1))
  do m=1,mme
     u2e(m)=u2e(m)-avgv
     x3e(m)=x2e(m)
     y3e(m)=y2e(m)
     z3e(m)=z2e(m)
     u3e(m)=u2e(m)
     w3e(m)=w2e(m)
     mue2(m) = mue(m)
     mue3(m) = mue(m)
  enddo

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

  !      return
end subroutine ldel
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine avge
  use gem_com
  use gem_equil
  implicit none
  INTEGER :: m,n,ip,i,j,k,j1,j2,j3,j4
  real :: bfldp,ter
  real :: wx0,wx1,wz0,wz1
  real :: r,qr,th,cost,sint,b,vfac,js,jv,dum,dum1,pidum,vel
  real :: xt,yt,zt,emac=10.,dvfac,dlamb,lamb,dely,g,myavewe,avwe,avwen
  integer :: isign,il,ie,bfcnt,mynobge
  real :: w(0:1,0:1,0:1,0:1)
  real :: mytotw(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1), &
       totw(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1), &
       wght(1:mmx)
  integer :: ng(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1), &
       myng(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1)
  complex :: cdum

  dvfac = emac/negrd
  dlamb = 2.0/nlgrd
  pidum = sqrt(pi)*pi2

  if(idg==1)write(*,*)'before 1 loop',myid
  mytotw = 0.
  myng = 0
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
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
     ter = wx0*t0e(i)+wx1*t0e(i+1)
     b=1.-tor+tor*bfldp
     vfac = 0.5*(emass*u3e(m)**2 + 2.*mue3(m)*b)
     vel = sqrt(2.*vfac/emass)
     lamb = u3e(m)/(vel+1.e-6)

     xt = x3e(m)
     i=int(xt/dx)
     i = min(i,im-1)
     j = int(y3e(m)/dy)
     j = min(j,jmx-1)
     ie = int((vfac/ter)/dvfac)
     il = int((lamb+1.0)/dlamb)
     il = min(il,nlgrd-1)
     if(ie.gt.(negrd-1))then
        goto 100
     end if

     mytotw(i,j,ie,il) = mytotw(i,j,ie,il)+w3e(m)

     myng(i,j,ie,il) = myng(i,j,ie,il)+1

100  continue
  end do

  bfcnt = imx*jmx*negrd*nlgrd
  if(idg==1)write(*,*)'before reduce'

  call MPI_ALLREDUCE(mytotw, totw, &
       bfcnt,MPI_REAL8,       &
       MPI_SUM,GRID_COMM,ierr)

  call MPI_ALLREDUCE(myng, ng, &
       bfcnt,MPI_INTEGER,       &
       MPI_SUM,GRID_COMM,ierr)

  if(idg==1)write(*,*)'before 2 loop',myid
  myavewe = 0.
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
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
     ter = wx0*t0e(i)+wx1*t0e(i+1)
     b=1.-tor+tor*bfldp
     vfac = 0.5*(emass*u3e(m)**2 + 2.*mue3(m)*b)
     vel = sqrt(2.*vfac/emass)
     lamb = u3e(m)/(vel+1.e-6)

     xt = x3e(m)
     i=int(xt/dx)
     i = min(i,im-1)
     j = int(y3e(m)/dy)
     j = min(j,jmx-1)
     ie = int((vfac/ter)/dvfac)
     il = int((lamb+1.0)/dlamb)
     il = min(il,nlgrd-1)
     if(ie.gt.(negrd-1))goto 200
     wght(m) = totw(i,j,ie,il)/ng(i,j,ie,il)
200  continue
     myavewe = myavewe+abs(wght(m))
  end do
  call MPI_ALLREDUCE(myavewe, avwe, &
       1            ,MPI_REAL8,       &
       MPI_SUM,MPI_COMM_WORLD,ierr)
  avwe = avwe/tmm(1)

  do m=1,mme
     w3e(m) = wght(m)
     w2e(m) = w3e(m)
  end do
  if(myid==0)write(*,*)'avwe',avwe

  !      return
end subroutine avge
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine avgi(ns)
  use gem_com
  use gem_equil
  implicit none
  INTEGER :: ns,m,n,ip,i,j,k,j1,j2,j3,j4
  real :: bfldp,ter
  real :: wx0,wx1,wz0,wz1
  real :: r,qr,th,cost,sint,b,vfac,js,jv,dum,dum1,pidum,vel
  real :: xt,yt,zt,emac=10.,dvfac,dlamb,lamb,dely,g,myavewe,avwi,avwen
  integer :: isign,il,ie,bfcnt,mynobge
  real :: w(0:1,0:1,0:1,0:1)
  real :: mytotw(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1), &
       totw(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1), &
       wght(1:mmx)
  integer :: ng(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1), &
       myng(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1)
  complex :: cdum

  dvfac = emac/negrd
  dlamb = 2.0/nlgrd
  pidum = sqrt(pi)*pi2

  if(idg==1)write(*,*)'before 1 loop',myid
  mytotw = 0.
  myng = 0
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
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
     ter = wx0*t0s(ns,i)+wx1*t0s(ns,i+1)
     b=1.-tor+tor*bfldp
     vfac = 0.5*(mims(ns)*u3(ns,m)**2 + 2.*mu(ns,m)*b)
     vel = sqrt(2.*vfac/mims(ns))
     lamb = u3(ns,m)/(vel+1.e-6)

     xt = x3(ns,m)
     i=int(xt/dx)
     i = min(i,im-1)
     j = int(y3(ns,m)/dy)
     j = min(j,jmx-1)
     ie = int((vfac/ter)/dvfac)
     il = int((lamb+1.0)/dlamb)
     il = min(il,nlgrd-1)
     if(ie.gt.(negrd-1))then
        goto 100
     end if

     mytotw(i,j,ie,il) = mytotw(i,j,ie,il)+w3(ns,m)

     myng(i,j,ie,il) = myng(i,j,ie,il)+1

100  continue
  end do

  bfcnt = imx*jmx*negrd*nlgrd
  if(idg==1)write(*,*)'before reduce'

  call MPI_ALLREDUCE(mytotw, totw, &
       bfcnt,MPI_REAL8,       &
       MPI_SUM,GRID_COMM,ierr)

  call MPI_ALLREDUCE(myng, ng, &
       bfcnt,MPI_INTEGER,       &
       MPI_SUM,GRID_COMM,ierr)

  if(idg==1)write(*,*)'before 2 loop',myid
  myavewe = 0.
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
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
     ter = wx0*t0s(ns,i)+wx1*t0s(ns,i+1)
     b=1.-tor+tor*bfldp
     vfac = 0.5*(mims(ns)*u3(ns,m)**2 + 2.*mu(ns,m)*b)
     vel = sqrt(2.*vfac/mims(ns))
     lamb = u3(ns,m)/(vel+1.e-6)

     xt = x3(ns,m)
     i=int(xt/dx)
     i = min(i,im-1)
     j = int(y3(ns,m)/dy)
     j = min(j,jmx-1)
     ie = int((vfac/ter)/dvfac)
     il = int((lamb+1.0)/dlamb)
     il = min(il,nlgrd-1)
     if(ie.gt.(negrd-1))goto 200
     wght(m) = totw(i,j,ie,il)/ng(i,j,ie,il)
200  continue
     myavewe = myavewe+abs(wght(m))
  end do
  call MPI_ALLREDUCE(myavewe, avwi, &
       1                ,MPI_REAL8,       &
       MPI_SUM,MPI_COMM_WORLD,ierr)
  avwi = avwi/tmm(1)

  do m=1,mme
     w3(ns,m) = wght(m)
     w2(ns,m) = w3(ns,m)
  end do
  !      if(myid==0)write(*,*)'avwi',avwi

  !      return
end subroutine avgi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rstpe
  use gem_com
  use gem_equil
  implicit none
  INTEGER :: m,n,ip,i,j,k,j1,j2,j3,j4,j5,nobge
  real :: wx0,wx1,wz0,wz1,bfldp,ter,gdum
  real :: favx=1,favy=1,favz=1,fave=0,favl=1
  real :: wx(0:1),wy(0:1),wz(0:1),we(0:1),wp(0:1),vte
  real :: r,qr,th,jacp,b,vfac,js,jv,dum,dum1,dum2,pidum,vel,jfnp
  real :: xt,yt,zt,emac=10.,dvfac,dlamb,lamb,dely,g,myavewe,avwe
  integer :: isign,il,ie,bfcnt,mynobge
  real :: w(0:1,0:1,0:1,0:1,0:1)
  real :: mytotw(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1), &
       totw(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1)
  real :: wght(1:mmx),y(1:mmx)
  real :: myg(0:imx,0:jmx,0:1,0:negrd,0:nlgrd), &
       tg(0:imx,0:jmx,0:1,0:negrd,0:nlgrd)
  integer :: ng(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1), &
       myng(0:imx-1,0:jmx-1,0:negrd-1,0:nlgrd-1)
  real ::    h(0:imx,0:jmx,0:1,0:negrd,0:nlgrd), &
       h1(0:imx,0:jmx,0:1,0:negrd,0:nlgrd), &
       lbfs(0:imx,0:jmx,0:negrd,0:nlgrd), &
       rbfs(0:imx,0:jmx,0:negrd,0:nlgrd), &
       lbfr(0:imx,0:jmx,0:negrd,0:nlgrd), &
       rbfr(0:imx,0:jmx,0:negrd,0:nlgrd)
  real :: myncell(0:imx-1,0:jmx-1),ncell(0:imx-1,0:jmx-1)  
  real :: mydwcell(0:imx-1,0:jmx-1),dwcell(0:imx-1,0:jmx-1)  
  real :: mydwpcell(0:imx-1,0:jmx-1),dwpcell(0:imx-1,0:jmx-1)  
  real :: mydwecell(0:imx-1,0:jmx-1),dwecell(0:imx-1,0:jmx-1)  
  real :: myavwcell(0:imx-1,0:jmx-1),avwcell(0:imx-1,0:jmx-1)  
  real :: dwn(0:imx-1,0:jmx-1),dwe(0:imx-1,0:jmx-1)  
  real :: dwp(0:imx-1,0:jmx-1)

  dvfac = emac/negrd
  dlamb = 2.0/nlgrd
  pidum = sqrt(pi)*pi2
  vte = sqrt(amie*t0e(nr/2))

  if(idg==1)write(*,*)'before 1 loop',myid
  h = 0.
  mytotw = 0.
  mynobge = 0
  myg = 0.
  myng = 0
  myncell = 0.
  myavwcell = 0.

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
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
     jfnp = wz0*jfn(k)+wz1*jfn(k+1)
     jacp = wx0*wz0*jacob(i,k)+wx0*wz1*jacob(i,k+1) &
          +wx1*wz0*jacob(i+1,k)+wx1*wz1*jacob(i+1,k+1)
     jacp = jacp*jfnp
     ter = wx0*t0e(i)+wx1*t0e(i+1)
     b=1.-tor+tor*bfldp
     vfac = 0.5*(emass*u3e(m)**2 + 2.*mue3(m)*b)/ter
     vel = sqrt(2.*vfac*ter/emass)
     lamb = u3e(m)/(vel+1.e-6)

     xt = x3e(m)
     i=int(xt/dx)
     i = min(i,im-1)
     yt = y3e(m)
     j = int(y3e(m)/dy)
     j = min(j,jmx-1)
     myncell(i,j) = myncell(i,j)+1
     myavwcell(i,j) = myavwcell(i,j)+abs(w3e(m))
     k=0
     wx(0)=1.-favx+(float(i+1)-xt/dx)*favx
     wx(1)=1.-wx(0)
     wy(0)=1.-favy+(float(j+1)-yt/dy)*favy
     wy(1)=1.-wy(0)
     wz(0)=1.-favz+(float(gclr*kcnt+k+1)-z3e(m)/dz)*favz
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
     jv = sqrt(2.)*vte**3*sqrt(vfac+1.e-3)
     dum = totvol/tmm(1)/(dx*dy*dz*dvfac*dlamb*js*jv)

     myng(i,j,ie,il) = myng(i,j,ie,il)+1
     mytotw(i,j,ie,il) = mytotw(i,j,ie,il)+w3e(m)

     do j1 = 0,1
        do j2 = 0,1
           do j3 = 0,1
              do j4 = 0,1
                 do j5 = 0,1
                    h(i+j1,j+j2,j3,ie+j4,il+j5) =  &
                         h(i+j1,j+j2,j3,ie+j4,il+j5) + &
                         dum*w(j1,j2,j3,j4,j5)*w3e(m)

                    myg(i+j1,j+j2,j3,ie+j4,il+j5) =  &
                         myg(i+j1,j+j2,j3,ie+j4,il+j5) + &
                         dum*w(j1,j2,j3,j4,j5)
                 end do
              end do
           end do
        end do
     end do
100  continue
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
  do m=1,mme
     wght(m) = w3e(m)
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
     bfldp = wx0*wz0*bfld(i,k)+wx0*wz1*bfld(i,k+1) &
          +wx1*wz0*bfld(i+1,k)+wx1*wz1*bfld(i+1,k+1)
     jfnp = wz0*jfn(k)+wz1*jfn(k+1)
     jacp = wx0*wz0*jacob(i,k)+wx0*wz1*jacob(i,k+1) &
          +wx1*wz0*jacob(i+1,k)+wx1*wz1*jacob(i+1,k+1)
     jacp = jacp*jfnp
     ter = wx0*t0e(i)+wx1*t0e(i+1)
     b=1.-tor+tor*bfldp
     vfac = 0.5*(emass*u3e(m)**2 + 2.*mue3(m)*b)/ter
     vel = sqrt(2.*vfac*ter/emass)
     lamb = u3e(m)/(vel+1.e-6)

     xt = x3e(m)
     i=int(xt/dx)
     i = min(i,im-1)
     yt = y3e(m)
     j = int(y3e(m)/dy)
     j = min(j,jmx-1)
     k=0
     wx(0)=1.-favx+(float(i+1)-xt/dx)*favx
     wx(1)=1.-wx(0)
     wy(0)=1.-favy+(float(j+1)-yt/dy)*favy
     wy(1)=1.-wy(0)
     wz(0)=1.-favz+(float(gclr*kcnt+k+1)-z3e(m)/dz)*favz
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

200  continue
     myavewe = myavewe+abs(wght(m))
     mydwcell(i,j) = mydwcell(i,j)+eprs*(w3e(m) -wght(m))
     mydwecell(i,j) = mydwecell(i,j)+eprs*(w3e(m) -wght(m))*vfac
     mydwpcell(i,j) = mydwpcell(i,j)+eprs*(w3e(m) -wght(m))*u3e(m)
     w3e(m) = w3e(m)*(1.-eprs)+wght(m)*eprs
     w2e(m) = w3e(m)
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

  dum = 0.
  dum1 = 0.
  dum2 = 0.
  avwe = 0.
  do i = 0,imx-1
     do j = 0,jmx-1
        dum = dum+dwcell(i,j)
        dum1 = dum1+dwecell(i,j)
        dum2 = dum2+dwpcell(i,j)
        avwe = avwe+avwcell(i,j)
     end do
  end do

  adwn = abs(dum)/avwe
  adwe = abs(dum1)/(avwe*1)
  adwp = abs(dum2)/(avwe*sqrt(amie))

  !      if(myid==0)write(*,*)'nobge, avwe= ',nobge,avwe

  !      return
end subroutine rstpe
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gam(u,v)
  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none
  !     
  complex :: u(0:imx-1),v(0:imx-1)
  INTEGER :: n,i,j,k,k1,i1

  call ccfft('z',0,kmx,1.0,tmpz,coefz,workz,0)

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
        call ccfft('z',1,kmx,1.0,tmpz,coefz,workz,0)
     end if

     call MPI_BCAST(tmpz,kmx,MPI_DOUBLE_COMPLEX,master, &
          tube_comm,ierr)
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
     do k = 0,kmx-1
        k1 = k
        if(k>kmx/2)k1=kmx-k
        if(k1>0)tmpz(k) = 0.
     end do
     call ccfft('z',-1,kmx,1.0,tmpz,coefz,workz,0)
     tmpz = tmpz/kmx

     v(i1) = tmpz(gclr)
  end do

  !      return
end subroutine gam
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine gkpsL(nstep,ip)   
  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none
  INTEGER :: ns
  real :: lbfr(0:imx,0:jmx)
  real :: lbfs(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
  real,DIMENSION(1:5) :: b,ter,gam0,gam1
  real :: b2,delyz,th,bf,gamphi
  real :: kx,ky,wx0,wz0,wx1,wz1
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
  real :: r,grp,gxdgyp,grdgtp,gthp
  real,dimension(:),allocatable :: akx,aky
  real,dimension(:,:),allocatable :: akx2
  real,dimension(:,:,:),allocatable :: formphi
  real,dimension(:,:,:),allocatable :: formfe
  real :: sgnx,sgny,sz,myfe
  INTEGER :: i,j,k,i1,k1,l,m,n,ifirst,nstep,ip
  INTEGER :: l1,m1,myk
  complex :: temp3d(0:imx-1,0:jmx-1,0:1),tempdv(0:imx-1,0:jmx-1,0:1)
  complex :: aphik(0:imx-1),myaphik(0:imx-1)
  real :: myrmsphi
  real :: myaph(0:imx),aph(0:imx)
  real :: myden(0:imx,0:jmx,0:1),rho1(0:imx,0:jmx,0:1)

  save formphi,formfe,ifirst,akx,aky,akx2

  !     form factors....
  if (ifirst.ne.-99) then
     ALLOCATE( akx(0:imx-1),aky(0:jmx-1))
     ALLOCATE( formphi(0:imx-1,0:jmx-1,0:kmx),akx2(0:imx-1,0:kmx))
     ALLOCATE( formfe(0:imx-1,0:jmx-1,0:kmx))

     do l=0,im-1
        do m=0,jm-1
           do n = 0,km
              r = lr0
              i1 = int((r-rin)/dr)
              i1 = min(i1,nr-1)
              wx0 = (rin+(i1+1)*dr-r)/dr
              wx1 = 1.-wx0
              do ns = 1, nsm
                 ter(ns) = wx0*t0s(ns,i1)+wx1*t0s(ns,i1+1)
              enddo
              k = n
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
              radiusp = wx0*wz0*radius(i1,k)+wx0*wz1*radius(i1,k+1) &
                   +wx1*wz0*radius(i1+1,k)+wx1*wz1*radius(i1+1,k+1) 

              if(l.ge.(im/2+1)) then
                 l1=im-l
                 sgnx=-1.
              else
                 l1=l
                 sgnx=1.
              endif
              if(m.ge.(jm/2+1)) then
                 m1=jm-m
                 sgny=-1.
              else
                 m1=m
                 sgny=1.
              endif
              akx(l) = sgnx*2.*pi*float(l1)/lx
              aky(m) = sgny*2.*pi*float(m1)/ly
              ky=sgny*2.*pi*float(m1)/ly
              kx=sgnx*2.*pi*float(l1)/lx

              bf=bfldp
              b2=(xshape*xshape*kx*kx + yshape*yshape*ky*ky)
              b2=b2**c4*(1-onemd)
              do ns = 1, nsm
                 b(ns)=mims(ns)*(kx*kx*grp**2 + ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                      +2*dydrp*lr0/q0*qhatp*grdgtp+(lr0/q0/radiusp)**2)+2.*kx*ky*gxdgyp)/(bf*bf) &
                      *ter(ns)/(q(ns)*q(ns))
                 call srcbes(b(ns),gam0(ns),gam1(ns))
              enddo

              !   formfactor in gkps
              formfe(l,m,n) = 1.-gam0(1)
              if(b(1).lt.3.e-5)then
                 formphi(l,m,n) = 0.
              else if(b(1).gt.bcut)then
                 formphi(l,m,n) = 0.
              else
                 gamphi = 0.
                 do ns = 1, nsm
                    gamphi = gamphi + q(ns)*cn0s(ns)*gn0s(ns,l)/ter(ns)*(1-gam0(ns))
                 enddo
                 formphi(l,m,n)=exp(-b2)/(fradi*xn0e(nr/2)*cn0e/t0e(nr/2)+gamphi) &
                      /float(imx*jmx)
              end if
              if(l1.eq.0.and.m1.eq.0)formphi(l,m,n)=0.
              if(m1>0.and.m1<nlow)formphi(l,m,n)=0.
              if(abs(ky).gt.kycut)formphi(l,m,n) = 0.
              if(abs(kx).gt.kxcut)formphi(l,m,n) = 0.
              if(m1.ne.mlk.and.onemd.eq.1)formphi(l,m,n) = 0.
              !                  if(l1.ne.llk.and.onemd.eq.1)formphi(l,m,n) = 0.
              if(ky.eq.0.)then
                 if(l1.eq.0)then
                    akx2(l,n) = 0.
                 else
                    akx2(l,n) = exp(-b2)/(1.-gam0(1))/float(km)
                 end if
              end if
           enddo
        enddo
     enddo

     ifirst=-99

  endif

  !   now do field solve...

  !      phi = 0.
  !      return
  temp3d = 0.
  aphik = 0.
  myaphik = 0.

  myden = rho+den0apa+fradi*(phi/ntube*xn0e(nr/2)*cn0e/t0e(nr/2)-den0)
  call MPI_ALLREDUCE(myden(0:im,0:jm,0:1),  &
       rho1(0:im,0:jm,0:1),             &
       (imx+1)*(jmx+1)*2,MPI_REAL8,       &
       MPI_SUM,GRID_COMM,ierr)

  !  find rho(kx,ky)
  do k=0,mykm
     n = GCLR*kcnt+k
     do j=0,jm-1
        do i=0,im-1
           temp3d(i,j,k)=rho1(i,j,k)
        enddo
     enddo

     do i = 0,imx-1
        do j = 0,jmx-1
           tmpy(j) = temp3d(i,j,k)
        end do
        call ccfft('y',-1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3d(i,j,k) = tmpy(j)   !rho(ky,x)
        end do
     end do

     do j = 0,jmx-1
        do i = 0,imx-1
           tmpx(i) = temp3d(i,j,k)
        end do
        call ccfft('x',-1,imx,1.0,tmpx,coefx,workx,0)
        do i = 0,imx-1
           temp3d(i,j,k) = tmpx(i)   ! rho(kx,ky)
        end do
     end do
  enddo

  !   find aphik(k_x) from rho(kx,ky=0,z)
  if(iadi.eq.1)then
     do k = 0,mykm-1
        do i = 0,imx-1
           myk=GCLR*kcnt+k         
           myaphik(i) = myaphik(i)+temp3d(i,0,k)*akx2(i,myk)
        end do
     end do
     call MPI_ALLREDUCE(myaphik(0:imx-1),aphik(0:imx-1),  &
          imx,MPI_DOUBLE_COMPLEX,MPI_SUM,TUBE_COMM,ierr)
  end if

  !  from rho(kx,ky) to phi(kx,ky)
  do k=0,mykm
     myk=GCLR*kcnt+k         
     do j=1,jm-1
        do i=0,im-1
           temp3d(i,j,k)=temp3d(i,j,k)*formphi(i,j,myk)
        enddo
     enddo
     do i = 0,im-1
        temp3d(i,0,k) = (temp3d(i,0,k)+aphik(i))*formphi(i,0,myk)
     end do
  enddo

  !      if(myid.eq.0)write(*,*)ip,abs(temp3d(1,2,0)),abs(temp3d(5,6,0))
  !      call filtor(temp3d(0:imx-1,0:jmx-1,0:1))
  if(idg==1)write(*,*)'pass filtor', myid

  ! calculate fe(nstep)
  myfe=0.
  do k=0,mykm-1
     myk=GCLR*kcnt+k
     do j=0,jm-1
        do i=0,im-1
           myfe=myfe+formfe(i,j,myk)*abs(temp3d(i,j,k))**2 !cabs
        enddo
     enddo
  enddo
  call MPI_ALLREDUCE(myfe,fe(nstep),1,MPI_REAL8, &
       MPI_SUM,TUBE_COMM,ierr)     

  !  from phi(kx,ky) to phi(x,y)
  do k=0,mykm
     do i = 0,imx-1
        do j = 0,jmx-1
           tmpy(j) = temp3d(i,j,k)
        end do
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3d(i,j,k) = tmpy(j)  ! phi(x,y)
        end do
     end do
  end do

  !      call filtbl(temp3d(0:imx-1,0:jmx-1,0:1))  !need be very careful. What's een is not what's expected.

  do k=0,mykm
     n = GCLR*kcnt+k
     do j = 0,jmx-1
        do i = 0,imx-1
           tmpx(i) = temp3d(i,j,k)
        end do
        call ccfft('x',1,imx,1.0,tmpx,coefx,workx,0)

        do i = 0,imx-1
           temp3d(i,j,k) = tmpx(i)  ! phi(ky,x)
        end do
     end do
  end do

100 continue
  do i = 0,im-1
     do j = 0,jm-1
        do k = 0,mykm
           phi(i,j,k) = temp3d(i,j,k)
        end do
     end do
  end do

  !    x-y boundary points 
  call enfxy(phi(:,:,:))
  call enfz(phi(:,:,:))
  !      call filter(phi(:,:,:))

  if(izonal.eq.1)return
  do i=0,im
     myaph(i)=0.
     aph(i)=0.
  enddo
  do k=0,mykm-1
     do j=0,jm-1
        do i=0,im
           myaph(i)=myaph(i)+phi(i,j,k) 
        enddo
     enddo
  enddo

  call MPI_ALLREDUCE(myaph(0:imx),aph(0:imx),imx+1,  &
       MPI_REAL8,MPI_SUM,TUBE_COMM,ierr)
  aph = aph/(jmx*kmx)
  do i = 0,imx
     do j = 0,jmx
        do k = 0,mykm
           phi(i,j,k) = phi(i,j,k)-aph(i)
        end do
     end do
  end do
  !      return
end subroutine gkpsL

!      End of gkps....
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine ezampL(nstep,ip)   

  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none
  real :: lbfr(0:imx,0:jmx)
  real :: lbfs(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)

  real :: b,b2,gam0,gam1,th,delyz,bf,rkper
  real :: kx,ky,wx0,wz0,wx1,wz1,ter,terc
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
  real :: r,grp,gxdgyp,grdgtp,gthp
  real,DIMENSION(:),ALLOCATABLE :: akx,aky
  real,DIMENSION(:,:,:), ALLOCATABLE :: formapa
  real :: sgnx,sgny,sz,xfrac=1.0
  INTEGER :: i,j,k,i1,k1,l,m,n,ifirs,ip,nstep
  INTEGER :: l1,m1,myk
  complex :: temp3d(0:imx-1,0:jmx-1,0:1)
  real :: myrmsapa
  real :: myden(0:imx,0:jmx,0:1),jtmp(0:imx,0:jmx,0:1)

  save formapa,ifirs,akx,aky

  !     form factors....
  if (ifirs.ne.-99) then
     !     
     ALLOCATE(akx(0:imx-1),aky(0:jmx-1))
     ALLOCATE(formapa(0:imx-1,0:jmx-1,0:kmx))

     formapa = 0.
     do l=0,im-1
        do m=0,jm-1
           do n = 0,km
              r = lr0
              i1 = int((r-rin)/dr)
              i1 = min(i1,nr-1)
              wx0 = (rin+(i1+1)*dr-r)/dr
              wx1 = 1.-wx0
              ter = wx0*t0i(i1)+wx1*t0i(i1+1)
              terc = wx0*t0c(i1)+wx1*t0c(i1+1)

              k = n
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

              if(l.ge.(im/2+1)) then
                 l1=im-l
                 sgnx=-1.
              else
                 l1=l
                 sgnx=1.
              endif
              if(m.ge.(jm/2+1)) then
                 m1=jm-m
                 sgny=-1.
              else
                 m1=m
                 sgny=1.
              endif
              akx(l) = sgnx*2.*pi*float(l1)/lx
              aky(m) = sgny*2.*pi*float(m1)/ly
              ky=sgny*2.*pi*float(m1)/ly
              kx=sgnx*2.*pi*float(l1)/lx

              bf=bfldp
              b=(kx*kx*grp**2 + ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                   +2*dydrp*lr0/q0*qhatp*grdgtp)+2.*kx*ky*gxdgyp)
              b2=(xshape*xshape*kx*kx + yshape*yshape*ky*ky)
              b2=0. !b2**c4*(1-onemd)

              !                  call srcbes(b,gam0,gam1)
              !  formfactor in ezamp
              formapa(l,m,n)=exp(-b2)*beta/((b+xfrac*beta*amie*gn0e(l)*cn0e+beta*q(1)*q(1)/mims(1)*gn0s(1,l)*cn0s(1)*isiap+1.e-10) &
                   *float(imx*jmx))
              if(l1.eq.0.and.m1.eq.0)formapa(l,m,n)=0.
              if(m1>0.and.m1<nlow)formapa(l,m,n)=0.
              if(abs(ky).gt.kycut)formapa(l,m,n) = 0.
              if(abs(kx).gt.kxcut)formapa(l,m,n) = 0.
              if(m1.ne.mlk.and.onemd.eq.1)formapa(l,m,n) = 0.
              !                  if(l1.ne.llk.and.onemd.eq.1)formapa(l,m,n) = 0.
              !                  if(b.gt.bcut)formapa(l,m,n) = 0.
           enddo
        enddo
     enddo

     ifirs=-99

  endif

  !   now do field solve...

  temp3d = 0.
  myden = jion*ision-upar-amie*upa00
  call MPI_ALLREDUCE(myden(0:im,0:jm,0:1),  &
       jtmp(0:im,0:jm,0:1),             &
       (imx+1)*(jmx+1)*2,MPI_REAL8,       &
       MPI_SUM,GRID_COMM,ierr)

  !  find jtot(kx,ky)
  do k=0,mykm
     n = GCLR*kcnt+k
     do j=0,jm-1
        do i=0,im-1
           temp3d(i,j,k)=jtmp(i,j,k)+amie*apar(i,j,k)*gn0e(i)*cn0e
        enddo
     enddo

     do i = 0,imx-1
        do j = 0,jmx-1
           tmpy(j) = temp3d(i,j,k)
        end do
        call ccfft('y',-1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3d(i,j,k) = tmpy(j)   !jtot(ky,x)
        end do
     end do

     do j = 0,jmx-1
        do i = 0,imx-1
           tmpx(i) = temp3d(i,j,k)
        end do
        call ccfft('x',-1,imx,1.0,tmpx,coefx,workx,0)
        do i = 0,imx-1
           temp3d(i,j,k) = tmpx(i)   ! jtot(kx,ky)
        end do
     end do
  enddo

  !  from jtot(kx,ky) to apar(kx,ky)
  do k=0,mykm
     myk=GCLR*kcnt+k         
     do j=0,jm-1
        do i=0,im-1
           temp3d(i,j,k)=temp3d(i,j,k)*formapa(i,j,myk)
        enddo
     enddo
  enddo
  !      if(myid.eq.0)write(*,*)abs(temp3d(1,2,0)),abs(temp3d(5,6,0))
  !      call filtor(temp3d(0:imx-1,0:jmx-1,0:1))

  !  from apar(kx,ky) to apar(x,y)
  do k=0,mykm
     n = GCLR*kcnt+k
     do j = 0,jmx-1
        do i = 0,imx-1
           tmpx(i) = temp3d(i,j,k)
        end do
        call ccfft('x',1,imx,1.0,tmpx,coefx,workx,0)
        do i = 0,imx-1
           temp3d(i,j,k) = tmpx(i)  ! j(ky,x)
        end do
     end do

     do i = 0,imx-1
        do j = 0,jmx-1
           tmpy(j) = temp3d(i,j,k)
        end do
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3d(i,j,k) = tmpy(j)  ! apar(x,y)
        end do
     end do
  end do

  do i = 0,im-1
     do j = 0,jm-1
        do k = 0,mykm
           apar(i,j,k) = temp3d(i,j,k)
        end do
     end do
  end do

  !    x-y boundary points 
  call enfxy(apar(:,:,:))
  call enfz(apar(:,:,:))
  !      call filter(apar(:,:,:))

  !      return
end subroutine ezampL
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine dpdtL(ip)   

  use gem_com
  use gem_equil
  use gem_fft_wrapper

  implicit none

  INTEGER :: ns
  real :: lbfr(0:imx,0:jmx)
  real :: lbfs(0:imx,0:jmx)
  real :: rbfr(0:imx,0:jmx)
  real :: rbfs(0:imx,0:jmx)
  real,DIMENSION(1:5) :: b,ter,gam0,gam1
  real :: b2,delyz,th,bf,gamphi
  real :: kx,ky,wx0,wz0,wx1,wz1
  real :: dbdrp,dbdtp,grcgtp,bfldp,fp,radiusp,dydrp,qhatp,psipp,jfnp
  real :: r,grp,gxdgyp,grdgtp,gthp
  real,dimension(:),allocatable :: akx,aky
  real,dimension(:,:,:),allocatable :: formphi
  real,dimension(:,:),allocatable :: akx2
  real :: sgnx,sgny,sz,myfe
  INTEGER :: i,j,k,i1,k1,l,m,n,ifirst,nstep,ip
  INTEGER :: l1,m1,myk
  complex :: temp3d(0:imx-1,0:jmx-1,0:1)
  real :: myaph(0:imx),aph(0:imx)
  real :: myden(0:imx,0:jmx,0:1)

  save formphi,ifirst,akx,aky,akx2

  !     form factors....
  if (ifirst.ne.-99) then
     allocate(akx(0:imx-1),aky(0:jmx-1))
     allocate(formphi(0:imx-1,0:jmx-1,0:kmx),akx2(0:imx-1,0:kmx))

     do l=0,im-1
        do m=0,jm-1
           do n = 0,km
              r = lr0
              i1 = int((r-rin)/dr)
              i1 = min(i1,nr-1)
              wx0 = (rin+(i1+1)*dr-r)/dr
              wx1 = 1.-wx0
              do ns = 1, nsm
                 ter(ns) = wx0*t0s(ns,i1)+wx1*t0s(ns,i1+1)
              enddo
              k = n
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

              if(l.ge.(im/2+1)) then
                 l1=im-l
                 sgnx=-1.
              else
                 l1=l
                 sgnx=1.
              endif
              if(m.ge.(jm/2+1)) then
                 m1=jm-m
                 sgny=-1.
              else
                 m1=m
                 sgny=1.
              endif
              akx(l) = sgnx*2.*pi*float(l1)/lx
              aky(m) = sgny*2.*pi*float(m1)/ly
              ky=sgny*2.*pi*float(m1)/ly
              kx=sgnx*2.*pi*float(l1)/lx

              bf=bfldp
              b2=(xshape*xshape*kx*kx + yshape*yshape*ky*ky)
              b2=b2**c4*(1-onemd)
              do ns = 1, nsm
                 b(ns)=mims(ns)*(kx*kx*grp**2 + ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                      +2*dydrp*lr0/q0*qhatp*grdgtp)+2.*kx*ky*gxdgyp)/(bf*bf) &
                      *ter(ns)/(q(ns)*q(ns))
                 call srcbes(b(ns),gam0(ns),gam1(ns))
              enddo

              !   formfactor in dpdt
              if(b(1).lt.3.e-5)then
                 formphi(l,m,n) = 0.
              else if(b(1).gt.bcut)then
                 formphi(l,m,n) = 0.
              else
                 gamphi = 0.
                 do ns = 1, nsm
                    gamphi = gamphi + q(ns)*cn0s(ns)*gn0s(ns,l)/ter(ns)*(1-gam0(ns))
                 enddo
                 formphi(l,m,n)=exp(-b2)/(0.+gamphi) &
                      /float(imx*jmx)
              end if
              if(m1>0.and.m1<nlow)formphi(l,m,n)=0.
              if(abs(ky).gt.kycut)formphi(l,m,n) = 0.
              if(abs(kx).gt.kxcut)formphi(l,m,n) = 0.
              if(m1.ne.mlk.and.onemd.eq.1)formphi(l,m,n) = 0.
              !                  if(l1.ne.llk.and.onemd.eq.1)formphi(l,m,n) = 0.

           enddo
        enddo
     enddo

     ifirst=-99

  endif

  !   now do field solve...

  temp3d = 0.

  myden = drhodt
  call MPI_ALLREDUCE(myden(0:im,0:jm,0:1),  &
       drhodt(0:im,0:jm,0:1),             &
       (imx+1)*(jmx+1)*2,MPI_REAL8,       &
       MPI_SUM,GRID_COMM,ierr)

  !  find rho(kx,ky)
  do k=0,mykm
     n = GCLR*kcnt+k
     do j=0,jm-1
        do i=0,im-1
           temp3d(i,j,k)=drhodt(i,j,k)
        enddo
     enddo

     do i = 0,imx-1
        do j = 0,jmx-1
           tmpy(j) = temp3d(i,j,k)
        end do
        call ccfft('y',-1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3d(i,j,k) = tmpy(j)   !rho(ky,x)
        end do
     end do

     do j = 0,jmx-1
        do i = 0,imx-1
           tmpx(i) = temp3d(i,j,k)
        end do
        call ccfft('x',-1,imx,1.0,tmpx,coefx,workx,0)
        do i = 0,imx-1
           temp3d(i,j,k) = tmpx(i)   ! rho(kx,ky)
        end do
     end do
  enddo

  !  from rho(kx,ky) to phi(kx,ky)
  do k=0,mykm
     myk=GCLR*kcnt+k         
     do j=0,jm-1
        do i=0,im-1
           temp3d(i,j,k)=temp3d(i,j,k)*formphi(i,j,myk)
        enddo
     enddo
  enddo
  !      call filtor(temp3d(0:imx-1,0:jmx-1,0:1))
  !  from phi(kx,ky) to phi(x,y)
  do k=0,mykm
     n = GCLR*kcnt+k
     do j = 0,jmx-1
        do i = 0,imx-1
           tmpx(i) = temp3d(i,j,k)
        end do
        call ccfft('x',1,imx,1.0,tmpx,coefx,workx,0)

        do i = 0,imx-1
           temp3d(i,j,k) = tmpx(i)  ! phi(ky,x)
        end do
     end do

     do i = 0,imx-1
        do j = 0,jmx-1
           tmpy(j) = temp3d(i,j,k)
        end do
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3d(i,j,k) = tmpy(j)  ! phi(x,y)
        end do
     end do
  end do

100 continue
  do i = 0,im-1
     do j = 0,jm-1
        do k = 0,mykm
           dphidt(i,j,k) = temp3d(i,j,k)
        end do
     end do
  end do

  if(izonal.eq.1)goto 200
  do i=0,im
     myaph(i)=0.
     aph(i)=0.
  enddo
  do k=0,mykm-1
     do j=0,jm-1
        do i=0,im
           myaph(i)=myaph(i)+dphidt(i,j,k) 
        enddo
     enddo
  enddo
  !
  call MPI_ALLREDUCE(myaph(0:imx),aph(0:imx),imx+1, &
       MPI_REAL8, &
       MPI_SUM, &
       TUBE_COMM,ierr)
  aph = aph/(jmx*kmx)
  do i = 0,imx
     do j = 0,jmx
        do k = 0,mykm
           dphidt(i,j,k) = dphidt(i,j,k)-aph(i)
        end do
     end do
  end do

200 continue
  !    x-y boundary points 
  call enfxy(dphidt(:,:,:))
  call enfz(dphidt)
  !      if(idpbf==0)call filter(dphidt(:,:,:))
  if(idpbf==1)call flty(dphidt(:,:,:),1)

  !      return
end subroutine dpdtL

!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine ftcamp
  use gem_com
  use gem_fft_wrapper
  implicit none

  complex :: mycampf(0:mynf-1)
  real :: aomega(0:nfreq-1)
  integer :: i,j,k,thisf,nsize,ir
  real :: domega,om,dum1,dum2,dum3,dum4,dum5,x,gam,peak(0:6,1:5)
  complex :: sbuf(0:mynf-1),rbuf(0:nfreq-1)

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
     gam = log(dum2/dum1)/(nsize-400)
     do i = 0,nsize-1
        !         camp(:,i) = camp(:,i)/exp(gam*i)
        !         camp(:,i) = exp(IU*1.5e-3*dt*ifskp*i)*exp(gam*i)  !test FT effect
     end do

     do j = 0,nfreq-1
        om = aomega(j)
        campf(ir,j) = 0.
        do i = 0,nsize-1
           campf(ir,j) = campf(ir,j)+camp(ir,i)*exp(-IU*om*dt*ifskp*i)
        end do
     end do

     !find 5 peaks
     dum1 = 0.
     do i=0,nfreq-1
        x = abs(campf(ir,i))
        if(x>dum1)then
           dum1 = x
        end if
     end do
     peak(ir,1)=dum1
     dum2 = 0.
     do i=0,nfreq-1
        x = abs(campf(ir,i))
        if(x>dum2.and. x .ne. dum1)then
           dum2 = x
        end if
     end do
     peak(ir,2)=dum2
     dum3 = 0.
     do i=0,nfreq-1
        x = abs(campf(ir,i))
        if(x>dum3 .and. x.ne.dum1 .and. x.ne.dum2)then
           dum3 = x
        end if
     end do
     peak(ir,3)=dum3
     dum4 = 0.
     do i=0,nfreq-1
        x = abs(campf(ir,i))
        if(x>dum4 .and. x.ne.dum1 .and. x.ne.dum2 .and. x.ne.dum3)then
           dum4 = x
        end if
     end do
     peak(ir,4)=dum4
     dum5 = 0.
     do i=0,nfreq-1
        x = abs(campf(ir,i))
        if(x>dum5 .and. x.ne.dum1 .and. x.ne.dum2 .and. x.ne.dum3 .and. x.ne.dum4)then
           dum5 = x
        end if
     end do
     peak(ir,5)=dum5
  end do

  do k = 0,glst
     if(tclr==0 .and. gclr==k)then
        open(15, file='freq', status='unknown',position='append')
        do ir=0,6
           do i = 0,nfreq-1
              x = abs(campf(ir,i))
              if(x==peak(ir,1)) write(15,10)i,aomega(i),abs(campf(ir,i))**2
              if(x==peak(ir,2)) write(15,10)i,aomega(i),abs(campf(ir,i))**2
              if(x==peak(ir,3)) write(15,10)i,aomega(i),abs(campf(ir,i))**2
              if(x==peak(ir,4)) write(15,10)i,aomega(i),abs(campf(ir,i))**2
              if(x==peak(ir,5)) write(15,10)i,aomega(i),abs(campf(ir,i))**2
           end do
           write(15,*)'frequency, gclr, ir=', gclr,ir
           do i = 0,nfreq-1
              write(15,10)i,aomega(i),abs(campf(ir,i))**2
           end do
           write(15,*)'end ir=',ir
        end do
        close(15)
     end if
     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
10   format(1x,i6, 3(2x,e12.5))
  end do
  close(15)
  !      return
end subroutine ftcamp
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine colli(ip,n)

  use gem_com
  use gem_equil

  implicit none

  integer :: i,ip,k,m,n,ncol,icol
  real :: edum,vdum,dum,dum1,ptch,vti,r,qr,th,cost,b
  real :: h_x,h_coll,x,eps,dtcol,uold,hee,nue,ter
  real :: wx0,wx1,wz0,wz1

  ncol = 1
  if(ip.eq.1)dtcol = dt/ncol*2
  if(ip.eq.0)dtcol = dt/ncol
  if(rneui==0.0)return
  do k = 1,mm(1)
     r=x3(1,k)-0.5*lx+lr0

     m = int(z3(1,k)/delz)
     wz0 = ((m+1)*delz-z3(1,k))/delz
     wz1 = 1-wz0
     th = wz0*thfnz(m)+wz1*thfnz(m+1)
     i = int((r-rin)/dr)
     wx0 = (rin+(i+1)*dr-r)/dr
     wx1 = 1.-wx0
     m = int((th+pi)/dth)
     wz0 = (-pi+(m+1)*dth-th)/dth
     wz1 = 1.-wz0
     ter = wx0*t0i(i)+wx1*t0i(i+1)        
     b = wx0*wz0*bfld(i,m)+wx0*wz1*bfld(i,m+1) &
          +wx1*wz0*bfld(i+1,m)+wx1*wz1*bfld(i+1,m+1) 
     uold = u3(1,k)
     edum = b*mu(1,k)+0.5*mims(1)*u3(1,k)*u3(1,k)
     vdum = sqrt(2.*edum/mims(1))
     ptch = u3(1,k)/vdum
     eps = edum/ter
     x = sqrt(eps)
     h_x    = 4.0*x*x/(3.0*sqrt(pi))
     h_coll = h_x/sqrt(1.0+h_x**2)
     !         h_coll = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*DERF(x)
     ! collision frequency for experimental profiles
     hee = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*erf(x)
     nue=rneui
     dum1 = 1/(eps+0.1)**1.5  
     !         dum = dtcol*rneu*dum1
     dum = dtcol*nue*dum1
     !         if(x<0.3)dum=0.0
     do icol = 1,ncol
        ptch =ptch-dum*ptch+sqrt(dum*(1.-ptch*ptch)) &
             *sign(1.0,ran2(iseed)-0.5)
        ptch = min(ptch,0.999)
        ptch = max(ptch,-0.999)
     end do
     u3(1,k) = vdum*ptch
     mu(1,k) = 0.5*mims(1)*vdum*vdum*(1.-ptch*ptch)/b
  end do
  !      return
end subroutine colli
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine flty(u,isbl)   

  use gem_com
  use gem_equil
  use gem_fft_wrapper

  implicit none

  INTEGER :: i,j,k,k1,l,m,n,ifirst,nstep,ip,INFO,isbl,ism
  INTEGER :: l1,m1,myk,myj,ix,ikx,ierror
  real :: u(0:imx,0:jmx,0:1)
  complex :: temp3d(0:imx-1,0:jmx-1,0:1),v(0:imx-1),u1(0:imx-1)
  real :: filterx(0:imx-1,0:jmx-1)
  real :: kx,ky,sgnx,sgny,dum,dum1,b2,myuav(0:imx-1)
  real :: myjaca(0:imx-1),jaca(0:imx-1),uav(0:imx-1)

  do j = 0,jmx-1
     do i = 0,imx-1
        if(i.ge.(im/2+1)) then
           l1=im-i
           sgnx=-1.
        else
           l1=i
           sgnx=1.
        endif
        if(j.ge.(jm/2+1)) then
           m1=jm-j
           sgny=-1.
        else
           m1=j
           sgny=1.
        endif

        ky=sgny*2.0*pi*float(m1)/ly
        kx=sgnx*2.0*pi*float(l1)/lx
        b2 = xshape**2*kx**2+yshape**2*ky**2

        filterx(i,j) = 1.0/(imx*jmx)

        if(abs(ky)>kycut)filterx(i,j) = 0.0
        if(abs(kx)>kxcut)filterx(i,j) = 0.0
        if(onemd==1.and.m1.ne.mlk)filterx(i,j) = 0.0
        !            if(j==0.and.i<=nzcrt)filterx(i,j) = 0.0
        !            if(m1>0.and.m1<nlow)filterx(i,j) = 0.0 
     end do
  end do

  do k=0,mykm
     n = GCLR*kcnt+k
     do j=0,jm-1
        do i=0,im-1
           temp3d(i,j,k)=u(i,j,k)
        enddo
     enddo

     do i = 0,imx-1
        do j = 0,jmx-1
           tmpy(j) = temp3d(i,j,k)
        end do
        call ccfft('y',-1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3d(i,j,k) = tmpy(j)   !rho(ky,x)
        end do
     end do

     do j = 0,jmx-1
        do i = 0,imx-1
           tmpx(i) = temp3d(i,j,k)
        end do
        call ccfft('x',-1,imx,1.0,tmpx,coefx,workx,0)
        do i = 0,imx-1
           temp3d(i,j,k) = tmpx(i)   ! rho(kx,ky)
        end do
     end do
  enddo

  do k = 0,1
     do j = 0,jmx-1
        do i = 0,imx-1
           temp3d(i,j,k) = temp3d(i,j,k)*filterx(i,j)
        end do
     end do
  end do

  do k=0,mykm
     do j = 0,jmx-1
        do i = 0,imx-1
           tmpx(i) = temp3d(i,j,k)
        end do
        call ccfft('x',1,imx,1.0,tmpx,coefx,workx,0)

        do i = 0,imx-1
           temp3d(i,j,k) = tmpx(i)  ! phi(ky,x)
        end do
     end do
  end do

  if(isbl==1)call filtbl(temp3d(:,:,:))
  do k = 0,mykm
     do i = 0,imx-1
        do j = 0,jmx-1
           tmpy(j) = temp3d(i,j,k)
        end do
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3d(i,j,k) = tmpy(j)  ! phi(x,y)
        end do
     end do
  end do

100 continue
  do i = 0,im-1
     do j = 0,jm-1
        do k = 0,mykm
           u(i,j,k) = temp3d(i,j,k)
        end do
     end do
  end do

  !    x-y boundary points 
  call enfxy(u(:,:,:))
  call enfz(u(:,:,:))

  !      return
end subroutine flty

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine mdamp(u,n)

  !     calculate mode histories. calculate modes for u at timestep n
  !     and store in modehis(mode,n).

  use gem_com
  use gem_equil
  use gem_fft_wrapper

  implicit none

  real :: u(0:imx,0:jmx,0:1)
  INTEGER :: n,i,j,k,m
  complex :: tmp3d(0:imx,0:jmx,0:1),cdum

  call ccfft('z',0,kmx,1.0,tmpz,coefz,workz,0)

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
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
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

  !      return
end subroutine mdamp

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mdampa(u,n)

  !     calculate mode histories. calculate modes for u at timestep n
  !     and store in modehis(mode,n).

  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none
  !     
  real :: u(0:imx,0:jmx,0:1)
  INTEGER :: n,i,j,k,m
  complex :: tmp3d(0:imx,0:jmx,0:1),cdum

  call ccfft('z',0,kmx,1.0,tmpz,coefz,workz,0)

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
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
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

  !      return
end subroutine mdampa

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine mdampd(u,v)

  !     calculate mode histories. calculate modes for u at timestep n
  !     and store in modehis(mode,n).

  use gem_com
  use gem_equil
  use gem_fft_wrapper

  implicit none
  !     
  real :: u(0:imx,0:jmx,0:1),v(0:100)
  INTEGER :: n,i,j,k,m
  complex :: tmp3d(0:imx,0:jmx,0:1),cdum

  call ccfft('z',0,kmx,1.0,tmpz,coefz,workz,0)

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
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
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

end subroutine mdampd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zon(u,v)
  use gem_com
  use gem_fft_wrapper
  implicit none
  INTEGER :: i,j,k,k1,l,m,n
  real :: u(0:imx,0:jmx,0:1),v(0:imx)
  real :: myjac(0:imx),mydbr(0:imx),totjac(0:imx)
  real :: dum,dum1

  do i = 0,nxpp-1
     dum = 0.
     dum1 = 0.
     do j = 0,jm-1
        dum = dum+u(i,j,0)*jac(i,0)
        dum1 = dum1+jac(i,0)
     end do
     mydbr(i) = dum
     myjac(i) = dum1
  end do
  call MPI_ALLREDUCE(mydbr,v,nxpp,  &
       MPI_REAL8,MPI_SUM,           &
       tube_comm,ierr)
  call MPI_ALLREDUCE(myjac,totjac,nxpp,  &
       MPI_REAL8,MPI_SUM,           &
       tube_comm,ierr)
  v(0:imx-1) = v(0:imx-1)/totjac(0:imx-1)

  !      return
end subroutine zon
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine joule
  use gem_com
  use gem_fft_wrapper
  implicit none
  INTEGER :: i,j,k,k1,l,m,n
  real :: v(0:imx),upaneq0(0:imx,0:jmx,0:1),upaneq8(0:imx,0:jmx,0:1)
  real :: exneq1(0:imx,0:jmx,0:1),eyneq1(0:imx,0:jmx,0:1)
  real :: exneq2(0:imx,0:jmx,0:1),eyneq2(0:imx,0:jmx,0:1)
  real :: exneq3(0:imx,0:jmx,0:1),eyneq3(0:imx,0:jmx,0:1)
  real :: uptot(0:imx,0:jmx,0:1),dbyneq0(0:imx,0:jmx,0:1),dbxneq0(0:imx,0:jmx,0:1)
  real :: myjac(0:imx),mydbr(0:imx),totjac(0:imx)
  real :: dum,dum1,dum2

  uptot = upart+amie*upa0t
  !      call fltx(uptot,0,1)
  call neq8(1,ex,exneq1)
  call neq8(1,ey,eyneq1)
  call neq8(2,ex,exneq2)
  call neq8(2,ey,eyneq2)
  call neq8(3,ex,exneq3)
  call neq8(3,ey,eyneq3)
  !      call neq8(uptot,upaneq8)
  !      call zon(upart+amie*upa0t,v)
  dum = 0.
  dum1 = 0.
  dum2 = 0.
  do i = 0,nxpp-1
     do j = 0,jm-1
        dum = dum+((delby(i,j,0)*ggx(i,0))**2 + (delbx(i,j,0))**2*ggy2(i,0)-2*delby(i,j,0)*delbx(i,j,0)*ggxdgy(i,0)) * jac(i,0)
        dum1 = dum1+uptot(i,j,0)*bdgxcgy(i,0)*(delbx(i,j,0)*exneq1(i,j,0)+delby(i,j,0)*eyneq1(i,j,0)) * jac(i,0)
        !            dum2 = dum2+upaneq0(i,j,0)*ezneq0(i,j,0)*bdgrzn(i,0)*jac(i,0)
        dum2 = dum2+uptot(i,j,0)*bdgxcgy(i,0)*(delbx(i,j,0)*exneq2(i,j,0)+delby(i,j,0)*eyneq2(i,j,0)) * jac(i,0)
     end do
  end do
  call MPI_ALLREDUCE(dum,tot_field_e,1,  &
       MPI_REAL8,MPI_SUM,           &
       tube_comm,ierr)
  call MPI_ALLREDUCE(dum1,tot_joule,1,  &
       MPI_REAL8,MPI_SUM,           &
       tube_comm,ierr)

  call MPI_ALLREDUCE(dum2,tot_joule1,1,  &
       MPI_REAL8,MPI_SUM,           &
       tube_comm,ierr)

end subroutine joule
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine neq0(u,v)   
  use gem_com
  use gem_fft_wrapper
  implicit none
  INTEGER :: i,j,k,k1,l,m,n,ifirst,nstep
  real :: u(0:imx,0:jmx,0:1),v(0:imx,0:jmx,0:1)
  real :: dum,dum1

  do k = 0,1
     do i = 0,imx
        dum = 0.
        do j = 0,jmx-1
           dum = dum+u(i,j,k)
        end do
        dum = dum/jmx
        v(i,:,k) = dum
     end do
  end do

  call enfxy(v(:,:,:))
  call enfz(v)

end subroutine neq0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine neq8(n,u,v)
   
  use gem_com
  use gem_fft_wrapper

  implicit none

  INTEGER :: i,j,k,k1,l,m,n,ifirst,nstep
  complex :: temp3dxy(0:imx-1,0:jmx-1,0:1)
  real :: u(0:imx,0:jmx,0:1),v(0:imx,0:jmx,0:1)
  real :: dum,dum1

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
        call ccfft('y',-1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3dxy(i,j,k) = tmpy(j)   !rho(ky,x)
        end do
     end do
  enddo

  do j = 0,jmx-1
     if(j .ne. n .and. j .ne. jmx-n) temp3dxy(:,j,:) = 0.
  end do

  do k=0,mykm
     do i = 0,imx-1
        do j = 0,jmx-1
           tmpy(j) = temp3dxy(i,j,k)
        end do
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
        do j = 0,jmx-1
           temp3dxy(i,j,k) = tmpy(j)  ! phi(x,y)
        end do
     end do
  end do

  do i = 0,imx-1
     do j = 0,jm-1
        do k = 0,mykm
           v(i,j,k) = temp3dxy(i,j,k)/jmx
        end do
     end do
  end do

  call enfxy(v(:,:,:))
  call enfz(v)

end subroutine neq8
