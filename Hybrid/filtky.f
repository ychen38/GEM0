!       3D Flux Tube Toroidal Electromagnetic GK Delta-f Code
!   global variables...

         use gem_com
         use fft_wrapper
	implicit none
	integer :: n,i,j,k,ip

	call initialize
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
         end do
	 lasttm=MPI_WTIME()
	 tottm=lasttm-starttm
!	 write(*,*)'ps time=',pstm,'tot time=',tottm
	 call MPI_FINALIZE(ierr)
         end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hybinit

      use gem_com
      implicit none
      GCLR = int(MyId/ntube)
      GLST = numprocs/ntube-1
      TCLR = mod(MyId,ntube)

!***MPI variables

!      if (GCLR.eq.GLST) 
!     %     mykm=km-GLST*mykm
      mykm = 1
      rngbr = GCLR+1
      if (GCLR.eq.GLST)rngbr = 0
      
      lngbr = GCLR-1
      if (GCLR.eq.0) lngbr = GLST
!     
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine init
      
      use gem_com
      implicit none
      character*(62) dumchar
      INTEGER :: i,j,k,n,ns,idum
      INTEGER :: mm1,lr1
      REAL(8) :: mims1,tets1,q1,kappan,kappat,r,qr,th,cost

      IU=cmplx(0.,1.)
      open(115,file='gem.in')
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) dumchar
      read(115,*) imx,jmx,kmx,mmx,nmx,nsmx,modemx,ntube
      read(115,*) dumchar
      read(115,*) im,jm,km,lx,ly,lz,ncar,xgplot
      read(115,*) dumchar
      read(115,*) dt,nm,nsm,xshape,yshape,zshape
      read(115,*) dumchar
      read(115,*) kymode,iput,iget,ision,peritr,llk,mlk,onemd,izon
      read(115,*) dumchar
      read(115,*) nplot,xnplt,modem,contu,pskip,wmax,imovie
      read(115,*) dumchar
      read(115,*) cut,amp,tor,ishift,fradi,kxcut,kycut,bcut
      read(115,*) dumchar
      read(115,*) br0,lr0,q0,qp,width
      read(115,*) dumchar
      read(115,*) c1,c2,c3,c4,ifluid,isg,amie,rneu
      read(115,*) dumchar
      read(115,*) beta,nonlin,nonline,vwidth,vwidthe,vcut,isuni,idg
	call new_gem_com()
      if(iget.eq.0)then
         if(myid.eq.master)open(9, file='plot', status='unknown')
         if(myid.eq.master)open(11, file='flux', status='unknown')
         if(myid.eq.master)open(15, file='zonal_amp', status='unknown')
      end if
      if(iget.eq.1)then
         if(myid.eq.master)open(9, file='plot', & 
             status='unknown',position='append')
         if(myid.eq.master)open(11, file='flux', &
             status='unknown',position='append')
         if(myid.eq.master)open(15, file='zonal_amp', &
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

!      mykm=int(km/numprocs)
!      if (MyId.eq.Last) mykm=km-Last*mykm
!     write(*,*)'mykm = ', mykm
!     
!     begin reading species info, ns=1,nsm...
      if(nsm.le.0) write(*,*)'invalid nsm',nsm
      read(115,*) dumchar
      ptr(1)=1
      do  200 ns=1,nsm
         read(115,*) dumchar
         read(115,*) mm1,tets1,mims1,q1,lr1
         read(115,*) dumchar
         read(115,*) kappan,kappat,kapt(2)
         tmm(ns)=mm1
         mm(ns)=int(mm1/numprocs)
         if (MyId.eq.Last) mm(ns)=mm1-Last*mm(ns)
!     write(*,*)'in init  ',Myid,mm(ns)
         tets(ns)=tets1
         mims(ns)=mims1
         q(ns)=q1
         lr(ns)=lr1
         kapn(ns)=kappan
         kapt(ns)=kappat
         if (ns.lt.nsm) ptr(ns+1)=ptr(ns)+mm(ns)
 200  continue
      tets(2) = 1.
      mims(2) = 1./(amie-1)
      lr(2) = 1
      q(2) = -1.
      kapn(2) = kapn(1)
!      kapt(2) = 0.!kapt(1)

      if(iget.eq.1) amp=0.
      pi=4.0*atan(1.0)
      pi2 = pi*2.
!     totvol is the square for now...
      dx=lx/float(im)
      dy=ly/float(jm)
      dz=lz/float(km)
      totvol=lx*ly*lz
      n0=float(tmm(1))/totvol
      e0=lr0/q0/br0
!     
      do 10 i=0,im
         xg(i)=dx*float(i)
 10   continue
      do 12 j=0,jm
         yg(j)=dy*float(j)
 12   continue
      kcnt=1
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
      do i = 0,imx
         r = xg(i)-0.5*lx+lr0
         qr = q0+qp*(r-lr0)
         qfld(i) = qr
         do k = 0,mykm
            th = (zg(k)-0.5*lz)/(br0*q0)
            cost = cos(th)
            bfld(i,k) = 1.-lr0/br0*cost*tor
            rmaj(i,k) = br0+r*cost
            dzp(i,k) = dz*rmaj(i,k)*qfld(i)/br0/q0
            if(tor.lt.0.1)dzp(i,k)=dz
            cfx(i,k) = sin(th)*tor
            cfy(i,k) = lr0*qr/r/q0*cos(th)+lr0/q0*qp*th*sin(th)*tor
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
         write(9,*)'dt,beta= ',dt, beta
         write(9,*)'kapt,kapn,kpte= ', kappat, kappan,kapt(2)
         write(9,*)'amp = ',amp
         write(9,*)'peritr,ifluid= ',peritr,ifluid
         write(9,*)'tor,nonlin= ',tor,nonlin
         write(9,*)'isuni= ',isuni, 'amie= ',amie
         write(9,*)'kxcut,kycut,bcut= ',kxcut,kycut,bcut
         write(9,*)'fradi,isg= ',fradi,isg
         write(9,*)'llk,mlk,onemd =',llk,mlk,onemd
         write(9,*)'vwidth (i,e), vcut= ', vwidth,vwidthe,vcut
         write(9,*)'mm1= ',mm1
      end if

      if(myid.eq.master)then
         write(*,*)'dt,beta= ',dt, beta
         write(*,*)'kapt,kapn,kpte= ', kappat, kappan,kapt(2)
         write(*,*)'amp = ',amp
         write(*,*)'peritr,ifluid= ',peritr,ifluid
         write(*,*)'tor,nonlin= ',tor,nonlin
         write(*,*)'isuni= ',isuni, 'amie= ',amie
         write(*,*)'kxcut,kycut,bcut= ',kxcut,kycut,bcut
         write(*,*)'fradi,isg= ',fradi,isg
         write(*,*)'llk,mlk,onemd =',llk,mlk,onemd
         write(*,*)'vwidth (i,e), vcut= ', vwidth,vwidthe,vcut
         write(*,*)'mm1= ',mm1
      end if
      close(115)
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ppush(n)

      use gem_com
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,dum,vxdum
      INTEGER :: m,i,j,k,l,n
      INTEGER :: np_old,np_new
      REAL(8) :: rhog,vfac,kap,vpar,pidum
      REAL(8) :: b,th,r,b,enerb,cost,sint,qr,laps,sz
      REAL(8) :: xt,xs,yt,xdot,ydot,zdot,pzdot,edot,pzd0

      pidum = 1./(pi*2)**1.5*vwidth**3
      do m=1,mm(1)
         r=x2(m)-0.5*lx+lr0
         qr=q0+qp*(r-lr0)
         th=(z2(m)-0.5*lz)/(br0*q0)
         cost=cos(th)
         sint=sin(th)
         b=1.-lr0/br0*cost*tor

         rhog=sqrt(2.*b*mu(m)*mims(1))/(q(1)*b)
         sz=(z2(m)-0.5*lz)*e0*qp/q0*tor !e0=lr0/q0/br0


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
            xs=x2(m)+rwx(1,l)*rhog
            yt=y2(m)+(rwy(1,l)+sz*rwx(1,l))*rhog
!
!   particle can go out of bounds during gyroavg...
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)

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
         vfac = 0.5*(mims(1)*u2(m)**2 + 2.*mu(m)*b)
         kap = kapn(1) - (1.5-vfac)*kapt(1)

         vpar = u2(m)-q(1)/mims(1)*aparp*nonlin
         enerb=(mu(m)+mims(1)*vpar*vpar/b)/q(1)*tor
         vxdum = (eyp+vpar/b*delbxp)
         xdot = vxdum*nonlin -tor*enerb*sint/br0
         ydot = (-exp1+vpar/b*delbyp)*nonlin     &
             -tor*b*enerb*(cost/br0+e0*qp*th*sint)
         zdot =  vpar !-enerb*q0/r*cost
!     %        -vpar*q0*br0/lr0*(-delbyp+lr0/q0*qp*th*delbxp)*tor
!     %        -1./b*q0*br0/r*(exp1+qp*th*eyp)*tor

         pzd0 = tor*(-b*mu(m)/mims(1)*lr0*sint/(q0*br0*br0))
         pzdot = pzd0 + (q(1)/mims(1)*ezp  &
             +q(1)/mims(1)*(-xdot*delbyp+ydot*delbxp+zdot*dadzp) &
             +tor*(-mu(m)/mims(1)/b*(delbyp*sint/br0-delbxp/br0  &
                                     *(cost+lr0/q0*qp*th*sint)   &
                                     -dadzp*q0/r*cost) )         &
             -tor*vpar/b/br0*(exp1*sint+(cost+lr0/q0*qp*th*sint)*eyp  &
                               +q0*br0/r*cost*ezp))*nonlin
         edot = xdot*exp1+ydot*eyp+zdot*ezp                      &
             -q(1)/mims(1)*mu(m)*r/qr/br0/br0*sint*aparp*tor     &
             +q(1)*vpar*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)
!         edot = mu(m)*(-cost/br0*vxdum+r/br0*sint*
!     %               (zdot-vpar+enerb*q0/r*cost )/q0/br0)*tor
!     %        +mims(1)*u2(m)*(pzdot-pzd0)

         x3(m) = x1(m) + 2.*dt*xdot
         y3(m) = y1(m) + 2.*dt*ydot
         z3(m) = z1(m) + 2.*dt*zdot
         u3(m) = u1(m) + 2.*dt*pzdot

         dum = 1-w2(m)*nonlin*0.
         if(isuni.eq.1)dum = exp(-vfac)*pidum
         vxdum = eyp+vpar/b*delbxp
!     %        -(2*u2(m)*aparp+q(1)/mims(1)*aparp*aparp)/(b*b)/br0*sint
         w3(m)=w1(m) + 2.*dt*(vxdum*kap + edot)*dum
         
!         if(x3(m)>lx .or. x3(m)<0.)w3(m) = 0.

         u1(m)=u2(m)  + .25*( u3(m) - u1(m) )
         x1(m)=x2(m)  + .25*( x3(m) - x1(m) )
         y1(m)=y2(m)  + .25*( y3(m) - y1(m) )
         z1(m)=z2(m)  + .25*( z3(m) - z1(m) )
         w1(m)=w2(m)  + .25*( w3(m) - w1(m) )

         laps=anint((z3(m)/lz)-.5)*(1-peritr)
         r=x3(m)-0.5*lx+lr0
         qr=q0+qp*(r-lr0)
         y3(m)=dmod(y3(m)-laps*2*pi*qr*lr0/q0+80.*ly,ly)
         z3(m)=dmod(z3(m)+8.*lz,lz)
         x3(m)=dmod(x3(m)+8.*lx,lx)         

      enddo

      np_old=mm(1)
      call init_pmove(z3,np_old,lz,ierr)

      call pmove(x1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mu,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit


      call end_pmove(ierr)
      mm(1)=np_new

      return
      end

!-----------------------------------------------------------------------

      subroutine cpush(n)

      use gem_com
      implicit none
      INTEGER :: n
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1
      INTEGER :: m,i,j,k,ns,l
      INTEGER :: np_old,np_new
      REAL(8) :: vfac,vpar,vxdum,dum,xdot,ydot,zdot,pzdot,edot,pzd0
      REAL(8) :: rhog,xt,yt,zt,kap,xs,pidum
      REAL(8) :: b,th,r,b,enerb,cost,sint,qr,laps,sz
      REAL(8) :: myke,mypfl,myavewi
      REAL(8) :: myefl,mynos
      REAL(8) :: ketemp,pfltemp
      REAL(8) :: efltemp,nostemp
      REAL(8) :: sbuf(10),rbuf(10)

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
         qr=q0+qp*(r-lr0)
         th=(z3(m)-0.5*lz)/(br0*q0)
         cost=cos(th)
         sint=sin(th)
         b=1.-lr0/br0*cost*tor

         rhog=sqrt(2.*b*mu(m)*mims(1))/(q(1)*b)
         sz=(z3(m)-0.5*lz)*e0*qp/q0*tor !e0=lr0/q0/br0

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
            xs=x3(m)+rwx(1,l)*rhog
            yt=y3(m)+(rwy(1,l)+sz*rwx(1,l))*rhog
!   BOUNDARY
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)

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

         vfac = 0.5*(mims(1)*u3(m)**2 + 2.*mu(m)*b)
         kap = kapn(1) -(1.5-vfac)*kapt(1)

         vpar = u3(m)-q(1)/mims(1)*aparp*nonlin
         enerb=(mu(m)+mims(1)*vpar*vpar/b)/q(1)*tor
         vxdum = (eyp+vpar/b*delbxp)

         xdot = vxdum*nonlin-tor*enerb*sint/br0
         ydot = (-exp1+vpar/b*delbyp)*nonlin  &
                 -tor*b*enerb*(cost/br0+e0*qp*th*sint)
         zdot =  vpar !-enerb*q0/r*cost
!     %        -vpar*q0*br0/lr0*(-delbyp+lr0/q0*qp*th*delbxp)*tor
!     %        -1./b*q0*br0/r*(exp1+qp*th*eyp)*tor

         pzd0 = tor*(-b*mu(m)/mims(1)*lr0*sint/(q0*br0*br0))
         pzdot = pzd0 + (q(1)/mims(1)*ezp                         &
             +q(1)/mims(1)*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)  &
             +tor*(-mu(m)/mims(1)/b*(delbyp*sint/br0-delbxp/br0   &
                                     *(cost+lr0/q0*qp*th*sint)    &
                                     -dadzp*q0/r*cost) )          &
             -tor*vpar/b/br0*(exp1*sint+(cost+lr0/q0*qp*th*sint)*eyp   &
                               +q0*br0/r*cost*ezp))*nonlin

         edot = xdot*exp1+ydot*eyp+zdot*ezp                      &
             -q(1)/mims(1)*mu(m)*r/qr/br0/br0*sint*aparp*tor     &
             +q(1)*vpar*(-xdot*delbyp+ydot*delbxp+zdot*dadzp)
!         edot = mu(m)*(-cost/br0*vxdum+r/br0*sint*
!     %                 (zdot-vpar+enerb*q0/r*cost )/q0/br0)
!     %        +mims(1)*u2(m)*(pzdot-pzd0)

         x3(m) = x1(m) + 0.5*dt*xdot
         y3(m) = y1(m) + 0.5*dt*ydot
         z3(m) = z1(m) + 0.5*dt*zdot
         u3(m) = u1(m) + 0.5*dt*pzdot

         dum = 1-w3(m)*nonlin*0.
         if(isuni.eq.1)dum = exp(-vfac)*pidum
         vxdum = eyp+vpar/b*delbxp
!     %        -(2*u3(m)*aparp+q(1)/mims(1)*aparp*aparp)/(b*b)/br0*sint
         w3(m) = w1(m) + 0.5*dt*(vxdum*kap+edot)*dum 

         if(abs(w3(m)).gt.1..and.nonlin==1)then
            w3(m) = 0.
            w2(m) = 0.
            w1(m) = 0.
         end if

!         if(x3(m)>lx .or. x3(m)<0.)w3(m) = 0.
         laps=anint((z3(m)/lz)-.5)*(1-peritr)
         r=x3(m)-0.5*lx+lr0
         qr=q0+qp*(r-lr0)
         y3(m)=dmod(y3(m)-laps*2*pi*qr*lr0/q0+80.*ly,ly)
         z3(m)=dmod(z3(m)+8.*lz,lz)
         x3(m)=dmod(x3(m)+8.*lx,lx)         

!     particle diagnostics done here because info is available...
         mypfl=mypfl + w3(m)*(eyp+vpar*delbxp/b) 
         myefl=myefl + vfac*w3(m)*(eyp/b+vpar*delbxp/b)
         myke=myke + vfac*w3(m)
         mynos=mynos + w3(m)
         myavewi = myavewi+abs(w3(m))

!     xn becomes xn-1...
         u1(m)=u2(m)
         x1(m)=x2(m)
         y1(m)=y2(m)
         z1(m)=z2(m)
         w1(m)=w2(m)
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
      efl(1,n)=mims(1)/tets(1)*efltemp/( float(tmm(1)) )
      ke(1,n)=ketemp/( 2.*float(tmm(1))*mims(1) )
      np_old=mm(1) 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call init_pmove(z3,np_old,lz,ierr)
!     
      call pmove(x1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mu,np_old,np_new,ierr)
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
      implicit none
      INTEGER :: i,j,k,ip
      real(8) :: ux(0:imx,0:jmx,0:1),uy(0:imx,0:jmx,0:1)
      real(8) :: tmp(0:imx,0:jmx,0:1),uoverb(0:imx,0:jmx,0:1)

      call gradu(phi(ip,:,:,:),ux,uy)
      ex(ip,:,:,:) = -ux(:,:,:)
      ey(ip,:,:,:) = -uy(:,:,:)

      delbx = 0.
      delby = 0.
      if(ifluid.eq.1)then
         call gradu(apar(ip,:,:,:),ux,uy)
         delbx(ip,:,:,:) = uy(:,:,:)
         delby(ip,:,:,:) = -ux(:,:,:)
      end if

      call gradu(dene(ip,:,:,:),ux,uy)
      dnedx(:,:,:) = ux(:,:,:)
      dnedy(:,:,:) = uy(:,:,:)
      do i = 0,imx
         do j = 0,jmx
            do k = 0,1
               uoverb(i,j,k) = upar(ip,i,j,k)/bfld(i,k)
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
      real(8) :: mydene(0:imx,0:jmx,0:1),myupar(0:imx,0:jmx,0:1)

      rho=0.
      den=0.
      jpar = 0.
      myden = 0.
      myjpar = 0.
      ns=1
if(idg.eq.1)write(*,*)'enter ion grid1',mm(1)
      do m=1,mm(1)
         dv=float(lr(1))*(dx*dy*dz)
         r=x3(m)-0.5*lx+lr0
         qr=q0+qp*(r-lr0)
         th=(z3(m)-0.5*lz)/(br0*q0)
         cost=cos(th)
         sint=sin(th)
         b=1.-lr0/br0*cost*tor
         rhog=sqrt(2.*b*mu(m)*mims(1))/(q(1)*b)
         sz=(z3(m)-0.5*lz)*e0*qp/q0*tor
         vfac=0.5*(mims(1)*u3(m)**2 + 2.*mu(m)*b )
!         vpar = u3(m)
         wght=w3(m)/dv

         vpar = u3(m) !linearly correct

!    now do 1,2,4 point average, where lr is the no. of points...
         do 200 l=1,lr(1)
            xs=x3(m)+rwx(1,l)*rhog
            yt=y3(m)+(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)

            include "gridngp.h"
 200     continue
      enddo
if(idg.eq.1)write(*,*)myid,'pass ion grid1'
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!   enforce periodicity
      call enforce(myden(:,:,:))
      call enforce(myjpar)
      call filter(myden(:,:,:))
      call filter(myjpar(:,:,:))

      do 400 ns=1,nsm
         do 410 i=0,im
            do 420 j=0,jm
               do 430 k=0,mykm
                  myden(i,j,k)=q(ns)*myden(i,j,k)/n0
                  myjpar(i,j,k) = q(ns)*myjpar(i,j,k)/n0*ifluid
 430           continue
 420        continue
 410     continue
 400  continue

	call MPI_ALLREDUCE(myden(0:im,0:jm,0:1),  &
     		den(1,0:im,0:jm,0:1),             &
     		(imx+1)*(jmx+1)*2,MPI_REAL8,       &
     		MPI_SUM,GRID_COMM,ierr)
	call MPI_ALLREDUCE(myjpar(0:im,0:jm,0:1), &
     		jpar(0:im,0:jm,0:1),              &
     		(imx+1)*(jmx+1)*2,MPI_REAL8,       &
     		MPI_SUM,GRID_COMM,ierr)

      do 440 ns=1,nsm
         do 450 i=0,im
            do 460 j=0,jm
               do 470 k=0,mykm
                  rho(i,j,k)=rho(i,j,k)+den(ns,i,j,k)
 470           continue
 460        continue
 450     continue
 440  continue

! electrons density and current
      vte = sqrt(amie)
      dene(ip,:,:,:) = 0.
      upar(ip,:,:,:) = 0.
      mydene = 0.
      myupar = 0.
if(idg.eq.1)write(*,*)'enter electron grid1'
      do m=1,mm(2)
         dv=(dx*dy*dz)

         r=x3e(m)-0.5*lx+lr0
         qr=q0+qp*(r-lr0)
         th=(z3e(m)-0.5*lz)/(br0*q0)
         cost=cos(th)
         sint=sin(th)
         b=1.-lr0/br0*cost*tor
         rhog=0.
         sz=(z3e(m)-0.5*lz)*e0*qp/q0*tor
         vfac=0.5*(mims(2)*u3e(m)**2 + 2.*mue3(m)*b )

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
      call filter(mydene(:,:,:))
      call filter(myupar(:,:,:))

      do  i=0,im
            do  j=0,jm
               do  k=0,mykm
                  mydene(i,j,k)= mydene(i,j,k)/n0*ifluid
                  myupar(i,j,k) = myupar(i,j,k)/n0*ifluid
               end do
            end do
         end do
	call MPI_ALLREDUCE(mydene(0:im,0:jm,0:1),   &
     		dene(ip,0:im,0:jm,0:1),             &
     		(imx+1)*(jmx+1)*2,MPI_REAL8,         &
     		MPI_SUM,GRID_COMM,ierr)
	call MPI_ALLREDUCE(myupar(0:im,0:jm,0:1),   &
            	upar(ip,0:im,0:jm,0:1),             &
     		(imx+1)*(jmx+1)*2,MPI_REAL8,         &
    		MPI_SUM,GRID_COMM,ierr)

 999  continue

      do i = 0,im
         do j = 0,jm
            do k = 0,mykm
               rho(i,j,k) = ision*rho(i,j,k) + dene(ip,i,j,k)*q(2)
            enddo
         enddo
      enddo      
                 
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine modes2(u,modehis,n)

!     calculate mode histories. calculate modes for u at timestep n
!     and store in modehis(mode,n).

      use gem_com
      use fft_wrapper
      implicit none
!     
      REAL(8) :: u(0:1,0:imx,0:jmx,0:1)
      COMPLEX(8) :: modehis(modemx,0:nmx)
      COMPLEX(8) :: modebuf
      INTEGER :: n,i,j,k,l,m,ii

      INTEGER :: mode,j,jj,thek,oproc,ifirst

!     

      if(n.eq.0) return

      do 100 mode=1,modem
         oproc=int(nmode(mode)/kcnt*ntube)

         if (MyId.eq.oproc) then
            thek=0
            do j=0,jm-1
               do i=0,im-1
                  tmpx(i)=u(0,i,j,thek)
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
      implicit none
      INTEGER :: m,ns,iflag,n,i,j,k,ip
      character*70 fname
      character*3 holdmyid

      write(holdmyid,'(I3.3)') MyId
      fname=directory//'dump_'//holdmyid//'.b'
      if(iflag.eq.1) then
         open(139+MyId,file=fname,form='unformatted',status='old')
         ns=1
         read(139+MyId)mm(ns),ncurr,tcurr
         do 110 m=ptr(ns),ptr(ns)+mm(ns)-1
            read(139+MyId) x1(m),y1(m),z1(m),u1(m),w1(m),mu(m)
            read(139+MyId) x2(m),y2(m),z2(m),u2(m),w2(m)
            w2(m)=w2(m)/cut
            w1(m)=w1(m)/cut
            x3(m)=x2(m)
            y3(m)=y2(m)
            z3(m)=z2(m)
            u3(m)=u2(m)
            w3(m)=w2(m)
 110     continue
         read(139+MyId)mm(2)
         do  m=1,mm(2)
            read(139+MyId) x1e(m),y1e(m),z1e(m),u1e(m),w1e(m),mue1(m)
            read(139+MyId) x2e(m),y2e(m),z2e(m),u2e(m),w2e(m),mue2(m),ipass(m)
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
         ns=1
         write(139+MyId)mm(ns),n+1,tcurr-dt
         do 210 m=ptr(ns),ptr(ns)+mm(ns)-1
            write(139+MyId) x1(m),y1(m),z1(m),u1(m),w1(m),mu(m)
            write(139+MyId) x2(m),y2(m),z2(m),u2(m),w2(m)
 210     continue
         write(139+MyId)mm(2)
         do  m=1,mm(2)
            write(139+MyId) x1e(m),y1e(m),z1e(m),u1e(m),w1e(m),mue1(m)
            write(139+MyId) x2e(m),y2e(m),z2e(m),u2e(m),w2e(m),mue2(m),ipass(m)
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
      implicit none
      INTEGER :: i,j
      REAL(8) :: dely,aweight,r,qr

!         peritr = 0
      if (GCLR.eq.Master.and.peritr.eq.0) then
	 do i=0,im
            r=float(i)*dx-0.5*lx+lr0
            qr=q0+qp*(r-lr0)
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
            qr=q0+qp*(r-lr0)
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
      use fft_wrapper
      implicit none
      real(8) :: lbfr(0:imx,0:jmx)
      real(8) :: lbfs(0:imx,0:jmx)
      real(8) :: rbfr(0:imx,0:jmx)
      real(8) :: rbfs(0:imx,0:jmx)
      REAL(8) :: b,b2,gam0,gam1,delyz,th,bf
      REAL(8) :: kx,ky
      real(8),dimension(:),allocatable :: akx,aky
      COMPLEX(8),dimension(:,:,:),allocatable :: fac
      REAL(8),dimension(:,:),allocatable :: akx2
      REAL(8),dimension(:,:,:),allocatable :: formphi
      REAL(8),dimension(:,:,:),allocatable :: formfe
      complex(8),dimension(:,:,:),allocatable :: fmjpe
      REAL(8) :: sgnx,sgny,sz,myfe
      INTEGER :: i,j,k,k1,l,m,n,ifirst,nstep,ip,isdvj0=0
      INTEGER :: l1,m1,myk,izonal
      COMPLEX(8) :: temp3d(0:imx-1,0:jmx-1,0:1),tempdv(0:imx-1,0:jmx-1,0:1)
      COMPLEX(8) :: aphik(0:imx-1),myaphik(0:imx-1)
      real(8) :: myrmsphi
      real(8) :: myaph(0:imx),aph(0:imx)


      save formphi,fmjpe,formfe,ifirst,akx,aky,fac,akx2
      izonal = 1

!     form factors....
      if (ifirst.ne.-99) then
         ALLOCATE( akx(0:imx-1),aky(0:jmx-1))
         ALLOCATE( fac(0:imx-1,0:jmx-1,0:kmx))
         ALLOCATE( formphi(0:imx-1,0:jmx-1,0:kmx),akx2(0:imx-1,0:kmx))
         ALLOCATE( formfe(0:imx-1,0:jmx-1,0:kmx))
         allocate( fmjpe(0:imx-1,0:jmx-1,0:kmx))

         do 51 l=0,im-1
            do 52 m=0,jm-1
               do 53 n = 0,km
                  th=( float(n)*dz-0.5*lz)/(br0*q0)
                  delyz = lr0/q0*qp*lx*th*ishift
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
                  kx=sgnx*2.*pi*float(l1)/lx-(delyz*ky)/lx
                  fac(l,m,n) = exp(IU*(aky(m)*delyz)/lx*xg(l))

                  bf=1.-lr0/br0*cos(th)*tor
                  sz=(float(n)*dz - 0.5*lz)*e0*qp/q0*tor
                  b=mims(1)*(kx*kx + ky*ky*(1. + sz**2)+2.*kx*ky*sz)/(bf*bf)
                  b2=(xshape*xshape*kx*kx + yshape*yshape*ky*ky)
                  b2=0. !b2**c4*(1-onemd)
                  
                  call srcbes(b,gam0,gam1)
                  
!   formfactor in gkps
                  formfe(l,m,n) = 1.-gam0
                  fmjpe(l,m,n) = kapn(1)/bf*IU*ky*(gam0-1.)
                  if(b.lt.3.e-5.and.fradi.eq.0.)then
                     formphi(l,m,n) = 0.
                  else if(b.gt.bcut)then
                     formphi(l,m,n) = 0.
                  else
                     formphi(l,m,n)=exp(-b2)/(fradi+tets(1)*(1.-gam0)) &
                         /float(imx*jmx)
                  end if 
                  if(l1.eq.0.and.m1.eq.0)formphi(l,m,n)=0.
!                  if(m1.eq.1)formphi(l,m,n)=0.
                  if(abs(ky).gt.kycut)formphi(l,m,n) = 0.
                  if(abs(kx).gt.kxcut)formphi(l,m,n) = 0.
                  if(m1.ne.mlk.and.onemd.eq.1)formphi(l,m,n) = 0.
!                  if(l1.ne.llk.and.onemd.eq.1)formphi(l,m,n) = 0.
                  if(ky.eq.0.)then
                     if(l1.eq.0)then
                        akx2(l,n) = 0.
                     else
                        akx2(l,n) = exp(-b2)/(1.-gam0)/float(km)
                     end if
                  end if
 53            continue
 52         continue
 51      continue

         ifirst=-99

      endif

!   now do field solve...

!      phi = 0.
!      return
      temp3d = 0.
      aphik = 0.
      myaphik = 0.

!  find rho(kx,ky)
      do k=0,mykm
         n = GCLR*kcnt+k
         do j=0,jm-1
            do i=0,im-1
               temp3d(i,j,k)=rho(i,j,k)+fradi*(phi(ip,i,j,k)-den0(i,j,k))
            enddo
         enddo
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3d(i,j,k)
            end do
            call ccfft('y',-1,jmx,1.,tmpy,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3d(i,j,k) = tmpy(j)   !rho(ky,x)
            end do
         end do
      end do
      call filtor(temp3d(0:imx-1,0:jmx-1,0:1))
      do k=0,mykm
         do j = 0,jmx-1
            do i = 0,imx-1
               temp3d(i,j,k) = temp3d(i,j,k)/fac(i,j,n)
            end do
         end do

         do j = 0,jmx-1
            do i = 0,imx-1
               tmpx(i) = temp3d(i,j,k)
            end do
            call ccfft('x',-1,imx,1.,tmpx,tmpx,coefx,workx,0)
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
      if(idg==1)write(*,*)'pass filtor', myid

! from phi(kx,ky) to divj0(kx,ky)
      do k=0,mykm
         myk=GCLR*kcnt+k         
         do j=0,jm-1
            do i=0,im-1
	       tempdv(i,j,k)=temp3d(i,j,k)*fmjpe(i,j,myk)
            enddo
         enddo
      enddo

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
         n = GCLR*kcnt+k
         do j = 0,jmx-1
            do i = 0,imx-1
               tmpx(i) = temp3d(i,j,k)
            end do
            call ccfft('x',1,imx,1.,tmpx,tmpx,coefx,workx,0)

            do i = 0,imx-1
               temp3d(i,j,k) = tmpx(i)  ! phi(ky,x)
            end do
         end do

         do i = 0,imx-1
            do j = 0,jmx-1
               temp3d(i,j,k) = temp3d(i,j,k)*fac(i,j,n)
            end do
         end do

         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3d(i,j,k)
            end do
            call ccfft('y',1,jmx,1.,tmpy,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3d(i,j,k) = tmpy(j)  ! phi(x,y)
            end do
         end do
      end do

!  from divj0(kx,ky) to divj0(x,y)
      if(isdvj0.eq.1)then
         do k=0,mykm
            n = GCLR*kcnt+k
            do j = 0,jmx-1
               do i = 0,imx-1
                  tmpx(i) = tempdv(i,j,k)
               end do
               call ccfft('x',1,imx,1.,tmpx,tmpx,coefx,workx,0)

               do i = 0,imx-1
                  tempdv(i,j,k) = tmpx(i) ! divj0(ky,x)
               end do
            end do

            do i = 0,imx-1
               do j = 0,jmx-1
                  tempdv(i,j,k) = tempdv(i,j,k)*fac(i,j,n)
               end do
            end do

            do i = 0,imx-1
               do j = 0,jmx-1
                  tmpy(j) = tempdv(i,j,k)
               end do
               call ccfft('y',1,jmx,1.,tmpy,tmpy,coefy,worky,0)
               do j = 0,jmx-1
                  tempdv(i,j,k) = tmpy(j) ! divj0(x,y)
               end do
            end do
         end do
      end if
!      goto 100

 100  continue
      do i = 0,im-1
         do j = 0,jm-1
            do k = 0,mykm
               phi(ip,i,j,k) = temp3d(i,j,k)
               divj0(i,j,k) = tempdv(i,j,k)
            end do
         end do
      end do

!    x-y boundary points 
      call enfxy(phi(ip,:,:,:))
      call enfz(phi(ip,:,:,:))
      call enfxy(divj0)
      call enfz(divj0)
      call filter(phi(ip,:,:,:))

      myrmsphi=0.
      rmsphi(nstep)=0.
      do k=0,mykm-1
         do j=0,jm-1
            do i=0,im-1
               myrmsphi=myrmsphi+phi(ip,i,j,k)*phi(ip,i,j,k)
!     myrmsphi=myrmsphi+abs(phi(ip,i,j,k))
            enddo
         enddo
      enddo
      call MPI_ALLREDUCE(myrmsphi,rmsphi(nstep),1, &
          MPI_REAL8,                               &
          MPI_SUM,TUBE_COMM,ierr)
      rmsphi(nstep)=sqrt(rmsphi(nstep)/(im*jm*km))

      if(izonal.eq.1)return
      do i=0,im
         myaph(i)=0.
         aph(i)=0.
      enddo
      do k=0,mykm-1
         do j=0,jm-1
            do i=0,im
	       myaph(i)=myaph(i)+phi(ip,i,j,k) 
            enddo
         enddo
      enddo

      call MPI_ALLREDUCE(myaph(0:imx),aph(0:imx),imx+1,  &
           MPI_REAL8,MPI_SUM,TUBE_COMM,ierr)
      aph = aph/(jmx*kmx)
      do i = 0,imx
         do j = 0,jmx
            do k = 0,mykm
               phi(ip,i,j,k) = phi(ip,i,j,k)-aph(i)
            end do
         end do
      end do
      return
      end

!      End of gkps....
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine eqmo(ip)
      use gem_com
      implicit none
      real(8) :: lbfr(0:imx,0:jmx)
      real(8) :: lbfs(0:imx,0:jmx)
      real(8) :: rbfr(0:imx,0:jmx)
      real(8) :: rbfs(0:imx,0:jmx)
      integer :: i,j,k,ip,isdadt
      real(8) :: eta

      isdadt = 0
      ez(ip,:,:,:) = 0.
      dadz = 0.
      dpdz = 0.
!      return
      do i = 0,im
         do j = 0,jm
            do k = 1,mykm-1
               ez(ip,i,j,k) = (phi(ip,i,j,k-1)-phi(ip,i,j,k+1))/(2.*dz)
!                   -dadt(i,j,k)*isdadt
               dadz(i,j,k) = (apar(ip,i,j,k+1)-apar(ip,i,j,k-1))/(2.*dz)
            end do
         end do
      end do   
      
      rbfs = phi(ip,:,:,0)

      call MPI_SENDRECV(rbfs(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8,rngbr,204,                    &
          lbfr(0,0),(imx+1)*(jmx+1),              &
          MPI_REAL8,lngbr,204,                    &
          TUBE_COMM,stat,ierr)

      lbfs=phi(ip,:,:,1)

      call MPI_SENDRECV(lbfs(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8,lngbr,205,                    &
          rbfr(0,0),(imx+1)*(jmx+1),              &
          MPI_REAL8,rngbr,205,                    &
          TUBE_COMM,stat,ierr)                    
      call MPI_BARRIER(TUBE_COMM,ierr)             
      do i=0,im
         do j=0,jm
            ez(ip,i,j,0)=(weightp(i)*lbfr(i,jpl(i,j)) &
                +weightpn(i)*lbfr(i,jpn(i,j))        &
                -phi(ip,i,j,1))/(2.*dz)
!                   -dadt(i,j,0)*isdadt
            ez(ip,i,j,1)=( phi(ip,i,j,0)-(weightm(i)*rbfr(i,jmi(i,j)) &
                +weightmn(i)*rbfr(i,jmn(i,j))))/(2.*dz)
!                   -dadt(i,j,1)*isdadt
         enddo
      enddo
      
      if(gclr==0)ez(ip,:,:,0)=ez(ip,:,:,1)
      if(gclr==glst)ez(ip,:,:,1)=ez(ip,:,:,0)

      rbfs = apar(ip,:,:,0)

      call MPI_SENDRECV(rbfs(0,0),(imx+1)*(jmx+1), &
          MPI_REAL8,rngbr,206,                    &
          lbfr(0,0),(imx+1)*(jmx+1),              &
          MPI_REAL8,lngbr,206,                    &
          TUBE_COMM,stat,ierr)

      lbfs = apar(ip,:,:,1)

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
                -apar(ip,i,j,1))/(2.*dz)

            dadz(i,j,1)=-( apar(ip,i,j,0)-(weightm(i)*rbfr(i,jmi(i,j)) &
                +weightmn(i)*rbfr(i,jmn(i,j))))/(2.*dz)
         enddo
      enddo

      if(gclr==0)dadz(:,:,0)=dadz(:,:,1)
      if(gclr==glst)dadz(:,:,1)=dadz(:,:,0)

      dpdz(:,:,:) = -ez(ip,:,:,:)

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine spec(n)
      use gem_com
      implicit none
      integer :: i,j,k,l,m,n

      i = tcurr-dt
      if(myid.eq.master)then
         if(izon.eq.0)then
            write(*,10)i,rmsphi(n),rmsapa(n),fe(n),efl(1,n),avewe(n),avewi(n)
!         write(*,11)pfl(1,n),pfl(2,n),pfl(3,n),efl(1,n),efl(2,n)
         end if
 10      format(1x,i6,6(2x,e10.3))
 11      format(6x,5(2x,e12.5))
 12      format(1x,i6,5(2x,e12.5))
         write(9,10)i,rmsphi(n),rmsapa(n),fe(n),avewe(n),avewi(n)
         write(11,12)i,pfl(1,n),pfl(2,n),pfl(3,n),efl(1,n),efl(2,n)
!         write(*,*)ftrap
      end if   
      
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ezamp(nstep,ip)   

      use gem_com
      use fft_wrapper
      implicit none
      real(8) :: lbfr(0:imx,0:jmx)
      real(8) :: lbfs(0:imx,0:jmx)
      real(8) :: rbfr(0:imx,0:jmx)
      real(8) :: rbfs(0:imx,0:jmx)

      REAL(8) :: b,b2,gam0,gam1,th,delyz,bf,rkper
      REAL(8) :: kx,ky
      REAL(8),DIMENSION(:),ALLOCATABLE :: akx,aky
      COMPLEX(8),DIMENSION(:,:,:),ALLOCATABLE :: fac
      REAL(8),DIMENSION(:,:,:), ALLOCATABLE :: formapa
      REAL(8) :: sgnx,sgny,sz,sf,xfrac=1.0
      INTEGER :: i,j,k,l,m,n,ifirs,ip,nstep
      INTEGER :: l1,m1,myk
      COMPLEX(8) :: temp3d(0:imx-1,0:jmx-1,0:1)
      REAL(8) :: myrmsapa


      save formapa,ifirs,akx,aky,fac

!     form factors....
      if (ifirs.ne.-99) then
!     
      ALLOCATE(akx(0:imx-1),aky(0:jmx-1))
      ALLOCATE(fac(0:imx-1,0:jmx-1,0:kmx))
      ALLOCATE(formapa(0:imx-1,0:jmx-1,0:kmx))

         formapa = 0.
         do 51 l=0,im-1
            do 52 m=0,jm-1
               do 53 n = 0,km
                  th=( float(n)*dz-0.5*lz)/(br0*q0)
                  delyz = lr0/q0*qp*lx*th*ishift
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
                  kx=sgnx*2.*pi*float(l1)/lx-dmod(delyz*ky,pi2)/lx
                  fac(l,m,n) = exp(IU*dmod(aky(m)*delyz,pi2)/lx*xg(l))
                  sf = (2.*sin(kx/2)/(kx+1.e-10))**2 &
                      *(2.*sin(ky/2)/(ky+1.e-10))**2
                  bf=1.-lr0/br0*cos(th)*tor
                  sz=(float(n)*dz - 0.5*lz)*e0*qp/q0*tor
                  rkper = sqrt(kx*kx+ky*ky)
                  b=(kx*kx + ky*ky*(1. + sz**2) +2.*kx*ky*sz)/(bf*bf)
                  b2=(xshape*xshape*kx*kx + yshape*yshape*ky*ky)
                  b2=0. !b2**c4*(1-onemd)

!                  call srcbes(b,gam0,gam1)
!  formfactor in ezamp
                  formapa(l,m,n)=exp(-b2)*beta/((b+xfrac*beta*amie+1.e-10) &
                      *float(imx*jmx))
                  if(l1.eq.0.and.m1.eq.0)formapa(l,m,n)=0.
!                  if(m1.eq.1)formapa(l,m,n)=0.
                  if(abs(ky).gt.kycut)formapa(l,m,n) = 0.
                  if(abs(kx).gt.kxcut)formapa(l,m,n) = 0.
                  if(m1.ne.mlk.and.onemd.eq.1)formapa(l,m,n) = 0.
!                  if(l1.ne.llk.and.onemd.eq.1)formapa(l,m,n) = 0.
                  if(b.gt.bcut)formapa(l,m,n) = 0.
 53            continue
 52         continue
 51      continue

         ifirs=-99

      endif

!   now do field solve...

      temp3d = 0.

!  find jtot(kx,ky)
      do k=0,mykm
         n = GCLR*kcnt+k
         do j=0,jm-1
            do i=0,im-1
               temp3d(i,j,k)=jpar(i,j,k)-upar(ip,i,j,k)  &
                     +amie*(xfrac*apar(ip,i,j,k)-upa0(i,j,k))
            enddo
         enddo
         
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3d(i,j,k)
            end do
            call ccfft('y',-1,jmx,1.,tmpy,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3d(i,j,k) = tmpy(j)   !jtot(ky,x)
            end do
         end do
      end do

      call filtor(temp3d(0:imx-1,0:jmx-1,0:1))

      do k=0,mykm
         do j = 0,jmx-1
            do i = 0,imx-1
               temp3d(i,j,k) = temp3d(i,j,k)/fac(i,j,n)
            end do
         end do

         do j = 0,jmx-1
            do i = 0,imx-1
               tmpx(i) = temp3d(i,j,k)
            end do
            call ccfft('x',-1,imx,1.,tmpx,tmpx,coefx,workx,0)
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

!  from apar(kx,ky) to apar(x,y)
      do k=0,mykm
         n = GCLR*kcnt+k
         do j = 0,jmx-1
            do i = 0,imx-1
               tmpx(i) = temp3d(i,j,k)
            end do
            call ccfft('x',1,imx,1.,tmpx,tmpx,coefx,workx,0)
            do i = 0,imx-1
               temp3d(i,j,k) = tmpx(i)  ! j(ky,x)
            end do
         end do

         do i = 0,imx-1
            do j = 0,jmx-1
               temp3d(i,j,k) = temp3d(i,j,k)*fac(i,j,n)
            end do
         end do

         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3d(i,j,k)
            end do
            call ccfft('y',1,jmx,1.,tmpy,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3d(i,j,k) = tmpy(j)  ! apar(x,y)
            end do
         end do
      end do

      do i = 0,im-1
         do j = 0,jm-1
            do k = 0,mykm
                  apar(ip,i,j,k) = temp3d(i,j,k)
            end do
         end do
      end do

!    x-y boundary points 
      call enfxy(apar(ip,:,:,:))
      call enfz(apar(ip,:,:,:))
      call filter(apar(ip,:,:,:))

      myrmsapa=0.
      rmsapa(nstep)=0.
      do k=0,mykm-1
         do j=0,jm-1
            do i=0,im-1
               myrmsapa=myrmsapa+apar(ip,i,j,k)*apar(ip,i,j,k)
            enddo
         enddo
      enddo
      call MPI_ALLREDUCE(myrmsapa,rmsapa(nstep),1, &
          MPI_REAL8,                               &
          MPI_SUM,TUBE_COMM,ierr)
      rmsapa(nstep)=sqrt(rmsapa(nstep)/(im*jm*km))

      return
      end
!!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine filter(u)

      use gem_com
      implicit none
      integer :: i,j,k,l,m,n,ip,NMDZ,NNI,ikz,ndum
      parameter(NMDZ=2)
      real(8) :: u(0:imx,0:jmx,0:1)
      real(8) :: temp(0:imx,0:jmx,0:1)
      real(8) :: cc0(0:imx,0:jmx,0:NMDZ-1),ss0(0:imx,0:jmx,0:NMDZ-1)
      real(8) :: cc1(0:imx,0:jmx,0:NMDZ-1),ss1(0:imx,0:jmx,0:NMDZ-1)
      real(8) ::  rkz,DW1,DW2
      parameter(NNI = 0)
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
      return
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

      if(kzlook==0)then
         do i = 0,im
            do j = 0,jm
               cc1(i,j,0) = 1./float(km)*cc1(i,j,0)*0.
               ss1(i,j,0) = 1./float(km)*ss1(i,j,0)
            end do
         end do
      end if
      if(kzlook>0)then
         do i = 0,im
            do j = 0,jm
               cc1(i,j,0) = 2./float(km)*cc1(i,j,0)
               ss1(i,j,0) = 2./float(km)*ss1(i,j,0)
            end do
         end do
      end if

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

      do k = 0,mykm-1
         tmpz(k) = tmp3d(llk,mlk,k)
      end do

      if(GCLR.ne.master)then
!         if(n.eq.1)write(*,*)n,myid, tmpz(2)
         call MPI_SEND(tmpz(0),mykm,MPI_DOUBLE_COMPLEX,master, &
             gclr,tube_comm,ierr)
      end if

      if(gclr.eq.master)then
         do i = 1,GLST
            call MPI_RECV(tmpz(i*mykm),mykm,MPI_DOUBLE_COMPLEX,i, &
                i,tube_comm,stat,ierr)
         end do

!         if(n.eq.1)then
!            do k = 0,km-1
!               write(*,*)n,myid,k,tmpz(k)
!            end do
!         end if   
      end if

      if(GCLR.eq.master) then
!     FT in z....
         call ccfft('z',1,kmx,1.,tmpz,tmpz,coefz,workz,0)
      end if

      call MPI_BCAST(tmpz,kmx,MPI_DOUBLE_COMPLEX,master, &
          tube_comm,ierr)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      do i = 1,kmx
         yyamp(i) = abs(tmpz(i-1))  !cabs
         yyre(i) = real(tmpz(i-1))
         yyim(i) = aimag(tmpz(i-1))
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
      subroutine loadi

      use gem_com
      implicit none
      INTEGER :: m,idum,ns,m1
      INTEGER :: np_old,np_new
      REAL(8) :: vpar,vperp2,r,qr,th,b,cost
      REAL(8) :: avgv,myavgv,avgw,myavgw


      cnt=int(tmm(1)/numprocs)

!   write(*,*)'in loader cnt, mm(1)=',cnt,mm(1)

      myavgv=0.
      avgv=0.
      avgw = 0.
      myavgw = 0.

      do 160 m=1,mm(1)

!     load a slab of ions...

         x2(m)=lx*revers(MyId*cnt+m,2) !ran2(iseed)
         y2(m)=ly*revers(MyId*cnt+m,3) !ran2(iseed)
         z2(m)=lz*revers(MyId*cnt+m,5) !ran2(iseed)

         call parperp(vpar,vperp2,m,pi,cnt,MyId)
!   normalizations will be done in following loop...

         r=x2(m)-0.5*lx+lr0
         qr=q0+qp*(r-lr0)
         th=(z2(m)-0.5*lz)/(br0*q0)
         cost=cos(th)
         b=1.-lr0/br0*cost*tor

         u2(m)=vpar/sqrt(mims(1))
         mu(m)=0.5*vperp2/b
         myavgv=myavgv+u2(m)

!    LINEAR: perturb w(m) to get linear growth...
         w2(m)=2.*amp*(revers(MyId*cnt+m,13)-0.5) !(ran2(iseed) - 0.5 )
         if(izon.eq.1)w2(m)=2.*amp*sin(x2(m)/lx*2*pi)  
         myavgw=myavgw+w2(m)
 160  continue

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
         x1(m)=x2(m)
         y1(m)=y2(m)
         z1(m)=z2(m)
         x3(m)=x2(m)
         y3(m)=y2(m)
         z3(m)=z2(m)
         u1(m)=u2(m)
         u3(m)=u2(m)
!         w2(m) = w2(m)-myavgw
         w1(m)=w2(m)
         w3(m)=w2(m)
 180  continue

      np_old=mm(1)
      call init_pmove(z1,np_old,lz,ierr)
      if(idg.eq.1)write(*,*)'pass init_pmove'
!     
      call pmove(x1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mu,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
!     
      call end_pmove(ierr)
      mm(1)=np_new
!     write(*,*)MyId,mm(1)
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
            x1(m1)=x2(m1)
            y1(m1)=y2(m1)
            z1(m1)=z2(m1)
            u1(m1)=u2(m1)

 250     continue
 260  continue

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine loade

      use gem_com
      implicit none
      integer :: m,idum,i,j,k
      real(8) :: vratio,dv,xt,yt,phip
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,vte
      dv = dx*dy*dz

      vratio = vwidthe/vwidth
!      write(*,*)'vratio= ',vratio

!    species two loaded on top of ions, for now species one
!    is ions and species two is electrons, additional species will
!    require modification of loader to initialize particle array...

      index = 0
      ipass = 0
      do 250 m=1,mm(1)
         x2e(m) = x1(m)
         y2e(m) = y1(m)
         z2e(m) = z1(m)
         u2e(m) = sqrt(mims(1)*amie)*u1(m)*vratio
         mue2(m) = tets(2)*mu(m)*vratio**2
         xt = x2e(m)
         yt = y2e(m)
         i=int(xt/dx)
         j=int(yt/dy)
         k=int(z2e(m)/dz)-gclr*kcnt

         wx0=float(i+1)-xt/dx 
         wx1=1.-wx0
         wy0=float(j+1)-yt/dy
         wy1=1.-wy0
         wz0=float(gclr*kcnt+k+1)-z2e(m)/dz
         wz1=1.-wz0

         phip = wx0*wy0*wz0*phi(0,i,j,k) &
           + wx1*wy0*wz0*phi(0,i+1,j,k) &
           + wx0*wy1*wz0*phi(0,i,j+1,k) &
           + wx1*wy1*wz0*phi(0,i+1,j+1,k) & 
           + wx0*wy0*wz1*phi(0,i,j,k+1) &
           + wx1*wy0*wz1*phi(0,i+1,j,k+1) & 
           + wx0*wy1*wz1*phi(0,i,j+1,k+1) &
           + wx1*wy1*wz1*phi(0,i+1,j+1,k+1) 

         w2e(m) = 2.*amp*(ran2(iseed)-0.5)  
         if(izon.eq.1)w2e(m) = -isg*phip

         x3e(m) = x2e(m)
         y3e(m) = y2e(m)
         z3e(m) = z2e(m)
         u3e(m) = u2e(m)
         w3e(m) = w2e(m)
         mue3(m) = mue2(m)

         x1e(m) = x2e(m)
         y1e(m) = y2e(m)
         z1e(m) = z2e(m)
         u1e(m) = u2e(m)
         w1e(m) = w2e(m)
         mue1(m) = mue2(m)
         index(m) = myid*mm(1)+m
 250  continue

      mm(2) = mm(1)

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine enforce(u)
      use gem_com
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
            th = (zg(k)-0.5*lz)/(q0*br0)
            dely = lr0/q0*qp*lx*th*ishift
            ydum = yg(j)+dely
            ydum = dmod(ydum+80.*ly,ly)
            jj = int(ydum/dy)
            wy1 = float(jj+1)-ydum/dy
            u(0,j,k) = u(0,j,k)+u(im,jj,k)*wy1+u(im,jj+1,k)*(1.-wy1)
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
         th = (zg(k)-0.5*lz)/(q0*br0)
         dely = -lr0/q0*qp*lx*th*ishift
         do 712 j=0,jm
            ydum = dmod(yg(j)+dely+80.*ly,ly)
            jj = int(ydum/dy)
            wy1 = float(jj+1)-ydum/dy
            u(im,j,k)=u(0,jj,k)*wy1+u(0,jj+1,k)*(1-wy1)
 712     continue
 711  continue

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine gradu(u,ux,uy)
      use gem_com
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
            th = (zg(k)-0.5*lz)/(q0*br0)
            dely = lr0/q0*qp*lx*th*ishift
            ydum = dmod(yg(j)+dely+80.*ly,ly)
            jj = int(ydum/dy)
            wy1 = float(jj+1)-ydum/dy
            ul=u(im-1,jj,k)*wy1+u(im-1,jj+1,k)*(1-wy1)
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
            th = (zg(k)-0.5*lz)/(q0*br0)
            dely = lr0/q0*qp*lx*th*ishift
            ydum = dmod(yg(j)+dely+80.*ly,ly)
            jj = int(ydum/dy)
            wy1 = float(jj+1)-ydum/dy
            ul=u(im-1,jj,k)*wy1+u(im-1,jj+1,k)*(1-wy1)
            ux(0,j,k)=(u(1,j,k)-ul)/(2.*dx)
         enddo
      enddo

      call enfxy(ux)
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine grady(u,uy)
      use gem_com
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
      implicit none
      real(8) :: u(0:imx,0:jmx,0:1)
      real(8) :: lbfr(0:imx,0:jmx)
      real(8) :: lbfs(0:imx,0:jmx)
      real(8) :: rbfr(0:imx,0:jmx)
      real(8) :: rbfs(0:imx,0:jmx)
      integer :: i,j
      return
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
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dgdtp,dpdzp,dadzp,aparp
      REAL(8) :: dum,vxdum,dum1,eps,x,h_x,h_coll
      INTEGER :: m,i,j,k,l,n
      INTEGER :: np_old,np_new
      REAL(8) :: vfac,kap,sz,vpar,ppar,vpdum,pidum
      REAL(8) :: b,th,r,b,enerb,cost,sint,qr,laps
      REAL(8) :: xt,xs,yt,zt,xdot,ydot,zdot,pzdot,edot,pzd0,pzd1
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1

      pidum = 1./(pi*2)**1.5*(vwidthe)**3
      do m=1,mm(2)
         r=x2e(m)-0.5*lx+lr0
         qr=q0+qp*(r-lr0)
         th=(z2e(m)-0.5*lz)/(br0*q0)
         cost=cos(th)
         sint=sin(th)
         b=1.-lr0/br0*cost*tor

         sz=(z2e(m)-0.5*lz)*e0*qp/q0*tor !e0=lr0/q0/br0

         xt = x2e(m)
         yt = y2e(m)

         include 'ppushlie.h'

         vfac = 0.5*(mims(2)*u2e(m)**2 + 2.*mue2(m)*b)
         kap = kapn(2) - (1.5-vfac)*kapt(2)

         ppar = u2e(m)
         vpar = u2e(m)-q(2)/mims(2)*aparp*nonline
         enerb=(mue2(m)+mims(2)*vpar*vpar/b)/q(2)*tor
         vxdum = (eyp+vpar/b*delbxp)
         xdot = vxdum*nonline  &
             -enerb*sint/br0*tor
         ydot = (-exp1+vpar/b*delbyp)*nonline &
             -enerb*b*(cost/br0+e0*qp*th*sint)*tor
         zdot =  vpar   !-enerb*q0/r*cost
!             -vpar*q0*br0/lr0*(-delbyp+lr0/q0*qp*th*delbxp)*tor
!             -1./b*q0*br0/r*(exp1+qp*th*eyp)*tor

         pzd0 = (-b*mue2(m)/mims(2)*lr0*sint/(q0*br0*br0))*tor
         pzd1 = q(2)/mims(2)*ezp  &
             +q(2)/mims(2)*(-xdot*delbyp+ydot*delbxp+zdot*dadzp) &
             +(-mue2(m)/mims(2)/b*(delbyp*sint/br0-delbxp/br0 &
                                     *(cost+lr0/q0*qp*th*sint) &
                                     -dadzp*q0/r*cost) ) *tor &
             -tor*vpar/b/br0*(exp1*sint+(cost+lr0/q0*qp*th*sint)*eyp &
                               +q0*br0/r*cost*ezp)
         pzdot = pzd0+pzd1*nonline
                      
         vpdum = vpar
         edot = q(2)*(xdot*exp1+ydot*eyp+vpdum*ezp) &
                +q(2)*pzdot*aparp &
                +q(2)*vpdum*(-xdot*delbyp+ydot*delbxp+vpdum*dadzp)

         x3e(m) = x1e(m) + 2.*dte*xdot
         y3e(m) = y1e(m) + 2.*dte*ydot
         z3e(m) = z1e(m) + 2.*dte*zdot
         u3e(m) = u1e(m) + 2.*dte*pzdot
         mue3(m) = mue1(m)

         eps = b*mue2(m)+0.5*mims(2)*u2e(m)*u2e(m)
         x = sqrt(eps)
         h_x    = 4.0*x*x/(3.0*sqrt(pi))
         h_coll = h_x/sqrt(1.0+h_x**2)
!         h_coll = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*DERF(x)
         dum1 = 0.5*rneu/eps**1.5*(1+h_coll)
         if(x<0.316)dum1=0.0

         dum = 1-w2e(m)*nonline*0.
         if(isuni.eq.1)dum = exp(-vfac)*pidum
         vxdum = eyp+vpdum/b*delbxp
!             -(2*u2e(m)*aparp+q(2)/mims(2)*aparp*aparp)/(b*b)/br0*sint
         w3e(m)=w1e(m) + 2.*dte*(  &
             vxdum*kap + edot-2*dum1*ppar*aparp     &
             +isg*(-dgdtp-vpdum*dpdzp+xdot*exp1+ydot*eyp &
                   +phip*(vxdum*kap*nonline-xdot*kapt(2)) &
                   +edot*phip*nonline) &
             )*dum
!         if(z3e(m)>lz .or. z3e(m)<0.)w3e(m) = 0. 
         u1e(m)=u2e(m)  + .25*( u3e(m) - u1e(m) )
         x1e(m)=x2e(m)  + .25*( x3e(m) - x1e(m) )
         y1e(m)=y2e(m)  + .25*( y3e(m) - y1e(m) )
         z1e(m)=z2e(m)  + .25*( z3e(m) - z1e(m) )
         w1e(m)=w2e(m)  + .25*( w3e(m) - w1e(m) )
         mue1(m) = mue2(m)

         laps=anint((z3e(m)/lz)-.5)*(1-peritr)
         r=x3e(m)-0.5*lx+lr0
         qr=q0+qp*(r-lr0)
         y3e(m)=dmod(y3e(m)-laps*2*pi*qr*lr0/q0+80.*ly,ly)
         z3e(m)=dmod(z3e(m)+8.*lz,lz)
         x3e(m)=dmod(x3e(m)+8.*lx,lx)         

      end do
      np_old=mm(2)
      call init_pmove(z3e,np_old,lz,ierr)

      call pmove(x1e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y1e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z1e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u1e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w1e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mue1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mue2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mue3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call pmove(index,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ipass,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit


      call end_pmove(ierr)
      mm(2)=np_new

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine cint(n)
      use gem_com
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dgdtp,dpdzp,dadzp,aparp
      REAL(8) :: dum,vxdum,dum1,eps,x,h_x,h_coll
      INTEGER :: m,i,j,k,l,n
      INTEGER :: np_old,np_new
      REAL(8) :: vfac,kap,sz,vpar,ppar,vpdum,pidum
      REAL(8) :: b,th,r,b,enerb,cost,sint,qr,laps
      REAL(8) :: xt,xs,yt,zt,xdot,ydot,zdot,pzdot,edot,pzd0,pzd1
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1
      REAL(8) :: myke,mypfl,myptrp
      REAL(8) :: myefl,mynos,myavewe
      REAL(8) :: ketemp,pfltemp,ptrptmp
      REAL(8) :: efltemp,nostemp
      REAL(8) :: sbuf(10),rbuf(10)
      REAL(8) :: mytotn,mytrap,totn,ttrap

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

      pidum = 1./(pi*2)**1.5*(vwidthe)**3
      do m=1,mm(2)
         r=x3e(m)-0.5*lx+lr0
         qr=q0+qp*(r-lr0)
         th=(z3e(m)-0.5*lz)/(br0*q0)
         cost=cos(th)
         sint=sin(th)
         b=1.-lr0/br0*cost*tor

         sz=(z3e(m)-0.5*lz)*e0*qp/q0*tor !e0=lr0/q0/br0

         xt = x3e(m)
         yt = y3e(m)

         include 'cpushlie.h'

         vfac = 0.5*(mims(2)*u3e(m)**2 + 2.*mue3(m)*b)
         kap = kapn(2) - (1.5-vfac)*kapt(2)

         ppar = u3e(m)
         vpar = u3e(m)-q(2)/mims(2)*aparp*nonline
         enerb=(mue3(m)+mims(2)*vpar*vpar/b)/q(2)*tor
         vxdum = (eyp+vpar/b*delbxp)
         xdot = vxdum*nonline &
             -enerb*sint/br0*tor
         ydot = (-exp1+vpar/b*delbyp)*nonline  &
             -enerb*b*(cost/br0+e0*qp*th*sint)*tor
         zdot =  vpar   !-enerb*q0/r*cost
!             -vpar*q0*br0/lr0*(-delbyp+lr0/q0*qp*th*delbxp)*tor
!             -1./b*q0*br0/r*(exp1+qp*th*eyp)*tor

         pzd0 = (-b*mue3(m)/mims(2)*lr0*sint/(q0*br0*br0))*tor
         pzd1 = q(2)/mims(2)*ezp  &
             +q(2)/mims(2)*(-xdot*delbyp+ydot*delbxp+zdot*dadzp) &
             +(-mue3(m)/mims(2)/b*(delbyp*sint/br0-delbxp/br0  &
                                     *(cost+lr0/q0*qp*th*sint) &
                                     -dadzp*q0/r*cost) )*tor  &
             -tor*vpar/b/br0*(exp1*sint+(cost+lr0/q0*qp*th*sint)*eyp &
                               +q0*br0/r*cost*ezp)
         pzdot = pzd0+pzd1*nonline

         vpdum = vpar
         edot = q(2)*(xdot*exp1+ydot*eyp+vpdum*ezp) &
                +q(2)*pzdot*aparp  &
                +q(2)*vpdum*(-xdot*delbyp+ydot*delbxp+vpdum*dadzp)

         x3e(m) = x1e(m) + 0.5*dte*xdot
         y3e(m) = y1e(m) + 0.5*dte*ydot
         z3e(m) = z1e(m) + 0.5*dte*zdot
         u3e(m) = u1e(m) + 0.5*dte*pzdot
         mue3(m) = mue1(m)

         eps = b*mue3(m)+0.5*mims(2)*u3e(m)*u3e(m)
         x = sqrt(eps)
         h_x    = 4.0*x*x/(3.0*sqrt(pi))
         h_coll = h_x/sqrt(1.0+h_x**2)
!         h_coll = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*DERF(x)
         dum1 = 0.5*rneu/eps**1.5*(1+h_coll)
         if(x<0.316)dum1=0.0

         dum = 1-w3e(m)*nonline*0.
         if(isuni.eq.1)dum = exp(-vfac)*pidum
         vxdum = eyp+vpdum/b*delbxp
!             -(2*u3e(m)*aparp+q(2)/mims(2)*aparp*aparp)/(b*b)/br0*sint
         w3e(m)=w1e(m) + 0.5*dte*(  &
             vxdum*kap + edot  -2*dum1*ppar*aparp    & 
             +isg*(-dgdtp-vpdum*dpdzp+xdot*exp1+ydot*eyp &
                   +phip*(vxdum*kap*nonline-xdot*kapt(2))  &
                   +edot*phip*nonline)  &
             )*dum

         if(abs(w3e(m)).gt.1..and.nonline==1)then
            w3e(m) = 0.
            w2e(m) = 0.
            w1e(m) = 0.
         end if

!         if(z3e(m)>lz .or. z3e(m)<0.)w3e(m) = 0. 
         mytotn = mytotn+1
         if(isuni.eq.1)mytotn = mytotn+exp(-vfac)
         mytrap = mytrap+1-ipass(m)
         if(isuni.eq.1)mytrap = mytrap+exp(-vfac)*(1-ipass(m))

         laps=anint((z3e(m)/lz)-.5)*(1-peritr)
         r=x3e(m)-0.5*lx+lr0
         qr=q0+qp*(r-lr0)
         y3e(m)=dmod(y3e(m)-laps*2*pi*qr*lr0/q0+80.*ly,ly)
         z3e(m)=dmod(z3e(m)+8.*lz,lz)
         x3e(m)=dmod(x3e(m)+8.*lx,lx)         

!     particle diagnostics done here because info is available...
         mypfl=mypfl + w3e(m)*(eyp+vpar*delbxp/b) 
         myptrp=myptrp + w3e(m)*eyp*(1-ipass(m)) 
         myefl=myefl + vfac*w3e(m)*(eyp/b+vpar*delbxp/b)
         myke=myke + vfac*w3e(m)
         mynos=mynos + w3e(m)
         myavewe = myavewe+abs(w3e(m))

         if(abs(z3e(m)-z2e(m)).gt.lz/2)ipass(m)=1

!     xn becomes xn-1...
         u1e(m)=u2e(m)
         x1e(m)=x2e(m)
         y1e(m)=y2e(m)
         z1e(m)=z2e(m)
         w1e(m)=w2e(m)
!    xn+1 becomes xn...
         u2e(m)=u3e(m)
         x2e(m)=x3e(m)
         y2e(m)=y3e(m)
         z2e(m)=z3e(m)
         w2e(m)=w3e(m)

      enddo

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
!      ke(2,n)=ketemp/( 2.*float(tmm(1))*mims(1) )
      ftrap = ttrap/totn

      np_old=mm(2)
      call init_pmove(z3e,np_old,lz,ierr)

      call pmove(x1e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y1e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z1e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u1e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w1e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w2e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w3e,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mue1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mue2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mue3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit

      call pmove(index,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(ipass,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit


      call end_pmove(ierr)
      mm(2)=np_new

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine unild
      use gem_com
      implicit none
      INTEGER :: i,j,k,l,n,m,idum,idtmp,ns,m1
      INTEGER :: np_old,np_new
      REAL(8) :: vpar,vperp2,dum1,dum2,vperp2,vpar,r,qr,th,cost,b
      REAL(8) :: avgv,myavgv,avgw,myavgw

      cnt=int(tmm(1)/numprocs)

!   write(*,*)'in loader cnt, mm(1)=',cnt,mm(1)

      if(myid.eq.master)then
         do idtmp = 1,last
            myavgv=0.
            avgv=0.
            avgw = 0.
            myavgw = 0.

            do 160 m=1,mm(1)
               x2(m)=lx*ran2(iseed)
               y2(m)=ly*ran2(iseed)
               z2(m)=lz*ran2(iseed)

               dum1 = vwidth*(ran2(iseed)-0.5)
               dum2 = vwidth*(ran2(iseed)-0.5)
               vperp2 = dum1**2+dum2**2
               vpar = vwidth*(ran2(iseed)-0.5)

!   normalizations will be done in following loop...

               r=x2(m)-0.5*lx+lr0
               qr=q0+qp*(r-lr0)
               th=(z2(m)-0.5*lz)/(br0*q0)
               cost=cos(th)
               b=1.-lr0/br0*cost*tor

               u2(m)=vpar/sqrt(mims(1))
               mu(m)=0.5*vperp2/b
               myavgv=myavgv+u2(m)

!    LINEAR: perturb w(m) to get linear growth...
               w2(m)=2.*amp*(ran2(iseed) - 0.5 )
               myavgw=myavgw+w2(m)
 160        continue

            myavgw = myavgw/mm(1)

            call MPI_SEND(x2(1),mm(1),MPI_REAL8,idtmp,3,  &
                MPI_COMM_WORLD,ierr)
            call MPI_SEND(y2(1),mm(1),MPI_REAL8,idtmp,4,  &
                MPI_COMM_WORLD,ierr)
            call MPI_SEND(z2(1),mm(1),MPI_REAL8,idtmp,5,  &
                MPI_COMM_WORLD,ierr)
            call MPI_SEND(w2(1),mm(1),MPI_REAL8,idtmp,6,  &
                MPI_COMM_WORLD,ierr)
            call MPI_SEND(mu(1),mm(1),MPI_REAL8,idtmp,7,  &
                MPI_COMM_WORLD,ierr)
            call MPI_SEND(myavgv,1,MPI_REAL8,idtmp,8,     &
                MPI_COMM_WORLD,ierr)
            call MPI_SEND(u2(1),mm(1),MPI_REAL8,idtmp,9,  &
                MPI_COMM_WORLD,ierr)
         end do

! master worry about itself here
         myavgv=0.
         avgv=0.
         avgw = 0.
         myavgw = 0.

         do  m=1,mm(1)
            x2(m)=lx*ran2(iseed)
            y2(m)=ly*ran2(iseed)
            z2(m)=lz*ran2(iseed)

            dum1 = vwidth*(ran2(iseed)-0.5)
            dum2 = vwidth*(ran2(iseed)-0.5)
            vperp2 = dum1**2+dum2**2
            vpar = vwidth*(ran2(iseed)-0.5)

!   normalizations will be done in following loop...

            r=x2(m)-0.5*lx+lr0
            qr=q0+qp*(r-lr0)
            th=(z2(m)-0.5*lz)/(br0*q0)
            cost=cos(th)
            b=1.-lr0/br0*cost*tor

            u2(m)=vpar/sqrt(mims(1))
            mu(m)=0.5*vperp2/b
            myavgv=myavgv+u2(m)

!    LINEAR: perturb w(m) to get linear growth...
            w2(m)=2.*amp*(ran2(iseed) - 0.5 )
            myavgw=myavgw+w2(m)
         end do

         myavgw = myavgw/mm(1)

      end if

      if(myid.ne.master)then
         call MPI_RECV(x2(1),mm(1),MPI_REAL8,master,3,  &
             MPI_COMM_WORLD,stat,ierr)
         call MPI_RECV(y2(1),mm(1),MPI_REAL8,master,4,  &
             MPI_COMM_WORLD,stat,ierr)
         call MPI_RECV(z2(1),mm(1),MPI_REAL8,master,5,  &
             MPI_COMM_WORLD,stat,ierr)
         call MPI_RECV(w2(1),mm(1),MPI_REAL8,master,6,  &
             MPI_COMM_WORLD,stat,ierr)
         call MPI_RECV(mu(1),mm(1),MPI_REAL8,master,7,  &
             MPI_COMM_WORLD,stat,ierr)
         call MPI_RECV(myavgv,1,MPI_REAL8,master,8,    &
             MPI_COMM_WORLD,stat,ierr)
         call MPI_RECV(u2(1),mm(1),MPI_REAL8,master,9, &
             MPI_COMM_WORLD,stat,ierr)
      end if   
         
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!        write(*,*)'myid ', myid,x2(10),mu(20),u2(20),w2(20)
!    subtract off avg. u...
      call MPI_ALLREDUCE(myavgv,avgv,1,  &
          MPI_REAL8, &
          MPI_SUM,MPI_COMM_WORLD,ierr)
      if(idg.eq.1)write(*,*)'all reduce'
      avgv=avgv/float(tmm(1))
      do 180 m=1,mm(1)
         u2(m)=u2(m)-avgv
         x1(m)=x2(m)
         y1(m)=y2(m)
         z1(m)=z2(m)
         x3(m)=x2(m)
         y3(m)=y2(m)
         z3(m)=z2(m)
         u1(m)=u2(m)
         u3(m)=u2(m)
!         w2(m) = w2(m)-myavgw
         w1(m)=w2(m)
         w3(m)=w2(m)
 180  continue

      np_old=mm(1)
      call init_pmove(z1,np_old,lz,ierr)
      if(idg.eq.1)write(*,*)'pass init_pmove'
!     
      call pmove(x1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(x3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(y3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(z3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(u3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w1,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w2,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(w3,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call pmove(mu,np_old,np_new,ierr)
      if (ierr.ne.0) call ppexit
      call end_pmove(ierr)

      mm(1)=np_new

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine nadi(ip)
      use gem_com
      implicit none
      INTEGER :: i,j,k,l,m,n,ifirst,nstep,ip

      dgdt(:,:,:) = -q(2)*dphidt(:,:,:)
      negp = 0.

      return 
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine drdt(ip)
      use gem_com
      implicit none
      integer :: i,j,k,ip,isdvj0,isdndt
      real(8) :: lbfr(0:imx,0:jmx)
      real(8) :: lbfs(0:imx,0:jmx)
      real(8) :: rbfr(0:imx,0:jmx)
      real(8) :: rbfs(0:imx,0:jmx)
      real(8) :: dum
      real(8) :: djdx(0:imx,0:jmx,0:1),djdy(0:imx,0:jmx,0:1)

      call gradx(jpex-upex,djdx)
      call grady(jpey-upey,djdy)
!      djdx = 0.
!      djdy = 0.
      isdvj0 = 0
      isdndt = 1
!      upa0(:,:,:) = apar(ip,:,:,:)
      do i = 0,im
         do j = 0,jm
            rbfs(i,j)=jpar(i,j,0)-(upar(ip,i,j,0)+amie*upa0(i,j,0))
         end do   
      end do
      call MPI_SENDRECV(rbfs,(imx+1)*(jmx+1),  &
          MPI_REAL8,rngbr,404, &
          lbfr,(imx+1)*(jmx+1), &
          MPI_REAL8,lngbr,404, &
          tube_comm,stat,ierr)
      do i = 0,im
         do j = 0,jm
            lbfs(i,j)=jpar(i,j,1)-(upar(ip,i,j,1)+amie*upa0(i,j,1))
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
            drhodt(ip,i,j,0) = -(jpar(i,j,1)-upar(ip,i,j,1) &
                -amie*upa0(i,j,1)-dum)/(2.*dz) &
                -djdx(i,j,0)-djdy(i,j,0) &
                +(dnidt(i,j,0)-dnedt(i,j,0))*isdndt &
                -divj0(i,j,0)*isdvj0  &
                -2./br0*(ex(ip,i,j,0)*cfx(i,0) &
                         +ey(ip,i,j,0)*cfy(i,0))*isg &
              +nonline*(-ex(ip,i,j,0)*dnedy(i,j,0)+ey(ip,i,j,0)*dnedx(i,j,0) &
                +delbx(ip,i,j,0)*dupadx(i,j,0)+delby(ip,i,j,0)*dupady(i,j,0))
            dum = weightm(i)*rbfr(i,jmi(i,j)) &
                +weightmn(i)*rbfr(i,jmn(i,j))
            drhodt(ip,i,j,1) = -(dum-(jpar(i,j,0)-upar(ip,i,j,0)  &
                -amie*upa0(i,j,0)))/(2.*dz)  &
                -djdx(i,j,1)-djdy(i,j,1) &
                +(dnidt(i,j,1)-dnedt(i,j,1))*isdndt  &
                -divj0(i,j,1)*isdvj0  &
                -2./br0*(ex(ip,i,j,1)*cfx(i,1)  &
                         +ey(ip,i,j,1)*cfy(i,1))*isg &
              +nonline*(-ex(ip,i,j,1)*dnedy(i,j,1)+ey(ip,i,j,1)*dnedx(i,j,1) &
                +delbx(ip,i,j,1)*dupadx(i,j,1)+delby(ip,i,j,1)*dupady(i,j,1))
         end do
      end do
      call enfxy(drhodt(ip,:,:,:))
      call enfz(drhodt(ip,:,:,:))
!      call filter(drhodt(ip,:,:,:))
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dpdt(ip)   
      use gem_com
      use fft_wrapper
      implicit none
      real(8) :: lbfr(0:imx,0:jmx)
      real(8) :: lbfs(0:imx,0:jmx)
      real(8) :: rbfr(0:imx,0:jmx)
      real(8) :: rbfs(0:imx,0:jmx)
      REAL(8) :: b,b2,gam0,gam1,delyz,th,bf
      REAL(8) :: kx,ky
      REAL(8),dimension(:),allocatable :: akx,aky
      REAL(8),dimension(:,:,:),allocatable :: formphi
      REAL(8),dimension(:,:),allocatable :: akx2
      REAL(8) :: sgnx,sgny,sz,myfe
      INTEGER :: i,j,k,k1,l,m,n,ifirst,nstep,ip
      INTEGER :: l1,m1,myk,izonal
      COMPLEX(8) :: temp3d(0:imx-1,0:jmx-1,0:1)
      real(8) :: myaph(0:imx),aph(0:imx)


      save formphi,ifirst,akx,aky,akx2
      izonal = 1

!     form factors....
      if (ifirst.ne.-99) then
      allocate(akx(0:imx-1),aky(0:jmx-1))
      allocate(formphi(0:imx-1,0:jmx-1,0:kmx),akx2(0:imx-1,0:kmx))

         do 51 l=0,im-1
            do 52 m=0,jm-1
               do 53 n = 0,km
                  th=( float(n)*dz-0.5*lz)/(br0*q0)
                  delyz = lr0/q0*qp*lx*th*ishift
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
                  kx=sgnx*2.*pi*float(l1)/lx-(delyz*ky)/lx

                  bf=1.-lr0/br0*cos(th)*tor
                  sz=(float(n)*dz - 0.5*lz)*e0*qp/q0*tor
                  b=mims(1)*(kx*kx + ky*ky*(1. + sz**2)+2.*kx*ky*sz)/(bf*bf)
                  b2=(xshape*xshape*kx*kx + yshape*yshape*ky*ky)
                  b2=0. !b2**c4*(1-onemd)
                  
                  call srcbes(b,gam0,gam1)
                  
!   formfactor in dpdt
                  if(b.lt.3.e-5)then
                     formphi(l,m,n) = 0.
                  else if(b.gt.bcut)then
                     formphi(l,m,n) = 0.
                  else
                     formphi(l,m,n)=exp(-b2)/(0.+tets(1)*(1.-gam0)) &
                         /float(imx*jmx)
                  end if 
!                  if(m1.eq.1)formphi(l,m,n)=0.
                  if(abs(ky).gt.kycut)formphi(l,m,n) = 0.
                  if(abs(kx).gt.kxcut)formphi(l,m,n) = 0.
                  if(m1.ne.mlk.and.onemd.eq.1)formphi(l,m,n) = 0.
!                  if(l1.ne.llk.and.onemd.eq.1)formphi(l,m,n) = 0.

 53            continue
 52         continue
 51      continue

         ifirst=-99

      endif

!   now do field solve...

      temp3d = 0.

!  find rho(kx,ky)
      do k=0,mykm
         n = GCLR*kcnt+k
         do j=0,jm-1
            do i=0,im-1
               temp3d(i,j,k)=drhodt(ip,i,j,k)
            enddo
         enddo
         
         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3d(i,j,k)
            end do
            call ccfft('y',-1,jmx,1.,tmpy,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3d(i,j,k) = tmpy(j)   !rho(ky,x)
            end do
         end do
      end do

      call filtor(temp3d(0:imx-1,0:jmx-1,0:1))

      do k=0,mykm
         do j = 0,jmx-1
            do i = 0,imx-1
               tmpx(i) = temp3d(i,j,k)
            end do
            call ccfft('x',-1,imx,1.,tmpx,tmpx,coefx,workx,0)
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

!  from phi(kx,ky) to phi(x,y)
      do k=0,mykm
         n = GCLR*kcnt+k
         do j = 0,jmx-1
            do i = 0,imx-1
               tmpx(i) = temp3d(i,j,k)
            end do
            call ccfft('x',1,imx,1.,tmpx,tmpx,coefx,workx,0)

            do i = 0,imx-1
               temp3d(i,j,k) = tmpx(i)  ! phi(ky,x)
            end do
         end do

         do i = 0,imx-1
            do j = 0,jmx-1
               tmpy(j) = temp3d(i,j,k)
            end do
            call ccfft('y',1,jmx,1.,tmpy,tmpy,coefy,worky,0)
            do j = 0,jmx-1
               temp3d(i,j,k) = tmpy(j)  ! phi(x,y)
            end do
         end do
      end do

 100  continue
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

 200  continue
!    x-y boundary points 
      call enfxy(dphidt(:,:,:))
!      call enfz(dphidt)
      call filter(dphidt(:,:,:))

      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine jie(ip,n)

!    source quantities are are calculated: n_i
!    right now only ion quantitities are calculated...

      use gem_com
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: enerb,vxdum,dum,xdot,ydot,pidum
      INTEGER :: m,n,i,j,k,l,ns,ip,nonfi,nonfe
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,vte
      REAL(8) :: sz,wght,wght0,r,th,cost,sint,b,qr,dv,kap
      REAL(8) :: xt,yt,rhog,pidum,vpar,xs,dely,vfac
      REAL(8) :: lbfs(0:imx,0:jmx)
      REAL(8) :: rbfs(0:imx,0:jmx)
      REAL(8) :: lbfr(0:imx,0:jmx)
      REAL(8) :: rbfr(0:imx,0:jmx)
      real(8) :: myjpar(0:imx,0:jmx,0:1),myjpex(0:imx,0:jmx,0:1)
      real(8) :: myjpey(0:imx,0:jmx,0:1),myupar(0:imx,0:jmx,0:1)
      real(8) :: myupex(0:imx,0:jmx,0:1),myupey(0:imx,0:jmx,0:1)
      real(8) :: mydnidt(0:imx,0:jmx,0:1),mydnedt(0:imx,0:jmx,0:1)

      nonfi = 1 
      nonfe = 1 
      myjpar = 0.
      myjpex = 0.
      myjpey = 0.
      mydnidt = 0.
      ns=1

      pidum = 1./(pi*2)**1.5*(vwidth)**3
      if(isuni.eq.0)pidum = 1.

      do m=1,mm(1)
         dv=float(lr(1))*(dx*dy*dz)

         r=x3(m)-0.5*lx+lr0
         qr=q0+qp*(r-lr0)
         th=(z3(m)-0.5*lz)/(br0*q0)
         cost=cos(th)
         sint=sin(th)
         b=1.-lr0/br0*cost*tor
         rhog=sqrt(2.*b*mu(m)*mims(1))/(q(1)*b)
         sz=(z3(m)-0.5*lz)*e0*qp/q0*tor
         vfac = 0.5*(mims(1)*u3(m)**2 + 2.*mu(m)*b)
         kap = kapn(1) - (1.5-vfac)*kapt(1)
         vpar = u3(m)
         wght=w3(m)/dv
         wght0 = exp(-vfac)/dv
         if(isuni.eq.0)wght0 = 1./dv

         exp1=0.
         eyp=0.
         delbxp = 0.
         delbyp = 0.
         aparp = 0.
         do 100 l=1,lr(1)
            xs=x3(m)+rwx(1,l)*rhog
            yt=y3(m)+(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            i=int(xt/dx+0.5)
            j=int(yt/dy+0.5)
            k=int(z3(m)/dz+0.5)-gclr*kcnt
            exp1=exp1 + ex(ip,i,j,k)
            eyp=eyp + ey(ip,i,j,k)
            delbxp = delbxp+delbx(ip,i,j,k)
            delbyp = delbyp+delby(ip,i,j,k)
            aparp = aparp+apar(ip,i,j,k)
 100     continue

         exp1=exp1/4.
         eyp=eyp/4.
         delbxp=delbxp/4.
         delbyp=delbyp/4.
         aparp = aparp/4.

         vpar = u3(m)-q(1)/mims(1)*aparp*0.
         enerb=(mu(m)+mims(1)*vpar*vpar/b)/q(1)*tor
         vxdum = (eyp+vpar/b*delbxp)

         xdot = vxdum*nonlin  &
             -tor*enerb*sint/br0
         ydot = (-exp1+vpar/b*delbyp)*nonlin  &
             -tor*b*enerb*(cost/br0+e0*qp*th*sint)

!    now do 1,2,4 point average, where lr is the no. of points...
         do 200 l=1,lr(1)
            xs=x3(m)+rwx(1,l)*rhog
            yt=y3(m)+(rwy(1,l)+sz*rwx(1,l))*rhog
            xt=dmod(xs+8.*lx,lx)
            yt=dmod(yt+80.*ly,ly)
            i=int(xt/dx+0.5)
            j=int(yt/dy+0.5)
            k=int(z3(m)/dz+0.5)-gclr*kcnt
            myjpar(i,j,k) = myjpar(i,j,k)+wght*vpar
            myjpex(i,j,k) = myjpex(i,j,k)+wght*xdot
            myjpey(i,j,k) = myjpey(i,j,k)+wght*ydot
            mydnidt(i,j,k) = mydnidt(i,j,k)+vxdum*kap*wght0
 200     continue

      enddo

!   enforce periodicity
      call enforce(myjpar)
      call enforce(myjpex)
      call enforce(myjpey)
      call enforce(mydnidt)
      call filter(myjpar)
      call filter(mydnidt)

      do 410 i=0,im
         do 420 j=0,jm
            do 430 k=0,mykm
               myjpar(i,j,k) = q(ns)*myjpar(i,j,k)/n0
               myjpex(i,j,k) = q(ns)*myjpex(i,j,k)/n0
               myjpey(i,j,k) = q(ns)*myjpey(i,j,k)/n0
               mydnidt(i,j,k) = q(ns)*mydnidt(i,j,k)/n0*pidum
 430        continue
 420     continue
 410  continue
      call MPI_ALLREDUCE(myjpar(0:im,0:jm,0:1), &
          jpar(0:im,0:jm,0:1), &
          (imx+1)*(jmx+1)*2,MPI_REAL8, &
          MPI_SUM,GRID_COMM,ierr)
      call MPI_ALLREDUCE(myjpex(0:im,0:jm,0:1), &
          jpex(0:im,0:jm,0:1),  &
          (imx+1)*(jmx+1)*2,MPI_REAL8, &
          MPI_SUM,GRID_COMM,ierr)
      call MPI_ALLREDUCE(myjpey(0:im,0:jm,0:1), &
          jpey(0:im,0:jm,0:1),  &
          (imx+1)*(jmx+1)*2,MPI_REAL8, &
          MPI_SUM,GRID_COMM,ierr)     
      call MPI_ALLREDUCE(mydnidt(0:im,0:jm,0:1), &
          dnidt(0:im,0:jm,0:1),  &
          (imx+1)*(jmx+1)*2,MPI_REAL8, &
          MPI_SUM,GRID_COMM,ierr)     
      

! electrons current
      pidum = 1./(pi*2)**1.5*(vwidthe)**3
      if(isuni.eq.0)pidum = 1.

      vte = sqrt(amie)
      myupex = 0.
      myupey = 0.
      mydnedt = 0.
      do m=1,mm(2)
         dv=(dx*dy*dz)

         r=x3e(m)-0.5*lx+lr0
         qr=q0+qp*(r-lr0)
         th=(z3e(m)-0.5*lz)/(br0*q0)
         cost=cos(th)
         sint=sin(th)
         b=1.-lr0/br0*cost*tor
         rhog=0.
         sz=(z3e(m)-0.5*lz)*e0*qp/q0*tor
         vfac = 0.5*(mims(2)*u3e(m)**2 + 2.*mue3(m)*b)
         kap = kapn(2) - (1.5-vfac)*kapt(2)

         wght=w3e(m)/dv
         wght0=exp(-vfac)/dv
         if(isuni.eq.0)wght0=1./dv

         xt = x3e(m)
         yt = y3e(m)
         i=int(xt/dx)
         j=int(yt/dy)
         k=int(z3e(m)/dz)-gclr*kcnt
         wx0=float(i+1)-xt/dx 
         wx1=1.-wx0
         wy0=float(j+1)-yt/dy
         wy1=1.-wy0
         wz0=float(gclr*kcnt+k+1)-z3e(m)/dz
         wz1=1.-wz0

         eyp = wx0*wy0*wz0*ey(ip,i,j,k) + wx1*wy0*wz0*ey(ip,i+1,j,k) &
           + wx0*wy1*wz0*ey(ip,i,j+1,k) + wx1*wy1*wz0*ey(ip,i+1,j+1,k) +&
           wx0*wy0*wz1*ey(ip,i,j,k+1) + wx1*wy0*wz1*ey(ip,i+1,j,k+1) +&
           wx0*wy1*wz1*ey(ip,i,j+1,k+1) + wx1*wy1*wz1*ey(ip,i+1,j+1,k+1)
         delbxp = wx0*wy0*wz0*delbx(ip,i,j,k) &
           + wx1*wy0*wz0*delbx(ip,i+1,j,k) &
           + wx0*wy1*wz0*delbx(ip,i,j+1,k) &
           + wx1*wy1*wz0*delbx(ip,i+1,j+1,k) &
           + wx0*wy0*wz1*delbx(ip,i,j,k+1) &
           + wx1*wy0*wz1*delbx(ip,i+1,j,k+1) &
           + wx0*wy1*wz1*delbx(ip,i,j+1,k+1) &
           + wx1*wy1*wz1*delbx(ip,i+1,j+1,k+1)

         vpar = u3e(m)
         if(abs(vpar/vte).gt.vcut)wght = 0.
         if(abs(vpar/vte).gt.vcut)wght0 = 0.
         enerb=(mue3(m)+mims(2)*vpar*vpar/b)/q(2)*tor
         vxdum = (eyp+vpar/b*delbxp)
         xdot = -enerb*sint/br0*tor
         ydot = -enerb*b*(cost/br0+e0*qp*th*sint)*tor

         i=int(xt/dx)
         j=int(yt/dy)
         k=int(z3e(m)/dz)-gclr*kcnt
         wx0=float(i+1)-xt/dx 
         wx1=1.-wx0
         wy0=float(j+1)-yt/dy
         wy1=1.-wy0
         wz0=float(gclr*kcnt+k+1)-z3e(m)/dz
         wz1=1.-wz0

         myupex(i,j,k)      =myupex(i,j,k)       &
             +wght*wx0*wy0*wz0*xdot 
         myupex(i+1,j,k)    =myupex(i+1,j,k)   &
             +wght*wx1*wy0*wz0*xdot
         myupex(i,j+1,k)    =myupex(i,j+1,k)  &
             +wght*wx0*wy1*wz0*xdot
         myupex(i+1,j+1,k)  =myupex(i+1,j+1,k) &
             +wght*wx1*wy1*wz0*xdot
         myupex(i,j,k+1)    =myupex(i,j,k+1)  &
             +wght*wx0*wy0*wz1*xdot
         myupex(i+1,j,k+1)  =myupex(i+1,j,k+1) &
             +wght*wx1*wy0*wz1*xdot
         myupex(i,j+1,k+1)  =myupex(i,j+1,k+1) &
            +wght*wx0*wy1*wz1*xdot
         myupex(i+1,j+1,k+1)=myupex(i+1,j+1,k+1) &
             +wght*wx1*wy1*wz1*xdot

         myupey(i,j,k)      =myupey(i,j,k)   &
             +wght*wx0*wy0*wz0*ydot
         myupey(i+1,j,k)    =myupey(i+1,j,k)  &
             +wght*wx1*wy0*wz0*ydot
         myupey(i,j+1,k)    =myupey(i,j+1,k)  &
             +wght*wx0*wy1*wz0*ydot
         myupey(i+1,j+1,k)  =myupey(i+1,j+1,k) &
             +wght*wx1*wy1*wz0*ydot
         myupey(i,j,k+1)    =myupey(i,j,k+1)  &
             +wght*wx0*wy0*wz1*ydot
         myupey(i+1,j,k+1)  =myupey(i+1,j,k+1)  &
             +wght*wx1*wy0*wz1*ydot
         myupey(i,j+1,k+1)  =myupey(i,j+1,k+1)  &
             +wght*wx0*wy1*wz1*ydot
         myupey(i+1,j+1,k+1)=myupey(i+1,j+1,k+1)  &
             +wght*wx1*wy1*wz1*ydot

         wght0 = wght0*vxdum*kap
         mydnedt(i,j,k) =mydnedt(i,j,k)+wght0*wx0*wy0*wz0
         mydnedt(i+1,j,k) =mydnedt(i+1,j,k)+wght0*wx1*wy0*wz0
         mydnedt(i,j+1,k) =mydnedt(i,j+1,k)+wght0*wx0*wy1*wz0
         mydnedt(i+1,j+1,k) =mydnedt(i+1,j+1,k)+wght0*wx1*wy1*wz0
         mydnedt(i,j,k+1) =mydnedt(i,j,k+1)+wght0*wx0*wy0*wz1
         mydnedt(i+1,j,k+1) =mydnedt(i+1,j,k+1)+wght0*wx1*wy0*wz1
         mydnedt(i,j+1,k+1) =mydnedt(i,j+1,k+1)+wght0*wx0*wy1*wz1
         mydnedt(i+1,j+1,k+1) =mydnedt(i+1,j+1,k+1)+wght0*wx1*wy1*wz1
      enddo

!   enforce periodicity
      call enforce(myupex(:,:,:))
      call enforce(myupey(:,:,:))
      call enforce(mydnedt(:,:,:))
      call filter(mydnedt(:,:,:))

      do  i=0,im
         do  j=0,jm
            do  k=0,mykm
               myupex(i,j,k) = myupex(i,j,k)/n0
               myupey(i,j,k) = myupey(i,j,k)/n0
               mydnedt(i,j,k) = mydnedt(i,j,k)/n0*pidum
            end do
         end do
      end do
      call MPI_ALLREDUCE(myupex(0:im,0:jm,0:1),  &
          upex(0:im,0:jm,0:1),  &
          (imx+1)*(jmx+1)*2,MPI_REAL8, &
          MPI_SUM,GRID_COMM,ierr)
      call MPI_ALLREDUCE(myupey(0:im,0:jm,0:1), &
          upey(0:im,0:jm,0:1), &
          (imx+1)*(jmx+1)*2,MPI_REAL8, &
          MPI_SUM,GRID_COMM,ierr)
      call MPI_ALLREDUCE(mydnedt(0:im,0:jm,0:1), &
          dnedt(0:im,0:jm,0:1), &
          (imx+1)*(jmx+1)*2,MPI_REAL8, &
          MPI_SUM,GRID_COMM,ierr)

 999  continue
      return
      end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lorentz(ip,n)
      use gem_com
      implicit none
      integer :: ip,k,n,ncol,icol
      real(8) :: edum,vdum,dum,dum1,ptch,vte,r,qr,th,cost,b
      real(8) :: h_x,h_coll,x,eps,dtcol,uold

      ncol = 1
      if(ip.eq.1)dtcol = dt/ncol*2
      if(ip.eq.0)dtcol = dt/ncol
      vte = sqrt(amie)
      if(rneu==0.0)return
      do k = 1,mm(2)
         r=x3e(k)-0.5*lx+lr0
         qr=q0+qp*(r-lr0)
         th=(z3e(k)-0.5*lz)/(br0*q0)
         cost=cos(th)
!         sint=sin(th)
         b=1.-lr0/br0*cost*tor
         uold = u3e(k)
         edum = b*mue3(k)+0.5*mims(2)*u3e(k)*u3e(k)
         vdum = sqrt(2.*edum/mims(2))
         ptch = u3e(k)/vdum
!         if(myid.eq.master)write(*,*)mue(k),vdum,k,ptch
!         if((1.-ptch*ptch).lt.0.)write(*,*)mue3(k),u3e(k),vdum,ptch
         eps = edum
         x = sqrt(eps)
         h_x    = 4.0*x*x/(3.0*sqrt(pi))
         h_coll = h_x/sqrt(1.0+h_x**2)
!         h_coll = exp(-x*x)/(x*sqrt(pi))+(1.0-1.0/(2.0*x*x))*DERF(x)
         dum1 = 1/eps**1.5*(1+h_coll)
         dum = dtcol*rneu*dum1
         if(x<0.316)dum=0.0
         do icol = 1,ncol
            ptch =ptch-dum*ptch+sqrt(dum*(1.-ptch*ptch)) &
                *dsign(1.d0,ran2(iseed)-0.5)
            ptch = dmin1(ptch,0.999d0)
            ptch = dmax1(ptch,-0.999d0)
         end do
         u3e(k) = vdum*ptch
         mue3(k) = 0.5*mims(2)*vdum*vdum*(1.-ptch*ptch)/b
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
      use fft_wrapper
	implicit none

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
	call weight
!     initialize particle quantities...
         if( cut.eq.0.) cut=1.
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

         call ccfft('x',0,imx,0.,tmpx,tmpx,coefx,workx,0)
         call ccfft('y',0,jmx,0.,tmpy,tmpy,coefy,worky,0)
         call ccfft('z',0,kmx,0.,tmpz,tmpz,coefz,workz,0)

         ncurr = 1

	end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine loader_wrapper
         use gem_com
	implicit none

	integer :: n,i,j,k,ip

         if(isuni.eq.1)call unild
         if(isuni.eq.0)call loadi
         if(ifluid.eq.1)call loade
         if(izon.eq.1)then
            do i = 1,200
               call accumulate(0,0)
               call poisson(0,0)
               if(mod(i,20).eq.0)call outd(0)
               if(ifluid.eq.1)call loade
            end do
         end if
         if(idg.eq.1)write(*,*)'past loader'
end subroutine loader_wrapper
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine accumulate(n,ip)
         use gem_com
	implicit none

	integer :: n,i,j,k,ip
	call grid1(ip,n)
	if(idg.eq.1)write(*,*)'pass grid1'
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine accumulate
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine poisson(n,ip)
         use gem_com
	implicit none
	integer :: n,i,j,k,ip,it,iter=1
        do it = 1,iter
           call den0_phi(ip,n,it)
           call gkps(n,ip)
        end do
	if(idg.eq.1)write(*,*)'pass gkps'
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine poisson
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine ampere(n,ip)
         use gem_com
	implicit none
        integer :: iter=7
	integer :: n,i,j,k,ip
        if(ifluid==1.and.beta.gt.1.e-8)then
           do i = 1,iter
              call jpar0(ip,n,i)
              call ezamp(n,ip)
           end do
        end if
	if(idg.eq.1)write(*,*)'pass ezamp'
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine ampere
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine split_weight(n,ip)
         use gem_com
	implicit none

	integer :: n,i,j,k,ip
            if(isg.gt.0..and.ifluid.eq.1)then
               call jie(ip,n)
	if(idg.eq.1)write(*,*)'pass jie'
               call drdt(ip)
	if(idg.eq.1)write(*,*)'pass drdt'
               call dpdt(ip)
	if(idg.eq.1)write(*,*)'pass dpdt'
               call nadi(ip)
	if(idg.eq.1)write(*,*)'pass nadi'
            end if
	if(idg.eq.1)write(*,*)'pass split_weight'
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine split_weight
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine field(n,ip)
         use gem_com
	implicit none
	integer :: n,i,j,k,ip
	call grad(ip)
	call eqmo(ip)
	if(idg.eq.1)write(*,*)'pass grad and eqmo'
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)        
end subroutine field
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine push_wrapper(n,ip)
         use gem_com
	implicit none
	integer :: n,i,j,k,ip

            if(ip.eq.1)call ppush(n)
            if(ip.eq.0)call cpush(n)
            if(idg.eq.1)write(*,*)'pass ppush'

            if(ip.eq.1.and.ifluid==1)call pint
!            if(ip.eq.1.and.ifluid==1)call countw(n)
            if(ip.eq.1.and.ifluid==1)call countu(n)
            if(ip.eq.0.and.ifluid==1)call cint(n)
            if(idg.eq.1)write(*,*)'pass pint'
!            if(ip.eq.0.and.ifluid==1)call countw(n)
            if(ip.eq.0.and.ifluid==1)call countu(n)
            if(ifluid==1)call lorentz(ip,n)
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)        		
end subroutine push_wrapper
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine diagnose(n)
         use gem_com
	implicit none
	integer :: n,i,j,k,ip

        call modes2(phi,pmodehis,n)
        if(idg.eq.1)write(*,*)'pass modes2'  
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)          
        call yveck(phi(0,:,:,:),n)
        if(idg.eq.1)write(*,*)'pass yvec'    
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)        

end subroutine diagnose
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
subroutine reporter(n)
         use gem_com
	implicit none
	integer :: n,i,j,k,ip

        if(mod(n,xnplt).eq.0) then
           call spec(n)
        endif

        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!     save particle arrays for restart if iput=1...
!     do this before the code crashes due to graphics problems
        if((iput.eq.1).and.mod(n+1,100).eq.0)call restart(2,n)

!     periodically make output for plots
        call outd(n)
        if(idg.eq.1)write(*,*)'pass outd'

end subroutine reporter
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine jpar0(ip,n,it)

!    source quantities are are calculated: n_i
!    right now only ion quantitities are calculated...

      use gem_com
      implicit none
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: enerb,vxdum,dum,xdot,ydot,avede0
      INTEGER :: m,n,i,j,k,l,ns,ip,it
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,vte
      REAL(8) :: sz,wght,wght0,r,th,cost,sint,b,qr,dv
      REAL(8) :: xt,yt,zt,rhog,pidum,vpar,xs,dely,vfac,pidum
      REAL(8) :: lbfs(0:imx,0:jmx)
      REAL(8) :: rbfs(0:imx,0:jmx)
      REAL(8) :: lbfr(0:imx,0:jmx)
      REAL(8) :: rbfr(0:imx,0:jmx)
      real(8) :: myupa(0:imx,0:jmx,0:1)

      pidum = 1./(pi*2)**1.5*(vwidthe)**3
      if(isuni.eq.0)pidum = 1.
      ns=1

!    electrons current due to f_M (p_para)
      vte = sqrt(amie)
      myupa = 0.
      upa0(:,:,:) = 0.  !apar(ip,:,:,:)
!      return
      if(it.eq.1)then
         apar(ip,:,:,:) = 0.
         upa0(:,:,:) = 0.
         return
      end if
      do m=1,mm(2)
         dv=(dx*dy*dz)

         r=x3e(m)-0.5*lx+lr0
         qr=q0+qp*(r-lr0)
         th=(z3e(m)-0.5*lz)/(br0*q0)
         cost=cos(th)
         sint=sin(th)
         b=1.-lr0/br0*cost*tor
         rhog=0.
         sz=(z3e(m)-0.5*lz)*e0*qp/q0*tor
         vfac = 0.5*(mims(2)*u3e(m)**2 + 2.*mue3(m)*b )

         wght=w3e(m)/dv
         if(abs(w3e(m)).gt.0.5)wght = 0.
         wght0 = exp(-vfac)/dv
         if(isuni.eq.0)wght0 = 1./dv
         vpar = u3e(m) !linearly correct
         if(abs(vpar/vte).gt.vcut)wght = 0.
         if(abs(vpar/vte).gt.vcut)wght0 = 0.

         xt=x3e(m)
         yt=y3e(m)
         zt=z3e(m)
         i=int(xt/dx)
         j=int(yt/dy)
         k=int(z3e(m)/dz)-gclr*kcnt
         wx0=float(i+1)-xt/dx 
         wx1=1.-wx0
         wy0=float(j+1)-yt/dy
         wy1=1.-wy0
         wz0=float(gclr*kcnt+k+1)-z3e(m)/dz
         wz1=1.-wz0

!         aparp = 0.
!         do l = 1,8
!            aparp = aparp+apk(l)*exp(IU*pi2*   &
!                (lapa(l)*xt/lx+mapa(l)*yt/ly+napa(l)*zt/lz))
!         end do
         aparp = wx0*wy0*wz0*apar(ip,i,j,k)  &
             + wx1*wy0*wz0*apar(ip,i+1,j,k) &
             + wx0*wy1*wz0*apar(ip,i,j+1,k) &
             + wx1*wy1*wz0*apar(ip,i+1,j+1,k) &
             + wx0*wy0*wz1*apar(ip,i,j,k+1) &
             + wx1*wy0*wz1*apar(ip,i+1,j,k+1) &
             + wx0*wy1*wz1*apar(ip,i,j+1,k+1) &
             + wx1*wy1*wz1*apar(ip,i+1,j+1,k+1)

         wght0 = wght0*aparp*vpar*vpar/amie 

         myupa(i,j,k)      =myupa(i,j,k)  &
             +wght0*wx0*wy0*wz0
         myupa(i+1,j,k)    =myupa(i+1,j,k) &
             +wght0*wx1*wy0*wz0
         myupa(i,j+1,k)    =myupa(i,j+1,k)  &
             +wght0*wx0*wy1*wz0
         myupa(i+1,j+1,k)  =myupa(i+1,j+1,k) &
             +wght0*wx1*wy1*wz0
         myupa(i,j,k+1)    =myupa(i,j,k+1)  &
             +wght0*wx0*wy0*wz1
         myupa(i+1,j,k+1)  =myupa(i+1,j,k+1) &
             +wght0*wx1*wy0*wz1
         myupa(i,j+1,k+1)  =myupa(i,j+1,k+1) &
             +wght0*wx0*wy1*wz1
         myupa(i+1,j+1,k+1)=myupa(i+1,j+1,k+1) &
             +wght0*wx1*wy1*wz1
      enddo

!   enforce periodicity
      call enforce(myupa(:,:,:))
      call filter(myupa(:,:,:))

      do  i=0,im
         do  j=0,jm
            do  k=0,mykm
               myupa(i,j,k)= myupa(i,j,k)/n0*pidum*ifluid
            end do
         end do
      end do
      call MPI_ALLREDUCE(myupa(0:im,0:jm,0:1),  &
          upa0(0:im,0:jm,0:1),  &
          (imx+1)*(jmx+1)*2,MPI_REAL8,  &
          MPI_SUM,GRID_COMM,ierr) 

      upa0(:,:,:) = upa0(:,:,:)+ &
                   (isg*phi(ip,:,:,:)+dene(ip,:,:,:))*apar(ip,:,:,:)*nonline
 999  continue

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine filtor(u)   
      use gem_com
      use fft_wrapper
      implicit none
      complex(8) :: lbfs(0:imx-1,0:jmx-1)
      complex(8) :: rbfr(0:imx-1,0:jmx-1)
      complex(8) :: u(0:imx-1,0:jmx-1,0:1)
      REAL(8) :: qr
      REAL(8) :: formz(0:kmx-1),smfac(0:kmx-1),kx,ky
      real(8) :: aky(0:jmx-1),dely(0:imx-1),dklz(0:imx-1,0:jmx-1)
      INTEGER :: i,j,k,l,m,n,ind
      INTEGER :: myk,mynum,id
      COMPLEX(8) :: temp3d(0:imx-1,0:jmx-1,0:1)
      complex(8),dimension(:),allocatable :: holdu,rbuf,sbuf

!      write(*,*)'enter filtor'
      mynum = imx*jmx/numprocs      
      allocate(holdu(0:imx*jmx*kmx/numprocs-1),  &
               rbuf(0:mynum-1), &
               sbuf(0:mynum-1))
      do 51 i = 0,im-1
         qr = q0+qp*(xg(i)-0.5*lx)
         dely(i) = dmod(-pi2*lr0/q0*qr+80.*ly,ly)*tor*0.
 51   continue
      do 52 j = 0,jm-1
         if(j.ge.(jm/2+1)) then
            aky(j) = -2.*pi*float(jm-j)/ly
         else
            aky(j) = 2.*pi*float(j)/ly
         end if
 52   continue
      do 53 k = 0,km-1
         if(k.ge.(km/2+1)) then
            m=km-k
         else
            m=k
         endif
         formz(k) = 1./float(kmx)
         if(m.gt.4)formz(k) = 0.
!         if(m==0)formz(k) = 0.
 53   continue
      do 54 j = 0,jm-1
         do i = 0,im-1
            dklz(i,j) = dmod(aky(j)*dely(i)+80*pi*2,pi*2)
            if(dklz(i,j)>pi)dklz(i,j)=dklz(i,j)-pi2
         end do
 54   continue
      i = 3
      do 55 k = 0,km-1
         if(k<i) then
            smfac(k)=1.!float(k)/float(i)
         else if(k>km-i)then
            smfac(k)=1.!float(km-k)/float(i)
         else   
            smfac(k)=1.
         endif
 55   continue

!      return
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
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

       call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!fft and fft back holdu to do filtering
       do ind = 0,mynum-1
          j = (myid*mynum+ind)/im
          i = myid*mynum+ind-j*im
          do k = 0,km-1
             tmpz(k) = holdu(ind*km+k) &
                            *exp(-IU*dklz(i,j)*k*dz/lz)
          end do

          call ccfft('z',-1,kmx,1.,tmpz,tmpz,coefz,workz,0)

          do k = 0,km-1
             tmpz(k) = tmpz(k)*formz(k)
          end do

          call ccfft('z',+1,kmx,1.,tmpz,tmpz,coefz,workz,0)

          do k = 0,km-1
             holdu(ind*km+k) = tmpz(k)  &
                              *exp(IU*dklz(i,j)*k*dz/lz) &
                              *smfac(k)
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
                u(i,j,1) = u(i,j,0) !rbfr(i,j)*exp(IU*aky(j)*dely(i))
             end do
          end do
       end if
          
!       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       return
       end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine den0_phi(ip,n,it)

!    source quantities are are calculated: n_i
!    right now only ion quantitities are calculated...

      use gem_com
      REAL(8) :: phip,exp1,eyp,ezp,delbxp,delbyp,dpdzp,dadzp,aparp
      REAL(8) :: enerb,vxdum,dum,xdot,ydot,avede0
      INTEGER :: m,n,i,j,k,l,ns,ip,it
      REAL(8) :: wx0,wx1,wy0,wy1,wz0,wz1,vte
      REAL(8) :: sz,wght,wght0,r,th,cost,sint,b,qr,dv
      REAL(8) :: xt,yt,zt,rhog,pidum,vpar,xs,dely,vfac,pidum
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
         den0(:,:,:) = phi(ip,:,:,:) 
         return
      end if
      do m=1,mm(2)
         dv=(dx*dy*dz)

         r=x3e(m)-0.5*lx+lr0
         qr=q0+qp*(r-lr0)
         th=(z3e(m)-0.5*lz)/(br0*q0)
         cost=cos(th)
         sint=sin(th)
         b=1.-lr0/br0*cost*tor
         rhog=0.
         sz=(z3e(m)-0.5*lz)*e0*qp/q0*tor
         vfac = 0.5*(mims(2)*u3e(m)**2 + 2.*mue3(m)*b )

         wght=w3e(m)/dv
         wght0 = exp(-vfac)/dv
         vpar = u3e(m)
         if(isuni.eq.0)wght0 = 1./dv
         if(abs(vpar/vte).gt.vcut)wght = 0.
         if(abs(vpar/vte).gt.vcut)wght0 = 0.

         xt=x3e(m)
         yt=y3e(m)
         zt=z3e(m)
         i=int(xt/dx)
         j=int(yt/dy)
         k=int(z3e(m)/dz)-gclr*kcnt
         wx0=float(i+1)-xt/dx 
         wx1=1.-wx0
         wy0=float(j+1)-yt/dy
         wy1=1.-wy0
         wz0=float(gclr*kcnt+k+1)-z3e(m)/dz
         wz1=1.-wz0

         phip = 0.
         phip = wx0*wy0*wz0*phi(ip,i,j,k)  &
             + wx1*wy0*wz0*phi(ip,i+1,j,k) &
             + wx0*wy1*wz0*phi(ip,i,j+1,k) &
             + wx1*wy1*wz0*phi(ip,i+1,j+1,k) &
             + wx0*wy0*wz1*phi(ip,i,j,k+1) &
             + wx1*wy0*wz1*phi(ip,i+1,j,k+1) &
             + wx0*wy1*wz1*phi(ip,i,j+1,k+1) &
             + wx1*wy1*wz1*phi(ip,i+1,j+1,k+1)

         wght0 = wght0*phip

         myden0(i,j,k)      =myden0(i,j,k) &
             +wght0*wx0*wy0*wz0
         myden0(i+1,j,k)    =myden0(i+1,j,k) &
             +wght0*wx1*wy0*wz0
         myden0(i,j+1,k)    =myden0(i,j+1,k) &
             +wght0*wx0*wy1*wz0
         myden0(i+1,j+1,k)  =myden0(i+1,j+1,k) &
             +wght0*wx1*wy1*wz0
         myden0(i,j,k+1)    =myden0(i,j,k+1) &
             +wght0*wx0*wy0*wz1
         myden0(i+1,j,k+1)  =myden0(i+1,j,k+1) &
             +wght0*wx1*wy0*wz1
         myden0(i,j+1,k+1)  =myden0(i,j+1,k+1) &
             +wght0*wx0*wy1*wz1
         myden0(i+1,j+1,k+1)=myden0(i+1,j+1,k+1) &
             +wght0*wx1*wy1*wz1
      enddo

!   enforce periodicity
      call enforce(myden0(:,:,:))

      do  i=0,im
         do  j=0,jm
            do  k=0,mykm
               myden0(i,j,k)= myden0(i,j,k)/n0*pidum*ifluid
            end do
         end do
      end do
      call MPI_ALLREDUCE(myden0(0:im,0:jm,0:1), &
          den0(0:im,0:jm,0:1), &
          (imx+1)*(jmx+1)*2,MPI_REAL8, &
          MPI_SUM,GRID_COMM,ierr)
      return
      end
!!!!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine countw(n)
      use gem_com
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
      do m=1,mm(2)
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
      do m = 1,mm(2)
         i = int(((w3e(m))-minw)/dw)
         mynpbin(i) = mynpbin(i)+1
         mywpbin(i) = mywpbin(i)+abs(w3e(m))
         if(abs(w3e(m)).gt.4.)then
            w3e(m) = 0.
            w2e(m) = 0.
            w1e(m) = 0.
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine countu(n)
      use gem_com
      implicit none

      INTEGER :: n,nbin,m,i,j,k
      parameter(nbin=200)
      REAL(8) :: myavue,avue,mymaxw,myminw,maxw,minw,vte,vratio
      REAL(8) :: sbuf(10),rbuf(10)
      REAL(8),dimension(:),allocatable :: wghtmax,wghtmin
      REAL(8) :: du,vpar,vperp2,r,qr,th,b,cost,dum1,dum2
      integer :: npbin(0:nbin),mynpbin(0:nbin)
      real(8) :: wpbin(0:nbin),mywpbin(0:nbin)

      allocate (wghtmax(0:numprocs-1),wghtmin(0:numprocs-1))

      cnt=int(tmm(1)/numprocs)
      vratio = vwidthe/vwidth
      vte=sqrt(amie)
      sbuf(1:10) = 0.
      rbuf(1:10) = 0.
      myavue = 0.
      avue = 0.
      mymaxw = 0.
      myminw = 0.
      maxw = 0.
      minw = 0.
      goto 100
      do m=1,mm(2)
         if((u3e(m)).gt.mymaxw)mymaxw = (u3e(m))
         if((u3e(m)).lt.myminw)myminw = (u3e(m))
         myavue = myavue+abs(u3e(m))
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

      du = (maxw-minw)/nbin
      mynpbin = 0
      npbin = 0
      wpbin = 0.
      mywpbin = 0.
      
 100  do m = 1,mm(2)
         i = int(((u3e(m))-minw)/du)
         mynpbin(i) = mynpbin(i)+1
         mywpbin(i) = mywpbin(i)+abs(u3e(m))
         if(abs(u3e(m)/vte).gt.vcut.and.isuni.eq.0)then
            call parperp(vpar,vperp2,m,pi,cnt,MyId)
            r=x2e(m)-0.5*lx+lr0
            qr=q0+qp*(r-lr0)
            th=(z2e(m)-0.5*lz)/(br0*q0)
            cost=cos(th)
            b=1.-lr0/br0*cost*tor
            u2e(m)=vpar/sqrt(mims(1))*sqrt(mims(1)*amie)*vratio
            mue3(m)=tets(2)*0.5*vperp2/b*vratio**2
            mue2(m) = mue3(m)
            mue1(m) = mue3(m)
            w3e(m) = 0.
            u3e(m) = u2e(m)
            u1e(m) = u2e(m)
            w2e(m) = 0.
            w1e(m) = 0.
         end if
         if(abs(u3e(m)/vte).gt.5.and.isuni.eq.1)then
            dum1 = vwidth*(ran2(iseed)-0.5)
            dum2 = vwidth*(ran2(iseed)-0.5)
            vperp2 = dum1**2+dum2**2
            vpar = vwidth*(ran2(iseed)-0.5)
            r=x3e(m)-0.5*lx+lr0
            qr=q0+qp*(r-lr0)
            th=(z3e(m)-0.5*lz)/(br0*q0)
            cost=cos(th)
            b=1.-lr0/br0*cost*tor
            u2e(m) = sqrt(mims(1)*amie)*vpar/sqrt(mims(1))*vratio
            mue3(m) = tets(2)*0.5*vperp2/b*vratio**2
            mue2(m) = mue3(m)
            mue1(m) = mue3(m)
            w3e(m) = 0.
            u3e(m) = u2e(m)
            u1e(m) = u2e(m)
            w2e(m) = 0.
            w1e(m) = 0.
         end if
      end do
      return
      call MPI_ALLREDUCE(mynpbin(0:nbin),npbin(0:nbin),nbin+1, &
          MPI_integer, &
          MPI_SUM, &
          MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(mywpbin(0:nbin),wpbin(0:nbin),nbin+1, &
          MPI_real8, &
          MPI_SUM, &
          MPI_COMM_WORLD,ierr)

      call MPI_ALLREDUCE(myavue,avue,1, &
          MPI_real8, &
          MPI_SUM, &
          MPI_COMM_WORLD,ierr)
      avue = avue/tmm(1)
      wpbin = wpbin/tmm(1)
      if(myid.eq.0.and.mod(n,xnplt).eq.0)then
!         write(*,*)'maxw,minw,avue= ', maxw, minw,avue
         do i = 0,nbin
!            write(*,*)n,(minw+i*du)/avue,npbin(i),wpbin(i)
         end do
      end if

      return
      end
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine filtky(u)   
      use gem_com
      use fft_wrapper
      implicit none
      complex(8) :: u(0:imx-1,0:jmx-1,0:1)
      REAL(8) :: qr
      real(8) :: aky(0:jmx-1),dely(0:imx-1)
      INTEGER :: i,j,k,l,m,n,ind,NNI
      INTEGER :: myk,mynum,id
      COMPLEX(8) :: temp(0:imx-1,0:jmx-1,0:1)
      COMPLEX(8) :: pfacl(0:imx-1,0:jmx-1),pfacr(0:imx-1,0:jmx-1)
      real(8) ::  DW1,DW2
      parameter(NNI = 3)
      parameter(DW1 = 0.5)
      parameter(DW2 = -1./6.)
      COMPLEX(8) :: ttemp3d(0:1,0:imx-1,0:jmx-1)
      COMPLEX(8) :: htemp3d(0:1,0:imx-1,0:jmx-1)
      COMPLEX(8) :: lbfs(0:imx-1,0:jmx-1) 
      COMPLEX(8) :: lbfr(0:imx-1,0:jmx-1)
      COMPLEX(8) :: rbfs(0:imx-1,0:jmx-1)
      COMPLEX(8) :: rbfr(0:imx-1,0:jmx-1)

!      write(*,*)'enter filtky'
      do 51 i = 0,im-1
         qr = q0+qp*(xg(i)-0.5*lx)
         dely(i) = dmod(-pi2*lr0/q0*qr+80.*ly,ly)*tor
 51   continue
      do 52 j = 0,jm-1
         if(j.ge.(jm/2+1)) then
            aky(j) = -2.*pi*float(jm-j)/ly
         else
            aky(j) = 2.*pi*float(j)/ly
         end if
 52   continue

      pfacl(:,:)=1.
      pfacr(:,:)=1.
      if(gclr==0)then
         do  j = 0,jm-1
            do i = 0,im-1
               pfacl(i,j) = exp(-IU*aky(j)*dely(i))
               pfacr(i,j) = 1.
            end do
         end do
      end if
      if(gclr==glst)then
         do  j = 0,jm-1
            do i = 0,im-1
               pfacl(i,j) = 1.
               pfacr(i,j) = exp(IU*aky(j)*dely(i))
            end do
         end do
      end if
      
      do i = 0,im-1
         do j = 0,jm-1
            do k = 0,1
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

	 do n=1,NNI

	 htemp3d=ttemp3d

!  build buffers and send grid info to neighbors

	 do j=0,jm-1
	   do i=0,im-1
	     rbfs(i,j)=htemp3d(0,i,j)
	     lbfs(i,j)=htemp3d(0,i,j)
	   enddo
	 enddo

	 call MPI_SENDRECV(rbfs(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_DOUBLE_COMPLEX,rngbr,9, &
     	                   lbfr(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_DOUBLE_COMPLEX,lngbr,9,  &
     		           TUBE_COMM,stat,ierr)

	 call MPI_SENDRECV(lbfs(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_DOUBLE_COMPLEX,lngbr,8,  &
     	                   rbfr(0:imx-1,0:jmx-1),imx*jmx,  &
     		MPI_DOUBLE_COMPLEX,rngbr,8,  &
     		           TUBE_COMM,stat,ierr)
!	 write(*,*)'past buff'

	 do j=0,jm-1
	   do i=0,im-1
	     ttemp3d(0,i,j)=( htemp3d(0,i,j)   &
     	       +DW1*(pfacl(i,j)*lbfr(i,j) + pfacr(i,j)*rbfr(i,j)))/ (1.+2*DW1)
	   enddo
	 enddo

	 htemp3d=ttemp3d

!  build buffers and send grid info to neighbors

	 do j=0,jm-1
	   do i=0,im-1
	     rbfs(i,j)=htemp3d(0,i,j)
	     lbfs(i,j)=htemp3d(0,i,j)
	   enddo
	 enddo

	 call MPI_SENDRECV(rbfs(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_DOUBLE_COMPLEX,rngbr,9, &
     	                   lbfr(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_DOUBLE_COMPLEX,lngbr,9, &
     		           TUBE_COMM,stat,ierr) 

	 call MPI_SENDRECV(lbfs(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_DOUBLE_COMPLEX,lngbr,8,  &
     	                   rbfr(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_DOUBLE_COMPLEX,rngbr,8, &
     		           TUBE_COMM,stat,ierr)
!	 write(*,*)'past buff'

	 do j=0,jm-1
	   do i=0,im-1
	     ttemp3d(0,i,j)=(  htemp3d(0,i,j) &
     	       +DW2*(pfacl(i,j)*lbfr(i,j) + pfacr(i,j)*rbfr(i,j)))/ (1.+2*DW2)
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
     		MPI_DOUBLE_COMPLEX,lngbr,10, &
     	                   rbfr(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_DOUBLE_COMPLEX,rngbr,10, &
     		           TUBE_COMM,stat,ierr)

      do i = 0,im-1
         do j = 0,jm-1
	    ttemp3d(mykm,i,j) = rbfr(i,j)*pfacr(i,j)
         end do
      end do


      do k=0,mykm
         do j=0,jm-1
            do i=0,im-1
	       temp(i,j,k)=ttemp3d(k,i,j)
            enddo
         enddo
      enddo

 100  continue
      do i = 0,im-1
         do j = 0,jm-1
            do k = 0,mykm
               u(i,j,k) = temp(i,j,k)
            end do
         end do
      end do

       return
       end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
