subroutine gkps_adiabatic_electron(nstep,ip)
  use gem_com
  use gem_equil
  use gem_fft_wrapper
  implicit none

  ! input variables
  integer :: nstep,ip

  ! global variables for matrix implement
  integer :: ifirst
  integer,dimension(:,:,:,:),allocatable :: ipiv,ipiv_zonal
  real,dimension(:,:,:),allocatable :: formphi
  complex,dimension(:,:,:,:),allocatable :: mx,mx_zonal

  ! local variables for matrix implement
  integer :: i,j,k,l,m,n,i1,k1,m1,ns,ix,ikx,icount,INFO
  real :: r,th,wx0,wx1,wz0,wz1,kx1,kx2,ky,gam0,gam1,sgny, &
          bfldp,dydrp,qhatp,grp,gthp,gxdgyp,grdgtp
  real,dimension(1:5) :: b1,b2,ter
  real :: gamb1(1:5,0:imx-1,0:jcnt-1,0:imx-1,0:1),gamb2(1:5,0:imx-1,0:jcnt-1,0:imx-1,0:1), &
          jacob_theta(imx-1,0:1),jacob_theta_tmp(imx-1,0:1)
  complex :: mx0(imx-1,imx-1,0:jcnt-1,0:1),mx_tmp(imx-1,imx-1,0:jcnt-1,0:1)
  character(len=70) fname
  character(len=5) holdmyid

  ! the local variables for field solver
  integer :: myj
  real :: u(0:imx,0:jmx,0:1),u_zonal(0:imx,0:jmx,0:1), &
          rho_zonal_tmp(0:imx,0:1),rho_zonal(0:imx,0:1),jacob_rho_tmp(0:imx,0:1),jacob_rho(0:imx,0:1)
  complex :: v(0:imx-1,0:jcnt-1,0:1),v_zonal(0:imx-1,0:jcnt-1,0:1),sl(1:imx-1,0:jcnt-1,0:1),temp3dxy(0:imx-1,0:jmx-1,0:1)
 
  ! save the global variables
  save ifirst,mx,mx_zonal,ipiv,ipiv_zonal,formphi

  if(idg==1)write(*,*)'enter gkps'

  write(holdmyid,'(I5.5)') MyId
  fname='./matrix/'//'mx_phi_'//holdmyid

  ! form factors....
  if(ifirst.ne.-99)then
    allocate(mx(imx-1,imx-1,0:jcnt-1,0:1),mx_zonal(imx-1,imx-1,0:jcnt-1,0:1))
    allocate(ipiv(imx-1,imx-1,0:jcnt-1,0:1),ipiv_zonal(imx-1,imx-1,0:jcnt-1,0:1))
    allocate(formphi(0:imx-1,0:jcnt-1,0:1))

    mx=0.0
    mx_zonal=0.0
    mx0=0.0
    mx_tmp=0.0

    ! restart the matrix
    if(iget.eq.1)then
      open(10000+MyId,file=fname,form='unformatted',status='old')
      read(10000+MyId)mx,mx_zonal,ipiv,ipiv_zonal
      close(10000+myid)
      goto 200
    end if

    ! construct the matrix for poisson solver
    do l=0,im-1
       do m=0,jcnt-1
          j=tclr*jcnt+m
          do i=0,im-1
             r=lr0-lx/2+i*dx
             i1=int((r-rin)/dr)
             i1=min(i1,nr-1)
             wx0=(rin+(i1+1)*dr-r)/dr
             wx1=1.-wx0

             do ns=1,nsm
                ter(ns)=wx0*t0s(ns,i1)+wx1*t0s(ns,i1+1)
             enddo

             do n=0,1
                k=gclr*kcnt+n
                k1=int(k*dz/delz)
                k1=min(k1,ntheta-1)
                wz0=((k1+1)*delz-dz*k)/delz
                wz1=1-wz0
                th=wz0*thfnz(k1)+wz1*thfnz(k1+1)
                k=int((th+pi)/dth)
                k=min(k,ntheta-1)
                wz0=(-pi+(k+1)*dth-th)/dth
                wz1=1.-wz0

                bfldp=wx0*wz0*bfld(i1,k)+wx0*wz1*bfld(i1,k+1) &
                      +wx1*wz0*bfld(i1+1,k)+wx1*wz1*bfld(i1+1,k+1)
                dydrp=wx0*wz0*dydr(i1,k)+wx0*wz1*dydr(i1,k+1) &
                      +wx1*wz0*dydr(i1+1,k)+wx1*wz1*dydr(i1+1,k+1)
                qhatp=wx0*wz0*qhat(i1,k)+wx0*wz1*qhat(i1,k+1) &
                      +wx1*wz0*qhat(i1+1,k)+wx1*wz1*qhat(i1+1,k+1)
                grp=wx0*wz0*gr(i1,k)+wx0*wz1*gr(i1,k+1) &
                    +wx1*wz0*gr(i1+1,k)+wx1*wz1*gr(i1+1,k+1)
                gthp=wx0*wz0*gth(i1,k)+wx0*wz1*gth(i1,k+1) &
                     +wx1*wz0*gth(i1+1,k)+wx1*wz1*gth(i1+1,k+1)
                gxdgyp=wx0*wz0*gxdgy(i1,k)+wx0*wz1*gxdgy(i1,k+1) &
                       +wx1*wz0*gxdgy(i1+1,k)+wx1*wz1*gxdgy(i1+1,k+1)
                grdgtp=wx0*wz0*grdgt(i1,k)+wx0*wz1*grdgt(i1,k+1) &
                       +wx1*wz0*grdgt(i1+1,k)+wx1*wz1*grdgt(i1+1,k+1)
                 
                m1=mstart+int((float(m)+1.0)/2)
                if(m==0)m1=0
                sgny=isgnft(m)

                ky=sgny*2.*pi*float(m1)/ly
                kx1=pi*float(l)/lx
                kx2=-pi*float(l)/lx
                 
                ! the b+ and b- value
                do ns=1,nsm
                   b1(ns)=mims(ns)*(kx1*kx1*grp**2 + &
                          ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                          +2*dydrp*lr0/q0*qhatp*grdgtp) &
                          +2*kx1*ky*gxdgyp)/(bfldp*bfldp)*ter(ns)/(q(ns)*q(ns))

                   b2(ns)=mims(ns)*(kx2*kx2*grp**2 + &
                          ky*ky*(dydrp**2*grp**2+(lr0/q0*qhatp*gthp)**2 &
                          +2*dydrp*lr0/q0*qhatp*grdgtp) &
                          +2*kx2*ky*gxdgyp)/(bfldp*bfldp)*ter(ns)/(q(ns)*q(ns))
                   ! construct the Gamma(b+) and Gamma(b-)
                   call srcbes(b1(ns),gam0,gam1)
                   gamb1(ns,l,m,i,n)=gam0
                   call srcbes(b2(ns),gam0,gam1)
                   gamb2(ns,l,m,i,n)=gam0
                enddo

                ! formfactor in gkps
                formphi(l,m,n) = 1./jmx
                if(abs(ky)>kycut)formphi(l,m,n) = 0.
                if(m1.ne.mlk.and.onemd==1)formphi(l,m,n) = 0.

             enddo
          enddo
       enddo
    enddo

    ! construct the matrix MX for Poisson sovler
    do k=0,1
       do j=0,jcnt-1
          do i=1,imx-1
             r=lr0-lx/2+i*dx
             i1=int((r-rin)/dr)
             i1=min(i1,nr-1)
             wx0=(rin+(i1+1)*dr-r)/dr
             wx1=1.-wx0
             do ns=1,nsm
                ter(ns)=wx0*t0s(ns,i1)+wx1*t0s(ns,i1+1)
             enddo
             do ix=1,imx-1
                mx(i,ix,j,k)=0.0
                mx0(i,ix,j,k)=0.0
                ! the adiabatic response of electron, the option is fradi
                if(i==ix)mx(i,ix,j,k) = fradi*gn0e(i)*cn0e/gt0e(i)
                do ikx=0,imx-1
                   do ns=1,nsm
                      ! the matrix MX with adiabatic electron
                      mx(i,ix,j,k)=mx(i,ix,j,k)+q(ns)*sin(ix*ikx*pi/imx)* &
                                   ((1-gamb1(ns,ikx,j,i,k))*exp(IU*ikx*i*pi/imx)- &
                                   (1-gamb2(ns,ikx,j,i,k))*exp(-IU*ikx*i*pi/imx)) &
                                   /ter(ns)*cn0s(ns)*gn0s(ns,i)/(IU*imx)
                      ! the matrix MX0 without adiabatic electron for the calcaulation of MX_zonal jycheng
                      mx0(i,ix,j,k)=mx0(i,ix,j,k)+q(ns)*sin(ix*ikx*pi/imx)* &
                                    ((1-gamb1(ns,ikx,j,i,k))*exp(IU*ikx*i*pi/imx)- &
                                    (1-gamb2(ns,ikx,j,i,k))*exp(-IU*ikx*i*pi/imx)) &
                                    /ter(ns)*cn0s(ns)*gn0s(ns,i)/(IU*imx)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    ! begin to calculate the MX_zonal jycheng
 
    mx_zonal=0.0
    do i=1,imx-1
       do k=0,1
          do ix=1,imx-1
             ! only consider the ky=0 part
             mx_zonal(i,ix,0,k)=mx0(i,ix,0,k)*jac(i,k)
             jacob_theta(i,k)=jac(i,k)
          enddo
       enddo
    enddo

     ! flux average of mx0
     mx_tmp=0.0
     icount=(imx-1)*(imx-1)*jcnt*2        
     call MPI_ALLREDUCE(mx_zonal,mx_tmp,icount,MPI_DOUBLE_COMPLEX,MPI_SUM,TUBE_COMM,ierr)
     mx_zonal=mx_tmp
     jacob_theta_tmp=0.0
     icount=(imx-1)*2
     call MPI_ALLREDUCE(jacob_theta,jacob_theta_tmp,icount,MPI_REAL8,MPI_SUM,TUBE_COMM,ierr)
     jacob_theta=jacob_theta_tmp
     do i=1,imx-1
        do k=0,1
           do ix=1,imx-1
              do j=0,jcnt-1
                 mx_zonal(i,ix,j,k)=mx_zonal(i,ix,j,k)/jacob_theta(i,k)
              enddo
           enddo
        enddo
     enddo


     mx_zonal(:,:,0,:)=cmplx(real(mx_zonal(:,:,0,:)),0.0)
     mx(:,:,0,:)=cmplx(real(mx(:,:,0,:)),0.0)

     do k=0,1
        do j=0,jcnt-1
           call ZGETRF(imx-1,imx-1,mx(:,:,j,k),imx-1,ipiv(:,:,j,k),INFO )
           call ZGETRF(imx-1,imx-1,mx_zonal(:,:,j,k),imx-1,ipiv_zonal(:,:,j,k),INFO )
        enddo
     enddo
     
     ! save the matrix for restart
     if(iget.eq.0) then
        open(10000+MyId,file=fname,form='unformatted',status='unknown')
        write(10000+MyId)mx,mx_zonal,ipiv,ipiv_zonal
        close(10000+myid)
        goto 200
     end if
     
     ! write the mx,mx_zonal with ky=0,k=0
     if(gclr==kmx/2.and.tclr==0.and.nstep==0)then
       open(20,file="mx",status='unknown')
       j=0
       k=0
       do i=1,imx-1
          do ix=1,imx-1
             if(abs(i-ix)<40)write(20,10)i,ix,mx(i,ix,j,k),mx_zonal(ix,i,j,k)
          enddo
       enddo
10     format(1x,i5,1x,i5,2(2x,e12.5,2x,e12.5))
       close(20)
     end if

200  ifirst=-99
  endif

  if(idg==1)write(*,*)'pass form factors'

  !now do field solve...

  ! calcaulte the zonal part of density
  do i=0,im
     do k=0,1
        rho_zonal_tmp(i,k)=0.0
        jacob_rho_tmp(i,k)=0.0
        do j=0,jm-1
           rho_zonal_tmp(i,k)=rho_zonal_tmp(i,k)+rho(i,j,k)*jac(i,k)
           jacob_rho_tmp(i,k)=jacob_rho_tmp(i,k)+jac(i,k)
        enddo
     enddo
  enddo

  rho_zonal=0.0
  icount=(im+1)*2
  call MPI_ALLREDUCE(rho_zonal_tmp,rho_zonal,icount,MPI_REAL8,MPI_SUM,mpi_comm_world,ierr)
  jacob_rho=0.0
  icount=(im+1)*2
  call MPI_ALLREDUCE(jacob_rho_tmp,jacob_rho,icount,MPI_REAL8,MPI_SUM,TUBE_COMM,ierr)
  do i=0,im
     do k=0,1
        rho_zonal(i,k)=rho_zonal(i,k)/jacob_rho(i,k)
     enddo
  enddo
 
  ! different source term for zonal and non-zonal term
  do i=0,nxpp
     do j=0,jm
        do k=0,1
           u(i,j,k)=rho(i,j,k)
           ! expand the zonal source term into x and y dimension
           u_zonal(i,j,k)=rho_zonal(i,k)
        enddo
     enddo
  enddo


  ! expand the source term into sin tranformation
  call dcmpy(u(0:imx-1,0:jmx-1,0:1),v)
  call dcmpy(u_zonal(0:imx-1,0:jmx-1,0:1),v_zonal)

  temp3dxy=0.0
  ! solve the zonal-phi
  do k=0,1
     j=0
     myj=jft(j)
     sl(1:imx-1,j,k)=v_zonal(1:imx-1,j,k)
     ! solve the MX_zonal 
     call ZGETRS('N',imx-1,1,mx_zonal(:,:,j,k),imx-1,ipiv_zonal(:,:,j,k), &
             sl(:,j,k),imx-1,INFO)
     temp3dxy(1:imx-1,myj,k)=sl(1:imx-1,j,k)
     temp3dxy(0,myj,k)=0.
  enddo

  !  from rho(kx,ky) to phi(kx,ky)
  do k = 0,1
     do j = 0,jmx-1
        do i = 0,imx-1
           temp3dxy(i,j,k)=temp3dxy(i,j,k)/jmx !*formphi(i,j,k)
        enddo
     enddo
  enddo

  !  from phi(kx,ky) to phi(x,y)
  do k=0,mykm
     do i=0,imx-1
        do j=0,jmx-1
          tmpy(j)=temp3dxy(i,j,k)
        enddo
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
        do j=0,jmx-1
           temp3dxy(i,j,k)=tmpy(j)  ! phi(x,y)
        enddo
     enddo
  enddo
 
  ! here phi is zonal-phi
  do i=0,nxpp-1
     do j=0,jm-1
        do k=0,mykm
           phi(i,j,k)=temp3dxy(i,j,k)
        enddo
     enddo
  enddo
  call enfz(phi(:,:,:))

  ! for the solver for ky=0 with ky!=0 part together, the difference is their source term, 
  ! for ky=0 term the source term is \delta n_i (ky=0) +<\phi>, for ky!=0 term, the source term is \delta n_i (ky!=0)
  ! the matrix is same, i.e. np-\phi (no matter k_y=0 or k_y!=0)
 
  ! the ky=0 source term
  do i=0,nxpp
     do j=0,jm
        do k=0,1
           u_zonal(i,j,k)=rho(i,j,k)+gn0e(i)*cn0e/gt0e(i)*phi(i,j,k)/real(ntube)
        enddo
     enddo
  enddo
  call dcmpy(u_zonal(0:imx-1,0:jmx-1,0:1),v_zonal)

  temp3dxy=0.
  do k=0,1
     do j=0,jcnt-1
        myj=jft(j)
        ! the ky!=0 source term
        sl(1:imx-1,j,k)=v(1:imx-1,j,k)
        ! the ky=0 source term
        if(myj==jft(0))sl(1:imx,j,k)=v_zonal(1:imx-1,j,k)
        call ZGETRS('N',imx-1,1,mx(:,:,j,k),imx-1,ipiv(:,:,j,k), &
             sl(:,j,k),imx-1,INFO)
        temp3dxy(1:imx-1,myj,k)=sl(1:imx-1,j,k)
        temp3dxy(0,myj,k)=0.
     enddo
  enddo

  !  from rho(kx,ky) to phi(kx,ky)
  do k=0,1
     do j=0,jmx-1
        do i=0,imx-1
           temp3dxy(i,j,k)=temp3dxy(i,j,k)/jmx !*formphi(i,j,k)
        enddo
     enddo
  enddo

  !  from phi(kx,ky) to phi(x,y)
  do k=0,mykm
     do i=0,imx-1
        do j=0,jmx-1
           tmpy(j)=temp3dxy(i,j,k)
        enddo
        call ccfft('y',1,jmx,1.0,tmpy,coefy,worky,0)
        do j=0,jmx-1
           temp3dxy(i,j,k)=tmpy(j)  ! phi(x,y)
        enddo
     enddo
  enddo

  do i=0,nxpp-1
     do j=0,jm-1
        do k=0,mykm
           phi(i,j,k)=temp3dxy(i,j,k)
        enddo
     enddo
  enddo

  !    x-y boundary points 
  call enfxy(phi(:,:,:))
  call enfz(phi(:,:,:))
  if(idg==1)write(*,*)'pass enfz', myid

end subroutine gkps_adiabatic_electron

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

