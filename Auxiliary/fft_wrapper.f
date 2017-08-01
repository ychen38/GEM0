module fft_wrapper
use gem_com
implicit none
complex(8), dimension(:,:), allocatable :: exfac,eyfac,ezfac
real(8), dimension(:,:), allocatable :: sxfac

contains

subroutine fft_init
   allocate (exfac(0:imx-1,0:imx-1),eyfac(0:jmx-1,0:jmx-1),ezfac(0:kmx-1,0:kmx-1))
   allocate(sxfac(0:imx-1,0:imx-1))               
end subroutine fft_init

subroutine ccfft(c,isign,n,scale,x,table,work,isys)
   character :: c
   real(8) :: scale
   integer :: isign,n,isys      
   integer :: i,j,k
   complex(8),dimension(0:) :: x
   real(8),dimension(:) :: table
   real(8),dimension(:) :: work
   complex(8) :: cdum,cxarr(0:imx-1),cyarr(0:jmx-1),czarr(0:kmx-1)   



!initialize
   if(isign==0)then
      if(c=='x')then
         do i = 0,imx-1
            do j = 0,imx-1
               exfac(i,j) = exp(IU*pi2*i*j/imx)
            end do
         end do
      end if
      if(c=='y')then
         do i = 0,jmx-1
            do j = 0,jmx-1
               eyfac(i,j) = exp(IU*pi2*i*j/jmx)
            end do
         end do
      end if
      if(c=='z')then
         do i = 0,kmx-1
            do j = 0,kmx-1
               ezfac(i,j) = exp(IU*pi2*i*j/kmx)
            end do
         end do
      end if
   end if

   if(isign==1)then
      if(c=='x')then
         do i = 0,imx-1
            cdum = 0.
            do j = 0,imx-1
               cdum = cdum+x(j)*exfac(i,j)
            end do
            cxarr(i) = cdum
         end do
         x(0:imx-1) = cxarr(0:imx-1)
! call zfftb(n,x,coefxp)
      end if
      if(c=='y')then
         do i = 0,jmx-1
            cdum = 0.
            do j = 0,jmx-1
               cdum = cdum+x(j)*eyfac(i,j)
            end do
            cyarr(i) = cdum
         end do
         x(0:jmx-1) = cyarr(0:jmx-1)
! call zfftb(n,x,coefyp)
      end if
      if(c=='z')then
         do i = 0,kmx-1
            cdum = 0.
            do j = 0,kmx-1
               cdum = cdum+x(j)*ezfac(i,j)
            end do
            czarr(i) = cdum
         end do
         x(0:kmx-1) = czarr(0:kmx-1)
! call zfftb(n,x,coefzp)
      end if
   end if

   if(isign==-1)then
      if(c=='x')then
         do i = 0,imx-1
            cdum = 0.
            do j = 0,imx-1
               cdum = cdum+x(j)/exfac(i,j)
            end do
            cxarr(i) = cdum
         end do
         x(0:imx-1) = cxarr(0:imx-1)
! call zfftb(n,x,coefxp)
      end if
      if(c=='y')then
         do i = 0,jmx-1
            cdum = 0.
            do j = 0,jmx-1
               cdum = cdum+x(j)/eyfac(i,j)
            end do
            cyarr(i) = cdum
         end do
         x(0:jmx-1) = cyarr(0:jmx-1)
! call zfftb(n,x,coefyp)
      end if
      if(c=='z')then
         do i = 0,kmx-1
            cdum = 0.
            do j = 0,kmx-1
               cdum = cdum+x(j)/ezfac(i,j)
            end do
            czarr(i) = cdum
         end do
         x(0:kmx-1) = czarr(0:kmx-1)
! call zfftb(n,x,coefzp)
      end if
   end if

   return
end subroutine ccfft

subroutine dsinf(init,x,inc1x,inc2x,inc1y,inc2y,n,m,scale,aux1,naux1,aux2,naux2)
   integer :: init,inc1x,inc2x,inc1y,inc2y,n,m,naux1,naux2,n2,i,j
   real(8) :: scale,dum
   real(8),dimension(:) :: x,aux1,aux2
   real(8) :: dxarr(0:imx-1)

!   n2=n/2-1
!initialize
   if(init/=0)then
      do i = 0,imx-1
         do j = 0,imx-1
            sxfac(i,j) = dsin(pi*i*j/imx)
         end do
      end do
! call dsinti(n2,wsave)
   end if

   if(init==0)then
      do i = 0,imx-1
         dum = 0.
         do j = 0,imx-1
            dum = dum+x(j)*sxfac(i,j)
         end do
         dxarr(i) = dum
      end do
      x(0:imx-1) = dxarr(0:imx-1)
   end if

   return
end subroutine dsinf

end module fft_wrapper

