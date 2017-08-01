program fft

parameter(imx=32,jmx=128)
COMPLEX(8) :: x(0:imx-1,0:jmx-1),y(0:jmx-1,0:imx-1)
COMPLEX(8) :: sum,b(0:imx-1,0:jmx-1),IU=cmplx(0.,1.)
  real(8) :: cfxp(20000),cfyp(20000)
  real(8) :: cfxn(20000),cfyn(20000)
  real(8) :: wrkx(20000),wrky(20000)

 call dcft(1,x,1,0,y,1,0,imx,1,+1, 1.d0,cfxp,20000,wrkx,20000)
 call dcft(1,x,1,0,y,1,0,jmx,1,+1, 1.d0,cfyp,20000,wrky,20000)

  do i = 0,imx-1
      do j = 0,jmx-1
         x(i,j) = cmplx(exp(i/(1.0+imx)*cos(j+0.6)),sin(2.*j))
      end do
  end do

  do i = 0,imx-1
      do j = 0,jmx-1
         write(*,*)i,j,x(i,j)
      end do
  end do
  stop

  pi = atan(1.0)*4    
  do i = 0,imx-1
      do j = 0,jmx-1
         sum = 0.
         do l = 0,imx-1
            do m = 0,jmx-1
               sum = sum+x(l,m)*exp(-IU*pi*2.0*(i*l/float(imx)+j*m/float(jmx)))
            end do
         end do
         b(i,j) = sum
      end do
  end do
      
  do j = 0,jmx-1    
  call dcft(0,x(0,j),1,0,x(0,j),1,0,imx,1,+1, 1.d0,cfxp,20000,wrkx,20000)      
  end do
  y=transpose(x)
  do i = 0,imx-1
  call dcft(0,y(0,i),1,0,y(0,i),1,0,jmx,1,+1, 1.d0,cfyp,20000,wrky,20000)
  end do

  x=transpose(y)
  do i = 0,imx-1
      do j = 0,jmx-1
         write(*,*)i,j,x(i,j),b(i,j)
      end do
  end do
stop
end
