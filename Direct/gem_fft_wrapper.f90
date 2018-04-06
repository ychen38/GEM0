module gem_fft_wrapper

  implicit none

  real :: coefxp(20000),coefyp(20000),coefzp(20000)
  real :: coefxn(20000),coefyn(20000),coefzn(20000)
  real :: workxx(20000),workyy(20000),workzz(20000)
  real :: wsave(50000)

contains

  subroutine ccfft(c,isign,n,scale,x,table,work,isys)
    character :: c
    real :: scale
    integer :: isign,n,isys      
    complex,dimension(0:) :: x
    real,dimension(:) :: table
    real,dimension(:) :: work

    !initialize
    if(isign==0)then
       if(c=='x')then
          ! call dcft(1,x,1,0,y,1,0,n,1,+1, 1.d0,coefxn,20000,workxx,20000)
          call zffti(n,coefxn)
          ! call dcft(1,x,1,0,y,1,0,n,1,-1, 1.d0,coefxp,20000,workxx,20000)
          call zffti(n,coefxp)
       end if
       if(c=='y')then
          ! call dcft(1,x,1,0,y,1,0,n,1,+1, 1.d0,coefyn,20000,workyy,20000)
          call zffti(n,coefyn)
          ! call dcft(1,x,1,0,y,1,0,n,1,-1, 1.d0,coefyp,20000,workyy,20000)
          call zffti(n,coefyp)
       end if
       if(c=='z')then
          ! call dcft(1,x,1,0,y,1,0,n,1,+1, 1.d0,coefzn,20000,workzz,20000)
          call zffti(n,coefzn)
          ! call dcft(1,x,1,0,y,1,0,n,1,-1, 1.d0,coefzp,20000,workzz,20000)
          call zffti(n,coefzp)
       end if
    end if

    if(isign==1)then
       if(c=='x')then
          ! call dcft(0,x,1,0,y,1,0,n,1,-isign, 1.d0,coefxp,20000,workxx,20000)
          call zfftb(n,x,coefxp)
       end if
       if(c=='y')then
          ! call dcft(0,x,1,0,y,1,0,n,1,-isign, 1.d0,coefyp,20000,workyy,20000)
          call zfftb(n,x,coefyp)
       end if
       if(c=='z')then
          ! call dcft(0,x,1,0,y,1,0,n,1,-isign, 1.d0,coefzp,20000,workzz,20000)
          call zfftb(n,x,coefzp)
       end if
    end if

    if(isign==-1)then
       if(c=='x')then
          ! call dcft(0,x,1,0,y,1,0,n,1,-isign, 1.d0,coefxn,20000,workxx,20000)
          call zfftf(n,x,coefxn)
       end if
       if(c=='y')then
          ! call dcft(0,x,1,0,y,1,0,n,1,-isign, 1.d0,coefyn,20000,workyy,20000)
          call zfftf(n,x,coefyn)
       end if
       if(c=='z')then
          ! call dcft(0,x,1,0,y,1,0,n,1,-isign, 1.d0,coefzn,20000,workzz,20000)
          call zfftf(n,x,coefzn)
       end if
    end if

    return
  end subroutine ccfft

  subroutine dsinf(init,x,inc1x,inc2x,inc1y,inc2y,n,m,scale,aux1,naux1,aux2,naux2)
    integer :: init,inc1x,inc2x,inc1y,inc2y,n,m,naux1,naux2,n2,i
    real :: scale
    real,dimension(:) :: x,aux1,aux2

    n2=n/2-1
    !initialize
    if(init/=0)then
       call dsinti(n2,wsave)
    end if

    if(init==0)then
       call dsint(n2,x(2),wsave)
       do i=2,n2+1
          x(i)=0.5d0*x(i)
       end do
       x(1)=0.d0
    end if

    return
  end subroutine dsinf

end module gem_fft_wrapper

