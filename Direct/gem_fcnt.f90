!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!   Functions, etc.
!   (things that don't change)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!v
real function revers(num,n)
  integer :: num,n,inum,iquot,irem
  real :: rev,power
  !    function to reverse the digits of num in base n.

  rev = 0.d0
  inum = num
  power = 1.d0

11 continue
  iquot = int(inum/n)
  irem = inum - n*iquot
  power = power/n
  rev = rev + irem*power
  inum = iquot
  if(inum.gt.0) goto 11

  revers = rev
  return
end function revers

!-----------------------------------------------------------------------

subroutine srcbes(biz,gam0,gam1)

  REAL :: t1,t2,biz,gam0,gam1 

  !.....Calculates gamma nought and gamma 1. (Abramowitz and Stegun).
446 if (biz.gt.3.75d0) go to 148
  t1=(biz/3.75d0)**2
  t2=exp(-biz)
  gam0=t2*((((((.0045813d0*t1+.0360768d0)*t1+.2659732d0)*t1+ &
       1.2067492d0)*t1+3.0899424d0)*t1+3.5156229d0)*t1+1.d0)
  gam1=t2*biz*((((((.00032411d0*t1+.00301532d0)*t1+.02658733d0) &
       *t1+.15084934d0)*t1+.51498869d0)*t1+.87890594d0)*t1+.5d0)
  go to 149
148 t2=1.d0/sqrt(biz)
  t1=3.75d0/biz
  gam0=t2*((((((((.00392377d0*t1-.01647633d0)*t1+.02635537d0) &
       *t1-.02057706d0)*t1+.00916281d0)*t1-.00157565d0)*t1+ &
       .00225319d0)*t1+.01328592d0)*t1+.39894228d0)
  gam1=t2*((((((((-.00420059d0*t1+.01787654d0)*t1-.02895312d0) &
       *t1+.02282967d0)*t1-.01031555d0)*t1+.00163801d0)*t1- &
       .00362018d0)*t1-.03988024d0)*t1+.39894228d0)
149 continue
  return
end subroutine srcbes

!-----------------------------------------------------------------------
