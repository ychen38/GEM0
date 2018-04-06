      real function erf(y)
      implicit real(a-h,o-z)
      real y
      real p(5),q(3),p1(8),q1(7),p2(5),q2(4)
!
!                                  coefficients for 0.0 .le. y .lt.
!                                  .477
!
                     p(1)=-.44422647396874d0
                        p(2)=10.731707253648d0
                        p(3)=15.915606197771d0
                        p(4)=374.81624081284d0
                        p(5)=2.5612422994823d-02
                     q(1)=17.903143558843d0
                        q(2)=124.82892031581d0
                        q(3)=332.17224470532d0
!
!                                  coefficients for .477 .le. y .le.
!                                  4.0
!
                     p1(1)=7.2117582508831d0
                        p1(2)=43.162227222057d0
                        p1(3)=152.98928504694d0
                        p1(4)=339.32081673434d0
                        p1(5)=451.91895371187d0
                        p1(6)=300.45926102016d0
                        p1(7)=-1.3686485738272d-07
                        p1(8)=.56419551747897d0
                     q1(1)=77.000152935229d0
                        q1(2)=277.58544474399d0
                        q1(3)=638.98026446563d0
                        q1(4)=931.35409485061d0
                        q1(5)=790.95092532790d0
                        q1(6)=300.45926095698d0
                        q1(7)=12.782727319629d0
!
!                                 coefficients for 4.0 .lt. y
!
                     p2(1)=-.22695659353969d0
                        p2(2)=-4.9473091062325d-02
                        p2(3)=-2.9961070770354d-03
                        p2(4)=-2.2319245973418d-02
                        p2(5)=-2.7866130860965d-01
                     q2(1)=1.0516751070679d0
                        q2(2)=.19130892610783d0
                        q2(3)=1.0620923052847d-02
                        q2(4)=1.9873320181714d0
!
!                                 constants
!
                     xmin=1.0d-8
                     xlarge=5.6875d0
                     xbig=25.90625d0
                     ssq1h=.70710678118655d0
                     ssqpi=.56418958354776d0
!
      kret = 0
    5 x = y
   10 isw = 1
      if (x.ge.0.0d0) go to 15
      isw = -1
      x = -x
   15 if (x.lt..477d0) go to 20
      if (x.le.4.0d0) go to 40
      if (kret*isw.gt.0) go to 50
      if (x.lt.xlarge) go to 55
      res = 1.d0
      if (kret.eq.0) go to 75
      res = 2.d0
      go to 80
!                                  abs(y) .lt. .477, evaluate
!                                  approximation for erf
   20 if (x.lt.xmin) go to 30
      xsq = x*x
      xnum = p(5)
      do 25 i=1,4
         xnum = xnum*xsq+p(i)
   25 continue
      xden = ((q(1)+xsq)*xsq+q(2))*xsq+q(3)
      res = x*xnum/xden
      go to 35
   30 res = x*p(4)/q(3)
   35 if (isw.eq.-1) res = -res
      if (kret.eq.0) go to 90
      res = 1.0d0-res
      go to 80
!                                  .477 .le. abs(y) .le. 4.0
!                                  evaluate approximation for erfc
   40 xsq = x*x
      xnum = p1(7)*x+p1(8)
      xden = x+q1(7)
      do 45 i=1,6
         xnum = xnum*x+p1(i)
         xden = xden*x+q1(i)
   45 continue
      res = xnum/xden
      go to 65
!                                  4.0 .lt. abs(y), evaluate
!                                  approximation for erfc
   50 if (x.gt.xbig) go to 85
   55 xsq = x*x
      xi = 1.0d0/xsq
      xnum = p2(4)*xi+p2(5)
      xden = xi+q2(4)
      do 60 i=1,3
         xnum = xnum*xi+p2(i)
         xden = xden*xi+q2(i)
   60 continue
      res = (ssqpi+xi*xnum/xden)/x
   65 res = res*exp(-xsq)
      if (kret.eq.0) go to 70
      if (isw.eq.-1) res = 2.0d0-res
      go to 80
   70 res = 1.0d0-res
   75 if (isw.eq.-1) res = -res
      go to 90
!                                  prepare to return
   80 if (kret.eq.2) res = res*0.5d0
      go to 90
   85 res = 0.0d0
   90 erf = res
      return
      end
