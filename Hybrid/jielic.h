                i=int(xt/dx)
                j=int(yt/dy)
                k=0 !int(z3e(m)/dz)-gclr*kcnt
     
                wx0=float(i+1)-xt/dx 
                wx1=1.-wx0
                wy0=float(j+1)-yt/dy
                wy1=1.-wy0
                wz0=float(gclr*kcnt+k+1)-zc3(m)/dz
                wz1=1.-wz0

                myjpex(i,j,k)      =myjpex(i,j,k)   &
     				+wght*xdot*wx0*wy0*wz0
                myjpex(i+1,j,k)    =myjpex(i+1,j,k) &
     				+wght*xdot*wx1*wy0*wz0
                myjpex(i,j+1,k)    =myjpex(i,j+1,k) &
     				+wght*xdot*wx0*wy1*wz0
                myjpex(i+1,j+1,k)  =myjpex(i+1,j+1,k) &
     				+wght*xdot*wx1*wy1*wz0
                myjpex(i,j,k+1)    =myjpex(i,j,k+1)  &
     				+wght*xdot*wx0*wy0*wz1
                myjpex(i+1,j,k+1)  =myjpex(i+1,j,k+1) &
     				+wght*xdot*wx1*wy0*wz1
                myjpex(i,j+1,k+1)  =myjpex(i,j+1,k+1) &
     				+wght*xdot*wx0*wy1*wz1
                myjpex(i+1,j+1,k+1)=myjpex(i+1,j+1,k+1) &
     				+wght*xdot*wx1*wy1*wz1

                myjpey(i,j,k)      =myjpey(i,j,k)   &
     				+wght*ydot*wx0*wy0*wz0
                myjpey(i+1,j,k)    =myjpey(i+1,j,k) &
     				+wght*ydot*wx1*wy0*wz0
                myjpey(i,j+1,k)    =myjpey(i,j+1,k) &
     				+wght*ydot*wx0*wy1*wz0
                myjpey(i+1,j+1,k)  =myjpey(i+1,j+1,k) &
     				+wght*ydot*wx1*wy1*wz0
                myjpey(i,j,k+1)    =myjpey(i,j,k+1)  &
     				+wght*ydot*wx0*wy0*wz1
                myjpey(i+1,j,k+1)  =myjpey(i+1,j,k+1) &
     				+wght*ydot*wx1*wy0*wz1
                myjpey(i,j+1,k+1)  =myjpey(i,j+1,k+1) &
     				+wght*ydot*wx0*wy1*wz1
                myjpey(i+1,j+1,k+1)=myjpey (i+1,j+1,k+1) &
     				+wght*ydot*wx1*wy1*wz1


                myjpar(i,j,k)      =myjpar(i,j,k) &
     				+wght*zdot*wx0*wy0*wz0
                myjpar(i+1,j,k)    =myjpar(i+1,j,k) &
     				+wght*zdot*wx1*wy0*wz0
                myjpar(i,j+1,k)    =myjpar(i,j+1,k)  &
     				+wght*zdot*wx0*wy1*wz0
                myjpar(i+1,j+1,k)  =myjpar(i+1,j+1,k) &
     				+wght*zdot*wx1*wy1*wz0
                myjpar(i,j,k+1)    =myjpar(i,j,k+1)  &
     				+wght*zdot*wx0*wy0*wz1
                myjpar(i+1,j,k+1)  =myjpar(i+1,j,k+1) &
     				+wght*zdot*wx1*wy0*wz1
                myjpar(i,j+1,k+1)  =myjpar(i,j+1,k+1) &
     				+wght*zdot*wx0*wy1*wz1
                myjpar(i+1,j+1,k+1)=myjpar(i+1,j+1,k+1) &
     				+wght*zdot*wx1*wy1*wz1

                mydnidt(i,j,k)      =mydnidt(i,j,k)   &
     				+wght1*wx0*wy0*wz0
                mydnidt(i+1,j,k)    =mydnidt(i+1,j,k) &
     				+wght1*wx1*wy0*wz0
                mydnidt(i,j+1,k)    =mydnidt(i,j+1,k) &
     				+wght1*wx0*wy1*wz0
                mydnidt(i+1,j+1,k)  =mydnidt(i+1,j+1,k) &
     				+wght1*wx1*wy1*wz0
                mydnidt(i,j,k+1)    =mydnidt(i,j,k+1)  &
     				+wght1*wx0*wy0*wz1
                mydnidt(i+1,j,k+1)  =mydnidt(i+1,j,k+1) &
     				+wght1*wx1*wy0*wz1
                mydnidt(i,j+1,k+1)  =mydnidt(i,j+1,k+1) &
     				+wght1*wx0*wy1*wz1
                mydnidt(i+1,j+1,k+1)=mydnidt(i+1,j+1,k+1) &
     				+wght1*wx1*wy1*wz1
