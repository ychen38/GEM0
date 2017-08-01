                i=int(xt/dx)
                j=int(yt/dy)
                k=0 !int(z3e(m)/dz)-gclr*kcnt
     
                wx0=float(i+1)-xt/dx 
                wx1=1.-wx0
                wy0=float(j+1)-yt/dy
                wy1=1.-wy0
                wz0=float(gclr*kcnt+k+1)-zc3(m)/dz
                wz1=1.-wz0
     
                myden(i,j,k)      =myden(i,j,k)   &
     				+wght*wx0*wy0*wz0
                myden(i+1,j,k)    =myden(i+1,j,k) &
     				+wght*wx1*wy0*wz0
                myden(i,j+1,k)    =myden(i,j+1,k) &
     				+wght*wx0*wy1*wz0
                myden(i+1,j+1,k)  =myden(i+1,j+1,k) &
     				+wght*wx1*wy1*wz0
                myden(i,j,k+1)    =myden(i,j,k+1)  &
     				+wght*wx0*wy0*wz1
                myden(i+1,j,k+1)  =myden(i+1,j,k+1) &
     				+wght*wx1*wy0*wz1
                myden(i,j+1,k+1)  =myden(i,j+1,k+1) &
     				+wght*wx0*wy1*wz1
                myden(i+1,j+1,k+1)=myden(i+1,j+1,k+1) &
     				+wght*wx1*wy1*wz1


                myjpar(i,j,k)      =myjpar(i,j,k) &
     				+wght*wx0*wy0*wz0*vpar
                myjpar(i+1,j,k)    =myjpar(i+1,j,k) &
     				+wght*wx1*wy0*wz0*vpar
                myjpar(i,j+1,k)    =myjpar(i,j+1,k)  &
     				+wght*wx0*wy1*wz0*vpar
                myjpar(i+1,j+1,k)  =myjpar(i+1,j+1,k) &
     				+wght*wx1*wy1*wz0*vpar
                myjpar(i,j,k+1)    =myjpar(i,j,k+1)  &
     				+wght*wx0*wy0*wz1*vpar
                myjpar(i+1,j,k+1)  =myjpar(i+1,j,k+1) &
     				+wght*wx1*wy0*wz1*vpar
                myjpar(i,j+1,k+1)  =myjpar(i,j+1,k+1) &
     				+wght*wx0*wy1*wz1*vpar
                myjpar(i+1,j+1,k+1)=myjpar(i+1,j+1,k+1) &
     				+wght*wx1*wy1*wz1*vpar

