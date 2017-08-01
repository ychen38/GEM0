                i=int(xt/dx)
                j=int(yt/dy)
                k=0 !int(z3e(m)/dz)-gclr*kcnt
     
                wx0=float(i+1)-xt/dx 
                wx1=1.-wx0
                wy0=float(j+1)-yt/dy
                wy1=1.-wy0
                wz0=float(gclr*kcnt+k+1)-zh3(m)/dz
                wz1=1.-wz0
     
                myhden(i,j,k)      =myhden(i,j,k)   &
     				+wght*wx0*wy0*wz0
                myhden(i+1,j,k)    =myhden(i+1,j,k) &
     				+wght*wx1*wy0*wz0
                myhden(i,j+1,k)    =myhden(i,j+1,k) &
     				+wght*wx0*wy1*wz0
                myhden(i+1,j+1,k)  =myhden(i+1,j+1,k) &
     				+wght*wx1*wy1*wz0
                myhden(i,j,k+1)    =myhden(i,j,k+1)  &
     				+wght*wx0*wy0*wz1
                myhden(i+1,j,k+1)  =myhden(i+1,j,k+1) &
     				+wght*wx1*wy0*wz1
                myhden(i,j+1,k+1)  =myhden(i,j+1,k+1) &
     				+wght*wx0*wy1*wz1
                myhden(i+1,j+1,k+1)=myhden(i+1,j+1,k+1) &
     				+wght*wx1*wy1*wz1


                myhpar(i,j,k)      =myhpar(i,j,k) &
     				+wght*wx0*wy0*wz0*vpar
                myhpar(i+1,j,k)    =myhpar(i+1,j,k) &
     				+wght*wx1*wy0*wz0*vpar
                myhpar(i,j+1,k)    =myhpar(i,j+1,k)  &
     				+wght*wx0*wy1*wz0*vpar
                myhpar(i+1,j+1,k)  =myhpar(i+1,j+1,k) &
     				+wght*wx1*wy1*wz0*vpar
                myhpar(i,j,k+1)    =myhpar(i,j,k+1)  &
     				+wght*wx0*wy0*wz1*vpar
                myhpar(i+1,j,k+1)  =myhpar(i+1,j,k+1) &
     				+wght*wx1*wy0*wz1*vpar
                myhpar(i,j+1,k+1)  =myhpar(i,j+1,k+1) &
     				+wght*wx0*wy1*wz1*vpar
                myhpar(i+1,j+1,k+1)=myhpar(i+1,j+1,k+1) &
     				+wght*wx1*wy1*wz1*vpar

                myhden0(i,j,k)      =myhden0(i,j,k)   &
     				+wght0*wx0*wy0*wz0
                myhden0(i+1,j,k)    =myhden0(i+1,j,k) &
     				+wght0*wx1*wy0*wz0
                myhden0(i,j+1,k)    =myhden0(i,j+1,k) &
     				+wght0*wx0*wy1*wz0
                myhden0(i+1,j+1,k)  =myhden0(i+1,j+1,k) &
     				+wght0*wx1*wy1*wz0
                myhden0(i,j,k+1)    =myhden0(i,j,k+1)  &
     				+wght0*wx0*wy0*wz1
                myhden0(i+1,j,k+1)  =myhden0(i+1,j,k+1) &
     				+wght0*wx1*wy0*wz1
                myhden0(i,j+1,k+1)  =myhden0(i,j+1,k+1) &
     				+wght0*wx0*wy1*wz1
                myhden0(i+1,j+1,k+1)=myhden0(i+1,j+1,k+1) &
     				+wght0*wx1*wy1*wz1
