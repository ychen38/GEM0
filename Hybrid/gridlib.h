                i=int(xt/dx)
                j=int(yt/dy)
                k=0 !int(z3e(m)/dz)-gclr*kcnt
     
                wx0=float(i+1)-xt/dx 
                wx1=1.-wx0
                wy0=float(j+1)-yt/dy
                wy1=1.-wy0
                wz0=float(gclr*kcnt+k+1)-zb3(m)/dz
                wz1=1.-wz0
     
                mybden(i,j,k)      =mybden(i,j,k)   &
     				+wght*wx0*wy0*wz0
                mybden(i+1,j,k)    =mybden(i+1,j,k) &
     				+wght*wx1*wy0*wz0
                mybden(i,j+1,k)    =mybden(i,j+1,k) &
     				+wght*wx0*wy1*wz0
                mybden(i+1,j+1,k)  =mybden(i+1,j+1,k) &
     				+wght*wx1*wy1*wz0
                mybden(i,j,k+1)    =mybden(i,j,k+1)  &
     				+wght*wx0*wy0*wz1
                mybden(i+1,j,k+1)  =mybden(i+1,j,k+1) &
     				+wght*wx1*wy0*wz1
                mybden(i,j+1,k+1)  =mybden(i,j+1,k+1) &
     				+wght*wx0*wy1*wz1
                mybden(i+1,j+1,k+1)=mybden(i+1,j+1,k+1) &
     				+wght*wx1*wy1*wz1


                mybpar(i,j,k)      =mybpar(i,j,k) &
     				+wght*wx0*wy0*wz0*vpar
                mybpar(i+1,j,k)    =mybpar(i+1,j,k) &
     				+wght*wx1*wy0*wz0*vpar
                mybpar(i,j+1,k)    =mybpar(i,j+1,k)  &
     				+wght*wx0*wy1*wz0*vpar
                mybpar(i+1,j+1,k)  =mybpar(i+1,j+1,k) &
     				+wght*wx1*wy1*wz0*vpar
                mybpar(i,j,k+1)    =mybpar(i,j,k+1)  &
     				+wght*wx0*wy0*wz1*vpar
                mybpar(i+1,j,k+1)  =mybpar(i+1,j,k+1) &
     				+wght*wx1*wy0*wz1*vpar
                mybpar(i,j+1,k+1)  =mybpar(i,j+1,k+1) &
     				+wght*wx0*wy1*wz1*vpar
                mybpar(i+1,j+1,k+1)=mybpar(i+1,j+1,k+1) &
     				+wght*wx1*wy1*wz1*vpar

