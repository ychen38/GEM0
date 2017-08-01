                i=int(xt/dx)
                j=int(yt/dy)
                k=0 !int(z3e(m)/dz)-gclr*kcnt
     
                wx0=float(i+1)-xt/dx 
                wx1=1.-wx0
                wy0=float(j+1)-yt/dy
                wy1=1.-wy0
                wz0=float(gclr*kcnt+k+1)-z3e(m)/dz
                wz1=1.-wz0
     
                mydene(i,j,k)      =mydene(i,j,k)   &
     				+wght2*wx0*wy0*wz0
                mydene(i+1,j,k)    =mydene(i+1,j,k) &
     				+wght2*wx1*wy0*wz0
                mydene(i,j+1,k)    =mydene(i,j+1,k) &
     				+wght2*wx0*wy1*wz0
                mydene(i+1,j+1,k)  =mydene(i+1,j+1,k) &
     				+wght2*wx1*wy1*wz0
                mydene(i,j,k+1)    =mydene(i,j,k+1)  &
     				+wght2*wx0*wy0*wz1
                mydene(i+1,j,k+1)  =mydene(i+1,j,k+1) &
     				+wght2*wx1*wy0*wz1
                mydene(i,j+1,k+1)  =mydene(i,j+1,k+1) &
     				+wght2*wx0*wy1*wz1
                mydene(i+1,j+1,k+1)=mydene(i+1,j+1,k+1) &
     				+wght2*wx1*wy1*wz1


                myupar(i,j,k)      =myupar(i,j,k) &
     				+wght3*wx0*wy0*wz0
                myupar(i+1,j,k)    =myupar(i+1,j,k) &
     				+wght3*wx1*wy0*wz0
                myupar(i,j+1,k)    =myupar(i,j+1,k)  &
     				+wght3*wx0*wy1*wz0
                myupar(i+1,j+1,k)  =myupar(i+1,j+1,k) &
     				+wght3*wx1*wy1*wz0
                myupar(i,j,k+1)    =myupar(i,j,k+1)  &
     				+wght3*wx0*wy0*wz1
                myupar(i+1,j,k+1)  =myupar(i+1,j,k+1) &
     				+wght3*wx1*wy0*wz1
                myupar(i,j+1,k+1)  =myupar(i,j+1,k+1) &
     				+wght3*wx0*wy1*wz1
                myupar(i+1,j+1,k+1)=myupar(i+1,j+1,k+1) &
     				+wght3*wx1*wy1*wz1

                mydltpe(i,j,k)      =mydltpe(i,j,k)   &
     				+wght*wx0*wy0*wz0
                mydltpe(i+1,j,k)    =mydltpe(i+1,j,k) &
     				+wght*wx1*wy0*wz0
                mydltpe(i,j+1,k)    =mydltpe(i,j+1,k) &
     				+wght*wx0*wy1*wz0
                mydltpe(i+1,j+1,k)  =mydltpe(i+1,j+1,k) &
     				+wght*wx1*wy1*wz0
                mydltpe(i,j,k+1)    =mydltpe(i,j,k+1)  &
     				+wght*wx0*wy0*wz1
                mydltpe(i+1,j,k+1)  =mydltpe(i+1,j,k+1) &
     				+wght*wx1*wy0*wz1
                mydltpe(i,j+1,k+1)  =mydltpe(i,j+1,k+1) &
     				+wght*wx0*wy1*wz1
                mydltpe(i+1,j+1,k+1)=mydltpe(i+1,j+1,k+1) &
     				+wght*wx1*wy1*wz1


                mydltpa(i,j,k)      =mydltpa(i,j,k) &
     				+wght1*wx0*wy0*wz0
                mydltpa(i+1,j,k)    =mydltpa(i+1,j,k) &
     				+wght1*wx1*wy0*wz0
                mydltpa(i,j+1,k)    =mydltpa(i,j+1,k)  &
     				+wght1*wx0*wy1*wz0
                mydltpa(i+1,j+1,k)  =mydltpa(i+1,j+1,k) &
     				+wght1*wx1*wy1*wz0
                mydltpa(i,j,k+1)    =mydltpa(i,j,k+1)  &
     				+wght1*wx0*wy0*wz1
                mydltpa(i+1,j,k+1)  =mydltpa(i+1,j,k+1) &
     				+wght1*wx1*wy0*wz1
                mydltpa(i,j+1,k+1)  =mydltpa(i,j+1,k+1) &
     				+wght1*wx0*wy1*wz1
                mydltpa(i+1,j+1,k+1)=mydltpa(i+1,j+1,k+1) &
     				+wght1*wx1*wy1*wz1

