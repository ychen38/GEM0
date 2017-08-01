                i=int(xt/dx)
                j=int(yt/dy)
                k=0

                wx0=float(i+1)-xt/dx 
                wx1=1.-wx0
                wy0=float(j+1)-yt/dy
                wy1=1.-wy0
                wz0=float(gclr*kcnt+k+1)-zb3(m)/dz
                wz1=1.-wz0

            exp1=exp1 + wx0*wy0*wz0*ex(i,j,k) + wx1*wy0*wz0*ex(i+1,j,k) &
            + wx0*wy1*wz0*ex(i,j+1,k) + wx1*wy1*wz0*ex(i+1,j+1,k) + &
            wx0*wy0*wz1*ex(i,j,k+1) + wx1*wy0*wz1*ex(i+1,j,k+1) + &
            wx0*wy1*wz1*ex(i,j+1,k+1) + wx1*wy1*wz1*ex(i+1,j+1,k+1)

            eyp=eyp + wx0*wy0*wz0*ey(i,j,k) + wx1*wy0*wz0*ey(i+1,j,k) &
            + wx0*wy1*wz0*ey(i,j+1,k) + wx1*wy1*wz0*ey(i+1,j+1,k) + &
            wx0*wy0*wz1*ey(i,j,k+1) + wx1*wy0*wz1*ey(i+1,j,k+1) + &
            wx0*wy1*wz1*ey(i,j+1,k+1) + wx1*wy1*wz1*ey(i+1,j+1,k+1)

            ezp =ezp + wx0*wy0*wz0*ez(i,j,k) + wx1*wy0*wz0*ez(i+1,j,k) &
            + wx0*wy1*wz0*ez(i,j+1,k) + wx1*wy1*wz0*ez(i+1,j+1,k) + &
            wx0*wy0*wz1*ez(i,j,k+1) + wx1*wy0*wz1*ez(i+1,j,k+1) + &
            wx0*wy1*wz1*ez(i,j+1,k+1) + wx1*wy1*wz1*ez(i+1,j+1,k+1)

            delbxp =delbxp + wx0*wy0*wz0*delbx(i,j,k)  &
            + wx1*wy0*wz0*delbx(i+1,j,k) &
            + wx0*wy1*wz0*delbx(i,j+1,k) &
            + wx1*wy1*wz0*delbx(i+1,j+1,k) & 
            + wx0*wy0*wz1*delbx(i,j,k+1) &
            + wx1*wy0*wz1*delbx(i+1,j,k+1) & 
            + wx0*wy1*wz1*delbx(i,j+1,k+1) &
            + wx1*wy1*wz1*delbx(i+1,j+1,k+1)

            delbyp =delbyp + wx0*wy0*wz0*delby(i,j,k) &
            + wx1*wy0*wz0*delby(i+1,j,k) &
            + wx0*wy1*wz0*delby(i,j+1,k)  &
            + wx1*wy1*wz0*delby(i+1,j+1,k)  &
            + wx0*wy0*wz1*delby(i,j,k+1)  &
            + wx1*wy0*wz1*delby(i+1,j,k+1)  &
            + wx0*wy1*wz1*delby(i,j+1,k+1)  &
            + wx1*wy1*wz1*delby(i+1,j+1,k+1)


