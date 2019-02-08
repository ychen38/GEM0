                i=int(xt/dx)
                j=int(yt/dy)
                k=0 !int(z3e(m)/dz)-gclr*kcnt

                wx0=float(i+1)-xt/dx 
                wx1=1.-wx0
                wy0=float(j+1)-yt/dy
                wy1=1.-wy0
                wz0=float(gclr*kcnt+k+1)-z2e(m)/dz
                wz1=1.-wz0
                x000=wx0*wy0*wz0
                x001=wx0*wy0*wz1
                x010=wx0*wy1*wz0
                x011=wx0*wy1*wz1
                x100=wx1*wy0*wz0
                x101=wx1*wy0*wz1
                x110=wx1*wy1*wz0
                x111=wx1*wy1*wz1
            exp1 = x000*ex(i,j,k) + x100*ex(i+1,j,k)  &
           + x010*ex(i,j+1,k) + x110*ex(i+1,j+1,k) + &
           x001*ex(i,j,k+1) + x101*ex(i+1,j,k+1) +   &
           x011*ex(i,j+1,k+1) + x111*ex(i+1,j+1,k+1)

            eyp = x000*ey(i,j,k) + x100*ey(i+1,j,k)  &
           + x010*ey(i,j+1,k) + x110*ey(i+1,j+1,k) + &
           x001*ey(i,j,k+1) + x101*ey(i+1,j,k+1) +   &
           x011*ey(i,j+1,k+1) + x111*ey(i+1,j+1,k+1)

            ezp = x000*ez(i,j,k) + x100*ez(i+1,j,k)  &
           + x010*ez(i,j+1,k) + x110*ez(i+1,j+1,k) + &
           x001*ez(i,j,k+1) + x101*ez(i+1,j,k+1) +   &
           x011*ez(i,j+1,k+1) + x111*ez(i+1,j+1,k+1)

            delbxp = x000*delbx(i,j,k) + x100*delbx(i+1,j,k)  &
           + x010*delbx(i,j+1,k) + x110*delbx(i+1,j+1,k) + &
           x001*delbx(i,j,k+1) + x101*delbx(i+1,j,k+1) +   &
           x011*delbx(i,j+1,k+1) + x111*delbx(i+1,j+1,k+1)

            delbyp = x000*delby(i,j,k) + x100*delby(i+1,j,k)  &
           + x010*delby(i,j+1,k) + x110*delby(i+1,j+1,k) + &
           x001*delby(i,j,k+1) + x101*delby(i+1,j,k+1) +   &
           x011*delby(i,j+1,k+1) + x111*delby(i+1,j+1,k+1)

            dgdtp = x000*dphidt(i,j,k) + x100*dphidt(i+1,j,k)  &
           + x010*dphidt(i,j+1,k) + x110*dphidt(i+1,j+1,k) + &
           x001*dphidt(i,j,k+1) + x101*dphidt(i+1,j,k+1) +   &
           x011*dphidt(i,j+1,k+1) + x111*dphidt(i+1,j+1,k+1)

            dpdzp = x000*dpdz(i,j,k) + x100*dpdz(i+1,j,k)  &
           + x010*dpdz(i,j+1,k) + x110*dpdz(i+1,j+1,k) + &
           x001*dpdz(i,j,k+1) + x101*dpdz(i+1,j,k+1) +   &
           x011*dpdz(i,j+1,k+1) + x111*dpdz(i+1,j+1,k+1)

            dadzp = x000*dadz(i,j,k) + x100*dadz(i+1,j,k)  &
           + x010*dadz(i,j+1,k) + x110*dadz(i+1,j+1,k) + &
           x001*dadz(i,j,k+1) + x101*dadz(i+1,j,k+1) +   &
           x011*dadz(i,j+1,k+1) + x111*dadz(i+1,j+1,k+1)

            aparp = x000*apar(i,j,k) + x100*apar(i+1,j,k)  &
           + x010*apar(i,j+1,k) + x110*apar(i+1,j+1,k) + &
           x001*apar(i,j,k+1) + x101*apar(i+1,j,k+1) +   &
           x011*apar(i,j+1,k+1) + x111*apar(i+1,j+1,k+1)

            phip = x000*phi(i,j,k) + x100*phi(i+1,j,k)  &
           + x010*phi(i,j+1,k) + x110*phi(i+1,j+1,k) + &
           x001*phi(i,j,k+1) + x101*phi(i+1,j,k+1) +   &
           x011*phi(i,j+1,k+1) + x111*phi(i+1,j+1,k+1)

