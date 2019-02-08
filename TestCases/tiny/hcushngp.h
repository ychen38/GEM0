                i=int(xt/dx+0.5)
                j=int(yt/dy+0.5)
                k=int(zh3(m)/dz+0.5)-gclr*kcnt

            exp1=exp1 + ex(i,j,k)
            eyp=eyp + ey(i,j,k)
            ezp =ezp + ez(i,j,k)
            delbxp = delbxp+delbx(i,j,k)
            delbyp = delbyp+delby(i,j,k)
            dpdzp = dpdzp+dpdz(i,j,k)     
            dadzp = dadzp+dadz(i,j,k)     
            aparp = aparp+apar(i,j,k)
