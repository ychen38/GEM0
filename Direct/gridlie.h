                i=int(xt/dx)
                j=int(yt/dy)
                k=0 !int(z3e(m)/dz)-gclr*kcnt
     
                wx0=float(i+1)-xt/dx 
                wx1=1.-wx0
                wy0=float(j+1)-yt/dy
                wy1=1.-wy0
                wz0=float(gclr*kcnt+k+1)-z3e(m)/dz
                wz1=1.-wz0
     
                mydene(i,j,k)      =mydene(i,j,k)+wght*w000(m)
		mydene(i+1,j,k)    =mydene(i+1,j,k)+wght*w100(m)
		mydene(i,j+1,k)    =mydene(i,j+1,k)+wght*w010(m)
		mydene(i+1,j+1,k)  =mydene(i+1,j+1,k)+wght*w110(m)
		mydene(i,j,k+1)    =mydene(i,j,k+1)+wght*w001(m)
		mydene(i+1,j,k+1)  =mydene(i+1,j,k+1)+wght*w101(m)
		mydene(i,j+1,k+1)  =mydene(i,j+1,k+1)+wght*w011(m)
		mydene(i+1,j+1,k+1)=mydene(i+1,j+1,k+1)+wght*w111(m)

                myupar(i,j,k)      =myupar(i,j,k)+wght*vpar*w000(m)
		myupar(i+1,j,k)    =myupar(i+1,j,k)+wght*vpar*w100(m)
		myupar(i,j+1,k)    =myupar(i,j+1,k)+wght*vpar*w010(m)
		myupar(i+1,j+1,k)  =myupar(i+1,j+1,k)+wght*vpar*w110(m)
		myupar(i,j,k+1)    =myupar(i,j,k+1)+wght*vpar*w001(m)
		myupar(i+1,j,k+1)  =myupar(i+1,j,k+1)+wght*vpar*w101(m)
		myupar(i,j+1,k+1)  =myupar(i,j+1,k+1)+wght*vpar*w011(m)
		myupar(i+1,j+1,k+1)=myupar(i+1,j+1,k+1)+wght*vpar*w111(m)

