                i=int(xt/dx+0.5)
                j=int(yt/dy+0.5)
                k=int(z3(m)/dz+0.5)-gclr*kcnt
     
                myden(i,j,k) = myden(i,j,k) + wght
		myjpar(i,j,k) = myjpar(i,j,k)+wght*vpar
                mydti(i,j,k) = mydti(i,j,k)+wght*vfac
