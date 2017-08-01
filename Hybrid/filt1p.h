	 do n=1,NNI

	 htemp3d=ttemp3d

!  build buffers and send grid info to neighbors

	 do j=0,jm-1
	   do i=0,im-1
	     rbfs(i,j)=htemp3d(mykm-1,i,j)
	     lbfs(i,j)=htemp3d(0,i,j)
	   enddo
	 enddo

	 call MPI_SENDRECV(rbfs(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_REAL8,rngbr,9, &
     	                   lbfr(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_REAL8,lngbr,9,  &
     		           TUBE_COMM,stat,ierr)

	 call MPI_SENDRECV(lbfs(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_REAL8,lngbr,8,  &
     	                   rbfr(0:imx-1,0:jmx-1),imx*jmx,  &
     		MPI_REAL8,rngbr,8,  &
     		           TUBE_COMM,stat,ierr)
!	 write(*,*)'past buff'

	 do j=0,jm-1
	   do i=0,im-1
	     ttemp3d(0,i,j)=( htemp3d(0,i,j)   &
     	       + DW1*(  &
     	       + weightp(i)*lbfr(i,jpl(i,j)) &
     	       + weightpn(i)*lbfr(i,jpn(i,j)) &
     	       + weightm(i)*rbfr(i,jmi(i,j)) &
     	       + weightmn(i)*rbfr(i,jmn(i,j)) ) &
     	        ) &
     	       / (1.+2*DW1)
	   enddo
	 enddo

	 htemp3d=ttemp3d

!  build buffers and send grid info to neighbors

	 do j=0,jm-1
	   do i=0,im-1
	     rbfs(i,j)=htemp3d(mykm-1,i,j)
	     lbfs(i,j)=htemp3d(0,i,j)
	   enddo
	 enddo

	 call MPI_SENDRECV(rbfs(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_REAL8,rngbr,9, &
     	                   lbfr(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_REAL8,lngbr,9, &
     		           TUBE_COMM,stat,ierr) 

	 call MPI_SENDRECV(lbfs(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_REAL8,lngbr,8,  &
     	                   rbfr(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_REAL8,rngbr,8, &
     		           TUBE_COMM,stat,ierr)
!	 write(*,*)'past buff'

	 do j=0,jm-1
	   do i=0,im-1
	     ttemp3d(0,i,j)=(  htemp3d(0,i,j) &
     	       + DW2*( nflag*htemp3d(1,i,j) &
     	       + weightp(i)*lbfr(i,jpl(i,j)) &
     	       + weightpn(i)*lbfr(i,jpn(i,j)) &
     	       + weightm(i)*rbfr(i,jmi(i,j)) &
     	       + weightmn(i)*rbfr(i,jmn(i,j)) ) &
     	        ) &
     	       / (1.+2*DW2)
	   enddo
	 enddo

	 enddo

!  assign ttemp3d(mykm,:,:)

      do i = 0,im-1
        do j = 0,jm-1
	   lbfs(i,j) = ttemp3d(0,i,j)
        end do
      end do

	 call MPI_SENDRECV(lbfs(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_REAL8,lngbr,10, &
     	                   rbfr(0:imx-1,0:jmx-1),imx*jmx, &
     		MPI_REAL8,rngbr,10, &
     		           TUBE_COMM,stat,ierr)

      do i = 0,im-1
         do j = 0,jm-1
	    ttemp3d(mykm,i,j) = weightm(i)*rbfr(i,jmi(i,j)) &
     	       + weightmn(i)*rbfr(i,jmn(i,j)) 
         end do
      end do
