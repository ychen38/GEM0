!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
MODULE gem_pputil
!
!  use fft_wrapper
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ppinit_mpi,ppinit_decomp,ppexit, init_pmove,end_pmove,pmove,guard
  PUBLIC :: ppsum, ppmax, ppmin
  PUBLIC :: pptransp, ppcfft2
  PUBLIC :: timera, disp
!
  INTEGER, SAVE :: me, nvp,npp,GCLR,TCLR
  INTEGER, SAVE :: pmove_tag=0
  INTEGER, SAVE :: TUBE_COMM,GRID_COMM
  REAL, DIMENSION(:), ALLOCATABLE, SAVE :: s_buf, r_buf
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: s_counts, s_displ
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: r_counts, r_displ
  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ipsend, iphole
!
  INTERFACE disp
     MODULE PROCEDURE dispi, dispr
     MODULE PROCEDURE disp2i, disp2r
  END INTERFACE
  INTERFACE ppsum
     MODULE PROCEDURE ppsum_r, ppsum_ra, ppsum_i, ppsum_ia
  END INTERFACE
  INTERFACE ppmax
     MODULE PROCEDURE ppmax_r, ppmax_ra, ppmax_i, ppmax_ia
  END INTERFACE
  INTERFACE ppmin
     MODULE PROCEDURE ppmin_r, ppmin_ra, ppmin_i, ppmin_ia
  END INTERFACE
  INTERFACE pptransp
     MODULE PROCEDURE pptransp_c, pptransp_r, pptransp_i
     MODULE PROCEDURE pptransp2_c, pptransp2_r, pptransp2_i
  END INTERFACE
  INTERFACE guard
     MODULE PROCEDURE guard2, guard3
  END INTERFACE
  INTERFACE PPCFFT2
     MODULE PROCEDURE ppcfft2_2d, ppcfft2_3d
  END INTERFACE
!
CONTAINS
!
!===========================================================================
  SUBROUTINE init_pmove(xp, np, lz, ierr)
    !
    use mpi
    !
    REAL, DIMENSION(:), INTENT(in) :: xp
    INTEGER, INTENT(in) :: np
    REAL, INTENT(in) :: lz
    INTEGER, INTENT(out) :: ierr
    !
    !  Local vars
    INTEGER :: nsize, ksize
    INTEGER :: i, ip, iz, ih, iwork
    REAL :: dzz, xt
    INTEGER, DIMENSION(0:nvp-1) :: isb
    !
    !----------------------------------------------------------------------
    !              0.   Allocate fixed size arrays
    !
    IF( .not. ALLOCATED(s_counts) ) ALLOCATE(s_counts(0:nvp-1))
    IF( .not. ALLOCATED(s_displ) ) ALLOCATE(s_displ(0:nvp-1))
    IF( .not. ALLOCATED(r_counts) ) ALLOCATE(r_counts(0:nvp-1))
    IF( .not. ALLOCATED(r_displ) ) ALLOCATE(r_displ(0:nvp-1))
    !
    !----------------------------------------------------------------------
    !              1.  Construct send buffer
    !
    dzz = lz / nvp
    s_counts = 0
    DO ip = 1,np
       xt = MODULO(xp(ip), lz)            !!! Assume periodicity
       iz = INT(xt/dzz)
       IF( iz .ne. GCLR )THEN
          s_counts(iz) = s_counts(iz)+1
       END IF
    END DO
    s_displ(0) = 0
    DO i=1,nvp-1
       s_displ(i) = s_displ(i-1) + s_counts(i-1)
    END DO
    !
    nsize = sum(s_counts)
    IF( .not. ALLOCATED(s_buf) ) THEN
       ksize=2*nsize         ! To prevent too much futur reallocations
       ALLOCATE(s_buf(1:ksize))
       ALLOCATE(ipsend(1:ksize))
       ALLOCATE(iphole(1:ksize))
    ELSE IF ( SIZE(s_buf) .LT. nsize ) THEN
       DEALLOCATE(s_buf)
       DEALLOCATE(ipsend)
       DEALLOCATE(iphole)
       ALLOCATE(s_buf(1:nsize))
       ALLOCATE(ipsend(1:nsize))
       ALLOCATE(iphole(1:nsize))
    END IF
    !
    !----------------------------------------------------------------------
    !              2.  Construct (sorted) pointers to holes
    !
    isb(0:nvp-1) = s_displ(0:nvp-1)
    ih = 0
    DO ip=1,np
       xt = MODULO(xp(ip), lz)            !!! Assume periodicity
       iz = INT(xt/dzz)
       IF( iz .ne. GCLR ) THEN
          isb(iz) = isb(iz)+1
          ipsend(isb(iz)) = ip
          ih = ih+1
          iphole(ih) = ip
       END IF
    END DO
    !
    !----------------------------------------------------------------------
    !              3.  Construct receive buffer
    !
    CALL MPI_ALLTOALL(s_counts, 1, MPI_INTEGER, &
         & r_counts, 1, MPI_INTEGER, TUBE_COMM, ierr)
    r_displ(0) = 0
    DO i=1,nvp-1
       r_displ(i) = r_displ(i-1) + r_counts(i-1)
    END DO
    !
    nsize = sum(r_counts)
    IF( .not. ALLOCATED(r_buf) ) THEN
       ksize=2*nsize         ! To prevent too much futur reallocations
       ALLOCATE(r_buf(1:ksize))
    ELSE IF ( SIZE(r_buf) .LT. nsize ) THEN
       DEALLOCATE(r_buf)
       ALLOCATE(r_buf(1:nsize))
    END IF
    !
    !  Check for part. array overflow
    ierr = 0
    nsize = np - sum(s_counts) + sum(r_counts)
    if( nsize .gt. size(xp) ) then
       write(*,*) 'PE', me, 'Particle array overflow'
       ierr = 1
    end if
    !    call ppsum(ierr)
    !
    !----------------------------------------------------------------------
    !              4.  Initialize tag
    !
    pmove_tag = 101
    !
  END SUBROUTINE init_pmove
!===========================================================================
  SUBROUTINE pmove(xp, np_old, np_new, ierr)
!
    use mpi
!
    REAL, DIMENSION(:), INTENT(inout) :: xp
    INTEGER, INTENT(in) :: np_old
    INTEGER, INTENT(out) :: np_new
    INTEGER, INTENT(out) :: ierr
!
!  Local vars
    INTEGER :: nsize
    INTEGER :: i, ip, iz, ih, isrc
    INTEGER :: nhole, mhole, nrrecv, nrsend, nptot_old, nptot
    INTEGER :: ind, count, tot_count, iwork
!
!  Local arrays
    INTEGER, DIMENSION(1:nvp) :: s_requ, r_requ, id_source
    INTEGER :: isrt, iend
    INTEGER :: status(MPI_STATUS_SIZE), arr_status(MPI_STATUS_SIZE,nvp)
!
!----------------------------------------------------------------------
!              1.  Fill send buffer
!
    DO i=0,nvp-1
       IF( s_counts(i) .GT. 0 ) THEN
          isrt = s_displ(i)+1
          iend = s_displ(i)+s_counts(i)
          s_buf(isrt:iend) = xp(ipsend(isrt:iend))
       END IF
    END DO
!----------------------------------------------------------------------
!              2.   Initiate non-blocking send/receive
!
    pmove_tag = pmove_tag+1
    nrrecv=0             !......... Start non-blocking receive
    DO i=0,nvp-1
       IF(r_counts(i) .GT. 0 ) THEN
          nrrecv=nrrecv+1
          id_source(nrrecv) = i
          isrt = r_displ(i)+1
          CALL MPI_IRECV(r_buf(isrt), r_counts(i), MPI_REAL8,&
               & i, pmove_tag, TUBE_COMM, r_requ(nrrecv), ierr)
       END IF
    END DO
    nrsend=0             !......... Start non-blocking SYNCHRONOUS send
    DO i=0,nvp-1
       IF(s_counts(i) .GT. 0 ) THEN
          nrsend=nrsend+1
          isrt = s_displ(i)+1
          CALL MPI_ISSEND(s_buf(isrt), s_counts(i), MPI_REAL8,&
               & i, pmove_tag, TUBE_COMM, s_requ(nrsend), ierr)
       END IF
    END DO
!
!----------------------------------------------------------------------
!              3.   Remove holes and compress part. arrays
!
    nhole = sum(s_counts)
    ip = np_old
    DO ih = nhole, 1, -1
       xp(iphole(ih)) = xp(ip)
       ip = ip-1
    END DO
    np_new = ip
!
!----------------------------------------------------------------------
!              4.   Store incoming part. to the part. arrays
!
    tot_count = 0
    DO i=1,nrrecv
       CALL MPI_WAITANY(nrrecv, r_requ, ind, status, ierr)
       isrc = id_source(ind)
       CALL MPI_GET_COUNT(status, MPI_REAL8, count, ierr)
       IF( count .ne. r_counts(isrc) ) THEN
          WRITE(*,*) 'PE',me, '  Counts mismatched from PE',isrc,&
               & count,  r_counts(isrc)
       END IF
       tot_count =  tot_count+count
    END DO
!
    IF( tot_count .GT. 0 ) THEN
       isrt = np_new + 1
       iend = np_new + tot_count
       xp(isrt:iend) = r_buf(1:tot_count)
       np_new = iend
    END IF
!
!----------------------------------------------------------------------
!              5.   Epilogue
!
!... Wait for any non-blocking comm. requests
    IF( nrsend.gt.0) CALL MPI_WAITALL(nrsend, s_requ, arr_status, ierr)
!
!... Check consistency
    ierr = 0
    CALL MPI_ALLREDUCE(np_old, nptot_old, 1, MPI_INTEGER,&
         & MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE(np_new, nptot, 1, MPI_INTEGER,&
         & MPI_SUM, MPI_COMM_WORLD, ierr)
    IF( nptot.ne.nptot_old ) THEN
       IF(me.eq.0) WRITE(*,*) 'PMOVE: mismatch in total numbers:',&
            & nptot_old, nptot
       ierr = 1
    END IF
!    call ppsum(ierr)
!
!----------------------------------------------------------------------!
  END SUBROUTINE pmove
!===========================================================================
  SUBROUTINE end_pmove(ierr)
!
    use mpi

    INTEGER, INTENT(OUT) :: ierr
!
!   Local vars
!----------------------------------------------------------------------!
  END SUBROUTINE end_pmove
!
!===========================================================================
  SUBROUTINE dispi(iarr,string)

    use mpi

    INTEGER, INTENT(in) :: iarr(:)
    character(len=*), INTENT(in) :: string
    INTEGER :: i, ierr
!
    DO i=0,nvp-1
       CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
       IF(i.eq.me) THEN
          WRITE(*,'(a,i2,2x,a,16(i6))') 'PE',me,string,iarr
       END IF
    END DO
!
  END SUBROUTINE dispi
!===========================================================================
  SUBROUTINE disp2i(arr,string)

    use mpi

    INTEGER, INTENT(in) :: arr(:,:)
    character(len=*), INTENT(in) :: string
    INTEGER :: i, j, ierr
!
    DO i=0,nvp-1
       CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
       IF(i.eq.me) THEN
          WRITE(*,'(a,i2,2x,a)') 'PE',me,string
          DO j=1,SIZE(arr,1)
             WRITE(*,'(16(i6))') arr(j,:)
          END DO
       END IF
    END DO
  END SUBROUTINE disp2i
!===========================================================================
  SUBROUTINE dispr(arr,string)

    use mpi

    REAL, INTENT(in) :: arr(:)
    character(len=*), INTENT(in) :: string
    INTEGER :: i, ierr
!
    DO i=0,nvp-1
       CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
       IF(i.eq.me) THEN
          WRITE(*,'(a,i2,2x,a,12(1pe10.2))') 'PE',me,string,arr
       END IF
    END DO
!
  END SUBROUTINE dispr
!===========================================================================
  SUBROUTINE disp2r(arr,string)
    use mpi
    REAL, INTENT(in) :: arr(:,:)
    character(len=*), INTENT(in) :: string
    INTEGER :: i, j, ierr
!
    DO i=0,nvp-1
       CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
       IF(i.eq.me) THEN
          WRITE(*,'(a,i2,2x,a)') 'PE',me,string
          DO j=1,SIZE(arr,1)
             WRITE(*,'(12(1pe10.2))') arr(j,:)
          END DO
       END IF
    END DO
!
  END SUBROUTINE disp2r
!======================================================================
!mpi initialization
subroutine ppinit_mpi(idproc,nproc)
  use mpi
  integer, intent(out) :: idproc,nproc
  integer :: ierr,npp

  call mpi_init(ierr)
  !set the size of the MPI in the mpi_comm_world MPI space
  call mpi_comm_size(mpi_comm_world,npp,ierr)
  !set the rank of the MPI in the mpi_comm_world MPI space
  call mpi_comm_rank(mpi_comm_world,me,ierr)
  nproc=npp
  idproc=me 
end subroutine ppinit_mpi
!particle decomp
subroutine ppinit_decomp(mype,nproc,ntube,com1,com2)
  use mpi
  integer, intent(in) :: mype,nproc,ntube
  integer, intent(out) :: com1,com2
  integer :: ierr

  !the rank of the particle decomposition
  gclr=int(mype/ntube)
  !the rank in the particle decomposition
  tclr=mod(mype,ntube)

  !set the particle MPI space
  call mpi_comm_split(mpi_comm_world,gclr,tclr,grid_comm,ierr)
  call mpi_comm_split(mpi_comm_world,tclr,gclr,tube_comm,ierr)

  com1=tube_comm
  com2=grid_comm
  nvp=nproc/ntube
end subroutine ppinit_decomp
!===========================================================================
  SUBROUTINE ppexit
    INTEGER :: ierr
    CALL MPI_FINALIZE(ierr)
    STOP
  END SUBROUTINE ppexit
!
!===========================================================================
!=========================GLOBAL SUM========================================
  SUBROUTINE ppsum_r(f)
    use mpi
    REAL, INTENT(INOUT) :: f
    REAL :: sum
    INTEGER :: ierr
!
    CALL MPI_ALLREDUCE(f, sum, 1, MPI_REAL8, MPI_SUM,&
         & MPI_COMM_WORLD, ierr)
    f=sum
  END SUBROUTINE ppsum_r
!===========================================================================
  SUBROUTINE ppsum_ra(f)
    use mpi
    REAL, DIMENSION (:), INTENT(INOUT) :: f
    REAL, DIMENSION (size(f)) :: sum
    INTEGER :: count, ierr
!
    count = size(f)
    CALL MPI_ALLREDUCE(f, sum, count, MPI_REAL8, MPI_SUM,&
         & MPI_COMM_WORLD, ierr)
    f=sum
  END SUBROUTINE ppsum_ra
!===========================================================================
  SUBROUTINE ppsum_i(f)
    use mpi
    INTEGER, INTENT(INOUT) :: f
    INTEGER :: sum
    INTEGER :: ierr
!
    CALL MPI_ALLREDUCE(f, sum, 1, MPI_INTEGER, MPI_SUM,&
         & MPI_COMM_WORLD, ierr)
    f=sum
  END SUBROUTINE ppsum_i
!===========================================================================
  SUBROUTINE ppsum_ia(f)
    use mpi
    INTEGER, DIMENSION (:), INTENT(INOUT) :: f
    INTEGER, DIMENSION (size(f)) :: sum
    INTEGER :: count, ierr
!
    count = size(f)
    CALL MPI_ALLREDUCE(f, sum, count, MPI_INTEGER, MPI_SUM,&
         & MPI_COMM_WORLD, ierr)
    f=sum
  END SUBROUTINE ppsum_ia
!===========================================================================
!=========================GLOBAL MAX========================================
  SUBROUTINE ppmax_r(f)
    use mpi
    REAL, INTENT(INOUT) :: f
    REAL :: max
    INTEGER :: ierr
!
    CALL MPI_ALLREDUCE(f, max, 1, MPI_REAL8, MPI_MAX,&
         & MPI_COMM_WORLD, ierr)
    f=max
  END SUBROUTINE ppmax_r
!===========================================================================
  SUBROUTINE ppmax_ra(f)
    use mpi
    REAL, DIMENSION (:), INTENT(INOUT) :: f
    REAL, DIMENSION (size(f)) :: max
    INTEGER :: count, ierr
!
    count = size(f)
    CALL MPI_ALLREDUCE(f, max, count, MPI_REAL8, MPI_MAX,&
         & MPI_COMM_WORLD, ierr)
    f=max
  END SUBROUTINE ppmax_ra
!===========================================================================
  SUBROUTINE ppmax_i(f)
    use mpi
    INTEGER, INTENT(INOUT) :: f
    INTEGER :: max
    INTEGER :: ierr
!
    CALL MPI_ALLREDUCE(f, max, 1, MPI_INTEGER, MPI_MAX,&
         & MPI_COMM_WORLD, ierr)
    f=max
  END SUBROUTINE ppmax_i
!===========================================================================
  SUBROUTINE ppmax_ia(f)
    use mpi
    INTEGER, DIMENSION (:), INTENT(INOUT) :: f
    INTEGER, DIMENSION (size(f)) :: max
    INTEGER :: count, ierr
!
    count = size(f)
    CALL MPI_ALLREDUCE(f, max, count, MPI_INTEGER, MPI_MAX,&
         & MPI_COMM_WORLD, ierr)
    f=max
  END SUBROUTINE ppmax_ia
!===========================================================================
!=========================GLOBAL MIN========================================
  SUBROUTINE ppmin_r(f)
    use mpi
    REAL, INTENT(INOUT) :: f
    REAL :: min
    INTEGER :: ierr
!
    CALL MPI_ALLREDUCE(f, min, 1, MPI_REAL8, MPI_MIN,&
         & MPI_COMM_WORLD, ierr)
    f=min
  END SUBROUTINE ppmin_r
!===========================================================================
  SUBROUTINE ppmin_ra(f)
    use mpi
    REAL, DIMENSION (:), INTENT(INOUT) :: f
    REAL, DIMENSION (size(f)) :: min
    INTEGER :: count, ierr
!
    count = size(f)
    CALL MPI_ALLREDUCE(f, min, count, MPI_REAL8, MPI_MIN,&
         & MPI_COMM_WORLD, ierr)
    f=min
  END SUBROUTINE ppmin_ra
!===========================================================================
  SUBROUTINE ppmin_i(f)
    use mpi
    INTEGER, INTENT(INOUT) :: f
    INTEGER :: min
    INTEGER :: ierr
!
    CALL MPI_ALLREDUCE(f, min, 1, MPI_INTEGER, MPI_MIN,&
         & MPI_COMM_WORLD, ierr)
    f=min
  END SUBROUTINE ppmin_i
!===========================================================================
  SUBROUTINE ppmin_ia(f)
    use mpi
    INTEGER, DIMENSION (:), INTENT(INOUT) :: f
    INTEGER, DIMENSION (size(f)) :: min
    INTEGER :: count, ierr
!
    count = size(f)
    CALL MPI_ALLREDUCE(f, min, count, MPI_INTEGER, MPI_MIN,&
         & MPI_COMM_WORLD, ierr)
    f=min
  END SUBROUTINE ppmin_ia
!=======================================================================
!==============================  PPTRANSP ==============================
  SUBROUTINE pptransp_c(a, b)
    use mpi
    complex, DIMENSION(:,:), INTENT(in) :: a
    complex, DIMENSION(:,:), INTENT(out) :: b
!
!   Local vars
    INTEGER :: nx, nxp, ny, nyp
    INTEGER :: i, j, isrt, iend, ierr
    INTEGER, DIMENSION(1:nvp) :: idsr
    complex, DIMENSION(1:SIZE(a,1)/nvp, 1:SIZE(a,2)) :: s_buf, r_buf
    INTEGER :: status(MPI_STATUS_SIZE)
!----------------------------------------------------------------------
!              0.   Check array dimensions
!
    nx = SIZE(a,1); nxp = nx/nvp
    nyp = SIZE(a,2); ny = nyp*nvp
    IF( ny .GT. SIZE(b,1) ) THEN
       WRITE(*,*) 'PPTRANSP: DIMENSION 1 of b too small'
       CALL ppexit
    END IF
    IF( nxp .GT. SIZE(b,2) ) THEN
       WRITE(*,*) 'PPTRANSP: DIMENSION 2 of b too small'
       CALL ppexit
    END IF
!
!   Determine send/receive proc. id
!
    DO i=1,nvp
       idsr(i) = iand(nvp-1, ieor(me, i-1))
    END DO
!----------------------------------------------------------------------
!              1.   Send/receive bloks
!
    DO i=1,nvp
       isrt = nxp*idsr(i) + 1
       iend = nxp*idsr(i) + nxp
       s_buf(1:nxp,1:nyp) = a(isrt:iend,1:nyp)
       CALL MPI_SENDRECV(s_buf, nxp*nyp, MPI_DOUBLE_COMPLEX, idsr(i), i,&
            & r_buf, nxp*nyp, MPI_DOUBLE_COMPLEX, idsr(i), i,&
            & MPI_COMM_WORLD, status, ierr)
       isrt = nyp*idsr(i) + 1
       iend = nyp*idsr(i) + nyp
       b(isrt:iend,1:nxp) = transpose(r_buf)
    END DO
  END SUBROUTINE pptransp_c
!=======================================================================
  SUBROUTINE pptransp_r(a, b)
    use mpi
    REAL, DIMENSION(:,:), INTENT(in) :: a
    REAL, DIMENSION(:,:), INTENT(out) :: b
!
!   Local vars
    INTEGER :: nx, nxp, ny, nyp
    INTEGER :: i, j, isrt, iend, ierr
    INTEGER, DIMENSION(1:nvp) :: idsr
    REAL, DIMENSION(1:SIZE(a,1)/nvp, 1:SIZE(a,2)) :: s_buf, r_buf
    INTEGER :: status(MPI_STATUS_SIZE)
!----------------------------------------------------------------------
!              0.   Check array dimensions
!
    nx = SIZE(a,1); nxp = nx/nvp
    nyp = SIZE(a,2); ny = nyp*nvp
    IF( ny .GT. SIZE(b,1) ) THEN
       WRITE(*,*) 'PPTRANSP: DIMENSION 1 of b too small'
       CALL ppexit
    END IF
    IF( nxp .GT. SIZE(b,2) ) THEN
       WRITE(*,*) 'PPTRANSP: DIMENSION 2 of b too small'
       CALL ppexit
    END IF
!
!   Determine send/receive proc. id
!
    DO i=1,nvp
       idsr(i) = iand(nvp-1, ieor(me, i-1))
    END DO
!----------------------------------------------------------------------
!              1.   Send/receive bloks
!
    DO i=1,nvp
       isrt = nxp*idsr(i) + 1
       iend = nxp*idsr(i) + nxp
       s_buf(1:nxp,1:nyp) = a(isrt:iend,1:nyp)
       CALL MPI_SENDRECV(s_buf, nxp*nyp, MPI_REAL8, idsr(i), i,&
            & r_buf, nxp*nyp, MPI_REAL8, idsr(i), i,&
            & MPI_COMM_WORLD, status, ierr)
       isrt = nyp*idsr(i) + 1
       iend = nyp*idsr(i) + nyp
       b(isrt:iend,1:nxp) = transpose(r_buf)
    END DO
  END SUBROUTINE pptransp_r
!=======================================================================
  SUBROUTINE pptransp_i(a, b)
    use mpi
    INTEGER, DIMENSION(:,:), INTENT(in) :: a
    INTEGER, DIMENSION(:,:), INTENT(out) :: b
!
!   Local vars
    INTEGER :: nx, nxp, ny, nyp
    INTEGER :: i, j, isrt, iend, ierr
    INTEGER, DIMENSION(1:nvp) :: idsr
    INTEGER, DIMENSION(1:SIZE(a,1)/nvp, 1:SIZE(a,2)) :: s_buf, r_buf
    INTEGER :: status(MPI_STATUS_SIZE)
!----------------------------------------------------------------------
!              0.   Check array dimensions
!
    nx = SIZE(a,1); nxp = nx/nvp
    nyp = SIZE(a,2); ny = nyp*nvp
    IF( ny .GT. SIZE(b,1) ) THEN
       WRITE(*,*) 'PPTRANSP: DIMENSION 1 of b too small'
       CALL ppexit
    END IF
    IF( nxp .GT. SIZE(b,2) ) THEN
       WRITE(*,*) 'PPTRANSP: DIMENSION 2 of b too small'
       CALL ppexit
    END IF
!
!   Determine send/receive proc. id
!
    DO i=1,nvp
       idsr(i) = iand(nvp-1, ieor(me, i-1))
    END DO
!----------------------------------------------------------------------
!              1.   Send/receive bloks
!
    DO i=1,nvp
       isrt = nxp*idsr(i) + 1
       iend = nxp*idsr(i) + nxp
       s_buf(1:nxp,1:nyp) = a(isrt:iend,1:nyp)
       CALL MPI_SENDRECV(s_buf, nxp*nyp, MPI_INTEGER, idsr(i), i,&
            & r_buf, nxp*nyp, MPI_INTEGER, idsr(i), i,&
            & MPI_COMM_WORLD, status, ierr)
       isrt = nyp*idsr(i) + 1
       iend = nyp*idsr(i) + nyp
       b(isrt:iend,1:nxp) = transpose(r_buf)
    END DO
  END SUBROUTINE pptransp_i
!===========================================================================
!=======================================================================
!==============================  PPTRANSP2==============================
!=======================================================================
  SUBROUTINE pptransp2_c(a, b)
    use mpi
    complex, DIMENSION(:,:,:), INTENT(in) :: a
    complex, DIMENSION(:,:,:), INTENT(out) :: b
!
!   Local vars
    INTEGER :: nx, nxp, ny, nyp, nz, nzp
    INTEGER :: j, iter, isrt, iend, ierr
    INTEGER, DIMENSION(1:nvp) :: idsr
    complex, DIMENSION(1:SIZE(a,1)/nvp, 1:SIZE(a,2), 1:SIZE(a,3)) :: s_buf, r_buf
    INTEGER :: status(MPI_STATUS_SIZE)
!----------------------------------------------------------------------
!              0.   Check array dimensions
!
    nx = SIZE(a,1); nxp = nx/nvp
    ny = SIZE(a,2); nyp = ny/nvp
    nzp = SIZE(a,3); nz = nzp*nvp
    IF( nz .GT. SIZE(b,1) ) THEN
       WRITE(*,*) 'PPTRANSP2: DIMENSION 1 of b too small'
       CALL ppexit
    END IF
    IF( ny .GT. SIZE(b,2) ) THEN
       WRITE(*,*) 'PPTRANSP2: DIMENSION 2 of b too small'
       CALL ppexit
    END IF
    IF( nxp .GT. SIZE(b,3) ) THEN
       WRITE(*,*) 'PPTRANSP2: DIMENSION 3 of b too small'
       CALL ppexit
    END IF
!
!   Determine send/receive proc. id
!
    DO iter=1,nvp
       idsr(iter) = iand(nvp-1, ieor(me, iter-1))
    END DO
!----------------------------------------------------------------------
!              1.   Send/receive bloks
!
    DO iter=1,nvp
       isrt = nxp*idsr(iter) + 1
       iend = nxp*idsr(iter) + nxp
       s_buf(1:nxp,1:ny,1:nzp) = a(isrt:iend,1:ny,1:nzp)
       CALL MPI_SENDRECV(s_buf, nxp*ny*nzp, MPI_DOUBLE_COMPLEX, idsr(iter), iter,&
            & r_buf, nxp*ny*nzp, MPI_DOUBLE_COMPLEX, idsr(iter), iter,&
            & MPI_COMM_WORLD, status, ierr)
       isrt = nzp*idsr(iter) + 1
       iend = nzp*idsr(iter) + nzp
       DO j=1,ny
          b(isrt:iend,j,1:nxp) = transpose(r_buf(1:nxp,j,1:nzp))
       END DO
    END DO
  END SUBROUTINE pptransp2_c
!=======================================================================
  SUBROUTINE pptransp2_r(a, b)
    use mpi
    REAL, DIMENSION(:,:,:), INTENT(in) :: a
    REAL, DIMENSION(:,:,:), INTENT(out) :: b
!
!   Local vars
    INTEGER :: nx, nxp, ny, nyp, nz, nzp
    INTEGER :: j, iter, isrt, iend, ierr
    INTEGER, DIMENSION(1:nvp) :: idsr
    REAL, DIMENSION(1:SIZE(a,1)/nvp, 1:SIZE(a,2), 1:SIZE(a,3)) :: s_buf, r_buf
    INTEGER :: status(MPI_STATUS_SIZE)
!----------------------------------------------------------------------
!              0.   Check array dimensions
!
    nx = SIZE(a,1); nxp = nx/nvp
    ny = SIZE(a,2); nyp = ny/nvp
    nzp = SIZE(a,3); nz = nzp*nvp
    IF( nz .GT. SIZE(b,1) ) THEN
       WRITE(*,*) 'PPTRANSP2: DIMENSION 1 of b too small'
       CALL ppexit
    END IF
    IF( ny .GT. SIZE(b,2) ) THEN
       WRITE(*,*) 'PPTRANSP2: DIMENSION 2 of b too small'
       CALL ppexit
    END IF
    IF( nxp .GT. SIZE(b,3) ) THEN
       WRITE(*,*) 'PPTRANSP2: DIMENSION 3 of b too small'
       CALL ppexit
    END IF
!
!   Determine send/receive proc. id
!
    DO iter=1,nvp
       idsr(iter) = iand(nvp-1, ieor(me, iter-1))
    END DO
!----------------------------------------------------------------------
!              1.   Send/receive bloks
!
    DO iter=1,nvp
       isrt = nxp*idsr(iter) + 1
       iend = nxp*idsr(iter) + nxp
       s_buf(1:nxp,1:ny,1:nzp) = a(isrt:iend,1:ny,1:nzp)
       CALL MPI_SENDRECV(s_buf, nxp*ny*nzp, MPI_REAL8, idsr(iter), iter,&
            & r_buf, nxp*ny*nzp, MPI_REAL8, idsr(iter), iter,&
            & MPI_COMM_WORLD, status, ierr)
       isrt = nzp*idsr(iter) + 1
       iend = nzp*idsr(iter) + nzp
       DO j=1,ny
          b(isrt:iend,j,1:nxp) = transpose(r_buf(1:nxp,j,1:nzp))
       END DO
    END DO
  END SUBROUTINE pptransp2_r
!=======================================================================
  SUBROUTINE pptransp2_i(a, b)
    use mpi
    INTEGER, DIMENSION(:,:,:), INTENT(in) :: a
    INTEGER, DIMENSION(:,:,:), INTENT(out) :: b
!
!   Local vars
    INTEGER :: nx, nxp, ny, nyp, nz, nzp
    INTEGER :: j, iter, isrt, iend, ierr
    INTEGER, DIMENSION(1:nvp) :: idsr
    INTEGER, DIMENSION(1:SIZE(a,1)/nvp, 1:SIZE(a,2), 1:SIZE(a,3)) :: s_buf, r_buf
    INTEGER :: status(MPI_STATUS_SIZE)
!----------------------------------------------------------------------
!              0.   Check array dimensions
!
    nx = SIZE(a,1); nxp = nx/nvp
    ny = SIZE(a,2); nyp = ny/nvp
    nzp = SIZE(a,3); nz = nzp*nvp
    IF( nz .GT. SIZE(b,1) ) THEN
       WRITE(*,*) 'PPTRANSP2: DIMENSION 1 of b too small'
       CALL ppexit
    END IF
    IF( ny .GT. SIZE(b,2) ) THEN
       WRITE(*,*) 'PPTRANSP2: DIMENSION 2 of b too small'
       CALL ppexit
    END IF
    IF( nxp .GT. SIZE(b,3) ) THEN
       WRITE(*,*) 'PPTRANSP2: DIMENSION 3 of b too small'
       CALL ppexit
    END IF
!
!   Determine send/receive proc. id
!
    DO iter=1,nvp
       idsr(iter) = iand(nvp-1, ieor(me, iter-1))
    END DO
!----------------------------------------------------------------------
!              1.   Send/receive bloks
!
    DO iter=1,nvp
       isrt = nxp*idsr(iter) + 1
       iend = nxp*idsr(iter) + nxp
       s_buf(1:nxp,1:ny,1:nzp) = a(isrt:iend,1:ny,1:nzp)
       CALL MPI_SENDRECV(s_buf, nxp*ny*nzp, MPI_INTEGER, idsr(iter), iter,&
            & r_buf, nxp*ny*nzp, MPI_INTEGER, idsr(iter), iter,&
            & MPI_COMM_WORLD, status, ierr)
       isrt = nzp*idsr(iter) + 1
       iend = nzp*idsr(iter) + nzp
       DO j=1,ny
          b(isrt:iend,j,1:nxp) = transpose(r_buf(1:nxp,j,1:nzp))
       END DO
    END DO
  END SUBROUTINE pptransp2_i
!===========================================================================
  SUBROUTINE timera(icntrl,string,eltime)
    use mpi
    INTEGER, INTENT(in) :: icntrl
    CHARACTER(len=*), INTENT(in) :: string
    REAL, OPTIONAL, INTENT (OUT) :: eltime
    REAL, SAVE :: startt=0.0
    REAL :: endt
    INTEGER :: ierr
!
    SELECT CASE (icntrl)
       CASE (-1)
          CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
          startt = MPI_WTIME()
          IF( PRESENT(eltime) ) eltime=startt
       CASE (1)
          CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
          endt = MPI_WTIME() - startt
          IF( PRESENT(eltime) ) THEN
             eltime=endt
          ELSE IF( me .EQ. 0 ) THEN
             write(*,'(a,a,1pe14.7)') string,' wall clock time = ', endt
          END IF
    END SELECT
  END SUBROUTINE timera
!
!
!===========================================================================
  SUBROUTINE guard2(f, nidbas, flag)
!
!  nidbas = 1 : Linear Spline
!  nidbas = 2 : Quadratic Spline
!  nidbas = 3 : Cubic Spline
!
!  flag=0 : Add guard cell data to field
!  flag=1 : Copy field to guard cells
!
    use mpi
    REAL, DIMENSION(:,:), INTENT(INOUT) :: f
    INTEGER, INTENT(IN) :: nidbas, flag
!
!  Local vars and arrays
    INTEGER :: left, right, ierr, tag=200, status(MPI_STATUS_SIZE)
    INTEGER :: n1, np
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: buffer
!----------------------------------------------------------------------
!
!!! Left and right node ( assume periodicity)
    left = MODULO(me-1, nvp)
    right = MODULO(me+1, nvp)
!
    n1 = SIZE(f,1)
    np = SIZE(f,2)
    IF( .NOT. ALLOCATED(buffer) ) ALLOCATE(buffer(n1))
!
    SELECT CASE(nidbas)
       CASE(1)
          IF(flag .EQ. 0) THEN
             CALL guard_lin_add
          ELSE
             CALL guard_lin_copy
          END IF
       CASE(2)
          IF(flag .EQ. 0) THEN
             CALL guard_quad_add
          ELSE
             CALL guard_quad_copy
          END IF
       CASE(3)
          IF(flag .EQ. 0) THEN
             CALL guard_cub_add
          ELSE
             CALL guard_cub_copy
          END IF
       CASE DEFAULT
          WRITE(*,*) 'NIDBAS = ',nidbas,' NOT IMPLEMENTED'
          CALL PPEXIT
    END SELECT
!
  CONTAINS
!...
    SUBROUTINE guard_lin_add
      tag = tag+1
      CALL MPI_SENDRECV(f(:,np), n1, MPI_REAL8, right, tag,&
           & buffer, n1, MPI_REAL8, left, tag,&
           & MPI_COMM_WORLD, status, ierr)
      f(:,1) = f(:,1) + buffer
      f(:,np) = 0.0
    END SUBROUTINE guard_lin_add
!...
    SUBROUTINE guard_quad_add
      tag = tag+1
      CALL MPI_SENDRECV(f(:,np), n1, MPI_REAL8, right, tag,&
           & buffer, n1, MPI_REAL8, left, tag,&
           & MPI_COMM_WORLD, status, ierr)
      f(:,2) = f(:,2) + buffer
      f(:,np) = 0.0
      tag = tag+1
      CALL MPI_SENDRECV(f(:,1), n1, MPI_REAL8, left, tag,&
           & buffer, n1, MPI_REAL8, right, tag,&
           & MPI_COMM_WORLD, status, ierr)
      f(:,np-1) = f(:,np-1) + buffer
      f(:,1) = 0.0
    END SUBROUTINE guard_quad_add
!...
    SUBROUTINE guard_cub_add
      tag = tag+1
      CALL MPI_SENDRECV(f(:,np), n1, MPI_REAL8, right, tag,&
           & buffer, n1, MPI_REAL8, left, tag,&
           & MPI_COMM_WORLD, status, ierr)
      f(:,3) = f(:,3) + buffer
      f(:,np) = 0.0
      tag = tag+1
      CALL MPI_SENDRECV(f(:,np-1), n1, MPI_REAL8, right, tag,&
           & buffer, n1, MPI_REAL8, left, tag,&
           & MPI_COMM_WORLD, status, ierr)
      f(:,2) = f(:,2) + buffer
      f(:,np-1) = 0.0
      tag = tag+1
      CALL MPI_SENDRECV(f(:,1), n1, MPI_REAL8, left, tag,&
           & buffer, n1, MPI_REAL8, right, tag,&
           & MPI_COMM_WORLD, status, ierr)
      f(:,np-2) = f(:,np-2) + buffer
      f(:,1) = 0.0
    END SUBROUTINE guard_cub_add
!...
    SUBROUTINE guard_lin_copy
      tag = tag+1
      CALL MPI_SENDRECV(f(:,1), n1, MPI_REAL8, left, tag,&
           & f(:,np), n1, MPI_REAL8, right, tag,&
           & MPI_COMM_WORLD, status, ierr)
    END SUBROUTINE guard_lin_copy
!...
    SUBROUTINE guard_quad_copy
      tag = tag+1
      CALL MPI_SENDRECV(f(:,2), n1, MPI_REAL8, left, tag,&
           & f(:,np), n1, MPI_REAL8, right, tag,&
           & MPI_COMM_WORLD, status, ierr)
      tag = tag+1
      CALL MPI_SENDRECV(f(:,np-1), n1, MPI_REAL8, right, tag,&
           & f(:,1), n1, MPI_REAL8, left, tag,&
           & MPI_COMM_WORLD, status, ierr)
    END SUBROUTINE guard_quad_copy
!...
    SUBROUTINE guard_cub_copy
      tag = tag+1
      CALL MPI_SENDRECV(f(:,2), n1, MPI_REAL8, left, tag,&
           & f(:,np-1), n1, MPI_REAL8, right, tag,&
           & MPI_COMM_WORLD, status, ierr)
      tag = tag+1
      CALL MPI_SENDRECV(f(:,3), n1, MPI_REAL8, left, tag,&
           & f(:,np), n1, MPI_REAL8, right, tag,&
           & MPI_COMM_WORLD, status, ierr)
      tag = tag+1
      CALL MPI_SENDRECV(f(:,np-2), n1, MPI_REAL8, right, tag,&
           & f(:,1), n1, MPI_REAL8, left, tag,&
           & MPI_COMM_WORLD, status, ierr)
    END SUBROUTINE guard_cub_copy
!...
  END SUBROUTINE guard2
!===========================================================================
  SUBROUTINE guard3(f, nidbas, flag)
!
!  nidbas = 1 : Linear Spline
!  nidbas = 2 : Quadratic Spline
!  nidbas = 3 : Cubic Spline
!
!  flag=0 : Add guard cell data to field
!  flag=1 : Copy field to guard cells
!
    use mpi
    REAL, DIMENSION(:,:,:), INTENT(INOUT) :: f
    INTEGER, INTENT(IN) :: nidbas, flag
!
!  Local vars and arrays
    INTEGER :: left, right, ierr, tag=200, status(MPI_STATUS_SIZE)
    INTEGER :: n1, n2, n1n2, np
    REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: buffer
!----------------------------------------------------------------------
!
!!! Left and right node ( assume periodicity)
    left = MODULO(me-1, nvp)
    right = MODULO(me+1, nvp)
!
    n1 = SIZE(f,1)
    n2 = SIZE(f,2)
    n1n2= n1*n2
    np = SIZE(f,3)
    IF( .NOT. ALLOCATED(buffer) ) ALLOCATE(buffer(n1,n2))
!
    SELECT CASE(nidbas)
       CASE(1)
          IF(flag .EQ. 0) THEN
             CALL guard_lin_add
          ELSE
             CALL guard_lin_copy
          END IF
       CASE(2)
          IF(flag .EQ. 0) THEN
             CALL guard_quad_add
          ELSE
             CALL guard_quad_copy
          END IF
       CASE(3)
          IF(flag .EQ. 0) THEN
             CALL guard_cub_add
          ELSE
             CALL guard_cub_copy
          END IF
       CASE DEFAULT
          WRITE(*,*) 'NIDBAS = ',nidbas,' NOT IMPLEMENTED'
          CALL PPEXIT
    END SELECT
!
  CONTAINS
!...
    SUBROUTINE guard_lin_add
      tag = tag+1
      CALL MPI_SENDRECV(f(:,:,np), n1n2, MPI_REAL8, right, tag,&
           & buffer, n1n2, MPI_REAL8, left, tag,&
           & MPI_COMM_WORLD, status, ierr)
      f(:,:,1) = f(:,:,1) + buffer
      f(:,:,np) = 0.0
    END SUBROUTINE guard_lin_add
!...
    SUBROUTINE guard_quad_add
      tag = tag+1
      CALL MPI_SENDRECV(f(:,:,np), n1n2, MPI_REAL8, right, tag,&
           & buffer, n1n2, MPI_REAL8, left, tag,&
           & MPI_COMM_WORLD, status, ierr)
      f(:,:,2) = f(:,:,2) + buffer
      f(:,:,np) = 0.0
      tag = tag+1
      CALL MPI_SENDRECV(f(:,:,1), n1n2, MPI_REAL8, left, tag,&
           & buffer, n1n2, MPI_REAL8, right, tag,&
           & MPI_COMM_WORLD, status, ierr)
      f(:,:,np-1) = f(:,:,np-1) + buffer
      f(:,:,1) = 0.0
    END SUBROUTINE guard_quad_add
!...
    SUBROUTINE guard_cub_add
      tag = tag+1
      CALL MPI_SENDRECV(f(:,:,np), n1n2, MPI_REAL8, right, tag,&
           & buffer, n1n2, MPI_REAL8, left, tag,&
           & MPI_COMM_WORLD, status, ierr)
      f(:,:,3) = f(:,:,3) + buffer
      f(:,:,np) = 0.0
      tag = tag+1
      CALL MPI_SENDRECV(f(:,:,np-1), n1n2, MPI_REAL8, right, tag,&
           & buffer, n1n2, MPI_REAL8, left, tag,&
           & MPI_COMM_WORLD, status, ierr)
      f(:,:,2) = f(:,:,2) + buffer
      f(:,:,np-1) = 0.0
      tag = tag+1
      CALL MPI_SENDRECV(f(:,:,1), n1n2, MPI_REAL8, left, tag,&
           & buffer, n1n2, MPI_REAL8, right, tag,&
           & MPI_COMM_WORLD, status, ierr)
      f(:,:,np-2) = f(:,:,np-2) + buffer
      f(:,:,1) = 0.0
    END SUBROUTINE guard_cub_add
!...
    SUBROUTINE guard_lin_copy
      tag = tag+1
      CALL MPI_SENDRECV(f(:,:,1), n1n2, MPI_REAL8, left, tag,&
           & f(:,:,np), n1n2, MPI_REAL8, right, tag,&
           & MPI_COMM_WORLD, status, ierr)
    END SUBROUTINE guard_lin_copy
!...
    SUBROUTINE guard_quad_copy
      tag = tag+1
      CALL MPI_SENDRECV(f(:,:,2), n1n2, MPI_REAL8, left, tag,&
           & f(:,:,np), n1n2, MPI_REAL8, right, tag,&
           & MPI_COMM_WORLD, status, ierr)
      tag = tag+1
      CALL MPI_SENDRECV(f(:,:,np-1), n1n2, MPI_REAL8, right, tag,&
           & f(:,:,1), n1n2, MPI_REAL8, left, tag,&
           & MPI_COMM_WORLD, status, ierr)
    END SUBROUTINE guard_quad_copy
!...
    SUBROUTINE guard_cub_copy
      tag = tag+1
      CALL MPI_SENDRECV(f(:,:,2), n1n2, MPI_REAL8, left, tag,&
           & f(:,:,np-1), n1n2, MPI_REAL8, right, tag,&
           & MPI_COMM_WORLD, status, ierr)
      tag = tag+1
      CALL MPI_SENDRECV(f(:,:,3), n1n2, MPI_REAL8, left, tag,&
           & f(:,:,np), n1n2, MPI_REAL8, right, tag,&
           & MPI_COMM_WORLD, status, ierr)
      tag = tag+1
      CALL MPI_SENDRECV(f(:,:,np-2), n1n2, MPI_REAL8, right, tag,&
           & f(:,:,1), n1n2, MPI_REAL8, left, tag,&
           & MPI_COMM_WORLD, status, ierr)
    END SUBROUTINE guard_cub_copy
!...
  END SUBROUTINE guard3
!===========================================================================
  SUBROUTINE PPCFFT2_2D(isign, f, g)
!
!   Perform a 2d fft: the first variable is the first array dimension
!   which is local while the second variable is the last dimension which
!   partitionned
!
    INTEGER, INTENT(IN) :: isign
    complex, DIMENSION (:,:), INTENT(INOUT) :: f
    complex, DIMENSION (:,:), INTENT(OUT) :: g
!
!  Local vars and arrays
    INTEGER :: nx, nyp, ny, nxp
    REAL :: tablex(2*SIZE(f,1)), tabley(2*SIZE(g,1))
    REAL :: workx(4*SIZE(f,1)), worky(4*SIZE(g,1))
    complex :: dummy
    INTEGER :: i, j
!
    nx = SIZE(f,1); nyp = SIZE(f,2)
    ny = SIZE(g,1); nxp = SIZE(g,2)
!
!    call ccfft(0, nx, 1.d0, dummy, dummy, tablex, workx, 0)
!    call ccfft(0, ny, 1.d0, dummy, dummy, tabley, worky, 0)
!
    do j=1,nyp
!       call ccfft(isign, nx, 1.d0, f(1,j), f(1,j), tablex, workx, 0)
    end do
    CALL pptransp(f, g)
    do i=1,nxp
!       call ccfft(isign, ny, 1.d0, g(1,i), g(1,i), tabley, worky, 0)
    end do
  END SUBROUTINE PPCFFT2_2D
!===========================================================================
  SUBROUTINE PPCFFT2_3D(isign, f, g)
!
!   Perform a 2d fft: the first variable is the first array dimension
!   which is local while the second variable is the last dimension which
!   partitionned
!
    INTEGER, INTENT(IN) :: isign
    complex, DIMENSION (:,:,:), INTENT(INOUT) :: f
    complex, DIMENSION (:,:,:), INTENT(OUT) :: g
!
!  Local vars and arrays
    INTEGER :: nx, nyp, ny, nxp, nz
    REAL :: tablex(2*SIZE(f,1)), tabley(2*SIZE(g,1))
    REAL :: workx(4*SIZE(f,1)), worky(4*SIZE(g,1))
    complex :: dummy
    INTEGER :: i, j, k
!
    nx = SIZE(f,1); nyp = SIZE(f,3)
    ny = SIZE(g,1); nxp = SIZE(g,3)
    nz = SIZE(f,2)
!
!    CALL ccfft(0, nx, 1.d0, dummy, dummy, tablex, workx, 0)
!    CALL ccfft(0, ny, 1.d0, dummy, dummy, tabley, worky, 0)
!
    DO j=1,nyp
       DO k=1,nz
!          CALL ccfft(isign, nx, 1.d0, f(1,k,j), f(1,k,j), tablex, workx, 0)
       END DO
    END DO
    CALL pptransp(f, g)
    DO i=1,nxp
       DO k=1,nz
!          CALL ccfft(isign, ny, 1.d0, g(1,k,i), g(1,k,i), tabley, worky, 0)
       END DO
    END DO
  END SUBROUTINE PPCFFT2_3D
!===========================================================================
END MODULE gem_pputil
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
