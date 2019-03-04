module gem_com

  !common data used for gem

  use mpi
  use gem_pputil

  implicit none

  INTERFACE
     real function revers(num,n)
     end function revers

     real function ran2(i)
     end function ran2

     real function en3(s)
       real :: s
     end function en3
  END INTERFACE

  integer :: imx,jmx,kmx,mmx,mmxe,nmx,nsmx,nsubd=8,&
       modemx,ntube,nxpp,ngdx=5,nb=6, &
       negrd=8,nlgrd=8

  character(len=70) outname
  REAL :: endtm,begtm,pstm
  REAL :: starttm,lasttm,tottm
  real :: aux1(50000),aux2(20000)
  real,dimension(:),allocatable :: workx,worky,workz
  complex,dimension(:),allocatable :: tmpx
  complex,dimension(:),allocatable :: tmpy
  complex,dimension(:),allocatable :: tmpz

  !          imx,jmx,kmx = max no. of grid pts in x,y,z
  !          mmx         = max no. of particles
  !          nmx         = max. no. of time steps
  !          nsmx        = max. no. of species (including tracer particles
  !                                             as a species)
  !          ntrmx       = max. no. of tracer particles
  !          modemx      = max. no. of mode history plots

  integer :: mme,mmb
  REAL, dimension(:,:),allocatable :: rwx,rwy
  INTEGER,dimension(:),allocatable :: mm,tmm,lr
  integer :: micell,mecell !jycheng
  integer :: nonlin1,nonlin2 !jycheng
  real :: mims3,q3 !jycheng
  REAL,dimension(:),allocatable :: tets,mims,q
  REAL,dimension(:),allocatable :: kapn, kapt
  INTEGER :: timestep,im,jm,km,mykm,iseed,nrst,nfreq,isft,mynf,ifskp,iphbf,iapbf,idpbf
  real,dimension(:),allocatable :: time
  REAL :: dx,dy,dz,pi,pi2,dt,dte,totvol,n0,n0e,tcurr,rmpp,rmaa,eprs
  REAL :: lx,ly,lz,xshape,yshape,zshape,pzcrit(5),pzcrite,encrit,tot_field_e,tot_joule,tot_joule1
  INTEGER :: nm,nsm,kcnt,jcnt,ncurr,llk,mlk,onemd,iflr,iorb
  integer :: izonal,adiabatic_electron,ineq0,iflut,nlow,ntor0,mstart
  REAL :: cut,amp,tor,amie,isg,rneu,rneui,emass,qel,mbeam,qbeam,teth,vexbsw,vparsw
  REAL :: c4,fradi,kxcut,kycut,bcut,ftrap,adwn,adwe,adwp,frmax
  INTEGER :: iput,iget,idg,kzlook,ision,isiap,peritr,iadi,ipred,icorr,jpred,jcorr
  REAL,DIMENSION(:,:),allocatable :: yyamp,yyre,yyim
  complex,dimension(:,:),allocatable :: camp,campf
  REAL :: br0,lr0,qp,width,e0,vwidth,vwidthe,vcut,vpp,vt0,yd0
  integer :: nonlin(5),nonline,ipara,isuni,ifluid,ishift,nopz,nopi(5),noen,nowe
  complex :: IU
  real,dimension(:),allocatable :: coefx,coefy,coefz
  complex,dimension(1:8) :: apk,ptk,dpdtk
  integer,dimension(1:8) :: lapa,mapa,napa
  real :: mrtio(0:1),aven,avptch
  integer :: icrs_sec,ipg,isphi
  integer,dimension(0:255) :: isgnft,jft

  !          im,jm,km = max no. of grid pts in x,y,z
  !          mm       = max no. of particles
  !          nm       = max. no. of time steps
  !          nsm      = max. no. of species (including the set of
  !                           tracer particles as a species)
  !          ntrm     = max. no. of tracer particles
  !          modem    = max. no. of mode history plots
  !          iput     = 1 save for restart into dump.b, =0 no save
  !          iget     = 1 restart from dump.b, =0 use loader

  !            field or grid quantities

  REAL,DIMENSION(:,:,:,:),allocatable :: den
  REAL,DIMENSION(:,:,:,:),allocatable :: dnidt,jpar,jpex,jpey,dti
  REAL,DIMENSION(:,:,:),allocatable :: rho,jion,jionx,jiony
  real,dimension(:,:,:),allocatable :: phi
  real,dimension(:,:,:),allocatable :: drhodt,dnedt,dphidt,drhoidt
  REAL,DIMENSION(:,:,:),allocatable :: ex
  REAL,DIMENSION(:,:,:),allocatable :: ey
  REAL,DIMENSION(:,:,:),allocatable :: ez
  REAL,DIMENSION(:,:,:),allocatable :: dpdz,dadz
  REAL,DIMENSION(:,:,:),allocatable :: delbx,delby
  REAL,DIMENSION(:),allocatable :: xg,yg,zg
  real,dimension(:,:,:),allocatable :: apar,dene
  real,dimension(:,:,:),allocatable :: upar,upart,delte
  real,dimension(:,:,:),allocatable :: upex,upey,upa0,den0,upazd,upa00,upa0t,den0apa
  real,dimension(:,:),allocatable :: cfx,cfy,jac,bmag,bdgxcgy,bdgrzn,ggxdgy,ggy2,ggx
  real,dimension(:),allocatable :: gn0e,gt0e,gt0i,avap 
  real,dimension(:,:),allocatable :: gn0s

  !          particle array declarations
  REAL,DIMENSION(:,:),allocatable :: mu,xii,pzi,eki,z0i,u0i
  REAL,DIMENSION(:,:),allocatable :: x2,y2,z2,u2
  REAL,DIMENSION(:,:),allocatable :: x3,y3,z3,u3
  REAL,DIMENSION(:,:),allocatable :: w2,w3

  REAL,DIMENSION(:),allocatable :: mue,xie,pze,eke,z0e,u0e
  REAL,DIMENSION(:),allocatable :: x2e,y2e,z2e,u2e,mue2
  REAL,DIMENSION(:),allocatable :: x3e,y3e,z3e,u3e,mue3
  REAL,DIMENSION(:),allocatable :: w2e,w3e
  real,dimension(:),allocatable :: ipass, index
  REAL,DIMENSION(:),allocatable :: w000,w001,w010,w011,w100,w101,w110,w111

  !              Various diagnostic arrays and scalars
  !    plotting constants

  INTEGER :: nplot,xnplt,imovie=1000000,nzcrt,npze,npzi,npzc,npzb
  REAL :: contu,wmax

  !    energy diagnostic arrays

  REAL,DIMENSION(:,:),allocatable :: ke
  REAL,DIMENSION(:),allocatable :: fe,te
  REAL,DIMENSION(:),allocatable :: rmsphi,rmsapa,avewe
  REAL,DIMENSION(:,:),allocatable :: nos,avewi

  !    flux diagnostics
  REAL,DIMENSION(:),allocatable :: vol
  REAL,DIMENSION(:,:),allocatable :: efle_es,efle_em,pfle_es,pfle_em
  REAL,DIMENSION(:,:,:),allocatable :: pfl_es,pfl_em,efl_es,efl_em
  REAL,DIMENSION(:,:),allocatable :: chii, chie, ddi
  REAL,DIMENSION(:),allocatable :: achii, achie, addi

  !   mode diagnositics
  INTEGER :: modem
  INTEGER,dimension(:),allocatable :: lmode,mmode,nmode
  complex,dimension(:,:),allocatable :: pmodehis
  real,dimension(:),allocatable :: mdhis,mdhisa,mdhisb,mdhisc,mdhisd
  complex,dimension(:,:),allocatable :: aparhis,phihis

  !   kr, ktheta spectrum plots
  REAL,DIMENSION(:,:),allocatable :: phik

  !     weighty variables
  INTEGER,dimension(:),allocatable :: deljp,deljm
  INTEGER,dimension(:,:),allocatable :: jpl
  INTEGER,dimension(:,:),allocatable :: jpn
  INTEGER,dimension(:,:),allocatable :: jmi
  INTEGER,dimension(:,:),allocatable :: jmn
  REAL,DIMENSION(:),allocatable :: weightp,weightm
  REAL,DIMENSION(:),allocatable :: weightpn,weightmn

  !blending variable
  complex,dimension(:,:,:,:),allocatable :: pol,pmtrx,pmtrxi
  complex,dimension(:,:),allocatable :: pfac

  !      MPI variables
  !  include '/usr/include/mpif.h'

  integer,parameter :: Master=0
  integer :: numprocs
  INTEGER :: MyId,Last,cnt,ierr
  INTEGER :: GRID_COMM,TUBE_COMM
  INTEGER :: GCLR,TCLR,GLST,TLST
  INTEGER :: stat(MPI_STATUS_SIZE)
  INTEGER :: lngbr,rngbr,idprv,idnxt

  character(len=*) directory
  parameter(directory='./dump/')

  character(len=*) outdir
  parameter(outdir='./out/')

  !real :: ran2,revers
  !integer :: mod
  !real :: amod
  save

contains

  subroutine new_gem_com()
    nxpp = imx !/ntube
    allocate(workx(4*imx),worky(4*jmx),workz(4*kmx))
    allocate(tmpx(0:imx-1))
    allocate(tmpy(0:jmx-1))
    allocate(tmpz(0:kmx-1))

    allocate(rwx(nsmx,4),rwy(nsmx,4))
    allocate(mm(nsmx),tmm(nsmx),lr(nsmx))
    allocate(tets(nsmx),mims(nsmx),q(nsmx))
    allocate(kapn(nsmx),kapt(nsmx))
    allocate(time(0:nmx))
    allocate(yyamp(jmx,0:4),yyre(jmx,0:4),yyim(jmx,0:4),camp(0:6,0:50000),campf(0:6,0:nfreq-1))
    allocate(aparhis(0:6,0:jcnt-1),phihis(0:6,0:jcnt-1))
    allocate(mdhis(0:100),mdhisa(0:100),mdhisb(0:100))
    allocate(mdhisc(0:100),mdhisd(0:100))
    allocate(coefx(100+8*imx),coefy(100+8*jmx),coefz(100+8*kmx))

    ALLOCATE( den(nsmx,0:nxpp,0:jmx,0:1),dti(nsmx,0:nxpp,0:jmx,0:1), &
         delte(0:nxpp,0:jmx,0:1))
    ALLOCATE( rho(0:nxpp,0:jmx,0:1),drhoidt(0:nxpp,0:jmx,0:1), &
         jion(0:nxpp,0:jmx,0:1),jionx(0:nxpp,0:jmx,0:1), &
         jiony(0:nxpp,0:jmx,0:1))
    allocate( phi(0:nxpp,0:jmx,0:1))
    allocate( drhodt(0:nxpp,0:jmx,0:1),dnedt(0:nxpp,0:jmx,0:1))
    allocate( dnidt(nsmx,0:nxpp,0:jmx,0:1),jpar(nsmx,0:nxpp,0:jmx,0:1),  &
         jpex(nsmx,0:nxpp,0:jmx,0:1),jpey(nsmx,0:nxpp,0:jmx,0:1))
    allocate( dphidt(0:nxpp,0:jmx,0:1))
    ALLOCATE( ex(0:nxpp,0:jmx,0:1)) 
    ALLOCATE( ey(0:nxpp,0:jmx,0:1)) 
    ALLOCATE( ez(0:nxpp,0:jmx,0:1))
    ALLOCATE( dpdz(0:nxpp,0:jmx,0:1),dadz(0:nxpp,0:jmx,0:1))
    ALLOCATE( delbx(0:nxpp,0:jmx,0:1),delby(0:nxpp,0:jmx,0:1))
    ALLOCATE( xg(0:nxpp),yg(0:jmx),zg(0:1))
    allocate( apar(0:nxpp,0:jmx,0:1),dene(0:nxpp,0:jmx,0:1))
    allocate( upar(0:nxpp,0:jmx,0:1),upart(0:nxpp,0:jmx,0:1),upex(0:nxpp,0:jmx,0:1), &
         upey(0:nxpp,0:jmx,0:1),upa0(0:nxpp,0:jmx,0:1), &
         den0(0:nxpp,0:jmx,0:1),upazd(0:nxpp,0:jmx,0:1),&
         upa00(0:nxpp,0:jmx,0:1),upa0t(0:nxpp,0:jmx,0:1),den0apa(0:nxpp,0:jmx,0:1))
    allocate( cfx(0:nxpp,0:1),cfy(0:nxpp,0:1),jac(0:nxpp,0:1))
    allocate( bmag(0:nxpp,0:1),bdgxcgy(0:nxpp,0:1),bdgrzn(0:nxpp,0:1))
    allocate( ggxdgy(0:nxpp,0:1),ggy2(0:nxpp,0:1),ggx(0:nxpp,0:1))
    allocate (gn0e(0:nxpp),gt0e(0:nxpp),gt0i(0:nxpp),avap(0:nxpp))
    allocate (gn0s(1:5,0:nxpp))
    !          particle array declarations
    allocate( mu(nsmx,1:mmx),xii(nsmx,1:mmx),pzi(nsmx,1:mmx), &
         eki(nsmx,1:mmx),z0i(nsmx,1:mmx),u0i(nsmx,1:mmx))
    allocate( x2(nsmx,1:mmx),y2(nsmx,1:mmx),z2(nsmx,1:mmx),u2(nsmx,1:mmx))
    allocate( x3(nsmx,1:mmx),y3(nsmx,1:mmx),z3(nsmx,1:mmx),u3(nsmx,1:mmx))
    allocate( w2(nsmx,1:mmx),w3(nsmx,1:mmx))

    allocate( mue(1:mmxe),xie(1:mmxe),pze(1:mmxe),eke(1:mmxe),z0e(1:mmxe),u0e(1:mmxe))
    allocate( x2e(1:mmxe),y2e(1:mmxe),z2e(1:mmxe),u2e(1:mmxe),mue2(1:mmxe))
    allocate( x3e(1:mmxe),y3e(1:mmxe),z3e(1:mmxe),u3e(1:mmxe),mue3(1:mmxe))
    allocate( w2e(1:mmxe),w3e(1:mmxe))
    allocate( ipass(1:mmxe), index(1:mmxe))
    allocate(w000(1:mmxe),w001(1:mmxe),w010(1:mmxe),w011(1:mmxe),&
         w100(1:mmxe),w101(1:mmxe),w110(1:mmxe),w111(1:mmxe))

    !              Various diagnostic arrays and scalars
    !    plotting constants

    !    energy diagnostic arrays

    ALLOCATE( ke(nsmx,0:nmx),fe(0:nmx),te(0:nmx))
    ALLOCATE( rmsphi(0:nmx),rmsapa(0:nmx),avewi(1:3,0:nmx),avewe(0:nmx))
    ALLOCATE( nos(nsmx,0:nmx))

    !    flux diagnostics
    ALLOCATE( vol(1:nsubd),efle_es(1:nsubd,0:nmx),pfle_es(1:nsubd,0:nmx), &
         pfl_es(1:3,1:nsubd,0:nmx),efl_es(1:3,1:nsubd,0:nmx), &
         pfle_em(1:nsubd,0:nmx),efle_em(1:nsubd,0:nmx), &
         pfl_em(1:3,1:nsubd,0:nmx),efl_em(1:3,1:nsubd,0:nmx))
    ALLOCATE( chii(1:nsubd,0:nmx),chie(1:nsubd,0:nmx),ddi(1:nsubd,0:nmx), &
         achii(1:nsubd),achie(1:nsubd),addi(1:nsubd))

    !   mode diagnositics
    allocate(lmode(modemx),mmode(modemx),nmode(modemx))
    allocate(pmodehis(modemx,0:nmx))

    !   kr, ktheta spectrum plots
    ALLOCATE( phik(imx,jmx))

    !     weighty variables
    ALLOCATE( deljp(0:nxpp),deljm(0:nxpp))
    ALLOCATE( jpl(0:nxpp,0:jmx))
    ALLOCATE( jpn(0:nxpp,0:jmx))
    ALLOCATE( jmi(0:nxpp,0:jmx))
    ALLOCATE( jmn(0:nxpp,0:jmx))
    ALLOCATE( weightp(0:nxpp),weightm(0:nxpp))
    ALLOCATE( weightpn(0:nxpp),weightmn(0:nxpp))

    !Blending variable
    ALLOCATE(pol(1:nb,0:imx-1,0:jmx-1,0:kmx),pfac(0:imx-1,0:jmx-1), &
         pmtrx(0:imx-1,0:jmx-1,1:nb,1:nb), &
         pmtrxi(0:imx-1,0:jmx-1,1:nb,1:nb))

  end subroutine new_gem_com

end module gem_com










