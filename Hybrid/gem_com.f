module gem_com
!common data used for gem
USE pputil
implicit none

INTERFACE
  real(8) function revers(num,n)
  end function revers

  real(8) function ran2(i)
  end function ran2

  real(8) function en3(s)
      real(8) :: s
  end function en3
END INTERFACE

integer :: imx,jmx,kmx,mmx,mmxe,nmx,nsmx,nsubd=8,&
                     modemx,ntube=8,nxpp,ngdx=5,nb=16, &
                     negrd=16,nlgrd=24,nsnap=100

	 character*70 outname
	 REAL(8) :: endtm,begtm,pstm
	 REAL(8) :: starttm,lasttm,tottm
         real(8) :: aux1(50000),aux2(20000)
         real(8),dimension(:),allocatable :: workx,worky,workz
         COMPLEX(8),dimension(:),allocatable :: tmpx
         COMPLEX(8),dimension(:),allocatable :: tmpy
         COMPLEX(8),dimension(:),allocatable :: tmpz

!          imx,jmx,kmx = max no. of grid pts in x,y,z
!          mmx         = max no. of particles
!          nmx         = max. no. of time steps
!          nsmx        = max. no. of species (including tracer particles
!                                             as a species)
!          ntrmx       = max. no. of tracer particles
!          modemx      = max. no. of mode history plots

REAL(8), dimension(:,:),allocatable :: rwx,rwy
INTEGER,dimension(:),allocatable :: mm,tmm,ptr,lr
REAL(8),dimension(:),allocatable :: tets,mims,q
REAL(8),dimension(:),allocatable :: kapn, kapt
INTEGER :: timestep,im,jm,km,mykm,iseed,nrst,icmprs,iez,nfreq,mynf,ifskp,ignion,ignden,ignapa
INTEGER :: mtrace,mmt,mme,nzcrt,ntor0,tmme
real(8),dimension(:),allocatable :: time
REAL(8) :: dx,dy,dz,pi,pi2,dt,dte,totvol,n0,n0e,tcurr,rmpp,rmaa,eprs,dtref,tload,tloadc
REAL(8) :: epsnx,epsny,epsnz,epsax,epsay,epsaz,epspx,epspz,etaohm
REAL(8) :: lx,ly,lz,xshape,yshape,zshape,kapti,kapte,pzcrit(5),encrit(5)
INTEGER :: nm,nsm,kcnt,jcnt,mstart,ncurr,llk,mlk,onemd,izon,kxz,iflr,iorb,iflrh,islow,irlk,nom
REAL(8) :: cut,amp,tor,amie,emass,qel,isg,rneu,gamion,gamele,gambm,gamdne,gamphi,gamapa,gamte,nugam
REAL(8) :: c1,c2,c3,c4,fradi,kxcut,kycut,bcut,ftrap,frmax,omlk(0:1000),iomlk(0:1000),oml,omu,delom
INTEGER :: iput,iget,idg,kzlook,ision,ish,isham,ishgk,isbgk,iscgk,isbam,iscam,peritr,iadi,isgkm,itrace,icncl,ispre,isdte,iske,ifrzt,idenek,ieqmo814
REAL(8),DIMENSION(:,:),allocatable :: yyamp,yyre,yyim
complex(8),dimension(:,:),allocatable :: camp,aparhis,phihis
complex(8),dimension(:),allocatable :: campf
complex(8),dimension(:,:,:),allocatable :: campft
REAL(8) :: br0,lr0,qp,width,e0,vwidth,vwidthe,vcut,vpp,vt0,yd0
integer :: nonlin,nonline,iflut,nonlinh,ipara,iparah,nonlint,iparat,isuni,ifluid,ishift,nopz,nopi(5),noen,nowe
COMPLEX(8) :: IU
real(8),dimension(:),allocatable :: coefx,coefy,coefz
complex(8),dimension(1:8) :: apk,ptk,dpdtk
integer,dimension(1:8) :: lapa,mapa,napa
real(8) :: mrtio(0:1),aven,avptch
integer :: icrs_sec,ipg,isphi,iwtdgstc,isft
integer,dimension(0:255) :: isgnft,jft
REAL(8) :: delpzt,pzttot,delekt,ekttot

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

REAL(8),DIMENSION(:,:,:,:),allocatable :: den
REAL(8),DIMENSION(:,:,:),allocatable :: rho
real(8),dimension(:,:,:),allocatable :: phi
REAL(8),DIMENSION(:,:,:),allocatable :: ex
REAL(8),DIMENSION(:,:,:),allocatable :: ey
REAL(8),DIMENSION(:,:,:),allocatable :: ez,ezs
REAL(8),DIMENSION(:,:,:),allocatable :: dnidt,dncdt,dnhdt,dnbdt,dphidt,drhodt
REAL(8),DIMENSION(:,:,:),allocatable :: dpdz,dadz
REAL(8),DIMENSION(:,:,:),allocatable :: delbx,delby
REAL(8),DIMENSION(:),allocatable :: xg,yg,zg,dfltz
real(8),dimension(:,:,:),allocatable :: apar,dene
real(8),dimension(:,:,:),allocatable :: upar,djpa
real(8),dimension(:,:,:),allocatable :: phis,denes,apars,upars,deltes
real(8),dimension(:,:,:),allocatable :: dene_p,apar_p,ddedt,phi_p
real(8),dimension(:,:,:),allocatable :: dti,dtc,dne,delte
real(8),dimension(:,:,:),allocatable :: jpar,jpex,jpey,cpar,cpex,cpey
real(8),dimension(:,:,:),allocatable :: hpar,hden,denh,hden0,hpex,hpey
real(8),dimension(:,:,:),allocatable :: reyn,reynh,maxw,maxwe,maxwh
real(8),dimension(:,:,:),allocatable :: reynix,reyniy,maxwix,maxwiy
real(8),dimension(:,:,:),allocatable :: reynhx,reynhy,maxwhx,maxwhy
real(8),dimension(:,:,:),allocatable :: reynbx,reynby,maxwbx,maxwby
real(8),dimension(:,:,:),allocatable :: reyncx,reyncy,maxwcx,maxwcy
real(8),dimension(:),allocatable :: reyn00,maxw00,drdt00,reynh00,maxwe00,maxwh00
real(8),dimension(:,:,:),allocatable :: bpar,bden,denb,bpex,bpey
real(8),dimension(:,:,:),allocatable :: dltpe,dltpa,denek,upark
real(8),dimension(:,:),allocatable :: cfx,cfy,jac,bmag,bdgxcgy,bdgrzn,ggr,gnuobx,gnuoby,gupae0
real(8),dimension(:,:,:),allocatable :: dnedx,dnedy,dupadx,dupady
real(8),dimension(:),allocatable :: gn0i,gt0i,gn0e,gt0e,gcpne,gcpte,dbr,dtr,denhr,hden0r
complex(8),dimension(:,:,:,:),allocatable :: phiom,dteom
real(8),dimension(:,:,:),allocatable :: phiti,dteti
real(8),dimension(:,:,:),allocatable :: ppex,ppey,ppazd,djedt,djdt

!          particle array declarations
REAL(8),DIMENSION(:),allocatable :: mu,xii,pzi
REAL(8),DIMENSION(:),allocatable :: x2,y2,z2,u2,eki,z0i
REAL(8),DIMENSION(:),allocatable :: x3,y3,z3,u3
REAL(8),DIMENSION(:),allocatable :: w2,w3

REAL(8),DIMENSION(:),allocatable :: muh,muh2,xih,vih,pzh,mui
REAL(8),DIMENSION(:),allocatable :: xh2,yh2,zh2,uh2,ekh,z0h
REAL(8),DIMENSION(:),allocatable :: xh3,yh3,zh3,uh3
REAL(8),DIMENSION(:),allocatable :: wh2,wh3,index,isrl

REAL(8),DIMENSION(:),allocatable :: mub,mub2,ipass
REAL(8),DIMENSION(:),allocatable :: xb2,yb2,zb2,ub2
REAL(8),DIMENSION(:),allocatable :: xb3,yb3,zb3,ub3
REAL(8),DIMENSION(:),allocatable :: wb2,wb3

REAL(8),DIMENSION(:),allocatable :: muc
REAL(8),DIMENSION(:),allocatable :: xc2,yc2,zc2,uc2
REAL(8),DIMENSION(:),allocatable :: xc3,yc3,zc3,uc3
REAL(8),DIMENSION(:),allocatable :: wc2,wc3

REAL(8),DIMENSION(:),allocatable :: mue,xie,pze,eke,z0e,u0e
REAL(8),DIMENSION(:),allocatable :: x2e,y2e,z2e,u2e,mue2
REAL(8),DIMENSION(:),allocatable :: x3e,y3e,z3e,u3e,mue3
REAL(8),DIMENSION(:),allocatable :: w2e,w3e
REAL(8),DIMENSION(:),allocatable :: w000,w001,w010,w011,w100,w101,w110,w111


!              Various diagnostic arrays and scalars
!    plotting constants

INTEGER :: nplot,xnplt,pskip,kymode,imovie,npzh,nenh
REAL(8) :: contu,wmax

!    energy diagnostic arrays

REAL(8),DIMENSION(:,:),allocatable :: ke
REAL(8),DIMENSION(:),allocatable :: fe,te
REAL(8),DIMENSION(:),allocatable :: rmsphi,rmsez,rmsapa,avewi,avewc,avewe,avewh,avewt,avewb
REAL(8),DIMENSION(:,:),allocatable :: nos

!    flux diagnostics
REAL(8),DIMENSION(:),allocatable :: vol
REAL(8),DIMENSION(:,:),allocatable :: efle,pfle
REAL(8),DIMENSION(:,:),allocatable :: pfl,efl

!   mode diagnositics
INTEGER :: modem
INTEGER,dimension(:),allocatable :: lmode,mmode,nmode,ntor
real(8),dimension(:),allocatable :: mdhis,mdhisa,mdhisb,mdhisc,mdhisd

COMPLEX(8),dimension(:,:),allocatable :: pmodehis

!   kr, ktheta spectrum plots
REAL(8),DIMENSION(:,:),allocatable :: phik

!	weighty variables
INTEGER,dimension(:),allocatable :: deljp,deljm
INTEGER,dimension(:,:),allocatable :: jpl
INTEGER,dimension(:,:),allocatable :: jpn
INTEGER,dimension(:,:),allocatable :: jmi
INTEGER,dimension(:,:),allocatable :: jmn
REAL(8),DIMENSION(:),allocatable :: weightp,weightm
REAL(8),DIMENSION(:),allocatable :: weightpn,weightmn

!blending variable
      complex(8),dimension(:,:,:,:),allocatable :: pol,pmtrx,pmtrxi
      complex(8),dimension(:,:),allocatable :: pfac

! 	MPI variables
!  include '/usr/include/mpif.h'
include 'mpif.h'
integer,parameter :: Master=0
integer :: numprocs
INTEGER :: MyId,Last,cnt,ierr
INTEGER :: GRID_COMM,TUBE_COMM
INTEGER :: GCLR,TCLR,GLST,TLST
INTEGER :: stat(MPI_STATUS_SIZE)
INTEGER :: lngbr,rngbr,idprv,idnxt

character * (*) directory
parameter(directory='./dump/')

character * (*) outdir
parameter(outdir='./out/')

!real(8) :: ran2,revers
!integer :: mod
!real(8) :: amod
save

contains
subroutine new_gem_com()
      nxpp = imx !/ntube
      allocate(workx(4*imx),worky(4*jmx),workz(4*kmx))
      allocate(tmpx(0:imx-1))
      allocate(tmpy(0:jmx-1))
      allocate(tmpz(0:kmx-1))

      allocate(rwx(nsmx,4),rwy(nsmx,4))
allocate(mm(nsmx),tmm(nsmx),ptr(nsmx),lr(nsmx))
allocate(tets(nsmx),mims(nsmx),q(nsmx))
allocate(kapn(nsmx),kapt(nsmx))
allocate(time(0:nmx))
allocate(yyamp(jmx,0:4),yyre(jmx,0:4),yyim(jmx,0:4),camp(0:6,0:50000),campf(0:nfreq-1),campft(0:6,0:nsnap-1,0:nfreq-1))
allocate(aparhis(0:6,0:jcnt-1),phihis(0:6,0:jcnt-1))
allocate(coefx(100+8*imx),coefy(100+8*jmx),coefz(100+8*kmx))

ALLOCATE( den(2,0:nxpp,0:jmx,0:1),dti(0:nxpp,0:jmx,0:1),dtc(0:nxpp,0:jmx,0:1), dne(0:nxpp,0:jmx,0:1), &
          delte(0:nxpp,0:jmx,0:1))
ALLOCATE( rho(0:nxpp,0:jmx,0:1))
allocate( phi(0:nxpp,0:jmx,0:1))
ALLOCATE( ex(0:nxpp,0:jmx,0:1)) 
ALLOCATE( ey(0:nxpp,0:jmx,0:1)) 
ALLOCATE( ez(0:nxpp,0:jmx,0:1),ezs(0:nxpp,0:jmx,0:1))
ALLOCATE( dpdz(0:nxpp,0:jmx,0:1),dadz(0:nxpp,0:jmx,0:1))
ALLOCATE( delbx(0:nxpp,0:jmx,0:1),delby(0:nxpp,0:jmx,0:1))
ALLOCATE( xg(0:nxpp),yg(0:jmx),zg(0:1),dfltz(0:1))
allocate( apar(0:nxpp,0:jmx,0:1),dene(0:nxpp,0:jmx,0:1))
allocate( ddedt(0:nxpp,0:jmx,0:1))
allocate( apar_p(0:nxpp,0:jmx,0:1),dene_p(0:nxpp,0:jmx,0:1),phi_p(0:nxpp,0:jmx,0:1))
allocate( upar(0:nxpp,0:jmx,0:1),djpa(0:nxpp,0:jmx,0:1),jpar(0:nxpp,0:jmx,0:1),denh(0:nxpp,0:jmx,0:1), &
          hpar(0:nxpp,0:jmx,0:1),hden(0:nxpp,0:jmx,0:1),hden0(0:nxpp,0:jmx,0:1))
allocate( upars(0:nxpp,0:jmx,0:1),phis(0:nxpp,0:jmx,0:1),&
          denes(0:nxpp,0:jmx,0:1),apars(0:nxpp,0:jmx,0:1),deltes(0:nxpp,0:jmx,0:1))
allocate( ppex(0:nxpp,0:jmx,0:1),ppey(0:nxpp,0:jmx,0:1),ppazd(0:nxpp,0:jmx,0:1), &
          djedt(0:nxpp,0:jmx,0:1),djdt(0:nxpp,0:jmx,0:1))
allocate( dltpe(0:nxpp,0:jmx,0:1),dltpa(0:nxpp,0:jmx,0:1),denek(0:nxpp,0:jmx,0:1),upark(0:nxpp,0:jmx,0:1))
allocate( cfx(0:nxpp,0:1),cfy(0:nxpp,0:1),jac(0:nxpp,0:1))
allocate( bmag(0:nxpp,0:1),bdgxcgy(0:nxpp,0:1),bdgrzn(0:nxpp,0:1))
allocate( dnedx(0:nxpp,0:jmx,0:1),dnedy(0:nxpp,0:jmx,0:1), &
          dupadx(0:nxpp,0:jmx,0:1),dupady(0:nxpp,0:jmx,0:1))
allocate(gn0i(0:nxpp),gn0e(0:nxpp),gt0i(0:nxpp),gt0e(0:nxpp),gcpne(0:nxpp),gcpte(0:nxpp)) 
allocate(dbr(0:nxpp),dtr(0:nxpp),ggr(0:nxpp,0:1),gnuobx(0:nxpp,0:1),gnuoby(0:nxpp,0:1),gupae0(0:nxpp,0:1)) 
allocate(denhr(0:nxpp),hden0r(0:nxpp))
allocate( phiom(0:nxpp,0:1,0:1,0:nom-1),phiti(0:nxpp,0:jmx,0:1),dteom(0:nxpp,0:1,0:1,0:nom-1),dteti(0:nxpp,0:jmx,0:1))
allocate( dnidt(0:nxpp,0:jmx,0:1),dncdt(0:nxpp,0:jmx,0:1),dnhdt(0:nxpp,0:jmx,0:1),dnbdt(0:nxpp,0:jmx,0:1), &
          dphidt(0:nxpp,0:jmx,0:1),drhodt(0:nxpp,0:jmx,0:1), cpar(0:nxpp,0:jmx,0:1),cpex(0:nxpp,0:jmx,0:1),cpey(0:nxpp,0:jmx,0:1), &
          jpex(0:nxpp,0:jmx,0:1),jpey(0:nxpp,0:jmx,0:1),hpex(0:nxpp,0:jmx,0:1),hpey(0:nxpp,0:jmx,0:1))
allocate( bpar(0:nxpp,0:jmx,0:1),bpex(0:nxpp,0:jmx,0:1),bpey(0:nxpp,0:jmx,0:1),bden(0:nxpp,0:jmx,0:1),denb(0:nxpp,0:jmx,0:1))
allocate( reynix(0:nxpp,0:jmx,0:1),reyniy(0:nxpp,0:jmx,0:1),maxwix(0:nxpp,0:jmx,0:1),maxwiy(0:nxpp,0:jmx,0:1))
allocate( reyn(0:nxpp,0:jmx,0:1),reynh(0:nxpp,0:jmx,0:1),maxw(0:nxpp,0:jmx,0:1),maxwe(0:nxpp,0:jmx,0:1),maxwh(0:nxpp,0:jmx,0:1))
allocate( reyn00(0:nxpp-1),maxw00(0:nxpp-1),drdt00(0:nxpp-1),reynh00(0:nxpp-1),maxwe00(0:nxpp-1),maxwh00(0:nxpp-1))

!          particle array declarations
allocate( mu(1:mmx),xii(1:mmx))
allocate( x2(1:mmx),y2(1:mmx),z2(1:mmx),u2(1:mmx),pzi(1:mmx),eki(1:mmx))
allocate( x3(1:mmx),y3(1:mmx),z3(1:mmx),u3(1:mmx),z0i(1:mmx))
allocate( w2(1:mmx),w3(1:mmx))

if(ishgk==1)then
allocate( muh(1:mmx),muh2(1:mmx),xih(1:mmx),vih(1:mmx),mui(1:mmx))
allocate( xh2(1:mmx),yh2(1:mmx),zh2(1:mmx),uh2(1:mmx),pzh(1:mmx),ekh(1:mmx))
allocate( xh3(1:mmx),yh3(1:mmx),zh3(1:mmx),uh3(1:mmx),z0h(1:mmx))
allocate( wh2(1:mmx),wh3(1:mmx),index(1:mmx),isrl(1:mmx))
allocate( reynhx(0:nxpp,0:jmx,0:1),reynhy(0:nxpp,0:jmx,0:1),maxwhx(0:nxpp,0:jmx,0:1),maxwhy(0:nxpp,0:jmx,0:1))
end if

if(isbgk==1)then
allocate( mub(1:mmx),mub2(1:mmx),ipass(1:mmx))
allocate( xb2(1:mmx),yb2(1:mmx),zb2(1:mmx),ub2(1:mmx))
allocate( xb3(1:mmx),yb3(1:mmx),zb3(1:mmx),ub3(1:mmx))
allocate( wb2(1:mmx),wb3(1:mmx))
allocate( reynbx(0:nxpp,0:jmx,0:1),reynby(0:nxpp,0:jmx,0:1),maxwbx(0:nxpp,0:jmx,0:1),maxwby(0:nxpp,0:jmx,0:1))
end if

if(iscgk==1)then
allocate( muc(1:mmx))
allocate( xc2(1:mmx),yc2(1:mmx),zc2(1:mmx),uc2(1:mmx))
allocate( xc3(1:mmx),yc3(1:mmx),zc3(1:mmx),uc3(1:mmx))
allocate( wc2(1:mmx),wc3(1:mmx))
allocate( reyncx(0:nxpp,0:jmx,0:1),reyncy(0:nxpp,0:jmx,0:1),maxwcx(0:nxpp,0:jmx,0:1),maxwcy(0:nxpp,0:jmx,0:1))
end if

!allocate( mut(1:mmx),mut2(1:mmx),xit(1:mmx),vit(1:mmx))
!allocate( xt2(1:mmx),yt2(1:mmx),zt2(1:mmx),ut2(1:mmx),pzt(1:mmx),pzti(1:mmx),ekt(1:mmx),ekti(1:mmx))
!allocate( xt3(1:mmx),yt3(1:mmx),zt3(1:mmx),ut3(1:mmx),z0t(1:mmx))
!allocate( wt2(1:mmx),wt3(1:mmx))

if(iske==1)then
allocate( mue(1:mmxe),xie(1:mmxe),pze(1:mmxe),eke(1:mmxe),z0e(1:mmxe),u0e(1:mmxe))
allocate( x2e(1:mmxe),y2e(1:mmxe),z2e(1:mmxe),u2e(1:mmxe),mue2(1:mmxe))
allocate( x3e(1:mmxe),y3e(1:mmxe),z3e(1:mmxe),u3e(1:mmxe),mue3(1:mmxe))
allocate( w2e(1:mmxe),w3e(1:mmxe))
allocate(w000(1:mmxe),w001(1:mmxe),w010(1:mmxe),w011(1:mmxe),&
         w100(1:mmxe),w101(1:mmxe),w110(1:mmxe),w111(1:mmxe))
end if
!              Various diagnostic arrays and scalars
!    plotting constants

!    energy diagnostic arrays

ALLOCATE( ke(nsmx,0:nmx),fe(0:nmx),te(0:nmx))
ALLOCATE( rmsphi(0:nmx),rmsez(0:nmx),rmsapa(0:nmx),avewi(0:nmx),avewc(0:nmx),avewe(0:nmx),avewh(0:nmx),avewt(0:nmx),avewb(0:nmx))
ALLOCATE( nos(nsmx,0:nmx))

!    flux diagnostics
ALLOCATE(vol(1:nsubd),efle(1:nsubd,0:nmx),pfle(1:nsubd,0:nmx), &
         pfl(nsmx+1,0:nmx),efl(nsmx,0:nmx))

!   mode diagnositics
allocate(lmode(modemx),mmode(modemx),nmode(modemx),ntor(0:jcnt-1))
allocate(pmodehis(modemx,0:nmx))
allocate(mdhis(0:100),mdhisa(0:100),mdhisb(0:100),mdhisc(0:100),mdhisd(0:100))

!   kr, ktheta spectrum plots
ALLOCATE( phik(imx,jmx))

!	weighty variables
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










