pro nstx
directory='./out'

numx=401L
numz=601L
nxz=numx*numz
npts=intarr(nxz)
x1=fltarr(nxz)
z1=fltarr(nxz)
phixz=fltarr(nxz)

k=0L
i=0L
n2=10000L
n1=5000L
nb=2*n2+2
nr=n1+1
n=nr*(n2+1)
x=fltarr(n1)
z=fltarr(n1)
phi=fltarr(n1)
r=fltarr(nr)
phimr=fltarr(25,nr)
phimi=fltarr(25,nr)
phima=fltarr(25,nr)
rplt=fltarr(25,nr)
mpol=intarr(25)
rmax=fltarr(25)
ymax=fltarr(25)

xb=fltarr(nb)
zb=fltarr(nb)

infl=directory+'/pol1d.out'

openr, 1, infl
for k=0L,n1-1L do begin
readf, 1, dum1,dum2
x(k) = dum1
phi(k) = dum2
endfor
close,1

infl=directory+'/polxz.out'

openr, 1, infl
for k=0L,nxz-1L do begin
readf, 1, dum1,dum2,dum3,dum4
npts(k) = dum1
x1(k) = dum2
z1(k) = dum3
phixz(k) = dum4
endfor
close,1

infl=directory+'/mpol.out'
openr, 1, infl
for m = 0,24 do begin
    for k=0,nr-1 do begin
        readf, 1, i,dum1,dum2,dum3,dum4
        mpol(m) = i
        r(k) = dum1
        phima(m,k) = dum2
        phimr(m,k) = dum3
        phimi(m,k) = dum4
    endfor
endfor
close,1

for k=0,nr-1 do begin
    print,mpol(15),k,r(k),phimr(15,k)
endfor

infl=directory+'/bd.out'

openr, 1, infl
for k=0,nb-1 do begin
readf, 1, dum1,dum2
xb(k) = dum1
zb(k) = dum2
endfor
close,1

xmin=min(xb)
xmax=max(xb)
zmin=min(zb)
zmax=max(zb)

set_plot, 'ps'
loadct,32

ldsq=0.5
xsq=[ldsq,-ldsq,-ldsq,ldsq,ldsq]
ysq=[ldsq,ldsq,-ldsq,-ldsq,ldsq]

ldtr=1
xtr=[0,-1,1]
ytr=[1,-0.5,-0.5]

A=findgen(31)*!pi*2/30
xcir=1.0*cos(A)
ycir=1.0*sin(A)

a=[-100000]
b=[-100000]

!p.multi=[0,1,1,0]
flnm='pol1d.ps'
device, file=flnm, /encapsulated, xsize=21, ysize=21,/color
xc=(xmax+xmin)/2
phimax=max(phi)
phimin=min(phi)
plot, x, phi,xrange=[0.,1.],xtitle='R',ytitle='phi',xstyle=1,ystyle=1,position=[0.2,0.2,0.8,0.8]
device, /close


flnm='pol.ps'
device, file=flnm, /encapsulated, xsize=21, ysize=21,/color
xc=(xmax+xmin)/2
zspan=(zmax-zmin)
plot, a, b,psym=2,xrange=[xc-zspan/2,xc+zspan/2],yrange=[zmin,zmax],xtitle='R',ytitle='Z',xstyle=1,ystyle=1,position=[0.2,0.2,0.8,0.8]
phimax=max(phixz)
phimin=min(phixz)

ldsq=0.3
xsq1=[ldsq,-ldsq,-ldsq,ldsq,ldsq]
ysq1=[ldsq,ldsq,-ldsq,-ldsq,ldsq]

USERSYM, xsq1,ysq1, thick=1,/fill

for k=0L,nxz-1L do begin
  i=k
  c=(phixz(i)-phimin)/(phimax-phimin)*75
   if(abs(npts(k)>0))then begin
      oplot,x1(i:i),z1(i:i),psym=8,color=c
   endif
endfor
oplot,xb,zb,psym=3
device, /close

;Specify the red component of each color:
RED = [0, 1, 1, 0, 0, 1]
;Specify the green component of each color:
GREEN = [0, 1, 0, 1, 0, 1]
;Specify the blue component of each color:
BLUE = [0, 1, 0, 0, 1, 0]
;Load the first six elements of the color table:
TVLCT, 255 * RED, 255 * GREEN, 255 * BLUE

hong=2
lu=3
lan=4
huang=5

flnm='mpol.ps'
!p.multi=[0,1,1,0]
    !P.font=1
    !P.Charthick=2.0
    !P.Charsize=2.
    !X.Charsize=1.0
    !Y.Charsize=1.0
    !P.Thick=3.
    !X.Thick=5.0
    !Y.Thick=5.0

device, file=flnm, /encapsulated, xsize=21, ysize=21,/color

rplt = phima

xmax=max(rplt)
xmin=min(rplt)

print,xmin,xmax

for m = 0,24 do begin
    rmax(m) = 0.
    ymax(m) = 0.
    rdum = 0.
    ydum = 0.
    ydum1 = 0.
    for k=0,nr-1 do begin
        if(abs(rplt(m,k)) ge ydum) then begin
            rdum = r(k)
            ydum = abs(rplt(m,k))
            ydum1 = rplt(m,k)
        endif
    endfor
    rmax(m) = rdum
    ymax(m) = ydum1
    print,m,mpol(m),rmax(m),ymax(m)
endfor


m1 = 0
mskip = 3
plot, r(0:n1-1), rplt(m1,0:n1-1), yrange=[xmin,xmax*1.05],xtitle='r/a',linestyle=0,color=0,position=[0.2,0.2,0.8,0.8]
oplot,r(0:n1-1), rplt(m1,0:n1-1), linestyle=1,color=0

for m = 1,24 do begin
    icolor = ((mpol(m)-6) MOD 4)+2
    oplot,r(0:n1-1), rplt(m,0:n1-1), linestyle=0,color=0;icolor
    oplot,r(0:n1-1), rplt(m,0:n1-1), linestyle=1,color=0;icolor
    if(mpol(m) le 9)then begin
        xyouts,rmax(m),ymax(m),alignment=0.5,String(mpol(m),format='(I1)'),font=1,charthick=2,charsize=1.5
    endif

    if(mpol(m) ge 10 and mpol(m) le 100)then begin
        xyouts,rmax(m),ymax(m),alignment=0.5,String(mpol(m),format='(I2)'),font=1,charthick=2,charsize=1.5
    endif
endfor

device,/close

return
end
