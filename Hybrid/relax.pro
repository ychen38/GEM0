pro relax
directory='./out'
nh = 2.0779e-2
cv = 7.8348
ncolumn=2
nrow=3
npgs = 10
nplts = ncolumn*nrow
nskip = 1
nbegin = (npgs-1)-nskip*(nplts-1)

;print,'npgs='
;read,npgs
rin=0.1
rout=0.9

imx=128L
s=fltarr(imx+1)
n0=fltarr(npgs,imx+1)
dn=fltarr(npgs,imx+1)

for k = 0,imx do begin
   s(k)=rin+(rout-rin)/float(imx)*float(k)
endfor

dumstring=''
infl=directory+'/denhr.out'
openr, 1, infl
readf,1,x1,x2,x3,x4
for i=0,npgs-1 do begin
    readf, 1,dumstring
    for k=0L,imx do begin
        readf, 1, dum
        dn(i,k) = dum
    endfor
endfor
close,1
infl=directory+'/hden0.out'
openr, 1, infl
readf,1,x1,x2,x3,x4
for i=0,npgs-1 do begin
    readf, 1,dumstring
    for k=0L,imx do begin
        readf, 1, dum
        n0(i,k) = dum
    endfor
endfor
close,1

nhxpp = fltarr(201)
sxpp = fltarr(201)
for k = 0,200 do begin
   sxpp(k)=rin+(rout-rin)/200.0*float(k)
endfor
infl='./xpp'
openr, 1, infl
readf,1,dumstring
for i=0,200 do begin
   readf, 1, x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11
   nhxpp(i) = x11
endfor
nhxpp = nhxpp*nh*cv
close,1

;set n0xpp on imx grids
n0xpp = fltarr(imx+1)
dr = (rout-rin)/200
dx = (rout-rin)/imx
for i = 0,imx do begin
   j = long((i*dx)/dr)
   n0xpp(i) = nhxpp(j)
endfor

set_plot, 'ps'
loadct,32
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

!p.multi=[0,ncolumn,nrow,0]
    !P.Charthick=5.0
    !P.Charsize=2.
    !X.Charsize=1.0
    !Y.Charsize=1.0
    !P.Thick=5.
    !X.Thick=5.0
    !Y.Thick=5.0
flnm='relax.ps'
device, file=flnm, /encapsulated, xsize=21, ysize=28,/color

phimax=max(n0xpp)
phimin=min(n0xpp)

idum=0
for k = 0,nplts-1 do begin
   k1 = nbegin+nskip*k
   plot, s(idum:imx-idum), n0xpp(idum:imx-idum),xrange=[0.,1.],xtitle='r/a',ytitle='nbeam',title=k1,xstyle=1,ystyle=1 ;,position=[0.2,0.2,0.8,0.8]
;   oplot, s(idum:imx-idum), n0(k1,idum:imx-idum)+dn(k1,idum:imx-idum),color=hong ;,position=[0.2,0.2,0.8,0.8]
   oplot, s(idum:imx-idum), n0xpp(idum:imx-idum)+dn(k1,idum:imx-idum),color=hong ;,position=[0.2,0.2,0.8,0.8]
endfor
device, /close

return
end
