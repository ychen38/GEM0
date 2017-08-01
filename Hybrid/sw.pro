pro sw
npgs = 7
nplts = 7
nskip = 1
mynf = 4
nbegin = (npgs-1)-nskip*(nplts-1)

k=0L
i=0L
j=0L
n=0L
nbinx=100
nbiny=32*mynf
NXP=nbinx
NYP=nbiny
jbgn = nbiny/2-NYP/2
jend = jbgn+NYP-1
dfpz=fltarr(nbinx,nbiny,npgs)
zplt = fltarr(nbinx,nbiny,npgs)
ytmp=fltarr(Nbiny)
z = fltarr(nxp,nyp)
y=fltarr(NYP)
x=fltarr(NXP)

infl='sweep'

xmin = 0.
xmax = 1.
ymin = -0.01
ymax = 0.01
for j = 0,NXP-1 do begin
   x(j) = xmin+(j+0.5)*(xmax-xmin)/NXP
endfor

for j = 0,NYP-1 do begin
   y(j) = ymin+(jbgn+j)*(ymax-ymin)/nbiny   ;+ymin+(j+0.5)*(ymax-ymin)/NYP
endfor

openr, 1, infl
for k=0,npgs-1 do begin
   for j = 0,nbiny-1 do begin
      for i = 0,nbinx-1 do begin
         readf, 1, x0,x1,x2,x3,x4,x5
         i1 = x0
         j1 = x1
         zplt(i1,j1,k) = x5
      endfor
   endfor
   for i = 0,nbinx-1 do begin
      for j = 0,nbiny-1 do begin
         ytmp(j) = zplt(i,j,k)
      endfor
      xdum = max(ytmp)
      for j = 0,nbiny-1 do begin
         zplt(i,j,k) = zplt(i,j,k)/xdum
      endfor
   endfor
endfor
close,1

set_plot, 'ps'
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

;loadct,32

ldsq=0.6
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


loadct, 32
print,!D.table_size

!p.multi = [0,2,4,0,0]
    !P.Charthick=1.0
    !P.Charsize=2.
    !X.Charsize=2.0
    !Y.Charsize=2.5
    !P.Thick=1.
    !X.Thick=1.0
    !Y.Thick=1.0

set_plot, 'ps'
flnm='sweep.ps'
device, file=flnm, /color, /encapsulated, $
      bits=8, xsize=21, ysize=27

for k = 0,nplts-1 do begin
   k1 = nbegin+nskip*k
   z(0:nxp-1,0:nyp-1)=zplt(0:NXP-1,jbgn:jend,k1)
   nx = nxp
   ny = nyp
;   z(0:nxp-1,0:nyp-1)=min_curve_surf(zplt(0:nbinx-1,0:nbiny-1,k1),NX=NXP,NY=NYP)

   rmax=max(z)
   rmin=min(z)
   nlvl=15
   cskip=5
   levs=dindgen(nlvl)
   rran=rmax-rmin
   levs=rran/nlvl*levs
   levs=levs+rmin

   loadct, 32
   contour, z, x,y(0:NYP-1),levels=levs, xrange=[xmin,xmax],yrange=[y(0),y(nyp-1)],$
         xstyle=1, ystyle=1, charsize=.8, xtitle='x',ytitle='y',$
         /fill, c_colors=indgen(nlvl)*cskip,title='';,position=[0.2,0.2,0.8,0.8]

;   contour, z, x,y(0:NYP-1),levels=levs, $
;         c_labels=labs, xstyle=1, ystyle=1, $
;         c_linestyle=(levs lt 0.0), /overplot

endfor  
device, /close
loadct,0

return
end
