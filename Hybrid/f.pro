pro f
npgs = 40
nplts = 12
nskip = 3
nbegin = (npgs-1)-nskip*(nplts-1)

k=0L
i=0L
j=0L
n=0L
nbinx=50
nbiny=100
NXP=nbinx
NYP=nbiny
dfpz=fltarr(nbinx,nbiny,npgs)
lmb=fltarr(nbinx,nbiny,npgs)
f0pz=fltarr(nbinx,nbiny,npgs)
f0v=fltarr(nbinx,nbiny,npgs)
dfv=fltarr(nbinx,nbiny,npgs)
numpz=lonarr(nbinx,nbiny,npgs)
zplt = fltarr(nbinx,nbiny,npgs)
z = fltarr(nxp,nyp)
y=fltarr(NYP,npgs)
x=fltarr(NXP)

infl='wtdfpz'

openr, 1, infl
for k=0,npgs-1 do begin
   readf,1,x1,x2,x3,x4
   xmin = x1
   xmax = x2
   pzmin = x3
   pzmax = x4

   for j = 0,NXP-1 do begin
      x(j) = xmin+(j+0.5)*(xmax-xmin)/NXP
   endfor

   for j = 0,NYP-1 do begin
      y(j,k) = pzmin+(j+0.5)*(pzmax-pzmin)/NYP
   endfor

   for j = 0,nbiny-1 do begin
      for i = 0,nbinx-1 do begin
         readf, 1, x0,x1,x2,x3,x4,x5,x6  ;,x2,x3,x4,x5,x6,x7,x8,x9
         i1 = x0
         j1 = x1
         numpz(i1,j1,k) = x2
         dfpz(i1,j1,k) = x3
         lmb(i1,j1,k) = x4
         f0v(i1,j1,k) = x5
         dfv(i1,j1,k) = x6
         zplt(i1,j1,k) = x3
         if(numpz(i1,j1) lt 1)then begin
;           zplt(i1,j1)=-1.0
         endif
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

!p.multi = [0,3,4,0,0]
set_plot, 'ps'
flnm='contour.ps'
device, file=flnm, /color, /encapsulated, $
      bits=8, xsize=21, ysize=27

for k = 0,nplts-1 do begin
   k1 = nbegin+nskip*k
   z(0:nxp-1,0:nyp-1)=zplt(0:NXP-1,0:NYP-1,k1)
   nx = nxp
   ny = nyp
;z=min_curve_surf(zplt,NX=NXP,NY=NYP)

   rmax=max(z)
   rmin=min(z)
   nlvl=15
   cskip=5
   levs=dindgen(nlvl)
   rran=rmax-rmin
   levs=rran/nlvl*levs
   levs=levs+rmin

   loadct, 32
   contour, z, x,y(0:NYP-1,k1),levels=levs, xrange=[xmin,xmax],yrange=[y(0,k1),y(nyp-1,k1)],$
         xstyle=1, ystyle=1, charsize=.8, xtitle='x',ytitle='y',$
         /fill, c_colors=indgen(nlvl)*cskip,title='';,position=[0.2,0.2,0.8,0.8]

;   contour, z, x,y(0:NYP-1,k1),levels=levs, $
;         c_labels=labs, xstyle=1, ystyle=1, $
;         c_linestyle=(levs lt 0.0), /overplot

   RED = [0, 1, 1, 0, 0, 1]
   GREEN = [0, 1, 0, 1, 0, 1]
   BLUE = [0, 1, 0, 0, 1, 0]
   TVLCT, 255 * RED, 255 * GREEN, 255 * BLUE
   hong=2
   lu=3
   lan=4
   huang=5

   USERSYM, xsq,ysq, thick=1,/fill
   for i=0,nbinx-1 do begin
      for j = 0,nbiny-1 do begin
         if (numpz(i,j,k1) le 1)then begin
            oplot,x(i:i),y(j:j,k1),psym=3,color=hong
         endif
         if (lmb(i,j,k1) lt 0)then begin
;            oplot,x(i:i),y(j:j,k1),psym=3,color=huang
         endif
      endfor
   endfor
endfor  
device, /close
loadct,0

!p.multi = [0,1,2,0,0]
set_plot, 'ps'
flnm='contour1.ps'
device, file=flnm, /color, /encapsulated, $
      bits=8, xsize=21, ysize=21

k1 = npgs-1 
z(0:nxp-1,0:nyp-1)=zplt(0:NXP-1,0:NYP-1,k1)
nx = nxp
ny = nyp
;z=min_curve_surf(zplt,NX=NXP,NY=NYP)

rmax=max(z)
rmin=min(z)
nlvl=15
cskip=5
levs=dindgen(nlvl)
rran=rmax-rmin
levs=rran/nlvl*levs
levs=levs+rmin

loadct, 32
contour, z, x,y(0:NYP-1,k1),levels=levs, xrange=[xmin,xmax],yrange=[y(0,k1),y(nyp-1,k1)],$
         xstyle=1, ystyle=1, charsize=.8, xtitle='x',ytitle='y',$
         /fill, c_colors=indgen(nlvl)*cskip,title='',position=[0.2,0.2,0.8,0.8]

contour, z, x,y(0:NYP-1,k1),levels=levs, $
         c_labels=labs, xstyle=1, ystyle=1, $
         c_linestyle=(levs lt 0.0), /overplot

RED = [0, 1, 1, 0, 0, 1]
GREEN = [0, 1, 0, 1, 0, 1]
BLUE = [0, 1, 0, 0, 1, 0]
TVLCT, 255 * RED, 255 * GREEN, 255 * BLUE
hong=2
lu=3
lan=4
huang=5

USERSYM, xsq,ysq, thick=1,/fill
for i=0,nbinx-1 do begin
   for j = 0,nbiny-1 do begin
      if (numpz(i,j,k1) le 1)then begin
         oplot,x(i:i),y(j:j,k1),psym=3,color=hong
      endif
   endfor
endfor

loadct,32
nlg = 50
lg=fltarr(nlg,nlg)
for i = 0,nlg-1 do begin
   for j = 0,nlg-1 do begin
      lg(i,j) = min(z)+(max(z)-min(z))/nlg*i
   endfor
endfor

contour, lg, levels=levs, xticks=0,yticks=0, $
         xtickname=replicate(' ',50), $
         ytickname=replicate(' ',50), $
         c_labels=labs, xstyle=4, ystyle=4, $ ;xstyle=4 suppresses entire axis
         /fill, c_colors=indgen(nlvl)*cskip,title='',position=[0.3,0.85,0.7,0.9]

xyouts,0,nlg+1.2,alignment=0.5,min(z)     ;used min(dfpz) first to find the minimum
xyouts,nlg,nlg+1.2,alignment=0.5,max(z)   ;use max(dfpz) first to find value

device, /close
loadct, 0

return
end
