pro gtc

k=1L

ngtc=200L
rgtc=fltarr(ngtc)
zgtc=fltarr(ngtc)
ngem=200L
rgem=fltarr(ngem)
zgem=fltarr(ngem)

rgtc85=fltarr(ngtc)
zgtc85=fltarr(ngtc)
rgem85=fltarr(ngem)
zgem85=fltarr(ngem)

rgtc2=fltarr(ngtc)
zgtc2=fltarr(ngtc)
rgem2=fltarr(ngem)
zgem2=fltarr(ngem)

infl='rho5.dat'

openr, 1, infl
for k=0L,ngtc-1L do begin
readf,1, x1,x2,x3
rgtc(k) = x2
zgtc(k) = x3
endfor
for k=0L,ngem-1L do begin
readf,1, x1,x2,x3
rgem(k) = x2
zgem(k) = x3
endfor

close, 1

infl='rho85.dat'

openr, 1, infl
for k=0L,ngtc-1L do begin
readf,1, x1,x2,x3
rgtc85(k) = x2
zgtc85(k) = x3
endfor
for k=0L,ngem-1L do begin
readf,1, x1,x2,x3
rgem85(k) = x2
zgem85(k) = x3
endfor

close, 1

infl='rho2.dat'

openr, 1, infl
for k=0L,ngtc-1L do begin
readf,1, x1,x2,x3
rgtc2(k) = x2
zgtc2(k) = x3
endfor
for k=0L,ngem-1L do begin
readf,1, x1,x2,x3
rgem2(k) = x2
zgem2(k) = x3
endfor

close, 1

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

!p.multi=[0,1,1,0]
    !P.Charthick=5.0
    !P.Charsize=2.
    !X.Charsize=1.0
    !Y.Charsize=1.0
    !P.Thick=5.
    !X.Thick=5.0
    !Y.Thick=5.0

flnm='gtc.ps'
device, file=flnm, /color,/encapsulated, xsize=21, ysize=21

xmax = max([max(rgtc),max(rgtc85)])
xmin = min([min(rgtc),min(rgtc85)])
xc=(xmax+xmin)/2
ymax = max([max(zgtc),max(zgtc85)])
ymin = min([min(zgtc),min(zgtc85)])
yspan=ymax-ymin

plot, rgtc,zgtc, title='', $
  xstyle=1, ystyle=1,xtitle=' ',$
  ytitle=' ',$ ;'e<!7u!3>/T',$
  xrange=[xc-yspan/2,xc+yspan/2], $
  yrange=[ymin,ymax],linestyle=0,$
  charsize=2,charthick=4,color=0,position=[0.3,0.3,0.9,0.9]
oplot,rgem,zgem,color=hong

oplot,rgtc85,zgtc85,color=0
oplot,rgem85,zgem85,color=hong

oplot,rgtc2,zgtc2,color=0
oplot,rgem2,zgem2,color=hong

device, /close

return
end
