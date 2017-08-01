pro sigmar
iplt = 4
n0=2000L
iphi=0
acs=546
omeg = 1./sqrt(0.00124)/sqrt(2)/(2.*3.33*1048.2)
print,'omeg= ', omeg
k=1L
rin = 0.1
rout = 0.9
;collision=0
infl='plot1'
t0 = fltarr(n0)
phi0 = fltarr(n0)

if iphi eq 1 then begin
openr, 1, infl
for k=0L,n0-1L do begin
readf,1, x1,x2,x3,x4,x5
t0(k) = x1/acs
phi0(k) = x2
endfor
close, 1
endif

infl='flux'
;t0 = fltarr(n0)
flx = fltarr(12,n0)

openr, 1, infl
for k=0L,n0-1L do begin
readf,1, x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12
t0(k) = x0/acs
flx(0,k) = x1
flx(1,k) = x2
flx(2,k) = x3
flx(3,k) = x4
flx(4,k) = x5
flx(5,k) = x6
flx(6,k) = x7
flx(7,k) = x8
flx(8,k) = x9
flx(9,k) = x10
flx(10,k) = x11
flx(11,k) = x12
endfor

close, 1

;average to smooth
nav = 2
phi0a = fltarr(n0)
flxa = fltarr(n0)

for k = 1,nav-1 do begin
    phi0a(k) = phi0(k)
    flxa(k) = flx(iplt,k)
endfor

for k = nav,n0-nav/2-1 do begin
   dum0=0
   for i = 1, nav do begin
       dum0 = dum0+phi0(k-nav/2+i)
   endfor
   phi0a(k) = dum0/nav

   dum0=0
   for i = 1, nav do begin
       dum0 = dum0+flx(iplt,k-nav/2+i)
   endfor
   flxa(k) = dum0/nav
endfor

for k = n0-nav/2,n0-1 do begin
    phi0a(k) = phi0a(n0-nav/2-1)
    flxa(k) = flxa(n0-nav/2-1)
endfor

for i = 0,11 do begin
   dum0=0
   nst = 1000L
   for k=nst,n0-1L do begin
      dum0 = dum0+flx(i,k)
   endfor
   dum0=dum0/(n0-nst)
   print,i,dum0
endfor
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

flnm='rmsphit.ps'
device, file=flnm, /color,/encapsulated;, xsize=21, ysize=21

if iphi eq 1 then begin
xmax = max([max(t0)])
ymax = max([max(phi0(0:n0-1))])
plot, t0(0:n0-nav/2),phi0a(0:n0-nav/2), title='', $
  xstyle=1, ystyle=1,xtitle='t (ms)',$
  ytitle='rms(phi)',$ ;'e<!7u!3>/T',$
  xrange=[0,xmax*1.2], $
  yrange=[0,ymax*1.01],linestyle=0,$
  charsize=2,charthick=4,color=0,position=[0.3,0.3,0.9,0.9]

device, /close
endif

flnm='flx.ps'
device, file=flnm, /color,/encapsulated;, xsize=21, ysize=21

xmax = max([max(t0)])
ymax = max([max(flxa(0:n0-1))])
plot, t0(0:n0-nav/2),flxa(0:n0-nav/2), title='', $
  xstyle=1, ystyle=1,xtitle='t (ms)',$
  ytitle='delta T_e',$ ;'e<!7u!3>/T',$
  xrange=[0,xmax*1.2], $
  yrange=[0,ymax*1.01],linestyle=0,$
  charsize=2,charthick=4,color=0,position=[0.3,0.3,0.9,0.9]

device, /close

return
end
