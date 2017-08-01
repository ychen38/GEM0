pro spectrum

numflt=20
k=1L
nfreq=100000L
jmx = 32
f=fltarr(nfreq)
pw=dblarr(nfreq)
pwa=dblarr(nfreq)
y = fltarr(nfreq)
t = fltarr(nfreq)

infl='spectrum'
dumchar=' '
openr, 1, infl
z=0.d
for k=0,nfreq-1 do begin
   readf,1,x1,x2,x3
   f(k) = x2
   pw(k) = x3
endfor
close, 1

infl='yyre'
dumchar=' '
openr, 1, infl
z=0.d
for k=0,nfreq-1 do begin
   readf,1,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10
   t(k) = x0
   y(k) = x3 ;sqrt(x3*x3+x4*x4)/jmx
endfor
close, 1

dw1 = 0.5
dw2 = -1.0/6.0
;average over neighboring frequencies
pwa = pw
for nflt = 1,numflt do begin
   for k=1,nfreq-2 do begin
      pwa(k) = (dw1*(pw(k+1)+pw(k-1))+pw(k))/(1+2*dw1)
   endfor
   pw = pwa
   for k=1,nfreq-2 do begin
      pwa(k) = (dw2*(pw(k+1)+pw(k-1))+pw(k))/(1+2*dw2)
   endfor
   pw = pwa
endfor

pw=pw/max(pw)

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

flnm='spectrum.ps'
device, file=flnm, /color, xsize=21, ysize=25,xoffset=0,yoffset=2 
;/encapsulated doesn't wrok with yoffset. offset and size in cm

plot, f(0:nfreq-1),pw(0:nfreq-1), $
      xstyle=1, ystyle=1,xtitle='f(KHz)', ytitle='power density', $
      xrange=[-120,0],yrange=[0,1],linestyle=0
device, /close

flnm='y.ps'
device, file=flnm, /color, xsize=21, ysize=25,xoffset=0,yoffset=2 
;/encapsulated doesn't wrok with yoffset. offset and size in cm

plot, t(0:nfreq-1),y(0:nfreq-1), xrange=[30000,35000],$
      xstyle=1, ystyle=1,xtitle='n',linestyle=0
device, /close

return
end
