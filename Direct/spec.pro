pro spec
ol = -.3e-2
ou = -ol

omega=1/sqrt(0.00124*2.)/1048
tstart = 2500.0
nstart = long(tstart/omega/5.0/10)
nend = nstart+300

IU = complex(0,1)
n = 3000L
time=fltarr(n)
phir=dblarr(7,n)
phii=dblarr(7,n)
aparr=dblarr(7,n)
apari=dblarr(7,n)
dt=dblarr(n)

k=1L
nr=7
nfreq=100L
f=fltarr(nfreq)
pw=dblarr(nr,nfreq)
pwa=dblarr(nr,nfreq)

infl='yyre2'
dumchar=' '
openr, 1, infl
z=float(0.0)
for k = 0,n-1 do begin
   for i = 0,6 do begin
      readf,1,z,idum,x1,x2,x3 ;,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17
      time(k) = z
      if(idum NE i)then print,'idum not equal to i"
      phir(i,k) =  x2
      phii(i,k) = x3

      readf,1,z,idum,x1,x2,x3 ;,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17
      if(idum NE i)then print,'idum not equal to i"
      aparr(i,k) =  x2
      apari(i,k) = x3
    endfor
endfor

for k = 0,n-2 do begin
   dt(k) = time(k+1)-time(k)
;   print,k,dt(k)
endfor
dt(n-1) = dt(n-2)
T = double(0.0)
T = time(n-1)-time(0)

;print,ol,ou
domeg = (ou-ol)/float(nfreq)
for i = 0,nfreq-1 do begin
   f(i) = ol+domeg*i
;   print,i,f(i)
endfor

for i = 0,6 do begin
   for j = 0,nfreq-1 do begin
      omeg = f(j)
      cdum = complex(0,0)
      for k = nstart,nend-1 do begin
         cdum = cdum+complex(phir(i,k),phii(i,k))*exp(-IU*omeg*time(k))*dt(k)
      endfor
      pw(i,j) = (abs(cdum/(time(nend-1)-time(nstart))))^2
   endfor
endfor

dw1 = 0.5
dw2 = -1.0/6.0
;average over neighboring frequencies
pwa = pw
for nflt = 0,2 do begin
   for i = 0,nr-1 do begin
      for k=1,nfreq-2 do begin
         pwa(i,k) = (dw1*(pw(i,k+1)+pw(i,k-1))+pw(i,k))/(1+2*dw1)
      endfor
   endfor
   pw = pwa
   for i = 0,nr-1 do begin
      for k=1,nfreq-2 do begin
         pwa(i,k) = (dw2*(pw(i,k+1)+pw(i,k-1))+pw(i,k))/(1+2*dw2)
      endfor
   endfor
   pw = pwa
endfor

pw=pw/max(pw)

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

!p.multi=[0,2,4,0]
    !P.Charthick=5.0
    !P.Charsize=2.
    !X.Charsize=1.0
    !Y.Charsize=1.0
    !P.Thick=5.
    !X.Thick=5.0
    !Y.Thick=5.0

flnm='spec.ps'
device, file=flnm, /color, xsize=21, ysize=25,xoffset=0,yoffset=2 
;/encapsulated doesn't wrok with yoffset. offset and size in cm

for i=0,nr-1 do begin
    plot, f(0:nfreq-1),pw(i,0:nfreq-1), title='ir='+String(i,format='(I1)'), $
      xstyle=1, ystyle=1,xtitle='omega', linestyle=0
endfor

device, /close

return
end
