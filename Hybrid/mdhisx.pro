pro mdhisx

omega=1/sqrt(0.008*2.5)/3000

k=1L
n=4746L
time=fltarr(n)
a0=fltarr(n)
a1=fltarr(n)
a2=fltarr(n)
a3=fltarr(n)
a4=fltarr(n)
a5=fltarr(n)
a6=fltarr(n)
a7=fltarr(n)
a8=fltarr(n)
a9=fltarr(n)
a10=fltarr(n)
a11=fltarr(n)
a12=fltarr(n)
a13=fltarr(n)
a14=fltarr(n)
a15=fltarr(n)
a16=fltarr(n)


infl='mdhis'
dumchar=' '
openr, 1, infl
z=0.d
for k = 0,n-1 do begin
   readf,1,t,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16
   time(k) = t*omega
   a0(k) = x0
   a1(k) = x1
   a2(k) = x2
   a3(k) = x3
   a4(k) = x4
   a5(k) = x5
   a6(k) = x6
   a7(k) = x7
   a8(k) = x8
   a9(k) = x9
   a10(k) = x10
   a11(k) = x11
   a12(k) = x12
   a13(k) = x13
   a14(k) = x14
   a15(k) = x15
   a16(k) = x16
endfor


close, 1



infl='mdhis-n14'
np=2600L
timep=fltarr(np)
a0p=fltarr(np)
a1p=fltarr(np)
dumchar=' '
openr, 1, infl
z=0.d
for k = 0,np-1 do begin
   readf,1,t,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16
   timep(k) = t*omega
   a0p(k) = x0
   a1p(k) = x1
endfor


close, 1

infl='mdhis-n18'
nq=3000L
timeq=fltarr(nq)
a0q=fltarr(nq)
a1q=fltarr(nq)
dumchar=' '
openr, 1, infl
z=0.d
for k = 0,nq-1 do begin
   readf,1,t,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16
   timeq(k) = t*omega
   a0q(k) = x0
   a1q(k) = x1
endfor


close, 1

set_plot, 'ps'

;dum=   [0   1    2    3    4    5    6    7    8    9    10   11   12   13   14   15    16    17] 
RED =   [0,  1,   1,   0,   0,   1,   .75, .75, 0,   0,   .75, .5,  .5,  0,   0,   .5,   0.2,  0.7]
GREEN = [0,  1,   0,   1,   0,   1,   .75,  0,  .75, 0,   .75, .5,  0,   .5,  0,   .5,   0.7,   0.2]
BLUE =  [0,  1,   0,   0,   1,   0,   .75,  0,  0,   .75,  0,  .5,  0,   0,   .5,   0,   1,   1]

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

flnm='mdhisx.ps'
device, file=flnm, /color, xsize=21, ysize=25,xoffset=0,yoffset=2 
;/encapsulated doesn't wrok with yoffset. offset and size in cm

plot, time, a0, title='', xtitle='!4x!l!3A!Nt', ytitle='rms(!4u!3)',linestyle=0, xrange=[0,500],yrange=[1e-6, 2e-1],/ylog
;oplot,time,a1,color=lan,linestyle=0

;oplot,timep-25,a1p,color=15,linestyle=0
;oplot,timeq+25,a1q,color=13,linestyle=0

oplot,time,a9,color=2,linestyle=0
oplot,time,a10,color=3,linestyle=0
oplot,time,a11,color=4,linestyle=0
oplot,time,a12,color=17,linestyle=0
oplot,time,a13,color=6,linestyle=0
oplot,time,a14,color=16,linestyle=0

oplot,time,a15,color=11,linestyle=0
oplot,time,a16,color=10,linestyle=0
;oplot,time,a11,color=11,linestyle=0
;oplot,time,a12,color=12,linestyle=0

;xyouts,10,0.256,'n=14 single',color=15
;xyouts,10,0.512,'n=18 single',color=13

xyouts,10,0.128,'n=0'
xyouts,10,0.064,'n=19',color=2
xyouts,10,0.032,'n=20',color=3
xyouts,10,0.016,'n=21',color=4
xyouts,10,0.008,'n=22',color=17
xyouts,10,0.004,'n=23',color=6
xyouts,10,0.002,'n=24',color=16
xyouts,10,0.001,'n=25',color=11
xyouts,10,0.0005,'n=26',color=10



device, /close



return
end
