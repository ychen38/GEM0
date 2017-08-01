pro mdhis

omega=1/sqrt(0.00124*2.)/1048
print,'Omeg_alf = ', omega

k=1L
n=5000L
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

;goto, jump1

infl='mdhis4'
np=1998L
timep=fltarr(np)
a0p=fltarr(np)
a1p=fltarr(np)
dumchar=' '
openr, 1, infl
z=0.d
for k = 0,np-1 do begin
   readf,1,t,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16
   timep(k) = t*omega+200
   a0p(k) = x0
   a1p(k) = x1
endfor


close, 1

;goto, jump2

infl='mdhis3'
nq=2250L
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


jump2:

set_plot, 'ps'

;dum=   [0   1    2    3    4    5     6     7    8    9    10   11   12   13   14   15    16    17] 
RED =   [0,  1,   1,   0,   0,   .9,   1.,  .2,   .4,  .5,   .2, .0,  .5,  0,   .2,   .5,   0.2,  0.7]
GREEN = [0,  1,   0,   1,   0,   .2,   .6,  .7,   .4,  .7,   .4, .3,  0,   .5,  .2,   .5,   0.7,   0.2]
BLUE =  [0,  0,   0,   0,   1,   .9,   .1,  1.,   .4,  .6,  .75,  .3,  0,   0,   .5,   0,   1,   1]

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

flnm='mdhis.ps'
device, file=flnm, /color, xsize=21, ysize=25,xoffset=0,yoffset=2 
;/encapsulated doesn't wrok with yoffset. offset and size in cm

plot, time, a0, title='', xtitle='!4x!l!3A!Nt', ytitle='rms(!4u!3)',linestyle=0, xrange=[0,time(n-1)],yrange=[1e-8, 2e-1],/ylog
;oplot,time,a1,color=lan,linestyle=0

;oplot,timep-25,a1p,color=15,linestyle=0
;oplot,timeq+25,a1q,color=13,linestyle=0

oplot,time,a1,color=3,linestyle=0
;oplot,time,a2,color=3,linestyle=0

;oplot,time,a3,color=4,linestyle=0
;oplot,time,a4,color=5,linestyle=0
;oplot,time,a5,color=6,linestyle=0
;oplot,time,a6,color=7,linestyle=0

;oplot,time,a7,color=8,linestyle=0
;oplot,time,a8,color=9,linestyle=0


;oplot,time,a8,color=11,linestyle=0
;oplot,time,a12,color=12,linestyle=0

;xyouts,10,0.256,'n=14 single',color=15
;xyouts,10,0.512,'n=18 single',color=13

xyouts,10,0.128,'n=0'
;xyouts,10,0.064,'n=3',color=2
xyouts,10,0.032,'n=4',color=3
;xyouts,10,0.016,'n=5',color=4
;xyouts,10,0.008,'n=6',color=5
;xyouts,10,0.004,'n=7',color=6
;xyouts,10,0.002,'n=8',color=7
;xyouts,10,0.001,'n=9',color=8
;xyouts,10,0.0005,'n=10',color=9



device, /close



return
end
