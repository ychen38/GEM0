pro reyn

npg=600L
k=1L
n=5000L
imx=128L
s=fltarr(imx)
rey=fltarr(imx)
reyh=fltarr(imx)
maxw=fltarr(imx)
maxwe=fltarr(imx)
tot=fltarr(imx)
rin=0.1
rout=0.7
ds=(rout-rin)/imx
for k = 0,imx-1 do begin
   s(k) = rin+k*ds
endfor

infl='stress'
dumchar=' '
openr, 1, infl
z=0.d
for k = 0,npg-1 do begin
   readf,1,t,rey
   readf,1,t,reyh
   readf,1,t,maxw
   readf,1,t,maxwe
   readf,1,t,tot
endfor

close, 1

for k = 0,imx-1 do begin
   print,rey(k),maxw(k),tot(k)
endfor

;goto, jump1


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

flnm='reyn.ps'
device, file=flnm, /color, xsize=21, ysize=25,xoffset=0,yoffset=2 
;/encapsulated doesn't wrok with yoffset. offset and size in cm

ymax=max([max(rey),max(reyh),max(maxw),max(maxwe)]) ;,max(tot)])
ymin=min([min(rey),min(reyh),min(maxw),min(maxwe)]) ;,min(tot)])

plot, s, rey, title='', xtitle='r/a', ytitle='',linestyle=0, xrange=[rin,rout],yrange=[ymin,ymax] 
;oplot,s,tot,color=lan,linestyle=0
;oplot,s,reyh,color=0,linestyle=1
oplot,s,maxw,color=lu,linestyle=0
;oplot,s,maxwe,color=lu,linestyle=1

oplot,s,maxw+rey,color=hong,linestyle=0



device, /close



return
end
