pro amp

mm=8000

t=fltarr(mm)
amp=fltarr(mm)
amp1=fltarr(mm)

infl='yyre'

openr, 1, infl
for k=0,mm-1 do begin
readf, 1,x1,x2,x3,x4
t(k) = x1
amp(k) = x2
amp1(k) = x3
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
!p.thick=2
flnm='amp.ps'
device, file=flnm, /encapsulated, /color,xsize=21, ysize=21

x1=min(amp1)
y1=max(amp1)
x2=min(amp)
y2=max(amp)
xmin=min([x1,x2])
xmax=max([y1,y2])
plot, t, amp1, title='', $
  xstyle=1, ystyle=1,xtitle='t',ytitle='amp',$
  linestyle=0,yrange=[xmin,xmax],$
  charsize=2,charthick=4,xticks=3
oplot,t,amp,color=hong
device, /close

return
end
