pro frq

k=1L
nr=7
nfreq=640L
f=fltarr(nr,nfreq)
pw=dblarr(nr,nfreq)
infl='freq'
dumchar=' '
openr, 1, infl
z=0.d
for i = 0,nr-1 do begin
    for k=0,4 do begin
        readf,1, dum,dum,dum
    endfor
    readf,1,dumchar ;frequency,ir=1
    print,dumchar
    for k=0,nfreq-1 do begin
        readf,1,x,y,z
        f(i,k) = y
        pw(i,k) = z
    endfor
    readf,1,dumchar ;camp(ir,t)
    print,dumchar
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

flnm='frq.ps'
device, file=flnm, /color, xsize=21, ysize=25,xoffset=0,yoffset=2 
;/encapsulated doesn't wrok with yoffset. offset and size in cm

for i=0,nr-1 do begin
    plot, f(i,0:nfreq-1),pw(i,0:nfreq-1), title='ir='+String(i,format='(I1)'), $
      xstyle=1, ystyle=1,xtitle='omega', linestyle=0
endfor

device, /close

return
end
