pro contour_xz,field,npg
;field='apa'
runcase=''
npage = 20
pages = 1*[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]+npg-npage-1
ind = lonarr(5)
infl='./out/out/'+field+'xz.out'
openr, 1, infl
readf, 1,ind
print, ind
im=ind(0)
jm=ind(1)
tm=ind(2)
nplot=ind(3)
dt=ind(4)
holdgrd=dblarr((im+1)*(jm+1)*tm)
grdq=dblarr(im+1,jm+1,tm)
readf, format='(e12.5)',1,holdgrd
close, 1


labs=[1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]

for k=0,tm-1 do begin
  for j=0,jm do begin
    for i=0,im do begin
      grdq(i,j,k)=holdgrd(i+j*(im+1)+k*(im+1)*(jm+1))
    endfor
  endfor
endfor

loadct, 32

!p.multi = [0,4,5,0,0]
set_plot, 'ps'
flnm=field+'xz'+runcase+'.ps'
device, file=flnm, /color, /encapsulated, $
      bits=8, xsize=21, ysize=26
for k = 0,npage-1 do begin
  rmax=max(grdq(0:im,0:jm,pages(k)))
  rmin=min(grdq(0:im,0:jm,pages(k)))
  levs=dindgen(15)
  rran=rmax-rmin
  levs=rran/15*levs
  levs=levs+rmin

  contour, grdq(0:im,0:jm,pages(k)), levels=levs, $
         xstyle=1, ystyle=1, charsize=.8, $
	 /fill, c_colors=indgen(15)*5,title=pages(k)

;  contour, grdq(0:im,0:jm,pages(k)), levels=levs, $
;         c_labels=labs, xstyle=1, ystyle=1, $
;	 c_linestyle=(levs lt 0.0), /overplot
endfor
device, /close
loadct, 0

return
end
