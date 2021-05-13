PRO aling_new, band,exptime,cube_d

chip = 1

;folder='im_dark/'
;folder='im_jitter_gains/'
folder='im_jitter_NOgains/'
;folder ='im_sky_ESOReflex/'

data='_NOgains'
exptime_A=2
py_pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/py_pruebas/'
indir = '/Users/amartinez/Desktop/PhD/HAWK/The_Brick/07_Cleancubes/054_'+band+'/dit_'+strn(exptime)+'/'+folder
indir_band='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/07_Cleancubes/054_'+band+'/dit_'+strn(exptime_A)+'/'+folder
mask_path='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/04_Makemask/054_'+band+'/dit_'+strn(exptime)+'/im/'

outdir='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/07.1_Reduce_aligned/054_'+band+'/dit_'+strn(exptime)+'/'+folder

;outdir=py_pruebas

wx=1000
wy=1000

for  chip=chip,4 do begin
s=500
sizex_new = 2048 + s
sizey_new = 2048 + s

total = fltarr(sizex_new,sizey_new,cube_d)
mask_cube=fltarr(sizex_new,sizey_new,cube_d);---------------MASK

;Cogemos las coordenadas de la primera imagen de DIT2 porque tiene CUMOFFSETX (e Y) =0
;Esas coordenada las vamos a colocar en la cabecera del cubo porque asi se alinea con WCS
;Debe haber una manera de mejor de poner todas las imagenes en un cubo alineadas cons WCS sin
;tener que robar las coordenas de otra imagen

;---------------------------------------
cube_band=readfits(indir_band+ 'im1'+'_NOgains_chip'+strn(chip)+'_dit'+strn(exptime_A)+'.fits',EXT=0,h0)
ra=sxpar(h0,'RA')
dec=sxpar(h0,'DEC')
;---------------------------------------

for i =1, cube_d do begin

new_ref = fltarr(sizex_new,sizey_new)
new_ref_mask=fltarr(sizex_new,sizey_new)

cube=readfits(indir+ 'im'+strn(i)+data+'_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=1,header_i)
cube0=readfits(indir+ 'im'+strn(i)+data+'_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=0,header_i0)
mask=readfits(mask_path + 'mask' + strn(chip) +'_dit'+strn(exptime)+ '.fits')


sx=size(cube)

x_off = strsplit(header_i0[601],'HIERARCH ESO SEQ CUMOFFSETX = ', ESCAPE = '/', /extract)
y_off = strsplit(header_i0[602],'HIERARCH ESO SEQ CUMOFFSETY = ', ESCAPE = '/', /extract) 

print,x_off
print,x_off

x_for_header=float(x_off[0])
y_for_header=float(y_off[0])

x_off = fix(x_off[0])
y_off = fix(y_off[0])

print,x_off
print,y_off

new_ref[s/2:sx[1]+s/2-1,s/2:sx[1]+s/2-1]=cube
new_ref_mask[s/2:sx[1]+s/2-1,s/2:sx[1]+s/2-1]=mask;---------------MASK


if i eq 1 then begin

header=readfits(indir+ 'im'+strn(i)+data+'_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=0,header0)
cube=readfits(indir+ 'im'+strn(i)+data+'_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=1,header1)
cube_canvas=readfits(indir+ 'im'+strn(i)+data+'_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=1,header_canvas)


small_ref = fltarr(wx,wy)
a = sizex_new/2-wx/2
b = sizey_new/2-wy/2
small_ref = cube[a+x_off:a+wx+x_off, b+y_off:b+wy+y_off]
small_ref[0:2,*]=1
small_ref[wx-2:wx,*]=1
small_ref[*,0:2]=1
small_ref[*,wy-2:wy]=1

x_CRPIX1=sxpar(header1,'CRPIX1')
y_CRPIX2=sxpar(header1,'CRPIX2')

sxaddpar, header_canvas, 'CRPIX1 ',x_CRPIX1+s/2
sxaddpar, header_canvas, 'CRPIX2 ',y_CRPIX2+s/2

total[s/2-x_off:sx[1]+s/2-1-x_off,s/2-y_off:sx[1]+s/2-1-y_off,i-1]=cube
mask_cube[s/2-x_off:sx[1]+s/2-1-x_off,s/2-y_off:sx[1]+s/2-1-y_off,i-1]=mask;---------------MASK

openw, outp, outdir+'xy_off_xy_alig_chip'+strn(chip)+'.txt', /get_lun
printf, outp, format='(7f13.3)', x_off, y_off, 0, 0
free_lun, outp


print, 'termina con imagen', i
endif else begin;#################################################################################
print, 'empieza con imagen', i

cabeza=readfits(indir+ 'im'+strn(i)+data+'_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=0,cabezai)
cube=readfits(indir+ 'im'+strn(i)+data+'_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=1,header_i)
mask=readfits(mask_path + 'mask' + strn(chip) +'_dit'+strn(exptime)+ '.fits');---------------MASK

       small_new = fltarr(wx,wy)
	   small_new = cube[a+x_off:a+wx+x_off, b+y_off:b+wy+y_off]
	   small_new[0:2,*]=1
	   small_new[wx-2:wx,*]=1
	   small_new[*,0:2]=1
	   small_new[*,wy-2:wy]=1
	   
	   correl_optimize, small_ref, small_new, x_off_s, y_off_s, MAGNIFICATION=4, /NUMPIX
	   print,'comparando con.....','im'+strn(i)+data+'_chip'+strn(chip)+'_dit'+strn(exptime)
       print, 'Offsets from correlation: ' + strn(x_off_s) + ', ' + strn(y_off_s)
       
       openw, outp, outdir+'xy_off_xy_alig_chip'+strn(chip)+'.txt', /get_lun, /APPEND
       printf, outp, format='(7f13.3)', x_off, y_off, x_off_s, y_off_s
       free_lun, outp

       total[s/2-x_off+x_off_s:sx[1]+s/2-1-x_off+x_off_s,s/2-y_off+y_off_s:sx[1]+s/2-1-y_off+y_off_s,i-1]=cube
       mask_cube[s/2-x_off+x_off_s:sx[1]+s/2-1-x_off+x_off_s,s/2-y_off+y_off_s:sx[1]+s/2-1-y_off+y_off_s,i-1]=mask;-----------------------MASK

endelse


endfor

;------------------------------------
sxaddpar, header_canvas, 'CRVAL1',ra
sxaddpar, header_canvas, 'CRVAL2',dec
        
sxaddpar, header0, 'RA',ra
sxaddpar, header0, 'DEC',dec
;--------------------------------------

writefits, outdir + 'cube_chip'+strn(chip)+'_canvas.fits', header, header0
writefits, outdir + 'cube_chip'+strn(chip)+'_canvas.fits', total,header_canvas,/app

writefits, outdir + 'mask_cube_chip'+strn(chip)+'_canvas.fits', header, header0
writefits, outdir + 'mask_cube_chip'+strn(chip)+'_canvas.fits', mask_cube,header_canvas,/app

print, '############## Done witn chip',chip,'##############'
print, '############## Done witn chip',chip,'##############'
print, '############## Done witn chip',chip,'##############'

endfor

end





