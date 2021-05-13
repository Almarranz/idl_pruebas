pro align_shift_images,band,exptime,cube_d


;chip = 1
;exptime=10
;band = 'H'
;cube_d=8
exptime_A=2
;folder='im_dark/'
folder='im_jitter_gains/'
;folder='im_jitter_NOgains/'
;folder ='im_sky_ESOReflex/'
py_pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/py_pruebas/'
indir = '/Users/amartinez/Desktop/PhD/HAWK/The_Brick/07_Cleancubes/054_'+band+'/dit_'+strn(exptime)+'/'+folder
indir_H='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/07_Cleancubes/054_H/dit_'+strn(exptime_A)+'/'+folder
;indir = '/Users/amartinez/Desktop/PhD/HAWK/The_Brick/06_Reduce/054_'+band+'/dit_'+strn(exptime)+'/'+folder

outdir='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/07.1_Reduce_aligned/054_'+band+'/dit_'+strn(exptime)+'/'+folder
;outdir=py_pruebas
mask_path='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/04_Makemask/054_'+band+'/dit_'+strn(exptime)+'/im/'
wx=1000
wy=1000

for  chip=1,4 do begin
;---------------------------------------
;cube_H=readfits(indir_H+ 'im1'+'_NOgains_chip'+strn(chip)+'_dit'+strn(exptime_A)+'.fits',EXT=0,h0)
;ra=sxpar(h0,'RA')
;dec=sxpar(h0,'DEC')

cube_H=readfits(indir_H+ 'idl_im1'+'_gains_chip'+strn(chip)+'_dit'+strn(exptime_A)+'.fits',EXT=0,h0)
ra=sxpar(h0,'RA')
dec=sxpar(h0,'DEC')
;---------------------------------------
s=0
sizex_new = 2048 + s
sizey_new = 2048 + s

borders=fltarr(sizex_new,sizex_new)
borders[*,*]=1
borders[*,0:15]=0
borders[*,sizex_new-15:sizex_new-1]=0
borders[0:15,*]=0
borders[sizex_new-15:sizex_new-1,*]=0

;writefits,py_pruebas+'borders.fits',borders



total = fltarr(sizex_new+500,sizey_new+500,cube_d)
mask_cube=fltarr(sizex_new+500,sizey_new+500,cube_d);---------------MASK
;for i =1, 5 do begin
for i =1, cube_d do begin


	;s=500
	s=0
	sizex_new = 2048 + s
	sizey_new = 2048 + s

	size_cubito= wx+s
	new_ref = fltarr(sizex_new+500,sizey_new+500)
	new_ref_mask=fltarr(sizex_new+500,sizey_new+500);---------------MASK
	;total = fltarr(sizex_new+500,sizey_new+500,5)

	;cube=readfits(indir+ 'im'+strn(i)+'_NOgains_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=0,header_i0)
	cube=readfits(indir+ 'idl_im'+strn(i)+'_gains_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=0,header_i0)
	mask=readfits(mask_path + 'mask' + strn(chip) +'_dit'+strn(exptime)+ '.fits')

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

	;sxaddpar, header1, 'CRPIX1 ',2163+s/2
	;sxaddpar, header1, 'CRPIX2 ',2164+s/2

	new_ref[500/2:sizex_new+500/2-1,500/2:sizex_new+500/2-1]=cube
	new_ref_mask[500/2:sizex_new+500/2-1,500/2:sizex_new+500/2-1]=mask;---------------MASK

	;new_ref=cube

if i eq 1 then begin
;total[s/2-x_off:sizex_new-s/2-1-x_off,s/2-y_off:sizex_new-s/2-1-y_off]=cube
; extract a small region from near image centre
	   ; that is used for final fine alignment
	   ;cubito_ref=fltarr(sizex_new,sizey_new)
       ;header=readfits(indir+ 'im'+strn(i)+'_NOgains_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=0,header0)
       ;cube=readfits(indir+ 'im'+strn(i)+'_NOgains_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=1,header1)
       
       cubito_ref=fltarr(sizex_new,sizey_new)
       header=readfits(indir+ 'idl_im'+strn(i)+'_gains_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=0,header0)
       cube=readfits(indir+ 'idl_im'+strn(i)+'_gains_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=1,header1)
       
       ;cube_canvas=readfits(indir+ 'im'+strn(i)+'_NOgains_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=1,header_canvas)
       ;cube_canvas=readfits(indir+ 'idl_im'+strn(i)+'_gains_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=1,header_canvas)
       
	   small_ref = fltarr(wx,wy)
	   a = sizex_new/2-wx/2
	   b = sizey_new/2-wy/2
	   small_ref = cube[a+x_off:a+wx+x_off, b+y_off:b+wy+y_off]
	   ;small_ref = cube[x_off:wx+x_off, y_off:wy+y_off]
	   ; set borders of small image to 1
	   small_ref[0:2,*]=1
	   small_ref[wx-2:wx,*]=1
	   small_ref[*,0:2]=1
	   small_ref[*,wy-2:wy]=1
	   
	   x_pix=sxpar(header1,'CRPIX1')
	   y_pix=sxpar(header1,'CRPIX2')
	   
	   sxaddpar, header1, 'CRPIX1 ',x_pix+s/2+x_for_header
       sxaddpar, header1, 'CRPIX2 ',y_pix+s/2+y_for_header
       
       sxaddpar, header_canvas, 'CRPIX1 ',x_pix+500/2
       sxaddpar, header_canvas, 'CRPIX2 ',y_pix+500/2
       
       ;cubito_ref[s/2:size_cubito-s/2,s/2:size_cubito-s/2]=small_ref
	   total[500/2-x_off:sizex_new+500/2-1-x_off,500/2-y_off:sizex_new+500/2-1-y_off,i-1]=cube
	   mask_cube[500/2-x_off:sizex_new+500/2-1-x_off,500/2-y_off:sizex_new+500/2-1-y_off,i-1]=mask;---------------MASK
	   
	   openw, outp, outdir+'xy_off_xy_alig_chip'+strn(chip)+'.txt', /get_lun, /APPEND
       printf, outp, format='(7f13.3)', x_off, y_off, 0, 0
       free_lun, outp

;writefits, py_pruebas + 'new_ref'+strn(i)+'.fits', new_ref;, /COMPRESS
;writefits, py_pruebas + 'cubit_ref'+strn(i)+'.fits',cubito_ref;, header, header0;, /COMPRESS   
;writefits, py_pruebas + 'small_ref'+strn(i)+'.fits',small_ref;, header, header0;, /COMPRESS
writefits, outdir + 'im'+strn(i)+'_chip'+strn(chip)+'.fits', header, header0;, /COMPRESS
writefits, outdir + 'im'+strn(i)+'_chip'+strn(chip)+'.fits', cube*borders, header1,/app;, /COMPRESS
;writefits, outdir + 'im'+strn(i)+'_chip'+strn(chip)+'_canvas.fits', header, header0
;writefits, outdir + 'im'+strn(i)+'_chip'+strn(chip)+'_canvas.fits', total,header_canvas,/app
print, 'termina con imagen', i
endif else begin
print, 'empieza con imagen', i
       
       ;cabeza=readfits(indir+ 'im'+strn(i)+'_NOgains_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=0,cabezai)
       ;cube=readfits(indir+ 'im'+strn(i)+'_NOgains_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=1,header_i)
       ;mask=readfits(mask_path + 'mask' + strn(chip) +'_dit'+strn(exptime)+ '.fits');---------------MASK
       ;cube_canvas=readfits(indir+ 'im'+strn(i)+'_NOgains_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=1,header_canvas)
       cube_canvas=readfits(indir+ 'idl_im'+strn(i)+'_gains_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=1,header_canvas)
       cabeza=readfits(indir+ 'idl_im'+strn(i)+'_gains_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=0,cabezai)
       cube=readfits(indir+ 'idl_im'+strn(i)+'_gains_chip'+strn(chip)+'_dit'+strn(exptime)+'.fits',EXTEN_NO=1,header_i)
       mask=readfits(mask_path + 'mask' + strn(chip) +'_dit'+strn(exptime)+ '.fits');---------------MASK
       
       
       
       small_new = fltarr(wx,wy)
	   ;a = sizex_new/2-wx/2
	   ;b = sizey_new/2-wy/2
	   ;small_new = new_ref[x_off:wx+x_off, y_off:wy+y_off]
	   small_new = cube[a+x_off:a+wx+x_off, b+y_off:b+wy+y_off]
	   ; set borders of small image to 1
	   small_new[0:2,*]=1
	   small_new[wx-2:wx,*]=1
	   small_new[*,0:2]=1
	   small_new[*,wy-2:wy]=1
	   
	   correl_optimize, small_ref, small_new, x_off_s, y_off_s, MAGNIFICATION=4, /NUMPIX
	   print,'comparando con.....','im'+strn(i)+'_NOgains_chip'+strn(chip)+'_dit'+strn(exptime)
       print, 'Offsets from correlation: ' + strn(x_off_s) + ', ' + strn(y_off_s)
       
       ;total = fltarr(sizex_new,sizey_new)
       ;cubito=fltarr(sizex_new,sizey_new)
       ;new_ref = fltarr(sizex_new,sizey_new)
       ;new_ref[s/2+x_off_s:sizex_new-s/2-1+x_off_s,s/2+y_off_s:sizex_new-s/2-1+y_off_s]=cube
       ;cubito[s/2+x_off_s:size_cubito-s/2+x_off_s,s/2+y_off_s:size_cubito-s/2+y_off_s]=small_new
       total[500/2-x_off+x_off_s:sizex_new+500/2-1-x_off+x_off_s,500/2-y_off+y_off_s:sizex_new+500/2-1-y_off+y_off_s,i-1]=cube
       mask_cube[500/2-x_off+x_off_s:sizex_new+500/2-1-x_off+x_off_s,500/2-y_off+y_off_s:sizex_new+500/2-1-y_off+y_off_s,i-1]=mask;-----------------------MASK
       sxaddpar, header_canvas, 'CRPIX1 ',x_pix+500/2
       sxaddpar, header_canvas, 'CRPIX2 ',y_pix+500/2
       
       openw, outp, outdir+'xy_off_xy_alig_chip'+strn(chip)+'.txt', /get_lun, /APPEND
       printf, outp, format='(7f13.3)', x_off, y_off, x_off_s, y_off_s
       free_lun, outp
       
       sxaddpar, header1, 'CRPIX1 ',x_pix+s/2+x_for_header+x_off_s; x_for_header aligns the pointings by wcs coordinates. All images 
       sxaddpar, header1, 'CRPIX2 ',y_pix+s/2+y_for_header+y_off_s; must have the header from the reference ima(the first one)                         
       
       cube_shift=shift(cube,x_off_s,y_off_s)
       mask_shift=shift(mask,x_off_s,y_off_s)
       small_new_shift=shift(small_new,x_off_s,y_off_s)
       ;writefits, py_pruebas + 'cubito_new'+strn(i)+'.fits',cubito;, header, header0;, /COMPRESS
       ;writefits, py_pruebas + 'small_new_shif'+strn(i)+'.fits',small_new_shift;, header, header0;, /COMPRESS;--------ESTOS SON LOS BUENOS PARA COMPARAR EL ALINGMENT
       ;writefits, py_pruebas + 'small_new'+strn(i)+'.fits',small_new;, header, header0;, /COMPRESS           ;--------ESTOS SON LOS BUENOS PARA COMPARAR EL ALINGMENT
       ;writefits, outdir + 'im'+strn(i)+'_chip'+strn(chip)+'.fits', cabeza, cabezai;, /COMPRESS
       ;writefits, outdir + 'im'+strn(i)+'_chip'+strn(chip)+'.fits', cube_shift, header_i,/app;, /COMPRESS
       writefits, outdir + 'im'+strn(i)+'_chip'+strn(chip)+'.fits', header, header0;, /COMPRESS
       writefits, outdir + 'im'+strn(i)+'_chip'+strn(chip)+'.fits', cube_shift*borders, header1,/app;, /COMPRESS
       
       writefits, outdir + 'mask'+strn(i)+'_chip'+strn(chip)+'.fits', header, header0;, /COMPRESS
       writefits, outdir + 'mask'+strn(i)+'_chip'+strn(chip)+'.fits', mask_shift*borders, header1,/app;, /COMPRESS
       
       ;writefits, outdir + 'im'+strn(i)+'_chip'+strn(chip)+'_canvas.fits', header, header0
       ;writefits, outdir + 'im'+strn(i)+'_chip'+strn(chip)+'_canvas.fits', total,header_canvas,/app
       
       ;------------------------------------
        sxaddpar, header_canvas, 'CRVAL1',ra
        sxaddpar, header_canvas, 'CRVAL2',dec
        
        sxaddpar, header0, 'RA',ra
        sxaddpar, header0, 'DEC',dec
       ;--------------------------------------
        
        
       ; cubeA... means that the aligment was done with respect to image 1 of the stack, and NOT with respect to
       ;image 1 of H DIT2, like cube... was.
       mascara='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/07.1_Reduce_aligned/054_'+band+'/dit_'+strn(exptime)+'/'+folder
       writefits, outdir + 'cubeA_chip'+strn(chip)+'_canvas.fits', header, header0
       writefits, outdir + 'cubeA_chip'+strn(chip)+'_canvas.fits', total,header_canvas,/app
       
       writefits, mascara + 'maskA_chip'+strn(chip)+'_canvas.fits', header, header0;---------------MASK
       writefits, mascara + 'maskA_chip'+strn(chip)+'_canvas.fits', mask_cube,header_canvas,/app;---------------MASK
       
       
       
endelse	





;writefits, py_pruebas + 'new_ref'+strn(i)+'.fits', header, header0;, /COMPRESS
;writefits, py_pruebas + 'new_ref'+strn(i)+'.fits', new_ref, header1,/app;, /COMPRESS



endfor

print,'############### Done with chip ',chip,'###################'
print,'############### Done with chip ',chip,'###################'
print,'############### Done with chip ',chip,'###################'
endfor

;writefits, py_pruebas + 'new_ref'+strn(i)+'.fits', header, header0;, /COMPRESS
;writefits, py_pruebas + 'new_ref'+strn(i)+'.fits', new_ref, header1,/app;, /COMPRESS

end
