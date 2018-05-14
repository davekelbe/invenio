FUNCTION batch_ica_generic, dir_upper, ICA=ica, PCA=pca
; batch ICA based on pre-chosen band sets and ROI's 
;  batch_ica_generic(directory, /ica, /pca)
;  Keywords /ica and /pca are optional. 
;  Only allowed to use one at a time
;  If none provided, defaults to ica 
; Example command
; print, batch_ica_generic('D:\ImageProcessing\OUTPUT\Processed-K0033\', /ica)
;   * Note trailing slash
;
; D:\ImageProcessing\OUTPUT\
;   Processed-K0033\
;     K0033_000003\
;       K0033_000003+envi     << This is where it will search for cube etc
;       K0033_000003+matlab
;     K0033_000004\
;     K0033_000005\ ... etc
;   Processed-K0033-best\
;     K0033_000003          << This is where the best images are stored
;     K0033_000004
;     K0033_000005 ... etc
;   Processed-K0033-jpg\    << This is where all of the jpgs ICAs are
;   Processed-K0033-tif\
;     K0033_000003\
;       K0033_000003+tif    << This is where the output is saved
;
;

last_character = strmid(dir_upper, strlen(dir_upper)-1,1)
print, last_character
if ~strmatch(last_character, '/') && ~strmatch(last_character, '\') then begin
  dir_upper = dir_upper + '/'
endif

if KEYWORD_SET(ica) then is_ica = 1 $ 
else if KEYWORD_SET(pca) then is_ica = 0 $ 
else is_ica = 1;

if is_ica then begin
  print, 'Running ICA'
endif else begin
  print, 'Running PCA'
endelse

free_all 
print_open_lun
; Restore all the base save files
envi, /restore_base_save_files
 
; NOT IMPLEMENTED is_double = 0; double-precision (1) or single-precision floating point (0)
 
; Initialize ENVI in the batch mode and send all errors and warnings to the file batch.txt
envi_batch_init, /no_status_window;, log_file='batch.txt'
; LUN 100 for batch.txt file 
ENVI_BATCH_STATUS_WINDOW, /OFF
; Contrast stretch of output images = 1% 
lowfrac = 0.005

; dir_uppper now function input 
;dir_upper = '/dirs/grad/djk2312/Deliver_DK_17-Jun-2015/'
dir_temp = dir_upper
;dir_source = '/Volumes/Corinth/MOTB/Flattened/'
;dir_source = '/dirs/grad/djk2312/Flattened/EMEL/'
info_slash = '/'
if (is_ica) then begin  
  dir_write = strmid(dir_upper,0,strlen(dir_upper)-1)+'-' + 'jpg_ica' + info_slash
endif else begin
  dir_write = strmid(dir_upper,0,strlen(dir_upper)-1)+'-' + 'jpg_pca' + info_slash
endelse
;folders = dialog_pickfile(path=dir_upper,/directory, /multiple_files)
;shelf_mark = 65
;shot_start = 1
;shot_end = 1

filepath_envitemp = dir_temp + 'temp.img'
set_tf = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
nbins = 5000; for histogram computation (98% scale)

if ~file_test(dir_write, /directory) then begin
   command = 'mkdir ' + dir_write;
   spawn, command 
endif 

last_character = strmid(dir_upper, strlen(dir_upper)-1,1)
if ~ strmatch(last_character, info_slash) then begin
    dir_upper = dir_upper + info_slash
endif

;SelectedFolders = dialog_pickfile(path=dir_upper, /DIRECTORY, /MULTIPLE_FILES, title='Please choose folios to process')
subdir = file_search(dir_upper, /test_directory, '*')
n_subdir = n_elements(subdir)

; Remove subdir to only upper level
is_valid = replicate(0,n_subdir,1)
for f = 0,n_subdir-1 do begin
  ix1 = strpos(subdir[f], 'matlab')
  ix2 = strpos(subdir[f], 'envi')
  ix3 = strpos(subdir[f], 'tiff')
  if (ix1 EQ -1) and (ix2 EQ -1) and (ix3 EQ -1) then begin
    is_valid[f] = 1
  endif
endfor

ix = where(is_valid eq 1);
subdir = subdir[ix]
n_subdir = n_elements(ix);

;for shot_seq = shot_start, shot_end do begin ;for each shot sequence 
for f = 0,n_subdir-1 do begin
  SelectedFolders = subdir[f] + info_slash
  SelectedFolders_noslash = strmid(SelectedFolders, 0, strlen(SelectedFolders)-1)
  ix_slash = strpos(SelectedFolders_noslash, info_slash, /reverse_search)
  cube = strmid(SelectedFolders_noslash,ix_slash+1, strlen(SelectedFolders_noslash)-ix_slash)
  ;cube =  string(shelf_mark, format='(I04)') + '_' + string(shot_seq, format='(I06)')
  print, 'Folio ' + cube
  dir_roi = dir_upper+cube+info_slash+cube+'+envi'+info_slash
  dir_matlab = dir_upper+cube+info_slash+cube+'+matlab'+info_slash
  dir_tif = dir_upper+cube+info_slash+cube+'+tiff'+info_slash

  ; Continue if no images for cube 
  tf_exist = file_test(dir_roi,/directory)
  if tf_exist EQ 0 then begin
    continue
  endif
  
   dir_out_up = dir_write + cube + info_slash
   if ~file_test(dir_out_up, /directory) then begin
      command = 'mkdir ' + dir_out_up;
      spawn, command 
   endif
   dir_out = dir_out_up;
;   if (is_ica) then begin
;    dir_out = dir_out_up + cube + "+pca_jpg" + info_slash
;   endif else begin
;    dir_out = dir_out_up + cube + "+ica_jpg" + info_slash
;   endelse 
   if ~file_test(dir_out, /directory) then begin
      command = 'mkdir ' + dir_out;
      spawn, command 
   endif 
   
  ; Read normalization values 
  ;filepath_summary = dir_matlab + cube + '_summary.txt'
  ; summary = read_ascii(filepath_summary)
  ; normalization_val = summary.FIELD1(3,*)
  ; FMT = 'A,F,F,D'
  ; READCOL, filepath_summary,F=FMT, wavelength,aperture, shutter,normalization

  ; Get the list of ROI's 
  wildcard = dir_roi+'*mask*.tif'
  filepath_roi = file_search(wildcard);
  n_roi = n_elements(filepath_roi)
  if (n_roi eq 1) then continue 
 
  ; Read cube text file 
  filepath_cubes = dir_upper + cube + info_slash + cube + '+matlab' + info_slash + cube + '_' + 'cube.txt'
  READSTRING, filepath_cubes, lines
  lines = transpose(lines);
  lines = [[lines],['']]
  n_lines = n_elements(lines);

  ; Find start of each new cube header in list 
  is_cube = strmatch(lines,'*Cube*')
  ; logical indexing
  ix = transpose(where(is_cube EQ 1))
  set_start = [[ix], [n_lines]]
  set_end = set_start - 1
  n_cubes = n_elements(set_start)-1
  ;set_end[n_cubes-1] = set_end[n_cubes-1] + 1
  set_start = transpose(set_start[0:n_cubes-1])
  set_end = transpose(set_end[1:n_cubes])
  set_num = set_end - set_start-1;
  n_set = n_elements(set_start)
  
  ; for each cube set  
  for s = 0,n_set-1 do begin 
    if ~set_tf[s] then continue
    print, format='(%"\tcube set %d")', s+1
    ; check if all rois are already produced 
;    all_exist = 1
;    for r = 0,n_roi-1 do begin 
;        filepath_cubetemp = dir_roi + 'temp' + 'c' + strtrim(s+1,2) + 'r' + strtrim(r+1,2) +'.img'
;        if ~file_test(filepath_cubetemp) then all_exist = 0
;    end
;    if all_exist then continue
        
    ;print, format='(%"\tn_bands %d")', set_num[s]
    ; load images 
    set = lines[set_start[s]]
    set_str = strmid(set, 4,2)
    bands = transpose(lines[set_start[s]+1:set_end[s]])
    
    ; for each band
    ; read and add to multiband image I 
    for b =0,set_num[s]-1 do begin
      band_no_space = strreplacef(bands[b], ' ', '')
      filepath_image = dir_tif  + cube + band_no_space + '.tif'
      im = 0; clear memory 
      im = read_tiff(filepath_image)
      im = double(im) 

      ;normalize image by calibration value 
      ;strfind, wavelength, band_no_space, INDEX=ix
      ;normval = double(normalization(ix))
      ;normval = normval[0]
      ;im = im[*,*]/normval
      
      if (b eq 0) then begin 
        imsize = size(im, /dimensions)
        dispx = floor(imsize[0]/2.5)
        dispy = floor(imsize[1]/2.5)
        I = 0; 
        I = make_array(imsize[0],imsize[1],set_num[s], /double)
      endif 
      ;i_scl = bytscl(i);    
      ;dispx = floor(imsize[0]/10)
      ;dispy = floor(imsize[1]/10)
      ;i_disp = congrid(i_scl,dispx,dispy)
      ;window, 0, xsize=dispx, ysize=dispy
      ;tv, i_disp
      I[*,*,b] = im  
      im = 0;
    endfor ; load each cube into image 
    ; Output is multiband image. Now make subset image temp from roi 
    ; for each roi 
      Iscl = congrid(I,dispx, dispy, set_num[s]); can this be moved down?
    
    for r = 0,n_roi-1 do begin 
      filepath_cubetemp = dir_roi + 'temp' + 'c' + strtrim(s+1,2) + 'r' + strtrim(r+1,2) +'.img'
      print, format='(%"\t\tROI %d")', r+1    
      
      ;print, format='(%"\t\tROI %d")', r+1
      roi = read_tiff(filepath_roi[r])
      roi_str = filepath_roi[r]
      roi_str = strmid(roi_str, strlen(roi_str)-6, 2)
      
      if ~file_test(filepath_cubetemp) then begin 
      
        ;roi_disp = congrid(roi,dispx,dispy)
        ;window, 0, xsize=dispx, ysize=dispy
        ;tvscl, roi_disp    
        ;
        ; Get number of elements for squareish roi 
        isroi = where(roi eq 0);
        dims_is = size(isroi, /dimensions)  
        if n_elements(isroi) lt 100 then continue
        dims_square = floor(sqrt(dims_is)) 
        dims_is2 = dims_square^2
        relsize = (dims_is2/double((imsize[0]*imsize[1])))*100
        ;print, format='(%"\t\tSize %f")', relsize
        
        ; make E (roi image data) a square image 
        E = 0; 
        for b = 0,set_num[s]-1 do begin
          if (b eq 0) then begin
            E = make_array(dims_square,dims_square,set_num[s], /uint)
          endif 
          ;filepath_image = dir_upper + cube + info_slash + 'tiff' + info_slash + cube + '+' + bands[b] + '.tif'
          temp = I[*,*,b];
          e_roi = 0;  
          e_roi = temp[isroi]
          e_roi = e_roi[0:dims_is2-1]
          e_roi = reform(e_roi,dims_square,dims_square)
          E[*,*,b] = e_roi
        endfor 
        envi_write_envi_file, E, out_name=filepath_cubetemp
        ;Save scaled image and remove full size image
        ;I = 0;
      endif
        
      ; Do ICA computation 
      if (is_ica) then begin
        filename_trans = dir_upper + cube + info_slash + cube + '+envi' + info_slash + cube + '_ica' + '_c' + set_str + '_m' + roi_str + '.trans'
      endif else begin
        filename_trans = dir_upper + cube + info_slash + cube + '+envi' + info_slash + cube + '_pca' + '_c' + set_str + '_m' + roi_str + '.csv'
      endelse
      if ~file_test(filename_trans) then begin
        ;print, format='(%"\t\t\tOpening envi file")'      
        ENVI_OPEN_FILE, filepath_cubetemp, r_fid=c_fid, /invisible, /no_interactive_query, /no_realize
        ; remove temp cube file 
        ; Get cube dimensions 
        ;print, format='(%"\t\t\tComputing statistics")'      
        envi_file_query, c_fid, dims=dims, nb=nb, bnames=bnames
        pos = indgen(nb)
        ;print, format='(%"\t\t\tPerforming ICA on subset")'  
        if (is_ica) then begin
          envi_doit, 'envi_ica_doit', $
          fid=c_fid, pos=pos, dims=dims,$
          out_trans_name = filename_trans, /invisible, /no_realize, /IN_MEMORY
        endif else begin
          ; Compute statistics
          ENVI_Doit, 'ENVI_Stats_Doit', $
             FID = c_fid, $
             DIMS = dims, $
             POS = Lindgen(nb), $
             MEAN = Avg, $
             EVAL = Eval, $
             EVEC = Evec, $
             COMP_FLAG = 5
            ; Perform PC rotation
          ENVI_Doit, 'PC_Rotate', $
             FID = c_fid, $
             DIMS = dims, $
             POS = Lindgen(nb), $
             MEAN = Avg, $
             EVAL = Eval, $
             EVEC = Evec, $
             OUT_DT = 4, $
             OUT_NB = nb, $
             /FORWARD, $
             /IN_MEMORY, $ 
             /NO_PLOT, $ 
             /INVISIBLE, $ 
             /NO_REALIZE        
        endelse
        ;Free lun 
        if c_fid lt 129 then free_lun, c_fid
        free_lun, 100
        command = 'rm ' + filepath_cubetemp
        spawn, command
        command = 'rm ' + filepath_cubetemp + '.hdr'
        spawn, command
        ;extract transform csv
        ;print, format='(%"\t\t\tExtracting transform")'
        if is_ica then begin
          envi_extract_trans_file_noenvi, filename_trans
        endif else begin 

          if (size(evec, /n_elements) eq 0) then begin
            foo = 1;
          endif else begin
            write_csv_data, evec, filename=filename_trans
          endelse
        endelse       
      endif
      
      ; Save jpg images 
      J = 0; clear memory 
      J = reform(Iscl,[dispx*dispy, set_num[s]])
      J = double(J)
      ;filename_csv = strmid(filename_trans,0,strlen(filename_trans)-5) + 'csv'
      ;csv = read_csv(filename_csv) 
      if is_ica then begin
      envi_check_save, /transform
      transform_read, filename_trans, forward_matrix=tx
      endif else begin
      ; stuff eigenvectors into tx directly 
      tx = evec
      endelse 
      ;write_csv_data, tx, filename=filename_csv
      ;print, format='(%"\t\t\tApplying ICA to full image")'
      J_tx = 0; clear memory 
      J_tx = tx##J   
      ICA = 0; clear memory 
      ICA = reform(J_tx,[dispx,dispy,set_num[s]])
      minval = min(J_tx, dimension=1)
      maxval = max(J_tx, dimension=1);
      ;dispx = floor(imsize[0]/2.5)
      ;dispy = floor(imsize[1]/2.5)
      
      for q = 0,set_num[s]-1 do begin ; for each band 
        ;print, format='(%"\t\t\t\tSaving output band %d")', q+1
        pdf = histogram(J_tx[*,q],min=minval[q], max=maxval[q],nbins=500,locations=locations)
        cdf = total(pdf, /cumulative) / (dispx*dispy)
        index_low = where(cdf LE lowfrac)
        index_hi = where(cdf GE 1-lowfrac)
        val_low = locations[index_low[n_elements(index_low)-1]]
        val_hi = locations[index_hi[0]]
        ;ICA_out = (ICA[*,*,q]-val_low)/(val_hi - val_low) * double(65535)
        ;ICA_out(where(ICA LT 0)) = 0
        ;ICA_out(where(ICA GT 65535)) = 65535
        ;ICA_out = uint(ICA_out) 
        
        ;ICA_out = 0; clear memory 
        ;ICA_out = cgscalevector(ICA[*,*,q], 0,65535, minvalue=val_low, maxvalue=val_hi, double=1)
        ;ICA_out = uint(ICA_out)
        ;ICA_jpg = 0; clear memory 
        ;ICA_jpg = congrid(ICA[*,*,q],dispx,dispy)
        ICA_jpg = byte(cgscalevector(ICA[*,*,q], 0,255, minvalue=val_low, maxvalue=val_hi, double=1))
                
        c_str = string(s+1, format='(I02)')
        r_str = string(r+1, format='(I02)')
        q_str = string(q+1, format='(I02)')
        
        ;filepath_tiff = dir_upper + cube + info_slash + 'tiff' + info_slash + cube + $
        ;  '_c' + c_Str + '_m' + r_str + '_b' + q_str +  '.tif'
        if is_ica then begin
          filepath_out = dir_out + cube + $
            '_ica' + '_c' + c_Str + '_m' + r_str + '_b' + q_str +  '.jpg'
        endif else begin
          filepath_out = dir_out + cube + $
            '_pca' + '_c' + c_Str + '_m' + r_str + '_b' + q_str +  '.jpg'          
        endelse
       ; write_tiff, filepath_tiff, ICA_out, /short
       if ~file_test(filepath_out) then begin
          write_jpeg, filepath_out, ICA_jpg, order=1, quality = 50
       endif
      endfor   ; band   
              
    endfor; for roi
  endfor; cube set 
  
  ; Copy true color to cyclone directory
;  dir_jpg = dir_write + cube + info_slash + 'jpg' + info_slash;
;  wildcard = dir_jpg +'*true*'
;  filepath_source = file_search(wildcard);
;  n_file = n_elements(filepath_source)
;  for f = 0,n_file-1 do begin
;    filename = file_basename(filepath_source[f])
;    target_filepath = dir_jpg + filename
;    command = 'cp ' + filepath_source[f] + ' ' + target_filepath
;    spawn, command
;  endfor
;  ; Repeat for KTK images
;  wildcard = dir_jpg +'*KTK*'
;  filepath_source = file_search(wildcard);
;  n_file = n_elements(filepath_source)
;  for f = 0,n_file-1 do begin
;    filename = file_basename(filepath_source[f])
;    target_filepath = dir_jpg + filename
;    command = 'cp ' + filepath_source[f] + ' ' + target_filepath
;    spawn, command
;  endfor 
endfor ;shot sequence  
RETURN, 1
end
