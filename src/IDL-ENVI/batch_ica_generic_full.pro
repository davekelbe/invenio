FUNCTION batch_ica_generic_full, dir_task
;   batch ICA based on pre-chosen band sets and ROI's 
;  
; Example command 
; print, batch_ica_generic_full('D:\ImageProcessing\OUTPUT\Processed-K0033+best\')
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
last_character = strmid(dir_task, strlen(dir_task)-1,1)
print, last_character
if ~strmatch(last_character, '/') && ~strmatch(last_character, '\') then begin
  dir_task = dir_task + '/'
endif


free_all 
print_open_lun
; Restore all the base save files
envi, /restore_base_save_files
is_double = 0; double-precision (1) or single-precision floating point (0)
 
; Initialize ENVI in the batch mode and send all errors and warnings to the file batch.txt
envi_batch_init, /no_status_window;, log_file='batch.txt'
; LUN 100 for batch.txt file 
ENVI_BATCH_STATUS_WINDOW, /OFF
; Contrast stretch of output images = 1% 
lowfrac = 0.005
info_slash = '\'

; dir_uppper now function input 
dir_task = dir_task;'/Volumes/NewZealand/Deliver_DJK_11-Jan-2016-best/'
;dir_task = dialog_pickfile(/DIRECTORY, title='Please choose folder of best to process')
dir_upper = strmid(dir_task,0, strlen(dir_task)-6)+info_slash
;dir_upper = '/dirs/grad/djk2312/Deliver_DK_17-Jun-2015/'
dir_temp = dir_upper
;dir_source = '/Users/Kelbe/Desktop/'
;dir_source = '/Volumes/NewZealand/Deliver-tiff_DJK_12-Jan-2016/'
;dir_task = strmid(dir_upper,0,strlen(dir_upper)-1)+'-best/'

;folders = dialog_pickfile(path=dir_upper,/directory, /multiple_files)
;shelf_mark = 65
;shot_start = 1
;shot_end = 1

filepath_envitemp = dir_temp + 'temp.img'
nbins = 5000; for histogram computation (98% scale)



;SelectedFolders = dialog_pickfile(path=dir_upper, /DIRECTORY, /MULTIPLE_FILES, title='Please choose folios to process')

; List files within task 
;shelf_mark = string(shelf_mark, format='(I04)')
subdir = file_search(dir_task, /test_directory, '*')
n_subdir = n_elements(subdir)
for s = 0,n_subdir-1 do begin
  leng_subdir = strlen(subdir[s])
  suffix = strmid(subdir[s], leng_subdir -3,3)
  if strcmp(suffix, 'jpg') then begin
    continue
  endif
  SelectedFolders = subdir[s] + info_slash
  SelectedFolders_noslash = strmid(SelectedFolders, 0, strlen(SelectedFolders)-1)
  ix_slash = strpos(SelectedFolders_noslash, info_slash, /reverse_search)
  cube = strmid(SelectedFolders_noslash,ix_slash+1, strlen(SelectedFolders_noslash)-ix_slash)
    ;cube =  shelf_mark + '_' + shot_seq
    print, 'Folio ' + cube
    dir_roi = dir_upper+cube+info_slash+cube+'+'+'envi'+info_slash
    dir_matlab = dir_upper+cube+info_slash+cube+'+'+'matlab'+info_slash
    dir_source = dir_upper+cube+info_slash+cube+'+'+'tiff'+info_slash
 
    
   
   ; Read normalization values 
   ;filepath_summary = dir_matlab + cube + '_summary.txt'
   ;summary = read_ascii(filepath_summary)
   ;normalization_val = summary.FIELD1(3,*)
   ;FMT = 'A,F,F,D'
   ;READCOL, filepath_summary,F=FMT, wavelength,aperture, shutter,normalization
   
  ; Read cube text file 
  filepath_cubes = dir_upper + cube + info_slash + cube + '+' + 'matlab' + info_slash + cube + '_' + 'cube.txt'
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
  
    filename_task = file_search(SelectedFolders, '*')
    n_f = n_elements(filename_task)
    for f = 0,n_f-1 do begin
      ; Process each task image 
      ; 
      processing_type = strmid(filename_task[f], $
        strlen(filename_task[f])-19,3);
        if strcmp(processing_type, 'ica') then begin 
          is_ica = 1;
        endif else if strcmp(processing_type, 'pca') then begin 
            is_ica = 0; 
        endif 
            
        if is_ica then begin
        dir_write = strmid(dir_upper,0,strlen(dir_upper)-1)+'-ica-tif'+info_slash
      endif else begin
        dir_write = strmid(dir_upper,0,strlen(dir_upper)-1)+'-pca-tif'+info_slash
      endelse
      if ~file_test(dir_write, /directory) then begin
        command = 'mkdir ' + dir_write;
        spawn, command
      endif

      dir_out_up = dir_write + cube + info_slash
      if ~file_test(dir_out_up, /directory) then begin
        command = 'mkdir ' + dir_out_up;
        spawn, command
      endif
      dir_out = dir_write + cube + info_slash + cube+'+'+'tif' + info_slash;
      if ~file_test(dir_out, /directory) then begin
        command = 'mkdir ' + dir_out;
        spawn, command
      endif
      
      ; load images 
      filename_current = filename_task[f]
      leng_filename_current = strlen(filename_current)
      suffix = strmid(filename_current, leng_filename_current -4,4)
      if strcmp(suffix, '.jpg') then begin
      endif else begin
      continue
      endelse
      task_cube = strmid(filename_current, leng_filename_current - 14, 2)
      task_mask = strmid(filename_current, leng_filename_current - 10, 2)
      task_band = strmid(filename_current, leng_filename_current - 6, 2)
      task_cube_int = uint(task_cube)
      task_mask_int = uint(task_mask)
      task_band_int = uint(task_band)
      
      set = lines[set_start[task_cube_int-1]]
      set_str = strmid(set, 4,2)
      bands = transpose(lines[set_start[task_cube_int-1]+1:set_end[task_cube_int-1]])   
     
      n_bands = set_num[task_cube_int-1]
      
      ; for each band
      ; read and add to multiband image I 
      for b =0,n_bands-1 do begin
        band_no_space = strreplacef(bands[b], ' ', '')
        filepath_image = dir_source + cube + band_no_space + '.tif'
        im = 0; clear memory 
        im = read_tiff(filepath_image)
        if (is_double) then begin
          im = double(im)
        endif else begin
          im = FLOAT(im) 
        endelse
        ;normalize image by calibration value 
        ;strfind, wavelength, band_no_space, INDEX=ix
        ;normval = double(normalization(ix))
        ;normval = normval[0]
        ;im = im[*,*]/normval
      
        if (b eq 0) then begin 
          imsize = size(im, /dimensions)
          ;dispx = floor(imsize[0]/2.5)
          ;dispy = floor(imsize[1]/2.5)
          I = 0; 
          if (is_double) then begin
            I = make_array(imsize[0],imsize[1],set_num[task_cube_int-1], /double)
          endif else begin         
            I = make_array(imsize[0],imsize[1],set_num[task_cube_int-1], /float)
          endelse
        endif 
        ;i_scl = bytscl(i);    
        ;dispx = floor(imsize[0]/10)
        ;dispy = floor(imsize[1]/10)
        ;i_disp = congrid(i_scl,dispx,dispy)
        ;window, 0, xsize=dispx, ysize=dispy
        ;tv, i_disp
        I[*,*,b] = im  
    endfor ; load each cube into image 
    ; Output is multiband image.
    im = 0; Clear memory 
    
    if is_ica then begin 
      filename_trans = dir_upper + cube + info_slash +cube+'+'+ 'envi' +info_slash + cube + '_ica' + '_c' + task_cube + '_m' + task_mask + '.trans'
    endif else begin 
      filename_trans = dir_upper + cube + info_slash +cube+'+'+ 'envi' +info_slash + cube + '_pca' + '_c' + task_cube + '_m' + task_mask + '.csv'     
    endelse    ; Save jpg images 
    J = 0; clear memory 
    J = reform(I,[imsize[0]*imsize[1], set_num[task_cube_int-1]])
    I = 0; Clear memory 
    J = double(J)
    ;filename_csv = strmid(filename_trans,0,strlen(filename_trans)-5) + 'csv'
    ;csv = read_csv(filename_csv) 
    if is_ica then begin 
      envi_check_save, /transform
      transform_read, filename_trans, forward_matrix=tx
    endif else begin 
      tx_struct = read_csv(filename_trans)
      ; Now must convert to array from structure
      tx = MAKE_ARRAY(n_bands, n_bands);
      for b = 0, n_bands-1 do begin
        tx[*,b] = tx_struct.(b)
      end
    endelse    ;write_csv_data, tx, filename=filename_csv
    ;print, format='(%"\t\t\tApplying ICA to full image")'
    J_tx = 0; clear memory 
    J_tx = tx##J  
    J = 0; clear memory 
    ICA = 0; clear memory 
    ICA = reform(J_tx,[imsize[0],imsize[1],set_num[task_cube_int-1]])
    minval = min(J_tx, dimension=1)
    maxval = max(J_tx, dimension=1);
    q = task_band_int-1
    ICA_q = ICA[*,*,q];
    ICA = 0; clear memory 
    ;dispx = floor(imsize[0]/2.5)
    ;dispy = floor(imsize[1]/2.5)
   

    
    
   ; for q_ix = 0,n_bands-1 do begin ; for each band 
        J_tx_q = J_tx[*,q];
        J_tx = 0;
        ;print, format='(%"\t\t\t\tSaving output band %d")', q+1
        pdf = histogram(J_tx_q,min=minval[q], max=maxval[q],nbins=500,locations=locations)
        cdf = total(pdf, /cumulative) / (imsize[0]*imsize[1])
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
        if (is_double) then begin
          ICA_tif = uint(cgscalevector(ICA_q, 0,65535, minvalue=val_low, maxvalue=val_hi, double=1))
        endif else begin
          ICA_tif = uint(cgscalevector(ICA_q, 0,65535, minvalue=val_low, maxvalue=val_hi, double=0))
        endelse
        ;ICA_jpg = byte(cgscalevector(ICA[*,*,q], 0,255, minvalue=val_low, maxvalue=val_hi, double=1))
                
        c_str = task_cube
        r_str = task_mask
        q_str = task_band
        
        ;filepath_tiff = dir_upper + cube + info_slash + 'tiff' + info_slash + cube + $
        ;  '_c' + c_Str + '_m' + r_str + '_b' + q_str +  '.tif'
        filepath_jpg = dir_out + cube + $
          '_c' + c_Str + '_m' + r_str + '_b' + q_str +  '.jpg'
        filepath_tif = dir_out + cube + $
          '_c' + c_Str + '_m' + r_str + '_b' + q_str +  '.tif'
       ; write_jpeg, filepath_jpg, ICA_jpg, order=1, quality=50
        write_tiff, filepath_tif, ICA_tif, /short
        ICA_tif = 0;
     ; endfor   ; band  
        
    endfor
endfor
  

RETURN, 1
end
