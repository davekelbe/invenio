PRO envi_extract_trans_file_noenvi, trans_filename
; Batch Processes PCA Images from the Sinai Palimpsests Project
; SYNTAX: 
;   Type in command line: IDL> pca, filename_cube, filename_mask, pos
;             will process filename_cube with mask filename_mask at band 
;             positions pos
; 
; VARIABLES: --------------------------------------------------------------- 
; 
;   FILENAME_CUBE = filename for ENVI cube image 
;   FILENAME_MASK = filename for ENVI mask image
;   POS  = the band subset 
;    
; OPTIONS: ---------------------------------------------------------------
;  
;   MAC        = Set to 0u if working on a PC computer; 1u if working on Mac or Unix
;                This changes the slash between \ and / path names
;   SESSION    = Set to a string referencing date of image collection 
;                 e.g, '2012-11' or '2011-5' for the December and May sessions
;                 This adjusts the bandlists appropriately 
;   DATADIR    = Set to a string 
;
; VERSION: ---------------------------------------------------------------
;   1.0
;   
;   Email bugs, comments, etc. to:
;   Dave Kelbe
;   djk2312@cis.rit.edu
; ---------------------------------------------------------------
 ;%%%%%%%%%%%%%   End Comments Section     %%%%%%%%%%%%%%%%%%%%%%%;
; ---------------------------------------------------------------
; 
;trans_filename = '/Volumes/Philippi/Processed3/0057_000029/envi/0057_000029_cube_norm_b_002_005_012_020_021_roi02.trans'

; Restore all the base save files
;envi, /restore_base_save_files 
 
; Initialize ENVI in the batch mode and send all errors and warnings to the file batch.txt
;envi_batch_init, log_file='batch.txt'

; restore one of the ENVI .sav files that has a reader for this file:
envi_check_save, /transform

; Call the undocumented routine 'transform_read' 
transform_read, trans_filename, forward_matrix=f, inverse_matrix=i
;out_filename = '/Volumes/Ephesus/Processed/0054_000013/envi/0054_000013_cube_norm_b_001_009_017_021_022_mask01.csv'
out_filename = strmid(trans_filename,0,strlen(trans_filename)-5)+'csv'

if (size(f, /n_elements) eq 0) then begin
    foo = 1;
endif else begin
    write_csv_data, f, filename=out_filename
  endelse 


END