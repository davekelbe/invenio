pro free_all
;+
; ROUTINE:  free_all
;
; PURPOSE:  free up all logical units
;
; USEAGE:   free_all
;
; OUTPUT:   none
;
; DISCUSSION:
;           Produces a list of all open files.  For each item in the
;           list the user may enter one of the following reponses:
;
;           quit -- do not close this file, quit free_all
;           view -- do not close this file, view list of remaining open files
;           all  -- close this file and all remaining open files
;           yes  -- close this file
;           no   -- do not close this file
;
;           Note: only the first letter of the response is actually used in
;           free_all.
;
; SIDE EFFECTS:
;           closes open files
;  
; EXAMPLE:  
;         IDL> openw,/get_lun,lun1,'junk1'
;         IDL> openw,/get_lun,lun2,'junk2'
;         IDL> openw,/get_lun,lun3,'junk3'
;         IDL> openw,/get_lun,lun4,'junk4'
;         IDL> openw,/get_lun,lun5,'junk5'
;         IDL> free_all
;         Close junk1, logical unit 101 (quit,view/all/yes/no)?: y
;         Close junk2, logical unit 102 (quit,view/all/yes/no)?: y
;         Close junk3, logical unit 103 (quit,view/all/yes/no)?: a 
;         Close junk4, logical unit 104 
;         Close junk5, logical unit 105 
;         
;
; AUTHOR:   Paul Ricchiazzi                        17 Mar 98
;           Institute for Computational Earth System Science
;           University of California, Santa Barbara
;           paul@icess.ucsb.edu
;
; REVISIONS:
;
;-
;
    
ans=''
for i=1, 128 do begin
  a=fstat(i)
  if a.open eq 1 then begin
    
    free_lun,i

  endif
endfor
end