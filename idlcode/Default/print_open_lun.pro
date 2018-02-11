pro print_open_lun
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
;        
;
;  
; EXAMPLE:  
;         
;         
;
; AUTHOR:   Dave Kelbe
;
; REVISIONS:
;
;-
;
    
ans=''
for i=1, 128 do begin
  a=fstat(i)
  if a.open eq 1 then begin
    
    print, format='(%"LUN %d open")', i
    

  endif
endfor
end