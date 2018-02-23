function varargout = gui_select_roi2(varargin)
% GUI_SELECT_ROI2 MATLAB code for gui_select_roi2.fig
%      GUI_SELECT_ROI2, by itself, creates a new GUI_SELECT_ROI2 or raises the existing
%      singleton*.
%
%      H = GUI_SELECT_ROI2 returns the handle to a new GUI_SELECT_ROI2 or the handle to
%      the existing singleton*.
%
%      GUI_SELECT_ROI2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SELECT_ROI2.M with the given input arguments.
%
%      GUI_SELECT_ROI2('Property','Value',...) creates a new GUI_SELECT_ROI2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_select_roi2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_select_roi2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help gui_select_roi2

% Last Modified by GUIDE v2.5 28-Aug-2014 20:28:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_select_roi2_OpeningFcn, ...
    'gui_OutputFcn',  @gui_select_roi2_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before gui_select_roi2 is made visible.
function gui_select_roi2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_select_roi2 (see VARARGIN)

% Choose default command line output for gui_select_roi2 
 
handles.kill = false;

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%set(gcf, 'KeyPressFcn',{@dispkeyevent,hObject,handles});
handles.b = 1;
handles.I = varargin{1}; 
 
handles.nb = uint8(numel(handles.I));
handles.is3band = false(handles.nb,1);
for i = 1:handles.nb;
    [~,~,b] = size(handles.I{i});
    if b==3; 
        handles.is3band(i) = true;
    end
end
handles.target_dir = varargin{2};
handles.bandnames = varargin{3};
handles.aux = varargin{4};

set(handles.listbox,'String',handles.bandnames)
set(handles.listbox,'Value',handles.b);

set(handles.list_roitype,'String',{'Free','Poly','Rect'});
set(handles.list_roitype,'value',uint8(1));

handles.m = 0;
handles.mask = logical([]);

set(handles.hideroi,'value',0);
handles.roi_tf = false;

handles.mapexists = false;

handles.colormap = [...
    [255 255 255]; % pink
    [0 0 255]; % blue
    [0 255 0]; % green
    [255 0 0]; % red
    [255 255 0]; % yellow
    [0 255 255]; % cyan
    [255 0 255]; % white??
    [176 48 96]; % yellow
    [46 139 187]; % yellow
    [160 32 240]; % yellow
    [255 127 80]]; % yellow

%set(handles.edit_savedir,'String',handles.savedir);
%set(handles.pushrotate,'Value',0)
%handles.nrot = 0;
axes(handles.axes1);
guidata(hObject, handles);
%handles.xlim = get(handles.axes1,'xlim');
%handles.ylim = get(handles.axes1,'ylim');
%[handles.nr, handles.nc, handles.nb] = size(handles.I);
%scale = round(max(handles.nc/1440, handles.nr/900));
[handles.nr, handles.nc, ~] = size(handles.I{1});

handles.map = uint8(zeros(handles.nr, handles.nc,3));
handles.map1 = handles.map(:,:,1);
handles.map2 = handles.map(:,:,2);
handles.map3 = handles.map(:,:,3);


handles.himg = imshow(handles.I{handles.b},[0,255]);
colormap('gray')
axis image

%
%handles.tif_dir = sprintf('%stif%s', handles.target_dir(1:end-4), handles.target_dir(end));
filepath_parchment_mask = sprintf('%s%s_parchment_mask.tif', handles.aux.path_tiff_dir{handles.aux.m},handles.aux.m_name{handles.aux.m});
filepath_overtext_mask = sprintf('%s%s_overtext_mask.tif', handles.aux.path_tiff_dir{handles.aux.m},handles.aux.m_name{handles.aux.m}); 
filepath_undertext_mask = sprintf('%s%s_undertext_mask.tif', handles.aux.path_tiff_dir{handles.aux.m},handles.aux.m_name{handles.aux.m}); 
mask = imread(filepath_parchment_mask);
if numel(size(mask)) > 1
    mask = mask(:,:,1);
    imwrite(mask, filepath_parchment_mask);
end
if ~islogical(mask)
    mask = logical(mask./max(mask));
    imwrite(mask, filepath_parchment_mask);
end
mask = imresize(mask,0.4);
handles.m = handles.m + 1;
handles.mask = cat(3,handles.mask,mask);
% axes(handles.axes1);
% handles.map1(mask) = handles.colormap(handles.m,1);
% handles.map2(mask) = handles.colormap(handles.m,2);
% handles.map3(mask) = handles.colormap(handles.m,3);
% handles.map = cat(3, handles.map1,handles.map2,handles.map3);
% handles.ismap = any((handles.map),3);
% handles.alpha = zeros(handles.nr, handles.nc)+0.1;
% handles.alpha(~handles.ismap) = 0;
% hold on;
% handles.h_map = imshow(handles.map);
% set(handles.h_map,'AlphaData',handles.alpha);
% handles.mapexists = true;

% Update handles structure 
%guidata(hObject, handles);
%uiwait(handles.figure1);
%handles.map = uint8(zeros(handles.nr, handles.nc,3));
%handles.map1 = handles.map(:,:,1);
%handles.map2 = handles.map(:,:,2);
%handles.map3 = handles.map(:,:,3);
mask = imread(filepath_overtext_mask);
if numel(size(mask)) > 1
    mask = mask(:,:,1);
    imwrite(mask, filepath_parchment_mask);
end
if ~islogical(mask)
    mask = logical(mask./max(mask));
    imwrite(mask, filepath_overtext_mask);
end
mask = imresize(mask,0.4);
handles.m = handles.m + 1;
axes(handles.axes1);
handles.mask = cat(3,handles.mask,mask);
handles.map1(mask) = handles.colormap(handles.m,1);
handles.map2(mask) = handles.colormap(handles.m,2);
handles.map3(mask) = handles.colormap(handles.m,3);
handles.map = cat(3, handles.map1,handles.map2,handles.map3);
handles.ismap = any((handles.map),3);
handles.alpha = zeros(handles.nr, handles.nc)+0.25;
handles.alpha(~handles.ismap) = 0;
%hold on;
%handles.h_map = imshow(handles.map);
%set(handles.h_map,'AlphaData',handles.alpha);
%handles.mapexists = true;

%handles.map = uint8(zeros(handles.nr, handles.nc,3));
%handles.map1 = handles.map(:,:,1);
%handles.map2 = handles.map(:,:,2);
%handles.map3 = handles.map(:,:,3);
mask = imread(filepath_undertext_mask);
if numel(size(mask)) > 1
    mask = mask(:,:,1);
    imwrite(mask, filepath_parchment_mask);
end
if ~islogical(mask)
    mask = logical(mask./max(mask));
    imwrite(mask, filepath_undertext_mask);
end
mask = imresize(mask,0.4);
handles.m = handles.m + 1;
axes(handles.axes1);
handles.mask = cat(3,handles.mask,mask);
handles.map1(mask) = handles.colormap(handles.m,1);
handles.map2(mask) = handles.colormap(handles.m,2);
handles.map3(mask) = handles.colormap(handles.m,3);
handles.map = cat(3, handles.map1,handles.map2,handles.map3);
handles.ismap = any((handles.map),3);
handles.alpha = zeros(handles.nr, handles.nc)+0.25;
handles.alpha(~handles.ismap) = 0;
hold on;
handles.h_map = imshow(handles.map);
set(handles.h_map,'AlphaData',handles.alpha);
handles.mapexists = true;


% handles.mask = cat(3,handles.mask,mask);
% delete(h_roi);
% handles.map1(mask) = handles.colormap(handles.m,1);
% handles.map2(mask) = handles.colormap(handles.m,2);
% handles.map3(mask) = handles.colormap(handles.m,3);
handles.mapexists = true;
handles.roi_tf = true;
%}

handles.xlim = get(handles.axes1,'xlim');
handles.ylim = get(handles.axes1,'ylim');
handles.xlimorig = get(handles.axes1,'xlim');
handles.ylimorig = get(handles.axes1,'ylim');
hold on;
%handles.h_roi = plot(handles.xlim(1),handles.ylim(1));
%handles.h_map = plot(handles.xlim(1),handles.ylim(1));
guidata(hObject, handles);

% UIWAIT makes gui_select_roi2 wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function iskill = gui_select_roi2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)%
guidata(hObject, handles);
%delete(handles.figure1);
%iskill = true;%handles.kill;

% --- Executes on button press in pushsave.
function pushsave_Callback(hObject, eventdata, handles)
% hObject    handle to pushsave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Choose output directory
guidata(hObject, handles);

if ispc();
    info_slash = '\';
else 
    info_slash = '/';
end
% Save mask images
envi_dir = sprintf('%senvi%s',handles.target_dir(1:end-4), info_slash);
ix_slash = strfind(handles.target_dir, '/');
cube = handles.target_dir(ix_slash(end-2)+1:ix_slash(end-1)-1);
%cube = handles.target_dir(end-15:end-5);

cd(envi_dir);
D = dir('*mask*tif');
n_d = numel(D);
current_masks = zeros(n_d,1);
for d = 1:n_d;
    if strcmp(D(d).name(1),'.');
        continue
    end
    current_masks(d) = str2double(D(d).name(17:18));
end
max_mask = max(current_masks);
max_mask = max([max_mask; 0]);
 

% copy first three masks from tif folder 
filepath_source = sprintf('%s%s_parchment_mask.tif', ...
    handles.aux.path_tiff_dir{handles.aux.m}, handles.aux.m_name{handles.aux.m});
mask_parchment = imread(filepath_source);
filepath_target = sprintf('%s%s_mask01.tif', envi_dir, cube);
command = sprintf('cp %s %s', filepath_source, filepath_target);
[~,~] = system(command);
filepath_source = sprintf('%s%s_overtext_mask.tif', ...
    handles.aux.path_tiff_dir{handles.aux.m}, handles.aux.m_name{handles.aux.m});
filepath_target = sprintf('%s%s_mask02.tif', envi_dir, cube);
command = sprintf('cp %s %s', filepath_source, filepath_target);
[~,~] = system(command);
filepath_source = sprintf('%s%s_undertext_mask.tif', ...
    handles.aux.path_tiff_dir{handles.aux.m}, handles.aux.m_name{handles.aux.m});
filepath_target = sprintf('%s%s_mask03.tif', envi_dir, cube);
command = sprintf('cp %s %s', filepath_source, filepath_target);
[~,~] = system(command);

% Determine resize factor
length_orig = size(mask_parchment,2);
width_orig = size(mask_parchment,1);
dims_curr = size(handles.mask(:,:,1));
length_curr = max(dims_curr);
width_curr = min(dims_curr);
[~,ix] = sort(dims_curr);
reference_dims = [width_orig, length_orig];
reference_dims = reference_dims(ix);

resize_l = length_orig/length_curr;
resize_w = width_orig/width_curr; 
resize_lw = max(resize_l, resize_w);


n_m = size(handles.mask,3);
for m = 4:n_m;
    filepath_mask = sprintf('%s%s_mask%02.0f.tif',envi_dir,cube,max_mask + m);
    M = imresize(handles.mask(:,:,m),resize_lw);
    M = M(1:reference_dims(1),1:reference_dims(2));
    %M(~mask_parchment) = 0;
    imwrite(M,filepath_mask, 'tiff');
end

% Save alpha image RGB 
% Save alpha image
J1 = handles.I{1}(:,:,1);
J2 = handles.I{1}(:,:,2);
J3 = handles.I{1}(:,:,3);
handles.colormap01 = handles.colormap./255;

wildcard = 'W365O22';
%cd(handles.aux.path_tiff_dir{handles.aux.m});
is_valid = cellfun(@(x) contains(x,wildcard), handles.aux.m_wavelength{handles.aux.m});
if sum(is_valid) == 0
    is_valid = cellfun(@(x) contains(x,'WBUVO22'), handles.aux.m_wavelength{handles.aux.m});
end
ix_valid = find(is_valid);
ix_valid = ix_valid(1);
%wavelength = handles.aux.m_wavelength_file_new;
filepath_I = sprintf('%s',handles.aux.path_tiff_dir{handles.aux.m},handles.aux.m_wavelength_file_new{handles.aux.m}{ix_valid});
I1 = imresize(double(imread(filepath_I))./65535, size(J1));
I1 = I1./max(I1(:));
I2 = I1;
I3 = I1;

D = dir('*_F.tif');
w_wavelength = remove_hiddenfiles(D);
n_w = numel(w_wavelength);
for w = 1:n_w
    ix_delimiter = strfind(w_wavelength{w},options_delimiter_wavelength);
    w_wavelength{w} = w_wavelength{w}(ix_delimiter:end);
end


for m = 2:n_m;
    I1(handles.mask(:,:,m)) = .75*I1(handles.mask(:,:,m)) + 0.25*(handles.colormap01(m,1));
    I2(handles.mask(:,:,m)) = .75*I2(handles.mask(:,:,m)) + 0.25*(handles.colormap01(m,2));
    I3(handles.mask(:,:,m)) = .75*I3(handles.mask(:,:,m)) + 0.25*(handles.colormap01(m,3));

    J1(handles.mask(:,:,m)) = .75*J1(handles.mask(:,:,m)) + 0.25*(handles.colormap01(m,1));
    J2(handles.mask(:,:,m)) = .75*J2(handles.mask(:,:,m)) + 0.25*(handles.colormap01(m,2));
    J3(handles.mask(:,:,m)) = .75*J3(handles.mask(:,:,m)) + 0.25*(handles.colormap01(m,3));
end
J = cat(3,J1,J2,J3);
I = cat(3,I1,I2,I3);

D = dir('*map*');
n_d = numel(D);
current_maps = zeros(n_d,1);
for d = 1:n_d;
    if strcmp(D(d).name(1),'.');
        continue
    end
    current_maps(d) = str2double(D(d).name(16:17));
end
max_map = max([current_maps;0]);
 
filepath_map = sprintf('%s%s_mapTr%02.0f.tif',envi_dir,cube,max_map + 1);
imwrite(uint8(J*255),filepath_map, 'jpeg','quality',80);
filepath_map = sprintf('%s%s_mapTrGray%02.0f.tif',envi_dir,cube,max_map + 1);
imwrite(uint8(I*255),filepath_map, 'jpeg','quality',80);

% Save alpha image KTK 
% J1 = handles.I{2}(:,:,1);
% J2 = handles.I{2}(:,:,2);
% J3 = handles.I{2}(:,:,3);
% for m = 1:n_m;
%     J1(handles.mask(:,:,m)) = .75*J1(handles.mask(:,:,m)) + 0.25*(handles.colormap01(m,1));
%     J2(handles.mask(:,:,m)) = .75*J2(handles.mask(:,:,m)) + 0.25*(handles.colormap01(m,2));
%     J3(handles.mask(:,:,m)) = .75*J3(handles.mask(:,:,m)) + 0.25*(handles.colormap01(m,3));
% end
% J = cat(3,J1,J2,J3);
% 
% D = dir('*map*');
% n_d = numel(D);
% current_maps = zeros(n_d,1);
% for d = 1:n_d;
%     if strcmp(D(d).name(1),'.');
%         continue
%     end
%     current_maps(d) = str2double(D(d).name(16:17));
% end
% max_map = max([current_maps;0]);
%  
% filepath_map = sprintf('%s%s_mapPs%02.0f.tif',envi_dir,cube,max_map + 1);
% imwrite(uint8(J*255),filepath_map, 'jpeg','quality',80);


guidata(hObject, handles); 
%varargout = handles.info 
delete(handles.figure1);
 

% --- Executes on button press in zoom.
function zoom_Callback(hObject, eventdata, handles)
% hObject    handle to zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zoom
% if get(handles.zoom,'Value') == 1;
%     zoom on;
% end
% if get(handles.zoom,'Value') == 0;
%     zoom off;
%     handles.xlim = get(handles.axes1,'xlim');
%     handles.ylim = get(handles.axes1,'ylim');
%     
%     % Update handles structure 
%     guidata(hObject, handles);
% end
    zoom on;
    guidata(hObject, handles);

  
% --- Executes on button press in pushresetzoom. % Reset Zoom
function pushresetzoom_Callback(hObject, eventdata, handles)
% hObject    handle to pushresetzoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.axes1,'xlim',handles.xlimorig);
set(handles.axes1,'ylim',handles.ylimorig);
handles.xlim = handles.xlimorig;
handles.ylim = handles.ylimorig;
guidata(hObject, handles);
% axes(handles.axes1);
% [handles.nr, handles.nc, handles.nb] = size(handles.I);
% temp = reshape(handles.I, handles.nr*handles.nc,handles.nb)*...
%     handles.adjust;
% handles.icurrent = reshape(temp, handles.nr,handles.nc);
% imagesc(handles.icurrent);
% colormap('gray')
% axis image

% --- Executes on button press in pushcancel.
function pushcancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushcancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%status = true;
guidata(hObject, handles);

close all force
return
%exit 
%close(handles.figure1);
%delete(handles.figure1);


% --- Executes on button press in pushaddtoroi.
function pushaddtoroi_Callback(hObject, eventdata, handles)
% hObject    handle to pushaddtoroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.m==0;
    return
end
if handles.roi_tf;
    delete(handles.h_map);
end
switch get(handles.list_roitype,'value');
    case 1
        h_roi = imfreehand;
        mask = createMask(h_roi,handles.himg); 
    case 2
        h_roi = impoly;
        mask = createMask(h_roi,handles.himg);
    case 3
        h_roi = imrect;
        mask = createMask(h_roi,handles.himg);
end
handles.mask(:,:,handles.m) = or(mask, handles.mask(:,:,handles.m));
delete(h_roi);
handles.map1(handles.mask(:,:,handles.m)) = handles.colormap(handles.m,1);
handles.map2(handles.mask(:,:,handles.m)) = handles.colormap(handles.m,2);
handles.map3(handles.mask(:,:,handles.m)) = handles.colormap(handles.m,3);
if handles.roi_tf
    handles.map = cat(3, handles.map1,handles.map2,handles.map3);
    handles.ismap = any((handles.map),3);
    handles.alpha = zeros(handles.nr, handles.nc)+0.25;
    handles.alpha(~handles.ismap) = 0;
    hold on;
    handles.h_map = imshow(handles.map);
    set(handles.h_map,'AlphaData',handles.alpha);
end
 
% Update handles structure
guidata(hObject, handles);
uiwait(handles.figure1);

  
% --- Executes on button press in pushnewroi.
function pushnewroi_Callback(hObject, eventdata, handles)
% hObject    handle to pushnewroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.m>0&&handles.roi_tf;
    delete(handles.h_map);
end
handles.m = handles.m + 1;
axes(handles.axes1);
switch get(handles.list_roitype,'value');
    case 1 
        h_roi = imfreehand;
        mask = createMask(h_roi,handles.himg);
    case 2
        h_roi = impoly;
        mask = createMask(h_roi,handles.himg);
    case 3
        h_roi = imrect;
        mask = createMask(h_roi,handles.himg);
end
if sum(mask(:)) > 10
    handles.mask = cat(3,handles.mask,mask);
    delete(h_roi);
    handles.map1(mask) = handles.colormap(handles.m,1);
    handles.map2(mask) = handles.colormap(handles.m,2);
    handles.map3(mask) = handles.colormap(handles.m,3);
    if handles.roi_tf;
        handles.map = cat(3, handles.map1,handles.map2,handles.map3);
        handles.ismap = any((handles.map),3);
        handles.alpha = zeros(handles.nr, handles.nc)+0.25;
        handles.alpha(~handles.ismap) = 0;
        hold on;
        handles.h_map = imshow(handles.map);
        set(handles.h_map,'AlphaData',handles.alpha);
    end
    handles.mapexists = true;
else
    handles.h_map = imshow(handles.map);
    set(handles.h_map,'AlphaData',handles.alpha);
    handles.m = handles.m - 1;
end

% Update handles structure 
guidata(hObject, handles);
uiwait(handles.figure1);


% --- Executes on selection change in listbox.
function listbox_Callback(hObject, eventdata, handles)
% hObject    handle to listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 
% Hints: contents = cellstr(get(hObject,'String')) returns listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox

if get(handles.listbox, 'value') > handles.nb;
    set(handles.listbox, 'value', 1);
end
if get(handles.listbox, 'value') < 1;
    set(handles.listbox, 'value', handles.nb);
end
handles.b = get(handles.listbox,'value');

axes(handles.axes1);
guidata(hObject, handles);
delete(handles.himg);
handles.himg = imagesc(handles.I{handles.b},[0,1]);
set(handles.axes1,'xlim', handles.xlim);
set(handles.axes1,'ylim', handles.ylim);
colormap('gray')
if handles.roi_tf && handles.mapexists;
hold on;
handles.h_map = imshow(handles.map);
set(handles.axes1,'xlim', handles.xlim);
set(handles.axes1,'ylim', handles.ylim);
set(handles.h_map,'AlphaData',handles.alpha);
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in list_roitype.
function list_roitype_Callback(hObject, eventdata, handles)
% hObject    handle to list_roitype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns list_roitype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from list_roitype
handles.roitype = get(handles.list_roitype,'value');


% --- Executes during object creation, after setting all properties.
function list_roitype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to list_roitype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end 


% --- Executes on button press in showroi.
function showroi_Callback(hObject, eventdata, handles)
% hObject    handle to showroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%guidata(hObject, handles);
if ~handles.roi_tf
axes(handles.axes1);
handles.map = cat(3, handles.map1,handles.map2,handles.map3);
handles.ismap = any((handles.map),3);
handles.alpha = zeros(handles.nr, handles.nc)+0.25;
handles.alpha(~handles.ismap) = 0;
hold on;
handles.h_map = imshow(handles.map);
set(handles.h_map,'AlphaData',handles.alpha);
handles.roi_tf = true;
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in hideroi.
function hideroi_Callback(hObject, eventdata, handles)
% hObject    handle to hideroi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%handles.roi_tf = get(handles.hideroi,'value');
if handles.roi_tf;
    delete(handles.h_map);
end 
% Update handles structure
handles.roi_tf = false;
guidata(hObject, handles);

 
% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
switch eventdata.Key      
    case 'z'  
        zoom_Callback(hObject, eventdata, handles);
    case 'r';
        pushresetzoom_Callback(hObject, eventdata, handles);
    case 'n';
        pushnewroi_Callback(hObject, eventdata, handles);
    case 'a';
        pushaddtoroi_Callback(hObject, eventdata, handles);
    case 'h'; 
         hideroi_Callback(hObject, eventdata, handles);
    case 's';   
        showroi_Callback(hObject, eventdata, handles);
    case 'leftarrow'; 
         hideroi_Callback(hObject, eventdata, handles);
    case 'rightarrow'; 
        showroi_Callback(hObject, eventdata, handles);
    case 'x';
        pushcancel_Callback(hObject, eventdata, handles);
    case 'c';
        pushsave_Callback(hObject, eventdata, handles);  
    case 'o';
        zoomok_Callback(hObject, eventdata, handles); 
    case 'downarrow';
        handles.listval = get(handles.listbox,'value');
        handles.listval = min(handles.nb,handles.listval+1);
     %   set(handles.listbox,'value',handles.listval);
        guidata(hObject, handles);
        listbox_Callback(hObject, eventdata, handles)
    case 'uparrow';
        handles.listval = get(handles.listbox,'value');
        handles.listval = max(1,handles.listval-1);
    %    set(handles.listbox,'value',handles.listval);
        guidata(hObject, handles);
        listbox_Callback(hObject, eventdata, handles)
end
 

 
% --- Executes on button press in zoomok.
function zoomok_Callback(hObject, eventdata, handles)
% hObject    handle to zoomok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    zoom off;
    handles.xlim = get(handles.axes1,'xlim');
    handles.ylim = get(handles.axes1,'ylim');
    guidata(hObject, handles);
    %eventdata.Key = 'n';
    %figure1_KeyPressFcn(hObject, eventdata, handles)
    
