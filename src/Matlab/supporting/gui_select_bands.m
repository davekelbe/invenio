function varargout = gui_select_bands(varargin)
% GUI_SELECT_BANDS MATLAB code for gui_select_bands.fig
%      GUI_SELECT_BANDS, by itself, creates a new GUI_SELECT_BANDS or raises the existing
%      singleton*. 
%
%      H = GUI_SELECT_BANDS returns the handle to a new GUI_SELECT_BANDS or the handle to
%      the existing singleton*.
%
%      GUI_SELECT_BANDS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_SELECT_BANDS.M with the given input arguments.
%
%      GUI_SELECT_BANDS('Property','Value',...) creates a new GUI_SELECT_BANDS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_select_bands_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_select_bands_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_select_bands

% Last Modified by GUIDE v2.5 14-Oct-2015 13:18:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @gui_select_bands_OpeningFcn, ...
    'gui_OutputFcn',  @gui_select_bands_OutputFcn, ...
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


% --- Executes just before gui_select_bands is made visible.
function gui_select_bands_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_select_bands (see VARARGIN)

% Choose default command line output for gui_select_bands


handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%set(gcf, 'KeyPressFcn',{@dispkeyevent,hObject,handles});
handles.b = 1;
handles.b_prev = 1;
handles.I = varargin{1};
handles.mapexists = false; 
handles.action = 'null';
handles.cubeno = 1;
handles.out_bandnames = [];

handles.nb = uint8(numel(handles.I));
handles.is3band = false(handles.nb,1);
for i = 1:handles.nb;
    [~,~,b] = size(handles.I{i});
    if b==3;
        handles.is3band(i) = true;
    end
end
handles.target_dir = varargin{2};
handles.cube = varargin{3};
handles.bandnames = varargin{4};

set(handles.listbox,'String',handles.bandnames)
set(handles.listbox,'Value',handles.b);

%set(handles.list_roitype,'String',{'Free','Poly','Rect'});
%set(handles.list_roitype,'value',uint8(3));

handles.m = 0;
handles.mask = logical([]);

handles.roi_tf = 1;

hanldles.mapexists = false;

handles.colormap = [...
    [255 0 0]; % red
    [0 255 0]; % green
    [0 0 255]; % blue
    [255 255 0]; % yellow
    [0 255 255]; % cyan
    [255 0 255]; % yellow
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

handles.map = zeros(handles.nr, handles.nc,3);
handles.map1 = handles.map(:,:,1);
handles.map2 = handles.map(:,:,2);
handles.map3 = handles.map(:,:,3);


handles.himg = imshow(handles.I{handles.b},[0,255]);
colormap('gray')
axis image

 
handles.xlim = get(handles.axes1,'xlim');
handles.ylim = get(handles.axes1,'ylim');
handles.xlimorig = get(handles.axes1,'xlim');
handles.ylimorig = get(handles.axes1,'ylim');
hold on;
handles.h_roi = plot(handles.xlim(1),handles.ylim(1)); 
handles.h_map = plot(handles.xlim(1),handles.ylim(1));
guidata(hObject, handles);


% UIWAIT makes gui_select_bands wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = gui_select_bands_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)%
%guidata(hObject, handles);
%delete(handles.figure1); 
 
% --- Executes on button press in pushsave.
function pushsave_Callback(hObject, eventdata, handles)
% hObject    handle to pushsave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Choose output directory
guidata(hObject, handles);  
   

% Read/create file 
handles.cube_dir = sprintf('%scube/', handles.target_dir(1:end-16));
if ~exist(handles.cube_dir);
    command = sprintf('mkdir %s', handles.cube_dir);
    [~,~] = system(command);
end
filename_cubetxt = sprintf('%scube.txt',handles.cube_dir);
if ~exist(filename_cubetxt,'file');
    fid = fopen(filename_cubetxt,'a+');
    fprintf(fid,sprintf('%s','Header.ncubes = 0'));
    fclose all;
end

% Gather last set of bands not added to cell
contents = cellstr(get(handles.listbox,'String'));
is_valid = false(handles.nb,1);
for f = 1:handles.nb;
    if ~isempty(strfind(contents{f},'green'));
        is_valid(f) = true;
    end
end
out_bandnames = handles.bandnames(is_valid);

% Determine whether to include or exclude last added band set 
include = true;
if numel(out_bandnames) == numel(handles.out_bandnames{end});
    n_bands = numel(out_bandnames);
    for b = 1:n_bands;
        if strcmp(out_bandnames{b},handles.out_bandnames{b})
            include = false;
        end
    end
end

% Add to existing set 
n_cube = numel(handles.out_bandnames);
if include;
    handles.out_bandnames{n_cube + 1} = out_bandnames;
end

% Determine max number of cubes 
fid = fopen(filename_cubetxt,'r');
content = textscan(fid, '%s', 'delimiter','\n');
content = content{1};
header = content{1}; 
ix = strfind(header, '=');
n_cubes = str2double(header(ix+1:end));
fclose all;

% Check to see if band set already exists 
A = regexp( fileread(filename_cubetxt), '\n', 'split')';
n_A = numel(A);
C = false(n_A,1);
for a = 1:n_A;
    if ~isempty(strfind(A{a},'Cube'));
        C(a) = true;
    end
end
ix = 1:n_A;
set_start = [ix(C) n_A];
set_end = set_start-1;
set_start = set_start(1:end-1);
set_end = set_end(2:end);
set_num = set_end - set_start -1;
n_set = numel(set_start);
set = cell(n_set,1);
for s = 1:n_set; 
   for i = 1:set_num(s);
    set{s}{i} = A{set_start(s)+i};
   end
   set{s} = set{s}';
end
n_new = numel(handles.out_bandnames);
is_duplicate = false(n_new,1);
for n = 1:n_new;
    n_curr = numel(handles.out_bandnames{n});
    for s = 1:n_set;
        % only continue if same length 
        if n_curr == set_num(s);
            % only continue if each band is identical 
           % is_match = false(n_curr,1);
            for b = 1:n_curr;
                if ~strcmp(handles.out_bandnames{n}{b}, set{s}{b});
                    break
                end
            end
            if b == n_curr;
                is_duplicate(n) = true;
            end
        end
    end
end

handles.out_bandnames = handles.out_bandnames(~is_duplicate);
    

% Write text file
fid = fopen(filename_cubetxt,'a+');
n_newcubes = numel(handles.out_bandnames);
for c = 1:n_newcubes;
    fprintf(fid, '\n\nCube%02.0f\n',n_cubes + c);
    n_bands = numel(handles.out_bandnames{c});
    for b = 1:n_bands-1;
        fprintf(fid, '%s\n',handles.out_bandnames{c}{b});
    end
    % Print last one without new line
    fprintf(fid, '%s',handles.out_bandnames{c}{n_bands});
end
fclose all;

%Update header
A = regexp( fileread(filename_cubetxt), '\n', 'split');
A{1} = sprintf('Header.ncubes = %02.0f \n',n_cubes + n_newcubes);
fid = fopen(filename_cubetxt, 'w');
fprintf(fid, '%s\n', A{:});
fclose(fid); 

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
delete(handles.figure1);


% --- Executes on selection change in listbox.
function listbox_Callback(hObject, eventdata, handles)
% hObject    handle to listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox
handles.b = get(handles.listbox,'value');
if handles.b ~=handles.b_prev;
    axes(handles.axes1);
    guidata(hObject, handles); 
    %delete(handles.himg); 
    handles.himg = imagesc(handles.I{handles.b},[0,255]);
    set(handles.axes1,'xlim', handles.xlim);
    set(handles.axes1,'ylim', handles.ylim);
    colormap('gray')
    if handles.roi_tf && handles.mapexists;
        hold on;
        handles.h_map = imshow(handles.map);
        set(handles.h_map,'AlphaData',handles.alpha);
    end
    handles.b_prev = handles.b;
end
if strcmp(handles.action,'add'); 
contents = cellstr(get(handles.listbox,'String'));
NewText = contents{get(handles.listbox,'Value')};
if ~isempty(strfind(NewText, '>'));
    NewText = handles.bandnames{get(handles.listbox,'Value')};
end
NewColor = sprintf('<HTML><FONT color="%s">%s', 'green', NewText);
NewContents = contents;
NewContents{get(handles.listbox,'Value')} = NewColor;
set(handles.listbox, 'String', NewContents)
end
if strcmp(handles.action,'remove'); 
contents = cellstr(get(handles.listbox,'String'));
NewText = contents{get(handles.listbox,'Value')};
if ~isempty(strfind(NewText, '>'));
    NewText = handles.bandnames{get(handles.listbox,'Value')};
end
NewColor = sprintf('<HTML><FONT color="%s">%s', 'red', NewText);
NewContents = contents;
NewContents{get(handles.listbox,'Value')} = NewColor;
set(handles.listbox, 'String', NewContents)
end
%set(handles.listbox,'Value',[])
% Update handles structure
%handles.action = 'null'; 
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
    case 'x';
        pushcancel_Callback(hObject, eventdata, handles);
    case 'c';
        pushsave_Callback(hObject, eventdata, handles);
    case 'o';  
        zoomok_Callback(hObject, eventdata, handles);
    case 'downarrow';
        handles.action = 'null'; 
        listval = get(handles.listbox,'value');
        listval = min(handles.nb,listval+1);
        set(handles.listbox,'value',listval);
        guidata(hObject, handles); 
        listbox_Callback(hObject, eventdata, handles)
        handles.b = listval;
        handles.b_prev = listval;
        guidata(hObject, handles);
    case 'uparrow'; 
        handles.action = 'null'; 
        listval = get(handles.listbox,'value');
        listval = max(1,listval-1);
        set(handles.listbox,'value',listval);
        guidata(hObject, handles);
        listbox_Callback(hObject, eventdata, handles)
        handles.b = listval;
        handles.b_prev = listval;
        guidata(hObject, handles);
    case 'leftarrow';
        handles.action = 'remove';
        listbox_Callback(hObject, eventdata, handles)
        guidata(hObject, handles);
    case 'rightarrow';
        handles.action = 'add';
        listbox_Callback(hObject, eventdata, handles)
        guidata(hObject, handles);
    case 'b';
        newbands_Callback(hObject, eventdata, handles)
        handles.out_bandnames = get(handles.listbox,'UserData');
        guidata(hObject, handles);
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
eventdata.Key = 'downarrow';
figure1_KeyPressFcn(hObject, eventdata, handles)



% --- Executes on button press in newbands.
function newbands_Callback(hObject, eventdata, handles)
% hObject    handle to newbands (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Reset all to null 
% Save 'gree' bands
handles.action = 'null';
contents = cellstr(get(handles.listbox,'String'));
is_valid = false(handles.nb,1);
for f = 1:handles.nb;
    if ~isempty(strfind(contents{f},'green'));
        is_valid(f) = true;
    end
end
out_bandnames = handles.bandnames(is_valid);
handles.out_bandnames{handles.cubeno} = out_bandnames;
handles.cubeno = handles.cubeno + 1;
set(handles.listbox,'string', handles.bandnames);
set(handles.listbox,'UserData', handles.out_bandnames);
 
guidata(hObject, handles);

   
   

% --- Executes on button press in addbutton.
function addbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
       handles.action = 'add';
        listbox_Callback(hObject, eventdata, handles)
        guidata(hObject, handles);


% --- Executes on button press in removebutton.
function removebutton_Callback(hObject, eventdata, handles)
% hObject    handle to removebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        handles.action = 'remove';
        listbox_Callback(hObject, eventdata, handles)
        guidata(hObject, handles);

