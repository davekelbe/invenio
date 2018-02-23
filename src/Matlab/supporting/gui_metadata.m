function varargout = gui_metadata(varargin)
% GUI_METADATA MATLAB code for gui_metadata.fig
%      GUI_METADATA, by itself, creates a new GUI_METADATA or raises the existing
%      singleton*.
% 
%      H = GUI_METADATA returns the handle to a new GUI_METADATA or the handle to
%      the existing singleton*.
%
%      GUI_METADATA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_METADATA.M with the given input arguments.
%
%      GUI_METADATA('Property','Value',...) creates a new GUI_METADATA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_metadata_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_metadata_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_metadata

% Last Modified by GUIDE v2.5 22-Jul-2013 09:35:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_metadata_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_metadata_OutputFcn, ...
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


% --- Executes just before gui_metadata is made visible.
function gui_metadata_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_metadata (see VARARGIN)
 
% Choose default command line output for gui_metadata
handles.output = hObject; 

% Update handles structure
guidata(hObject, handles);

set(handles.edit_program, 'Max',4)
set(handles.edit_description, 'Max',40)    
set(handles.edit_parent, 'Max',10)
set(handles.edit_comments, 'Max',5)
    
%set(gcf, 'KeyPressFcn',{@dispkeyevent,hObject,handles});
 handles.I = varargin{1};
 if size(handles.I,3)==1;
     set(handles.checksave,'Enable','off');
 end
 handles.info = varargin{2};
 handles.savename = varargin{3};
 handles.savedir = varargin{4};
 handles.refim = varargin{5};
  set(handles.edit_manuscript, 'String',handles.info.ObjectName);
  set(handles.edit_shelf, 'String',handles.info.Source)
  set(handles.edit_description, 'String',handles.info.DAT_File_Processing)
  set(handles.edit_parent, 'String',handles.info.ID_Parent_File)
  set(handles.edit_comments, 'String',handles.info.DAT_Processing_Comments)
  set(handles.edit_contrast, 'String',handles.info.DAT_Type_of_Contrast_Adjustment)
  set(handles.edit_creator, 'String',handles.info.Creator)
  set(handles.edit_program, 'String',handles.info.DAT_Processing_Program)
  set(handles.edit_savename,'String',handles.savename);
  set(handles.edit_savedir,'String',handles.savedir);
  %set(handles.pushrotate,'Value',0)
  %handles.nrot = 0;
 axes(handles.axes1);
 guidata(hObject, handles); 
% 
 imshow(handles.I);
 axis image
 axes(handles.axes2);
 imagesc(handles.refim);
 colormap('gray');
 axis image, axis off
foo =1;
  

  
% UIWAIT makes gui_metadata wait for user response (see UIRESUME)
 uiwait(handles.figure1);
 
% --- Outputs from this function are returned to the command line.
function varargout = gui_metadata_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) 
guidata(hObject, handles);
varargout{1} = get(handles.edit_manuscript,'string');
varargout{2} = get(handles.edit_shelf,'string');
varargout{3} = get(handles.edit_description,'string');
varargout{4} = get(handles.edit_parent,'string');
varargout{5} = get(handles.edit_comments,'string');
varargout{6} = get(handles.edit_contrast,'string');
varargout{7} = get(handles.edit_creator,'string');
varargout{8} = get(handles.edit_program,'string');
varargout{9} = get(handles.edit_savedir,'string');
varargout{10} = get(handles.edit_savename,'string');
k = strfind(varargout{10}, '.');
if numel(k)>0;
    varargout{10} = varargout{10}(1:k-1);
end
varargout{11} = get(handles.checksave,'value');
varargout{12} = get(handles.edit_nrot,'value');


delete(handles.figure1);

% --- Executes on button press in pushsave.
function pushsave_Callback(hObject, eventdata, handles)
% hObject    handle to pushsave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Choose output directory 
guidata(hObject, handles);
handles.savename = get(handles.edit_savename,'string');
handles.savedir = get(handles.edit_savedir,'string');

%{
cd(handles.savedir);
[handles.savename, handles.savedir] = uiputfile(handles.savename,'Please choose output directory');
handles.savename = handles.savename;
set(handles.edit_savename,'String',handles.savename);
set(handles.edit_savedir,'String',handles.savedir);
guidata(hObject, handles);
%}
%varargout = handles.info

%save(sprintf('%s/temp_on_off.mat',pwd),'on_off');
%save(sprintf('%s/temp_mask_bounds.mat',pwd),'mask_bounds');
%save(sprintf('%s/temp_mask_image.mat',pwd),'mask_image');
close all




function edit_shelf_Callback(hObject, eventdata, handles)
% hObject    handle to edit_shelf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_shelf as text
%        str2double(get(hObject,'String')) returns contents of edit_shelf as a double


% --- Executes during object creation, after setting all properties.
function edit_shelf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_shelf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_description_Callback(hObject, eventdata, handles)
% hObject    handle to edit_description (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_description as text
%        str2double(get(hObject,'String')) returns contents of edit_description as a double


% --- Executes during object creation, after setting all properties.
function edit_description_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_description (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_program_Callback(hObject, eventdata, handles)
% hObject    handle to edit_program (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_program as text
%        str2double(get(hObject,'String')) returns contents of edit_program as a double


% --- Executes during object creation, after setting all properties.
function edit_program_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_program (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_software_Callback(hObject, eventdata, handles)
% hObject    handle to edit_software (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_software as text
%        str2double(get(hObject,'String')) returns contents of edit_software as a double


% --- Executes during object creation, after setting all properties.
function edit_software_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_software (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_contrast_Callback(hObject, eventdata, handles)
% hObject    handle to edit_contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_contrast as text
%        str2double(get(hObject,'String')) returns contents of edit_contrast as a double


% --- Executes during object creation, after setting all properties.
function edit_contrast_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_contrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_creator_Callback(hObject, eventdata, handles)
% hObject    handle to edit_creator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_creator as text
%        str2double(get(hObject,'String')) returns contents of edit_creator as a double


% --- Executes during object creation, after setting all properties.
function edit_creator_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_creator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

 

function edit_comments_Callback(hObject, eventdata, handles)
% hObject    handle to edit_comments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_comments as text
%        str2double(get(hObject,'String')) returns contents of edit_comments as a double


% --- Executes during object creation, after setting all properties.
function edit_comments_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_comments (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_manuscript_Callback(hObject, eventdata, handles)
% hObject    handle to edit_manuscript (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_manuscript as text
%        str2double(get(hObject,'String')) returns contents of edit_manuscript as a double


% --- Executes during object creation, after setting all properties.
function edit_manuscript_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_manuscript (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_parent_Callback(hObject, eventdata, handles)
% hObject    handle to edit_parent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_parent as text
%        str2double(get(hObject,'String')) returns contents of edit_parent as a double


% --- Executes during object creation, after setting all properties.
function edit_parent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_parent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject,'waitstatus'),'waiting')
    uiresume(hObject);
else
    delete(hObject);
end



function edit_savename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_savename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_savename as text
%        str2double(get(hObject,'String')) returns contents of edit_savename as a double


% --- Executes during object creation, after setting all properties.
function edit_savename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_savename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checksave.
function checksave_Callback(hObject, eventdata, handles)
% hObject    handle to checksave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checksave
TF = get(hObject,'Value');    

 
function edit_savedir_Callback(hObject, eventdata, handles)
% hObject    handle to edit_savedir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_savedir as text
%        str2double(get(hObject,'String')) returns contents of edit_savedir as a double


% --- Executes during object creation, after setting all properties.
function edit_savedir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_savedir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushrotate.
function pushrotate_Callback(hObject, eventdata, handles)
% hObject    handle to pushrotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
guidata(hObject, handles);
handles.I = imrotate(handles.I,-90);
 axes(handles.axes1);
 imshow(handles.I);
 axes(handles.axes2);
 current_nrot= get(handles.edit_nrot,'Value');
 set(handles.edit_nrot,'Value',current_nrot+1);
 guidata(hObject, handles); 
 


function edit_nrot_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nrot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nrot as text
%        str2double(get(hObject,'String')) returns contents of edit_nrot as a double


% --- Executes during object creation, after setting all properties.
function edit_nrot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nrot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider7_Callback(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider8_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider9_Callback(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider10_Callback(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider11_Callback(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
