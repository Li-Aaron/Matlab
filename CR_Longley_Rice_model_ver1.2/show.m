function varargout = show(varargin)
% SHOW MATLAB code for show.fig
%      SHOW, by itself, creates a new SHOW or raises the existing
%      singleton*.
%
%      H = SHOW returns the handle to a new SHOW or the handle to
%      the existing singleton*.
%
%      SHOW('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOW.M with the given input arguments.
%
%      SHOW('Property','Value',...) creates a new SHOW or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before show_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to show_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help show

% Last Modified by GUIDE v2.5 10-Dec-2012 20:51:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @show_OpeningFcn, ...
                   'gui_OutputFcn',  @show_OutputFcn, ...
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


% --- Executes just before show is made visible.
function show_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to show (see VARARGIN)

% Choose default command line output for show
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes show wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = show_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_tx_x_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tx_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tx_x as text
%        str2double(get(hObject,'String')) returns contents of edit_tx_x as a double


% --- Executes during object creation, after setting all properties.
function edit_tx_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tx_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_tx_y_Callback(hObject, eventdata, handles)
% hObject    handle to edit_tx_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_tx_y as text
%        str2double(get(hObject,'String')) returns contents of edit_tx_y as a double


% --- Executes during object creation, after setting all properties.
function edit_tx_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_tx_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_rx_x_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rx_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rx_x as text
%        str2double(get(hObject,'String')) returns contents of edit_rx_x as a double


% --- Executes during object creation, after setting all properties.
function edit_rx_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rx_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_rx_y_Callback(hObject, eventdata, handles)
% hObject    handle to edit_rx_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_rx_y as text
%        str2double(get(hObject,'String')) returns contents of edit_rx_y as a double


% --- Executes during object creation, after setting all properties.
function edit_rx_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_rx_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_surf.
function pushbutton_surf_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_surf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global terrain
[N,M] = size(terrain);
surf(handles.axes1,terrain,'FaceColor','interp','EdgeColor','interp');
% set(handles.axes4,'Visible','on');
% plot(handles.axes4,1,1);%为了清除colorbar属性 实在没办法了
% colorbar(handles.axes4);

set(handles.axes1,'NextPlot','add');
Tx_height = str2double(get(handles.edit_Txh,'string'));
Rx_height = str2double(get(handles.edit_Rxh,'string'));
Tx_x = str2double(get(handles.edit_tx_x,'string'));   
Tx_y = str2double(get(handles.edit_tx_y,'string'));
Rx_x = str2double(get(handles.edit_rx_x,'string'));
Rx_y = str2double(get(handles.edit_rx_y,'string'));
Tx_Realheight = terrain(Tx_y+1,Tx_x+1)+Tx_height;
Rx_Realheight = terrain(Rx_y+1,Rx_x+1)+Rx_height;
Txshowx = [Tx_x Tx_x];
Txshowy = [Tx_y Tx_y];
Txshowz = [terrain(Tx_y+1,Tx_x+1) Tx_Realheight];
Rxshowx = [Rx_x Rx_x];
Rxshowy = [Rx_y Rx_y];
Rxshowz = [terrain(Rx_y+1,Rx_x+1) Rx_Realheight];
Connshowx = [Tx_x Rx_x];
Connshowy = [Tx_y Rx_y];
Connshowz = [Tx_Realheight Rx_Realheight];
plot3(handles.axes1,Txshowx,Txshowy,Txshowz,'LineWidth',2,...
    'Color',[.8 .0 .8],'Marker','x','MarkerSize',15);
plot3(handles.axes1,Rxshowx,Rxshowy,Rxshowz,'LineWidth',2,...
    'Color',[.0 .8 .0],'Marker','x','MarkerSize',15);
plot3(handles.axes1,Connshowx,Connshowy,Connshowz,'LineWidth',2,'Color',...
    [.0 .0 .0],'Marker','.','MarkerSize',1);
set(handles.axes1,'NextPlot','replace');
minz = min(min(terrain));
maxz = max(max(max(terrain))*1.2+0.01,max(Tx_Realheight,Rx_Realheight)*1.2+0.01);
axis(handles.axes1,[1 M 1 N minz maxz]);
set(handles.axes1,'Xcolor','r','Ycolor','b','Zcolor','k');

% --- Executes on button press in pushbutton_contourf.
function pushbutton_contourf_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_contourf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global terrain 
[N,M] = size(terrain);
contourf(handles.axes1,terrain);
% set(handles.axes4,'Visible','on');
% plot(handles.axes4,1,1);%为了清除colorbar属性 实在没办法了
% colorbar(handles.axes4);

set(handles.axes1,'NextPlot','add');
Tx_x = str2double(get(handles.edit_tx_x,'string'));   
Tx_y = str2double(get(handles.edit_tx_y,'string'));
Rx_x = str2double(get(handles.edit_rx_x,'string'));
Rx_y = str2double(get(handles.edit_rx_y,'string'));
Showx = [Tx_x Rx_x];
Showy = [Tx_y Rx_y];
plot(handles.axes1,Showx,Showy,'LineWidth',2,'Color',[0 0 0],'Marker','x','MarkerSize',15);
set(handles.axes1,'NextPlot','replace');
axis(handles.axes1,[1 M 1 N]);

% --- Executes on button press in pushbutton_contour3.
function pushbutton_contour3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_contour3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global terrain
[N,M] = size(terrain);
contour3(handles.axes1,terrain,15);
% set(handles.axes4,'Visible','on');
% plot(handles.axes4,1,1);%为了清除colorbar属性 实在没办法了
% colorbar(handles.axes4);

set(handles.axes1,'NextPlot','add');
Tx_height = str2double(get(handles.edit_Txh,'string'));
Rx_height = str2double(get(handles.edit_Rxh,'string'));
Tx_x = str2double(get(handles.edit_tx_x,'string'));   
Tx_y = str2double(get(handles.edit_tx_y,'string'));
Rx_x = str2double(get(handles.edit_rx_x,'string'));
Rx_y = str2double(get(handles.edit_rx_y,'string'));
Tx_Realheight = terrain(Tx_y+1,Tx_x+1)+Tx_height;
Rx_Realheight = terrain(Rx_y+1,Rx_x+1)+Rx_height;
Connshowx = [Tx_x Rx_x];
Connshowy = [Tx_y Rx_y];
Connshowz = [Tx_Realheight Rx_Realheight];
plot3(handles.axes1,Connshowx,Connshowy,Connshowz,'LineWidth',2,...
    'Color',[.0 .0 .0],'Marker','x','MarkerSize',15);
set(handles.axes1,'NextPlot','replace');
minz = min(min(terrain));
maxz = max(max(max(terrain))*1.2+0.01,max(Tx_Realheight,Rx_Realheight)*1.2+0.01);
axis(handles.axes1,[1 M 1 N minz maxz]);
set(handles.axes1,'Xcolor','r','Ycolor','b','Zcolor','k');



% --- Executes on button press in pushbutton_cal.
function pushbutton_cal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global terrain terrain_tx_rx lenth_of_flat flat_terrain
tx_x = str2double(get(handles.edit_tx_x,'string'));   
tx_y = str2double(get(handles.edit_tx_y,'string'));
rx_x = str2double(get(handles.edit_rx_x,'string'));
rx_y = str2double(get(handles.edit_rx_y,'string'));
[ flat_terrain, terrain_tx_rx, lenth_of_flat ] = cal_flat_terrain(tx_x,...
    tx_y,rx_x,rx_y,terrain);
[N,M] = size(terrain_tx_rx);
range_array = linspace(0,lenth_of_flat/100,max(M,N));
plot(handles.axes2,range_array,flat_terrain);
axis(handles.axes2,[0 lenth_of_flat/100 min(flat_terrain)-abs(min(flat_terrain))*0.2 ...
    max(flat_terrain)+abs(max(flat_terrain))*0.2]);
xlabel(handles.axes2,'distance in Km');
ylabel(handles.axes2,'height in meter');





% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global V_azimuth V_elevation
V_azimuth = get(hObject,'Value');
V_elevation = get(handles.slider2,'Value');
set(handles.axes1,'View',[ V_elevation V_azimuth ]);
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
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global V_azimuth V_elevation
V_azimuth = get(handles.slider1,'Value');
V_elevation = get(hObject,'Value');
set(handles.axes1,'View',[ V_elevation V_azimuth ]);

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global terrain
i = get(hObject,'value');
load terrain10
switch i
    case 1
        terrain = bjtu_terrain';
    case 2
        terrain = gauss_terrain';
    case 3
        terrain = peaks_terrain';
    case 4
        terrain = blade_terrain';
    case 5
        terrain = sharp_terrain';
    case 6
        terrain = hill_terrain';
end

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_freqmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_freqmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_freqmin as text
%        str2double(get(hObject,'String')) returns contents of edit_freqmin as a double


% --- Executes during object creation, after setting all properties.
function edit_freqmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_freqmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_freqmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_freqmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_freqmax as text
%        str2double(get(hObject,'String')) returns contents of edit_freqmax as a double


% --- Executes during object creation, after setting all properties.
function edit_freqmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_freqmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_TLC_type.
function pop_TLC_type_Callback(hObject, eventdata, handles)
% hObject    handle to pop_TLC_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_TLC_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_TLC_type


% --- Executes during object creation, after setting all properties.
function pop_TLC_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_TLC_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on button press in pushbutton_TLCal.
function pushbutton_TLCal_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_TLCal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global lenth_of_flat flat_terrain
TLC_type = get(handles.pop_TLC_type,'value');
Freq_min = str2double(get(handles.edit_freqmin,'string'));
Freq_max = str2double(get(handles.edit_freqmax,'string'));
Tx_height = str2double(get(handles.edit_Txh,'string'));
Rx_height = str2double(get(handles.edit_Rxh,'string'));
TGC_type = get(handles.pop_TGC,'value');
d = lenth_of_flat / 100; % distance in map 1 stand for 10m; d -->km
% ITM Model must uses at d > 1km
if d < 1
    msg= 'Distance must larger than 1 km';
    hs = msgbox(msg,'Error','error');
    ht = findobj(hs, 'Type', 'text');
    set(ht, 'FontSize', 12, 'Unit', 'normal');
elseif Freq_min < 20
    msg= 'Minimum frequency must larger than 20 MHz  00000';
    hs = msgbox(msg,'Error','error');
    ht = findobj(hs, 'Type', 'text');
    set(ht, 'FontSize', 12, 'Unit', 'normal');
elseif Freq_max > 40000
    msg= 'Maximum frequency must less than 40 GHz  00000';
    hs = msgbox(msg,'Error','error');
    ht = findobj(hs, 'Type', 'text');
    set(ht, 'FontSize', 12, 'Unit', 'normal');
else

Freq = logspace(log10(Freq_min),log10(Freq_max),60);

Lbf = cal_Lbf(Freq,d);
[Acr,func,loc1,loc2] = cal_Acr(flat_terrain,10,Freq,Tx_height,Rx_height,TGC_type);
switch func
    case 1
        set(handles.text_funbox,'string','Two-Ray Optics');
    case 2
        set(handles.text_funbox,'string','Diffraction Attenuation');
    case 3
        set(handles.text_funbox,'string','Scatter Attenuation');
end

if TLC_type == 1
    Plt = Lbf;
elseif TLC_type == 2
    Plt = Acr;
elseif TLC_type == 3
    Plt = Lbf + Acr;
end
semilogx(handles.axes3,Freq,Plt,'marker','o','markersize',5);
axis(handles.axes3,[Freq_min Freq_max min(Plt)-10 ...
     max(Plt)+10]);
xlabel(handles.axes3,'Frequency in MHz');
ylabel(handles.axes3,'Transmission Loss in dB');
set(handles.axes3,'XGrid','on','YGrid','on');
% change the axes2 plot

Tx_Realheight = flat_terrain(1)+Tx_height;
Rx_Realheight = flat_terrain(size(flat_terrain,2))+Rx_height;
S1_position = loc1*lenth_of_flat/size(flat_terrain,2)/100;
S2_position = loc2*lenth_of_flat/size(flat_terrain,2)/100;

range_array = linspace(0,lenth_of_flat/100,size(flat_terrain,2));
plot(handles.axes2,range_array,flat_terrain);

set(handles.axes2,'NextPlot','add');
plot(handles.axes2,[0,lenth_of_flat/100],[Tx_Realheight,Rx_Realheight],...
    'LineWidth',1,'Color','r','Marker','x','MarkerSize',12);%line between antennas
plot(handles.axes2,[0,S1_position],[Tx_Realheight,flat_terrain(loc1)],...
    'LineWidth',1,'Color',[.3 .3 .3],'Marker','x','MarkerSize',9);%line between antenna and terrain obstruct
plot(handles.axes2,[S2_position,lenth_of_flat/100],[flat_terrain(loc2),Rx_Realheight],...
    'LineWidth',1,'Color',[.3 .3 .3],'Marker','x','MarkerSize',9);%line between antenna and terrain obstruct
set(handles.axes2,'NextPlot','replace');
miny = min(flat_terrain)-abs(min(flat_terrain))*0.2;
maxy = max(max(flat_terrain)+abs(max(flat_terrain))*0.2,...
    max(Tx_Realheight,Rx_Realheight)+abs(max(Tx_Realheight,Rx_Realheight))*0.2);
axis(handles.axes2,[0 lenth_of_flat/100 miny maxy]);
legend(handles.axes2,'terrain','antennas','A & T');
xlabel(handles.axes2,'distance in Km');
ylabel(handles.axes2,'height in meter');
end






function edit_Txh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Txh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Txh as text
%        str2double(get(hObject,'String')) returns contents of edit_Txh as a double


% --- Executes during object creation, after setting all properties.
function edit_Txh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Txh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_Rxh_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Rxh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Rxh as text
%        str2double(get(hObject,'String')) returns contents of edit_Rxh as a double


% --- Executes during object creation, after setting all properties.
function edit_Rxh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Rxh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_TGC.
function pop_TGC_Callback(hObject, eventdata, handles)
% hObject    handle to pop_TGC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_TGC contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_TGC


% --- Executes during object creation, after setting all properties.
function pop_TGC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_TGC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
