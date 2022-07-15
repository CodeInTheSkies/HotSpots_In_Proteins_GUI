%% Main Function
function varargout = HSpot_FSlider_GUI(varargin)
%***************************************************************************************
% This software accompanies the paper "Filter-Based Methodology for the Location of
% Hot Spots in Proteins and Exons in DNA" by Ramachandran et al., "IEEE Transactions
% on Biomedical Engineering", Volume: 59, Issue: 6, June 2012.
% 
% Copyright (c) 2012 by P. Ramachandran.  All rights reserved.
% 
% This software may be used as is by individual researchers to carry out research but
% citation to the above paper would be expected.  For any other use, the permission
% of the author, Dr. P. Ramachandran (Email: rpara26@gmail.com), must be obtained.
% This is a demonstration software that comes with no warranty expressed or implied.
%***************************************************************************************

% HSPOT_FSLIDER_GUI M-file for HSpot_FSlider_GUI.fig
%      HSPOT_FSLIDER_GUI, by itself, creates a new HSPOT_FSLIDER_GUI or raises the existing
%      singleton*.
%
%      H = HSPOT_FSLIDER_GUI returns the handle to a new HSPOT_FSLIDER_GUI or the handle to
%      the existing singleton*.
%
%      HSPOT_FSLIDER_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HSPOT_FSLIDER_GUI.M with the given input arguments.
%
%      HSPOT_FSLIDER_GUI('Property','Value',...) creates a new HSPOT_FSLIDER_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before HSpot_FSlider_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HSpot_FSlider_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HSpot_FSlider_GUI

% Last Modified by GUIDE v2.5 02-Oct-2012 19:01:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HSpot_FSlider_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @HSpot_FSlider_GUI_OutputFcn, ...
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

%% HSpot_FSlider_GUI Opening Function
% --- Executes just before HSpot_FSlider_GUI is made visible.
function HSpot_FSlider_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HSpot_FSlider_GUI (see VARARGIN)

% Choose default command line output for HSpot_FSlider_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes HSpot_FSlider_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
hMainGui     = getappdata(0,        'hMainGui');
HSGuiHandles = getappdata(hMainGui, 'HSGuiHandles');
[alp, xn, varcol_1_1] = anotchfreqvariation_for_FreqVarGUI(HSGuiHandles.w0, HSGuiHandles.OmegaL,...
    HSGuiHandles.OmegaU, 10);
setappdata(hMainGui, 'alp', alp);
setappdata(hMainGui, 'xn', xn);
setappdata(hMainGui, 'varcol_1_1', varcol_1_1);


%% HSpot_FSlider_GUI Output Function
% --- Outputs from this function are returned to the command line.
function varargout = HSpot_FSlider_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% Slider2 Callback
% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
currfrq = get(hObject,'Value');
set(handles.text4, 'string', num2str(currfrq));
updateFigs(currfrq);


%% Slider2 Create Function
% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

hMainGui = getappdata(0, 'hMainGui');
HSGuiHandles = getappdata(hMainGui, 'HSGuiHandles');

set(hObject, 'max',   HSGuiHandles.OmegaU);
set(hObject, 'min',   HSGuiHandles.OmegaL);
set(hObject, 'value', HSGuiHandles.w0);


%% Text1 Create Function
% --- Executes during object creation, after setting all properties.
function text1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hMainGui = getappdata(0, 'hMainGui');
HSGuiHandles = getappdata(hMainGui, 'HSGuiHandles');
set(hObject, 'string', num2str(HSGuiHandles.OmegaU));

%% Text2 Create Function
% --- Executes during object creation, after setting all properties.
function text2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hMainGui = getappdata(0, 'hMainGui');
HSGuiHandles = getappdata(hMainGui, 'HSGuiHandles');
set(hObject, 'string', num2str(0.5*(HSGuiHandles.OmegaL + HSGuiHandles.OmegaU)));

%% Text3 Create Function
% --- Executes during object creation, after setting all properties.
function text3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hMainGui = getappdata(0, 'hMainGui');
HSGuiHandles = getappdata(hMainGui, 'HSGuiHandles');
set(hObject, 'string', num2str(HSGuiHandles.OmegaL));

%% Text4 Create Function
% --- Executes during object creation, after setting all properties.
function text4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
hMainGui = getappdata(0, 'hMainGui');
HSGuiHandles = getappdata(hMainGui, 'HSGuiHandles');
set(hObject, 'string', num2str(HSGuiHandles.w0));


%% Opt. Anotch Filter Design with Stability Constraints imposed
% This function designs a second order antinotch digital filter by using an
% optimization technique, starting from a second order allpass-based
% antinotch filter. 
% Here we impose two types of constraints. One, the coefficient d0 must be
% greater than zero, and two, the filter must be stable. 
function [x0, xkp1, iter, ndel] = opt_antinotch_with_stability(w0,epsi,K,NN)

R = 0.95; d00 = R^2;
d10 = -cos(pi*w0)*(1+d00);
x0 = [d00; d10];
% Dividing the frequency interval [0,pi] 
% Here we divide the interval [0,pi] into several small sections for
% performing numerical integration
% 'eta' is a small constant, and the idea is to integrate in the following
% interval: [0, w0-eta] Union [w0+eta, pi]
eta = 1e-6;
divi = 1/NN;
dw = divi*pi; % dw is the infinitesimal size
omega = [0:divi:(w0-eta) (w0+eta):divi:1]*pi;
% Performing the optimization using QuadProg
xk = x0; iter = 0;
options = optimset('display','off');
while iter < K
    [H,f,A,b] = omgrad(xk, omega, epsi, dw);
    [del, fval, extflg, outpt] = quadprog(H,f,A,b,[],[],[],[],[],options);
    if extflg ~= 1
        outpt, extflg
    end
    xkp1 = xk + del;
    xk = xkp1;
    ndel = norm(del);
    iter = iter + 1;
end

%% SubFunction used by "opt_antinotch_with_stability" to compute gradient, H, f, A, and b 
%  for using in QuadProg
function [H,f,A,b] = omgrad(xk, omega, epsi, dw)
z = exp(1i*omega);
d0k = xk(1);  d1k = xk(2); % Separating the 2 components of xk
gk1 = (z.^-2 - 1).*(1 + d1k*z.^-1 + z.^-2).*0.5./(1 + d1k*z.^-1 + d0k*z.^-2).^2;
gk2 = (d0k - 1).*(1 - z.^-2).*(z.^-1).*0.5./(1 + d1k*z.^-1 + d0k*z.^-2).^2;
gk = [gk1; gk2];
H = 2*real((gk*sqrt(dw))*(gk'*sqrt(dw))); 
% The multiplication by '2' is to adjust for the '0.5' while passing to QuadProg
%
Hak = 0.5.*(1 - d0k).*(1 - z.^-2)./(1 + d1k*z.^-1 + d0k*z.^-2);
% 'Hak' is the original transfer function
Hak = Hak(:).';
f = 2*real((gk*sqrt(dw))*(Hak'*sqrt(dw))); % This '2' is part of the formula
% Constraints
Abar = [1 1; 1 -1; -1 0; 1 0];
tau = [1-epsi; 1-epsi; 1-epsi; -epsi];
A = -Abar;
b = tau + Abar*xk;


function [alp, xn, varcol_1_1] = anotchfreqvariation_for_FreqVarGUI(w0, ll, ul, NN)
%% Anotch Freq Variation for HSpot_FSlider_GUI
% This function helps analyze the variation of the antinotch filter coefficients
% with slight variations in the design frequency. 
% It has been specifically customized to work with this "HSpot_FSlider_GUI". 
% function [E,g,alp,varcol,xn] = anotchfreqvariation_for_FreqVarGUI(w0, NN)
if (w0<0.06)||(w0>0.94)
    disp('wo too close to an edge; aborting now; you make changes to the program if necessary');
    return
end
om_range = ll:(ul-ll)/(NN-1):ul;
% NOTE: Have NN above as an even number. If it is odd, then we have a
% problem since one of the entries in om_range becomes exactly equal to the
% middle frequency fn, which will give divide-by-zero error.

s = warning('off', 'optim:quadprog:SwitchToMedScale');

K = 50; N1 = 2000; epsi = 0.05;

% Getting the filter coefficients exactly at the frequency fn (w0)
[x0, xkp2, iter2, ndel] = opt_antinotch_with_stability(w0,epsi,K,N1);
xn = xkp2;

coun = 0;
hwai = waitbar(coun,['Filter number ' num2str(coun) ' completed.'],...
    'name','Please wait! Performing Optimization...');
Kk = length(om_range);
for w1 = om_range
    [x0, xkp2, iter2, ndel] = opt_antinotch_with_stability(w1,epsi,K,N1);
    coun = coun+1;
    varcol(:,coun) = xkp2;
    waitbar(coun/Kk, hwai,...
        ['Filter number ' num2str(coun) ' completed.']);      
end
close(hwai);

hh = figure; axes; hhax = gca;

hMainGui = getappdata(0, 'hMainGui');
setappdata(hMainGui, 'FrqResFigH', hh);

B2 = [1 0 -1];
A2 = [1 xn(2) xn(1)];
h2 = (1 - xn(1))/2;
[yst,xwt] = freqz(h2*B2,A2,2048); 
ys2 = abs(yst); 
xw = xwt/pi;
h12 = plot(hhax, xw, ys2);
set(h12,'YDataSource','ys2');
set(h12,'XDataSource','xw');
line([w0 w0], [0 1.1],'color','r')
title('Frequency response of optimized antinotch filter WITH stability constraints');
axis([0 1 0 1.1]), grid on,

% Computation of alp1 and alp2
dn = xn(2);
E = [om_range(:)-w0 (om_range(:)-w0).^2];
g1 = varcol(2,:)-dn;
g=g1(:);
alp = pinv(E)*g; % Although both pinv as well as '\' give least-squares soln., 
                 % pinv has the least error norm; see HELP
varcol_1_1 = varcol(1,1);
warning(s);

%% Function to Update Figs
% --- Executes when called (usually called by HSpot_FSlider_GUI), and updates the
% current hot-spot figures and the freq. response plot of the second-order BPN
% filter based on the current position of the HSpot_FSlider_GUI slider
function updateFigs(currfrq)

hMainGui = getappdata(0, 'hMainGui');
Ad       = getappdata(hMainGui);

% Updating the frequency response plot
w0 = Ad.HSGuiHandles.w0;
E1 = [currfrq-w0 (currfrq-w0)^2];
d1new = Ad.xn(2) + Ad.alp(1)*E1(1) + Ad.alp(2)*E1(2);
A2 = [1 d1new Ad.varcol_1_1]; 
B2 = [1 0 -1];
h2 = (1 - Ad.varcol_1_1)/2;
[yst,xwt] = freqz(h2*B2,A2,2048); 
ys2 = abs(yst); 
xw = xwt/pi;
refreshdata(Ad.FrqResFigH, 'caller');

% Updating the hot-spot figures
for q = 1:length(Ad.FigHndl)
    y_in = Ad.eiipseqsset{q};
    for zeu = 1:Ad.Npass_Opt
        y = filtfilt(h2*B2, A2, y_in);
        y_in = y;
    end
    ener = y.^2;
    ener = ener/max(ener);
    BPNupdtd_enrs{q} = ener;
    avgenrval = sum(ener)/length(ener);
    axshdl = findobj(Ad.FigHndl(q),'type','axes');
    namtitl = get(get(axshdl,'title'), 'string');
    wrkstt = namtitl{2,:};
    knn = strfind(wrkstt,'Average');
    wrkstt(knn:end) = '';
    wrkstt = [wrkstt 'Average Magnitude: ' num2str(avgenrval) '.'];
    namtitl{2,:} = wrkstt;
    set(get(axshdl,'title'),'string',namtitl);
    refreshdata(Ad.FigHndl(q), 'caller');
end
setappdata(hMainGui, 'BPNupdtd_enrs', BPNupdtd_enrs);


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hMainGui = getappdata(0, 'hMainGui');
Ad       = getappdata(hMainGui);

for q = 1:length(Ad.FigHndl)
    ener = Ad.BPNupdtd_enrs{q};
    protAccsn = Ad.BPNupdtd_AccNos{q};
    % Display the peak (hot-spot) locations in the command window in descending
    % order of peak magnitudes (strengths)
    [~, locs] = findpeaks(ener, 'sortstr', 'descend');
    disp(' ')
    disp('************ Tuning Intermediate Results ***********')
    disp('************ Optimized BPN Filter ******************')
    disp(['Latest Hot-spot locations for protein ID ' protAccsn])
    disp('in descending order of peak score (magnitudes):')
    disp(' ')
    disp(num2str(locs))
    disp(' ')
    disp('****************************************************')
end
