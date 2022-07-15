function varargout = HSpot_Detect_GUI(varargin)

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

% HSPOT_DETECT_GUI M-file for HSpot_Detect_GUI.fig
%      HSPOT_DETECT_GUI, by itself, creates a new HSPOT_DETECT_GUI or raises the existing
%      singleton*.
%
%      H = HSPOT_DETECT_GUI returns the handle to a new HSPOT_DETECT_GUI or the handle to
%      the existing singleton*.
%
%      HSPOT_DETECT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HSPOT_DETECT_GUI.M with the given input arguments.
%
%      HSPOT_DETECT_GUI('Property','Value',...) creates a new HSPOT_DETECT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Hot_Spot_GUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to HSpot_Detect_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help HSpot_Detect_GUI

% Last Modified by GUIDE v2.5 15-Sep-2012 09:04:37
% Also modified and customized by Parameswaran Ramachandran


%% --- Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @HSpot_Detect_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @HSpot_Detect_GUI_OutputFcn, ...
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


%% --- OpeningFcn -- Executes just before HSpot_Detect_GUI is made visible.
function HSpot_Detect_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to HSpot_Detect_GUI (see VARARGIN)

% Choose default command line output for HSpot_Detect_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Obtaining the name and path of this m-file so that the path can be added
% to the Matlab path
pcfn = mfilename('fullpath');
[nameofthisfile, pathofthisfile] = strtok(fliplr(pcfn), filesep);
nameofthisfile = fliplr(nameofthisfile);   pathofthisfile = fliplr(pathofthisfile);
if pathofthisfile(end) == filesep
    pathofthisfile = pathofthisfile(1:(length(pathofthisfile)-1));
end

% Checking if the path of this file is already in the Matlab search path.
% If it is already in path, then we leave everything as is, i.e., no
% changes to the current search path will be made. In case if it is already
% NOT in path, then we add it to the beginning of the current search path.
% However, when the GUI is closed, the added path will be automatically
% removed. Thus, effectively, there will be no permanent changes made to
% the Matlab search path. Also see the "CloseRequest_Fcn"
% callback at the very end of this m-file.
if ~isempty(strfind(lower(path), lower(pathofthisfile)))
    already_in_path = true;
else
    already_in_path = false;
    addpath(pathofthisfile, '-begin')
end

handles.already_in_path = already_in_path;
handles.pathofthisfile = pathofthisfile;

% Assigning default values for the optimization-based filter parameters
handles.R_default = 0.95;
handles.Eta_default = 1e-6;
handles.NN_default = 10000;
handles.K_default = 50;
handles.Tau_default = 0.03;
handles.ShowAll = true;

% Assigning the default values for the "no. of passes" fields
handles.Npass_C_default = 1;
handles.Npass_O_default = 1;

% Updating the handles structure
guidata(hObject, handles);

ini_dir = pwd;

% Populate listbox1
load_listbox1(ini_dir,handles)

set(handles.text59,'String','Ready!');

% Setting some application data to use for 'HSpot FSlider GUI'
setappdata(0 , 'hMainGui' , gcf);
% UIWAIT makes HSpot_Detect_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%% --- Outputs from this function are returned to the command line.
function varargout = HSpot_Detect_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
 varargout{1} = handles.output;


%% --- Pushbutton1 Callback - "Compute Consensus Spectrum"
% Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'seleiipsqs') && ~isempty(handles.seleiipsqs)
    set(handles.text59,'String','Busy...'); pause(0.001),
    [DFT_freqs, prod_DFT_magsq] = protcharfreq(handles.seleiipsqs);
    [B,IX] = sort(prod_DFT_magsq,'descend');
    handles.charfrq = DFT_freqs(IX(1));
    guidata(hObject, handles);
    set(handles.text14,'String',num2str(handles.charfrq));
    figure, plot(DFT_freqs, prod_DFT_magsq),
    title(['Product DFT (Squared Magnitude) | Determined char. freq. = '...
        num2str(handles.charfrq) '.']), xlabel('Normalized frequency'), ylabel('Magnitude')
    
    % Assigning values for the filter specs of the inverse-Chebyshev filter
    %
    handles.fstop1 = handles.charfrq - 8*1e-3;
    set(handles.edit1,'string',num2str(handles.fstop1));
    %
    handles.fpass1 = handles.charfrq - 3*1e-3;
    set(handles.edit2,'string',num2str(handles.fpass1));
    %
    handles.fpass2 = handles.charfrq + 3*1e-3;
    set(handles.edit3,'string',num2str(handles.fpass2));
    %
    handles.fstop2 = handles.charfrq + 8*1e-3;
    set(handles.edit4,'string',num2str(handles.fstop2));
    %
    handles.astop1 = 30;
    set(handles.edit5,'string',num2str(handles.astop1));
    %
    handles.apass = 1;
    set(handles.edit6,'string',num2str(handles.apass));
    %
    handles.astop2 = 30;
    set(handles.edit7,'string',num2str(handles.astop2));
    
    %Assigning values for the filter parameters of the optimization-based
    %filter
    handles.R    =  handles.R_default;
    set(handles.edit9,'string',num2str(handles.R));
    %
    handles.w0   =  handles.charfrq;
    set(handles.edit10,'string',num2str(handles.w0));
    %
    handles.Eta  =  handles.Eta_default;
    set(handles.edit11,'string',num2str(handles.Eta));
    %
    handles.K    =  handles.K_default;
    set(handles.edit12,'string',num2str(handles.K));
    %
    handles.Tau =  handles.Tau_default;
    set(handles.edit13,'string',num2str(handles.Tau));    
    %
    handles.NN   =  handles.NN_default;
    set(handles.edit14,'string',num2str(handles.NN));
    %
    handles.OmegaL =  handles.w0 - 0.05;
    set(handles.edit19,'string',num2str(handles.OmegaL));
    %
    handles.OmegaU   =  handles.w0 + 0.05;
    set(handles.edit18,'string',num2str(handles.OmegaU));
    
    % Assigning the default values to the "no. of passes" fields
    handles.Npass_Cheby = handles.Npass_C_default;
    set(handles.edit16,'string',num2str(handles.Npass_Cheby));
    %
    handles.Npass_Opt   = handles.Npass_O_default;
    set(handles.edit17,'string',num2str(handles.Npass_Opt));
    %
    % Updating the handles structure
    guidata(hObject, handles);
    %
    set(handles.text59,'String','Ready!');
else
    beep
    errordlg(['Select a valid FASTA file containing the protein sequences and "Load" it. '...
        'Read the protein sequences and convert them to EIIP sequences. '...
        'Then, you can see the consensus spectrum.'],'Bad Input','modal')
    return
end

function pushbutton2_Callback(hObject, eventdata, handles)
%% --- Pushbutton2 Callback - "Predict (Inverse-Chebyshev)"
% Executes on button press in pushbutton2.
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'fHd') && isfield(handles,'Npass_Cheby') &&...
        isnumeric(handles.Npass_Cheby) && ~isnan(handles.Npass_Cheby)
    set(handles.text59,'String','Predicting...'); pause(0.001),
    for q = 1:length(handles.seleiipsqs)
        eiipse = handles.seleiipsqs{q};
        plot_ttl = ['Hot spots of protein with accession number ' handles.seleiipacc_nos{q} '.'];
        Hd = handles.fHd;
        fcof = coeffs(Hd);
        [b,a] = sos2tf(fcof.SOSMatrix, prod(fcof.ScaleValues));
        y_in = eiipse;
        for zeu = 1:handles.Npass_Cheby
            y = filtfilt(b, a, y_in);
            y_in = y;
        end
        ener = y.^2;
        ener = ener/max(ener);
        figure, plot(ener),
        title({plot_ttl;...
            ['Inv.-Cheby. filter; Seq. length: ' num2str(length(eiipse)) '; ' ...
            'Avg. magnitude: ' num2str(sum(ener)/length(ener)) '.']}),
        xlabel('Amino acids'), ylabel('Normalized signal power at char. frequency')
        % Display the peak (hot-spot) locations in the command window in descending
        % order of peak magnitudes (strengths)
        [~, locs] = findpeaks(ener, 'sortstr', 'descend');
        disp(' ')         
        disp('*********** Inverse Chebyshev Filter ***************')
        disp(['Hot-spot locations for protein ID ' handles.seleiipacc_nos{q}])
        disp('in descending order of peak score (magnitudes):')
        disp(' ')
        disp(num2str(locs))
        disp(' ')        
        disp('****************************************************')    
    end
    guidata(hObject, handles);
    set(handles.text59,'String','Ready!');
else
    beep
    errordlg(['First, design the filter by clicking on "Design"!'],'Bad Input','modal')
    set(handles.text59,'String','Ready!');
    return    
end


%%  --- Pushbutton3 Callback - "Freq. response (Inverse-Chebyshev)"
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text59,'String','Busy...'); pause(0.001),
if isfield(handles,'fHd')
    h = fvtool(handles.fHd);
    set(h,'MagnitudeDisplay','Magnitude','DesignMask','off','name','Inverse-Chebyshev Filter');
end
set(handles.text59,'String','Ready!');


%% --- Pushbutton4 Callback - "Filter Info (Inverse-Chebyshev)"
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text59,'String','Busy...'); pause(0.001),
if isfield(handles,'fHd')
    sizesos = size(handles.fHd.sosMatrix);
    f_ord = 2*sizesos(1);
    if datenum(version('-date')) > 732892 %732892 represents the release date of 
        % Matlab version (R2006b), i.e., August 03, 2006
        default_info = info(handles.fHd, 'long');
    else
        default_info = info(handles.fHd);
    end
    filinfo = strvcat(default_info(1:2,:), ['Filter Order : ' num2str(f_ord)], default_info(3:end,:));
%     hm = msgbox(filinfo,'Current Filter Information');
    disp(' ')
    disp('****** CURRENT INVERSE-CHEBYSHEV FILTER INFORMATION ******')
    disp(' ')
    disp(default_info)
    disp(' ')
    disp('****** END OF FILTER INFORMATION ******')
end
set(handles.text59,'String','Ready!');


%%  --- Pushbutton5 Callback - "About the Parameters (Optimized)"
% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pinfo1 = ['"R" is the pole radius of the initial filter that is used as'...
    ' the starting point for the optimization. The default value assumed is '...
    num2str(handles.R_default) '.'];

pinfo2 = ['"w0" is the center frequency or, in other words, the notch frequency of interest'...
    ' forming the narrow passband of the filter. This corresponds to the period-3 frequency.'];

pinfo3 = ['"Eta" is a small constant. It is used to obtain the interval over which the'...
    ' discrete integration is performed. This interval is defined as the union of'...
    ' [0, w0-Eta] and [w0+Eta, pi]. The default value assumed is ' num2str(handles.Eta_default) '.'];

pinfo4 = ['"K" is the user-specified number of iterations of optimization. The default'...
    ' value assumed is ' num2str(handles.K_default) '.'];

pinfo5 = ['"Tau" is the stability margin. The stability of the filter designed'...
    ' is guaranteed if the poles are inside the unit circle. In order to ensure robust stability'...
    ' "Tau" imposes an additional constraint such that the poles will be forced to lie inside'...
    ' the circle with radius 1-Tau. The default value assumed for "Tau" is '...
    num2str(handles.Tau_default) '.'];

pinfo6 = ['"NN" is the number of divisions of the above interval of discrete integration.'...
    ' The larger this number, the greater would be the number of divisions, resulting in'...
    ' a denser sampling of the interval. The default value assumed is '...
    num2str(handles.NN_default) '.'];

pinfo7 = ['"Omega_L" is the lower bound for the range of "w0".'...
    ' It can be varied by the user, but then the optimization has to be performed again with the new value.'...
    ' The default value assumed is w0-0.05.'];

pinfo8 = ['"Omega_U" is the upper bound for the range of "w0".'...
    ' It can be varied by the user, but then the optimization has to be performed again with the new value.'...
    ' The default value assumed is w0+0.05.'];

pinfo9 = ['NOTE: Although "Omega_L" and "Omega_U" are allowed to be varied by the user, the least-squares'...
    ' polynomial model is most accurate if "Omega_L" and "Omega_U" are close to "w0", and are symmetrically placed'...
    ' on either side of "w0".'];

pinfo = strvcat(pinfo1,' ',' ',pinfo2,' ',' ',pinfo3,' ',' ',pinfo4,' ',' ',pinfo5,' ',' ',pinfo6,...
    ' ',' ',pinfo7,' ',' ',pinfo8,' ',' ',pinfo9);
msgbox(pinfo,'About the filter parameters','help','replace');


%%  --- Pushbutton6 Callback - "Filter Info (Optimized)"
% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'x0') && isfield(handles,'xkp1') &&  isfield(handles,'iter') &&  isfield(handles,'ndel')
    set(handles.text59,'String','Busy...'); pause(0.001),   
%     Assign the filter info from the handles structure
    xkp1 = handles.xkp1; 
    
    B2 = [1 0 -1];
    A2 = [1 xkp1(2) xkp1(1)];
    h2 = (1 - xkp1(1))/2;
    
    filinfo = strvcat(...
        'Optimized BPN Filter',...
        '======================',...
        ' ',...
        ['Gain constant : ' num2str(h2)],...
        ['Numerator coefficients : ' num2str(B2)],...
        ['Denominator coefficients : ' num2str(A2)]...
        );
    
    msgbox(filinfo,'Current Filter Information');
    set(handles.text59,'String','Ready!');
end


%%  --- Pushbutton7 Callback - "Freq. responses (Optimized)"
% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'x0') && isfield(handles,'xkp1') &&  isfield(handles,'iter') &&  isfield(handles,'ndel')
    set(handles.text59,'String','Busy...'); pause(0.001),
%     Assign the filter info from the handles structure
    xkp1 = handles.xkp1; iter = handles.iter; ndel = handles.ndel;
       
    B2 = [1 0 -1];
    A2 = [1 xkp1(2) xkp1(1)];
    h2 = (1 - xkp1(1))/2;
    
    hn2 = fvtool(h2*B2,A2);
    set(hn2,'MagnitudeDisplay','Magnitude','DesignMask','off','name',...
        'Optimization Technique - Optimized BPN Filter');
    set(handles.text59,'String','Ready!');
end


function pushbutton8_Callback(hObject, eventdata, handles)
%%  --- Pushbutton8 Callback - "Predict (Optimized)"
% --- Executes on button press in pushbutton8.
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.PrdctFlg = true;
guidata(hObject, handles);
hMainGui = getappdata(0, 'hMainGui');
if isfield(handles,'x0') && isfield(handles,'xkp1') &&...
        isfield(handles,'iter') && isfield(handles,'ndel') &&...
        isfield(handles,'Npass_Opt') && isnumeric(handles.Npass_Opt) && ~isnan(handles.Npass_Opt)
    set(handles.text59,'String','Predicting...'); pause(0.001),
    for q = 1:length(handles.seleiipsqs)
        eiipse = handles.seleiipsqs{q};
        plot_ttl = ['Hot spots of protein with accession number ' handles.seleiipacc_nos{q}];
        x0 = handles.x0; xkp1 = handles.xkp1; iter = handles.iter; ndel = handles.ndel;
        %
        B2 = [1 0 -1];
        A2 = [1 xkp1(2) xkp1(1)];
        h2 = (1 - xkp1(1))/2;
        %
%         Apply the filter to the signal
        y_in = eiipse;
        eiipseqsset{q} = eiipse;
        for zeu = 1:handles.Npass_Opt
            y = filtfilt(h2*B2, A2, y_in);
            y_in = y;
        end
        ener = y.^2;
        ener = ener/max(ener);
        FigHndl(q) = figure; 
        h12 = plot(ener);
        set(h12,'YDataSource','ener');
        avgenrval = sum(ener)/length(ener);
        title({plot_ttl;...
            ['Optimized BPN filter; Seq. Length: ' num2str(length(eiipse)) '; ' ...
            'Average Magnitude: ' num2str(avgenrval) '.']}),
        xlabel('Amino acids'), ylabel('Normalized signal power at characteristic frequency')
        % Display the peak (hot-spot) locations in the command window in descending
        % order of peak magnitudes (strengths)
        [~, locs] = findpeaks(ener, 'sortstr', 'descend');
        disp(' ')
        disp('************ Optimized BPN Filter ******************')
        disp(['Hot-spot locations for protein ID ' handles.seleiipacc_nos{q}])
        disp('in descending order of peak score (magnitudes):')
        disp(' ')        
        disp(num2str(locs))
        disp(' ')        
        disp('****************************************************')            
        %
        BPNupdtd_AccNos{q}   = handles.seleiipacc_nos{q};
        BPNupdtd_enrs{q}     = ener;
    end
    setappdata(hMainGui, 'eiipseqsset', eiipseqsset);
    setappdata(hMainGui, 'Npass_Opt', handles.Npass_Opt);
    setappdata(hMainGui, 'FigHndl', FigHndl);
    %
    setappdata(hMainGui, 'BPNupdtd_AccNos', BPNupdtd_AccNos);
    setappdata(hMainGui, 'BPNupdtd_enrs', BPNupdtd_enrs);
    %
    set(handles.text59,'String','Ready!');
else
    beep
    errordlg(['First, design the filter by clicking on "Design"!'],'Bad Input','modal')
    return
end


%%  --- Pushbutton9 Callback - "Show only '.txt' files"
% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ShowAll = false;
guidata(hObject, handles);
load_listbox1(pwd,handles)


%%  --- Pushbutton10 Callback - "Show all files"
% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ShowAll = true;
guidata(hObject, handles);
load_listbox1(pwd,handles)


function pushbutton11_Callback(hObject, eventdata, handles)
%%  --- Pushbutton11 Callback - "Design Filter - Inverse-Chebyshev"
% --- Executes on button press in pushbutton11.
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

everything_ok = inv_cheby_is_everything_fine(handles);

if everything_ok
    set(handles.text59,'String','Designing...'); pause(0.1);
    Hd = HSpot_DF_design_file(handles.charfrq,handles.fstop1,handles.fpass1,handles.fpass2,...
        handles.fstop2,handles.astop1,handles.apass,handles.astop2,handles.match);
    handles.fHd = Hd;
    guidata(hObject, handles);    
    set(handles.text59,'String','Ready!');
%    order(Hd)
end
        

%%  --- Pushbutton12 Callback - "More Info about No. of Passes - Inverse-Chebyshev"
% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pinfo = ['This field takes in the no. of passes through the filter. The default value is 1.'...
    ' The greater the no. of passes, the more the effect of filtering. This is equivalent to filtering'...
    ' with a higher order filter. For example, if the value in this field is 2, then the protein sequence'...
    ' is filtered twice with a second order filter, which brings in the effect of filtering with a filter'...
    ' of order 4.'];

msgbox(pinfo,'No. of Passes','help','replace');


%%  --- Pushbutton13 Callback - "Design Filter - Optimized"
% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

everything_ok = opt_is_everything_fine(handles);

if everything_ok
    set(handles.text59,'String','Designing...'); pause(0.001),
%     Switch off the "SwitchToMedScale" warning
    s = warning('off', 'optim:quadprog:SwitchToMedScale');

    [x0, xkp1, iter, ndel] = opt_anotch_with_stab_for_GUI(handles.R, handles.w0,...
        handles.Eta, handles.K, handles.Tau, handles.NN);
    
%     Reset all warnings to previous settings
    warning(s);
    
%     Assign the filter info to the handles structure
    handles.x0 = x0; 
    handles.xkp1 = xkp1; 
    handles.iter = iter; 
    handles.ndel = ndel;
    guidata(hObject, handles);
    set(handles.text59,'String','Ready!');
end

    
%%  --- Pushbutton14 Callback - "More Info about No. of Passes - Optimized"
% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pinfo = ['This field takes in the no. of passes through the filter. The default value is 1.'...
    ' The greater the no. of passes, the more the effect of filtering. This is equivalent to filtering'...
    ' with a higher order filter. For example, if the value in this field is 2, then the protein sequence'...
    ' is filtered twice with a second order filter, which brings in the effect of filtering with a filter'...
    ' of order 4.'];

msgbox(pinfo,'No. of Passes','help','replace');

%%  --- Pushbutton25 Callback - "Opens the HSpot_FSlider_GUI"
% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hMainGui = getappdata(0, 'hMainGui');

if isappdata(hMainGui, 'FigHndl')
    FHexst = true;
    FigHndl = getappdata(hMainGui,'FigHndl');
else
    FHexst = false;
end

if ~handles.PrdctFlg
    beep
    errordlg(['FSlider_GUI is useful only for manual tuning, and hence will work only '...
        'with the "Predict" button. It is not designed to work with the '...
        '"Predict & AutoTune" button. Please close all figure windows and re-invoke '...
        'FSlider_GUI after pressing "Predict".'],'Bad Input','modal')
    return
elseif ~FHexst || ~any(ishandle(FigHndl))
    beep
    errordlg(['First select the sequences of interest and predict the hot spots '...
        'using the default filter! Leave the hot-spot figures open. '...
        'Then, you can fine tune using the HSpot_FSlider_GUI.'],'Bad Input','modal')
    return
else
    setappdata(hMainGui, 'HSGuiHandles', handles);
    HSpot_FSlider_GUI;
end


function pushbutton26_Callback(hObject, eventdata, handles)
%% --- Pushbutton26 Callback - "Load the file containing the FASTA protein sequences"
% --- Executes on button press in pushbutton26.
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text59,'String','Busy...'); pause(0.001),
handles.mainFASTAfile = handles.mainCHOSENfile;
fid = fopen(handles.mainFASTAfile{1});
accs_nums = {}; cc = 0;
while feof(fid) == 0
    tline = fgetl(fid);
    if tline(1) == '>'
        cc = cc+1;
        kpu = strfind(tline,'|');
        if ~isempty(kpu)
            accs_nums{cc} = tline(2:kpu-1);            
        else
            beep
            errordlg('ERROR! Delimiting character | not found! Check FASTA file!','Bad Input','modal')
            set(handles.text59,'String','Ready!');
            return
        end
    end
end
status = fclose(fid);
if status ~= 0
    disp('ERROR! File did not close properly!!')
end
handles.accs_nums = accs_nums;
guidata(hObject, handles);
load_listbox4(hObject, eventdata, handles);
set(handles.text59,'String','Ready!');


function pushbutton27_Callback(hObject, eventdata, handles)
%% --- Pushbutton27 Callback - "Convert the selected character sequences into EIIP sequences"
% --- Executes on button press in pushbutton27.
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text59,'String','Busy...'); pause(0.001),
if isfield(handles,'charsqs')
    handles.eiipsqs = {};
    for jj = 1:length(handles.charsqs)
        handles.eiipsqs{jj} = protchareiip(handles.charsqs{jj});
    end
else
    beep
    errordlg('Error! No character sequences found in memory! Read them first!','Bad Input','modal')
    set(handles.text59,'String','Ready!');
    return
end
guidata(hObject, handles);
load_listbox5(hObject, eventdata, handles);
set(handles.text59,'String','Ready!');


function pushbutton28_Callback(hObject, eventdata, handles)
%% --- Pushbutton28 Callback - "Read the selected protein character sequences"
% --- Executes on button press in pushbutton28.
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text59,'String','Busy...'); pause(0.001),
fid = fopen(handles.mainFASTAfile{1});
handles.charsqs = {}; char_sq = []; cc = 1; beechmein = false; 
while feof(fid) == 0
    tline = fgetl(fid);
    if isempty(strfind(tline,'>')) && beechmein
        char_sq = strcat(char_sq, tline);
    end
    if ~isempty(strfind(tline,'>'))
        if beechmein
            handles.charsqs{cc} = char_sq;
            cc = cc+1;
            beechmein = false;
        end
        if length(handles.selacc_nos) >= cc
            if ~isempty(strfind(tline,handles.selacc_nos{cc}))
                char_sq = [];
                beechmein = true;
            end
        end
    end
end
if beechmein
    handles.charsqs{cc} = char_sq;
end
status = fclose(fid);
if status ~= 0
    disp('ERROR! File did not close properly!!')
end
guidata(hObject, handles);
set(handles.text59,'String','Ready!');


function pushbutton38_Callback(hObject, eventdata, handles)
%% --- Pushbutton38 Callback - "Predicts & AUTO_TUNE"
% --- Executes on button press in pushbutton38.
% hObject    handle to pushbutton38 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.PrdctFlg = false;
guidata(hObject, handles);
if isfield(handles,'x0') && isfield(handles,'xkp1') &&...
        isfield(handles,'iter') && isfield(handles,'ndel') &&...
        isfield(handles,'Npass_Opt') && isnumeric(handles.Npass_Opt) && ~isnan(handles.Npass_Opt)
    set(handles.text59,'String','Predicting... AutoTuning...'); pause(0.001),
    [alp, xn, varcol_1_1] = anotchfreqvariation_for_AutoTuning(handles.w0, handles.OmegaL, handles.OmegaU, 10);
    w0 = handles.w0;
    for q = 1:length(handles.seleiipsqs)
        eiipse = handles.seleiipsqs{q};
        plot_ttl = ['Auto-Tuned Hot spots of protein with accession number ' handles.seleiipacc_nos{q}];
        % AUTOMATIC TUNING
        % First, loop through all freqs. between OmegaL and OmegaU and compute the total power of
        % the signal for each case and store
        for currfrq = handles.OmegaL:0.001:handles.OmegaU
            E1 = [currfrq-w0 (currfrq-w0)^2];
            d1new = xn(2) + alp(1)*E1(1) + alp(2)*E1(2);
            A2 = [1 d1new varcol_1_1];
            B2 = [1 0 -1];
            h2 = (1 - varcol_1_1)/2;
            y = filtfilt(h2*B2, A2, eiipse);
            ener = y.^2;            
            % Compute the total power of the signal for tuning
            sOFMtcngHSs = sum(ener);
            decdcrit(uint64(((currfrq-handles.OmegaL)/0.001)+1),:) = [currfrq sOFMtcngHSs];
        end
        % THEN, CHOOSE THE FREQ. THAT BEST MATCHES THE KNOWN HOT SPOTS
        [Cvl, dcidx] = max(decdcrit(:,2));
        goodfrq = decdcrit(dcidx,1);
        % Finally, recompute the corresponding "ener" with the best freq. chosen
        E1 = [goodfrq-w0 (goodfrq-w0)^2];
        d1new = xn(2) + alp(1)*E1(1) + alp(2)*E1(2);
        A2 = [1 d1new varcol_1_1];
        B2 = [1 0 -1];
        h2 = (1 - varcol_1_1)/2;
        y = filtfilt(h2*B2, A2, eiipse);
        ener = y.^2;
        ener = ener/max(ener);
        % PLOTTING
        figure, plot(ener),
        title({plot_ttl;...
            ['Optimized (& Tuned) BPN filter; Seq. Length: ' num2str(length(eiipse)) ...
            '; Tuned freq. ' num2str(goodfrq)]}),
        xlabel('Amino acids'), ylabel('Normalized signal power at characteristic frequency')
         % END OF PLOTTING
        % Display the peak (hot-spot) locations in the command window in descending
        % order of peak magnitudes (strengths)
        [~, locs] = findpeaks(ener, 'sortstr', 'descend');
        disp(' ')
        disp('******* AUTO-TUNED Optimized BPN Filter ************')
        disp(['Hot-spot locations for protein ID ' handles.seleiipacc_nos{q}])
        disp('in descending order of peak score (magnitudes):')
        disp(' ')        
        disp(num2str(locs))
        disp(' ')        
        disp('****************************************************')                
    end
    set(handles.text59,'String','Ready!');
else
    beep
    errordlg(['First, design the filter by clicking on "Design"!'],'Bad Input','modal')
    return
end

        
%% --- Listbox1 Callback
% Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

index_selected = get(handles.listbox1,'Value');
file_list = get(handles.listbox1,'String');
get(handles.figure1,'SelectionType');
if strcmp(get(handles.figure1,'SelectionType'),'open')
    filename = file_list{index_selected};
    if  handles.is_dir1(index_selected)
        cd (filename)
        load_listbox1(pwd,handles)
    end
else
    filename = file_list(index_selected);
    sel_type = handles.is_dir1(index_selected);
    if all(~sel_type)
        handles.mainCHOSENfile = filename;
        guidata(hObject, handles);
    end
end


%% --- Load Listbox1
function load_listbox1(dir_path,handles)
cd (dir_path)

if handles.ShowAll
    dir_struct = dir(dir_path); % If ShowAll is enabled, the get the whole list
else
    aaa = dir(fullfile(dir_path, '*.'));
    bbb = dir(fullfile(dir_path, '*.txt'));
    dir_struct = [aaa; bbb]; % Else, get only the list of directories and '.txt' files
end

handles.file_names1 = {dir_struct.name};
handles.is_dir1 = [dir_struct.isdir];
guidata(handles.figure1,handles);
set(handles.listbox1,'String',handles.file_names1,...
	'Value',1)
set(handles.text2,'String',pwd)


%% --- Listbox1 CreateFcn
% Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Listbox4 Callback
% --- Executes on selection change in listbox4.
function listbox4_Callback(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox4
selacc_indx = get(handles.listbox4,'Value');
accs_nums = get(handles.listbox4,'String');
selacc_nos = accs_nums(selacc_indx);
handles.selacc_nos = selacc_nos;
guidata(hObject, handles);


%% --- Load Listbox4
function load_listbox4(hObject, eventdata, handles)
set(handles.listbox4,'String',handles.accs_nums,'Value',1)
listbox4_Callback(hObject, eventdata, handles);


%% --- Listbox4 CreatFcn
% --- Executes during object creation, after setting all properties.
function listbox4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Listbox5 Callback
% --- Executes on selection change in listbox5.
function listbox5_Callback(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox5
set(handles.text59,'String','Loading...');
seleacc_indx = get(handles.listbox5,'Value');
eaccs_nums = get(handles.listbox5,'String');
handles.seleiipacc_nos = eaccs_nums(seleacc_indx);
handles.seleiipsqs = handles.eiipsqs(seleacc_indx);
handles.selcharsqs = handles.charsqs(seleacc_indx);
guidata(hObject, handles);
set(handles.text59,'String','Ready!');


%% --- Load Listbox5
function load_listbox5(hObject, eventdata, handles)
set(handles.listbox5,'String',handles.selacc_nos,'Value',1)
listbox5_Callback(hObject, eventdata, handles)


%% --- Listbox5 CreatFcn
% --- Executes during object creation, after setting all properties.
function listbox5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit1 Callback - FStop1 Text Field
function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.fstop1 = user_entry;
guidata(hObject, handles);

%% --- Edit1 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit2 Callback - FPass1 Text Field
function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.fpass1 = user_entry;
guidata(hObject, handles);

%% --- Edit2 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit3 Callback - FPass2 Text Field
function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.fpass2 = user_entry;
guidata(hObject, handles);

%% --- Edit3 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit4 Callback - FStop2 Text Field
function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.fstop2 = user_entry;
guidata(hObject, handles);

%% --- Edit4 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit5 Callback - AStop1 Text Field
function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.astop1 = user_entry;
guidata(hObject, handles);

%% --- Edit5 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit6 Callback - APass Text Field
function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.apass = user_entry;
guidata(hObject, handles);

%% --- Edit6 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit7 Callback - AStop2 Text Field
function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.astop2 = user_entry;
guidata(hObject, handles);

%% --- Edit7 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit9 Callback - R (Pole Radius) Text Field
function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.R = user_entry;
guidata(hObject, handles);

%% --- Edit9 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit10 Callback - w0 (Center Freq.) Text Field
function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.w0 = user_entry;
guidata(hObject, handles);

%% --- Edit10 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit11 Callback - Eta (Small Constant) Text Field
function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.Eta = user_entry;
guidata(hObject, handles);

%% --- Edit11 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit12 Callback - K (No. of Iterations) Text Field
function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.K = user_entry;
guidata(hObject, handles);

%% --- Edit12 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit13 Callback - Tau (Stability Margin) Text Field
function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.Tau = user_entry;
guidata(hObject, handles);

%% --- Edit13 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit14 Callback - NN (No. of Divisions) Text Field
function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.NN = user_entry;
guidata(hObject, handles);

%% --- Edit14 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit16 Callback - (No. of Passes) - Inverse-Chebyshev
function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.Npass_Cheby = user_entry;
guidata(hObject, handles);


%% --- Edit16 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Edit17 Callback - (No. of Passes) - Optimized
function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double

user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.Npass_Opt = user_entry;
guidata(hObject, handles);


%% --- Edit17 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit18_Callback(hObject, eventdata, handles)
%% --- Edit18 Callback - OmegaU
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double
user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.OmegaU = user_entry;
guidata(hObject, handles);


%% --- Edit18 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%% --- Edit19 Callback - OmegaL
function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double
user_entry = str2double(get(hObject,'string'));
if isnan(user_entry)
    errordlg('You must enter a numeric value','Bad Input','modal')
end
handles.OmegaL = user_entry;
guidata(hObject, handles);

%% --- Edit19 CreateFcn
% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Popup Menu1 Callback
% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        popupmenu1

val = get(hObject,'Value');
switch val
    case 1
        sel = 'passband';
    case 2
        sel = 'stopband';
end
handles.match = sel;
guidata(hObject, handles);

    
%% --- Popup Menu1 CreateFcn
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

handles.match = 'passband';
guidata(hObject, handles);


function everything_ok = inv_cheby_is_everything_fine(handles)
%% --- Function to check if all inputs are proper for the Inverse-Chebyshev based Technique
everything_ok = false;
allfreqs_zero_one = false;  allfreqs_sorted = false;
all_set_to_go = false; allvals_correct = false;

% Checking whether all fields in the handles structure exist
allfields_exist = isfield(handles,'fstop1') && isfield(handles,'fpass1') &&...
    isfield(handles,'fpass2') && isfield(handles,'fstop2') && isfield(handles,'astop1') &&...
    isfield(handles,'apass') && isfield(handles,'astop2');

% Checking if all values have been correctly entered
if allfields_exist
    allvals_correct = ~isnan(handles.fstop1) && ~isnan(handles.fpass1) && ~isnan(handles.fpass2) &&...
        ~isnan(handles.fstop2) && ~isnan(handles.astop1) && ~isnan(handles.apass) &&...
        ~isnan(handles.astop2);
end

% If all freq. values have been entered then checking if they are all in
% the correct range (between 0 and 1) and in the right
% order(fstop1 < fpass1 < fpass2 < fstop2)
if allfields_exist
    freqvec = [handles.fstop1 handles.fpass1 handles.fpass2 handles.fstop2];
    allfreqs_zero_one = (freqvec > 0) & (freqvec < 1);
    allfreqs_zero_one = all(allfreqs_zero_one);
    allfreqs_sorted = issorted(freqvec);
end

all_set_to_go = allfields_exist && allvals_correct && allfreqs_zero_one && allfreqs_sorted;

if ~(all_set_to_go) && isfield(handles,'charfrq')
    beep
    errordlg(['Error in the input values! All values must be correctly entered and all frequencies '...
        'must be in the open interval (0, 1) (normalized) and must obey the order '...
        'FStop1 < FPass1 < FPass2 < FStop2.'],'Bad Input','modal')
    return
elseif ~(all_set_to_go) && ~isfield(handles,'charfrq')
    beep
    errordlg(['Compute consensus spectrum and also correct/enter input values! '...
        'All values must be correctly entered and all frequencies '...
        'must be in the open interval (0, 1) (normalized) and must obey the order '...
        'FStop1 < FPass1 < FPass2 < FStop2.'],'Bad Input','modal')
    return
elseif (all_set_to_go) && ~isfield(handles,'charfrq')
    beep
    errordlg('Compute consensus spectrum!','Bad Input','modal')
    return
else
    everything_ok = true;
end

function everything_ok = opt_is_everything_fine(handles)
%% --- Function to check if all inputs are proper for the Optimization based Technique
everything_ok = false;
all_set_to_go = false; allvals_correct = false;
R_ok = false;  w0_ok = false; Tau_ok = false;

% Checking whether all fields in the handles structure exist
allfields_exist = isfield(handles,'R') && isfield(handles,'w0') &&...
    isfield(handles,'Eta') && isfield(handles,'K') && isfield(handles,'Tau') &&...
    isfield(handles,'NN');

% Checking if all values have been correctly entered
if allfields_exist
    allvals_correct = ~isnan(handles.R) && ~isnan(handles.w0) && ~isnan(handles.Eta) &&...
        ~isnan(handles.K) && ~isnan(handles.Tau) && ~isnan(handles.NN);
end

% If all parameter values have been entered then checking if they are all in
% the correct range (between 0 and 1) and in the right
% order(fstop1 < fpass1 < fpass2 < fstop2)
if allfields_exist && allvals_correct
    R_ok = (handles.R > 0) & (handles.R < 1);
    w0_ok = (handles.w0 > 0) & (handles.w0 < 1);
    Tau_ok = (handles.Tau > 0) & (handles.Tau < 1);
end

all_set_to_go = allfields_exist && R_ok && w0_ok && Tau_ok;

if ~(all_set_to_go) && isfield(handles,'charfrq')
    beep
    errordlg(['Error in the input values! All values must be correctly entered and R, w0, and Tau '...
        'must be in the open interval (0, 1).'],'Bad Input','modal')
    return
elseif ~(all_set_to_go) && ~isfield(handles,'charfrq')
    beep
    errordlg(['Compute consensus spectrum and also correct/enter input values! '...
        'All values must be correctly entered and R,w0, and Tau '...
        'must be in the open interval (0, 1).'],'Bad Input','modal')
    return
elseif (all_set_to_go) && ~isfield(handles,'charfrq')
    beep
    errordlg('Compute consensus spectrum!','Bad Input','modal')
    return
else
    everything_ok = true;
end


function [DFT_freqs, prod_DFT_magsq] = protcharfreq(seleiipsqs) 
%% --- ProtCharfreq Function - Computes the Consensus Spectrum
% This M-File computes the consensus
% spectrum, calculated as per the concepts in Cosic'94 paper, i.e., by
% using the RRM Model. 
% The inputs are the protein EIIP sequences.
%
% Written by: Paramesh Ramachandran.
% Last Modified: March, 2009. 
% -------------------------------------

fft_number = 2048;   

% Removing the zero-frequency component and then calculating the DFTs of
% each protein numerical sequence participating in the product
for p = 1:length(seleiipsqs)
   prot_avg = sum(seleiipsqs{p})/length(seleiipsqs{p});
   zero_frq_remov = seleiipsqs{p} - prot_avg;
   fou_array{p} = fft(zero_frq_remov, fft_number);
end
  
% Calculating the product of the DFTs to determine the characteristic freq.
% i.e., calculating the consensus spectrum
prod_DFT = ones(1, fft_number);
for p = 1:length(seleiipsqs)
   prod_DFT = prod_DFT.*abs(fou_array{p});
end
prod_DFT_magsq = (prod_DFT(1:length(prod_DFT)/2)).^2;
DFT_freqs = linspace(0,1,fft_number/2);


function [num_seq] = protchareiip(char_seq)
%% --- ProtCharEIIP function - Converts character sequences into EIIP sequences
% This function receives a character protein sequence as input, assigns the corresponding numerical EIIP
% values for the amino acids, and returns the numerical
% sequence as output.
len = length(char_seq);
for i = 1:len
    switch upper(char_seq(i))
        case 'L'
            num_seq(i) = 0;
        case 'I'
            num_seq(i) = 0;
        case 'N'
            num_seq(i) = 0.0036;
        case 'G'
            num_seq(i) = 0.0050;
        case 'V'
            num_seq(i) = 0.0057;
        case 'E'
            num_seq(i) = 0.0058;
        case 'P'
            num_seq(i) = 0.0198;
        case 'H'
            num_seq(i) = 0.0242;
        case 'K'
            num_seq(i) = 0.0371;
        case 'A'
            num_seq(i) = 0.0373;
        case 'Y'
            num_seq(i) = 0.0516;
        case 'W'
            num_seq(i) = 0.0548;
        case 'Q'
            num_seq(i) = 0.0761;
        case 'M'
            num_seq(i) = 0.0823;
        case 'S'
            num_seq(i) = 0.0829;
        case 'C'
            num_seq(i) = 0.0829;
        case 'T'
            num_seq(i) = 0.0941;
        case 'F'
            num_seq(i) = 0.0946;
        case 'R'
            num_seq(i) = 0.0959;
        case 'D'
            num_seq(i) = 0.1263;
        otherwise
            error('ERROR! INVALID AMINO ACID FOUND! CHECK CHAR. STRING!')
            beep on, beep
    end
end
num_seq = num_seq(:)';


%% --- HSpot_DF_design_file 
% 'charf' ---> inputs the characteristic frequency (the center frequency of
% the filter)
%HSPOT_DF_DESIGN_FILE Returns a discrete-time filter object.
function Hd = HSpot_DF_design_file(charf,Fstop1,Fpass1,Fpass2,Fstop2,Astop1,Apass,Astop2,match)
%
% M-File generated by MATLAB(R) 7.3 and the Signal Processing Toolbox 6.6.
%
% Chebyshev Type II Bandpass filter designed using FDESIGN.BANDPASS.
% All frequency values are normalized to 1.
% IF HARD ASSIGNMENT IS DESIRED, UNCOMMENT AND USE THE LINES BELOW
% Fstop1 = charf - 8*1e-3;       % First Stopband Frequency
% Fpass1 = charf - 3*1e-3;       % First Passband Frequency
% Fpass2 = charf + 3*1e-3;       % Second Passband Frequency
% Fstop2 = charf + 8*1e-3;       % Second Stopband Frequency
% Astop1 = 30;                   % First Stopband Attenuation (dB)
% Apass  = 1;                    % Passband Ripple (dB)
% Astop2 = 30;                   % Second Stopband Attenuation (dB)
% match  = 'passband';           % Band to match exactly

% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, ...
                      Astop2);
Hd = design(h, 'cheby2', 'MatchExactly', match);

function [x0, xkp1, iter, ndel] = opt_anotch_with_stab_for_GUI(R,w0,eta,K,Tau,NN)
%% --- Opt. Filter Design File
% This m-file designs a second order antinotch digital filter by using an
% optimization technique, starting from a second order allpass-based antinotch
% filter. For the optimization, we impose two types of constraints. One, the
% coefficient d0 must be greater than zero, and two, the filter must be stable.
%
% Written by: Parameswaran Ramachandran; Last Modified in 2009.
%
x0 = [R^2; -cos(pi*w0)*(1+R^2)];
% Dividing the frequency interval [0,pi] 
% Here we divide the interval [0,pi] into several small sections for
% performing numerical integration
% 'eta' is a small constant, and the idea is to integrate in the following
% interval: [0, w0-eta] Union [w0+eta, pi]
divi = 1/NN;
dw = divi*pi; % dw is the infinitesimal size
omega = [0:divi:(w0-eta) (w0+eta):divi:1]*pi;
% Performing the optimization using QuadProg
xk = x0; iter = 0;
options = optimset('display','off');
hwai = waitbar(iter,'Please wait... Performing Optimization...');
while iter < K
    [H,f,A,b] = omgrad(xk, omega, Tau, dw);
    [del, fval, extflg, outpt] = quadprog(H,f,A,b,[],[],[],[],[],options);
    if extflg ~= 1
        outpt, extflg
    end
    xkp1 = xk + del;
    xk = xkp1;
    ndel = norm(del);
    iter = iter + 1;
    waitbar(iter/K)
end
close(hwai);


function [H,f,A,b] = omgrad(xk, omega, Tau, dw)
%% --- Function to compute gradient, H, f, A, and b for using in QuadProg
z = exp(1i*omega);
d0k = xk(1);  d1k = xk(2); % Separating the 2 components of xk
gk1 = (z.^-2 - 1).*(1 + d1k*z.^-1 + z.^-2).*0.5./(1 + d1k*z.^-1 + d0k*z.^-2).^2;
gk2 = (d0k - 1).*(1 - z.^-2).*(z.^-1).*0.5./(1 + d1k*z.^-1 + d0k*z.^-2).^2;
gk = [gk1; gk2];
H = 2*real((gk*sqrt(dw))*(gk'*sqrt(dw))); 
% The multiplication by '2' is to adjust for the '0.5' while passing to QuadProg
Hak = 0.5.*(1 - d0k).*(1 - z.^-2)./(1 + d1k*z.^-1 + d0k*z.^-2);
% 'Hak' is the original transfer function
Hak = Hak(:).';
f = 2*real((gk*sqrt(dw))*(Hak'*sqrt(dw))); % This '2' is part of the formula
% Constraints
Abar = [1 1; 1 -1; -1 0; 1 0];
tauu = [1-Tau; 1-Tau; 1-Tau; -Tau];
A = -Abar;
b = tauu + Abar*xk;


function [alp, xn, varcol_1_1] = anotchfreqvariation_for_AutoTuning(w0, ll, ul, NN)
%% --- Anotch Freq Variation for Auto Tuning
% This function helps analyze the variation of the antinotch filter coefficients
% with slight variations in the design frequency. 
% It has been specifically customized to perform Automatic Tuning. The
% original version can be found in "HSpot_FSlider_GUI.m".
if (w0<0.06)||(w0>0.94)
    disp('wo too close to an edge; aborting now; Make changes to the program if necessary!');
    return
end
om_range = ll:(ul-ll)/(NN-1):ul;
% NOTE: Have NN above as an even number. If it is odd, then we have a
% problem since one of the entries in om_range becomes exactly equal to the
% middle frequency fn, which will give a divide-by-zero error.

s = warning('off', 'optim:quadprog:SwitchToMedScale');
K = 50; N1 = 2000; epsi = 0.05;
% Getting the filter coefficients exactly at the frequency fn (w0)
[x0, xkp2, iter2, ndel] = opt_antinotch_with_stab_for_AutoTuning(w0,epsi,K,N1);
xn = xkp2;
coun = 0;
hwai = waitbar(coun,['Filter number ' num2str(coun) ' completed.'],...
    'name','Please wait! Performing Optimization...');
Kk = length(om_range);
for w1 = om_range
    [x0, xkp2, iter2, ndel] = opt_antinotch_with_stab_for_AutoTuning(w1,epsi,K,N1);
    coun = coun+1;
    varcol(:,coun) = xkp2;
    waitbar(coun/Kk, hwai,...
        ['Filter number ' num2str(coun) ' completed.']);      
end
close(hwai);
% Computation of alp1 and alp2
dn = xn(2);
E = [om_range(:)-w0 (om_range(:)-w0).^2];
g1 = varcol(2,:)-dn;
g=g1(:);
alp = pinv(E)*g; % Although both pinv as well as '\' give least-squares soln., 
                 % pinv has the least error norm; see HELP
varcol_1_1 = varcol(1,1);
warning(s);


%% --- Opt. Anotch Filter Design with Stability Constraints imposed
% This function designs a second order antinotch digital filter by using the
% optimization technique. 
% Here we impose two types of constraints. One, the coefficient d0 must be
% greater than zero, and two, the filter must be stable. 
function [x0, xkp1, iter, ndel] = opt_antinotch_with_stab_for_AutoTuning(w0,epsi,K,NN)

R = 0.95; d00 = R^2; d10 = -cos(pi*w0)*(1+d00); x0 = [d00; d10];
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
    [H,f,A,b] = omgrad_for_AutoTuning(xk, omega, epsi, dw);
    [del, fval, extflg, outpt] = quadprog(H,f,A,b,[],[],[],[],[],options);
    if extflg ~= 1
        outpt, extflg
    end
    xkp1 = xk + del;
    xk = xkp1;
    ndel = norm(del);
    iter = iter + 1;
end

%% --- SubFunction used by "opt_antinotch_with_stab_for_AutoTuning" to compute gradient, H, f, A, and b 
%  for using in QuadProg
function [H,f,A,b] = omgrad_for_AutoTuning(xk, omega, epsi, dw)
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


%% --- Figure1_CloseRequestFcn
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if ~handles.already_in_path && ~isempty(strfind(lower(path), lower(handles.pathofthisfile)))
    rmpath(handles.pathofthisfile);
end
delete(hObject);


