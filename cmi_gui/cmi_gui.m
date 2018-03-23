function varargout = cmi_gui(varargin)
% CMI_GUI MATLAB code for cmi_gui.fig
%      CMI_GUI, by itself, creates a new CMI_GUI or raises the existing
%      singleton*.
%
%      H = CMI_GUI returns the handle to a new CMI_GUI or the handle to
%      the existing singleton*.
%
%      CMI_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CMI_GUI.M with the given input arguments.
%
%      CMI_GUI('Property','Value',...) creates a new CMI_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cmi_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cmi_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cmi_gui

% Last Modified by GUIDE v2.5 23-Mar-2018 11:11:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cmi_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @cmi_gui_OutputFcn, ...
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


% --- Executes just before cmi_gui is made visible.
function cmi_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cmi_gui (see VARARGIN)

% Choose default command line output for cmi_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cmi_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global cmiObj C

add_files_to_path()

cmiObj = CMIclass();
C = {};

function add_files_to_path()
% Adds the folders in this directory to path for when we cd
addpath('.', 'helpers', 'serial_analysis_scripts')
addpath(genpath('cmi_R2015a-master'))

% --- Outputs from this function are returned to the command line.
function varargout = cmi_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in set_directory_btn.
function set_directory_btn_Callback(hObject, eventdata, handles)
% hObject    handle to set_directory_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
debug_print("Setting Directory")
dir_location = select_dir_fn();
cd(dir_location)
debug_print("Directory set!")

% --- Executes on button press in load_c_btn.
function load_c_btn_Callback(hObject, eventdata, handles)
% hObject    handle to load_c_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global C
debug_print("Importing C from CSV");
file_location = select_file_fn();
if file_location == 0
    % File select was cancelled
    return
end
C = read_csv(file_location);
debug_print("Finished importing C from CSV")


% --- Executes on button press in serial_analysis_btn.
function serial_analysis_btn_Callback(hObject, eventdata, handles)
% hObject    handle to serial_analysis_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global C cmiObj
if invalid_C(C)
    return
end

debug_print('Running Serial Analysis')
out = serial_analysis_v2(cmiObj, C);
debug_print('Finished Serial Analysis')


% --- Executes on button press in tprm_btn.
function tprm_btn_Callback(hObject, eventdata, handles)
% hObject    handle to tprm_btn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global C cmiObj
if invalid_C(C)
    return
end

debug_print('Running tPRM Script')

debug_print('Finished tPRM Script')