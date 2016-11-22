%%%%%%%%%%%%%%%%%%%%%%%%%% EXPERIMENT PARAMETERS (edit as necessary)

% display
ptres = []; % [1280 960 60 32]; % [1152 864 60 32]; % [1024 768 60 32];  % display resolution. [] means to use current display resolution.

% fixation dot
fixationinfo = {uint8([255 0 0; 0 0 0; 255 255 255]) 0.5};  % dot colors and alpha value
fixationsize = 4;          % dot size in pixels
meanchange = 3;            % dot changes occur with this average interval (in seconds)
changeplusminus = 2;       % plus or minus this amount (in seconds)

% trigger
triggerkey = '5';          % stimulus starts when this key is detected
tfun = @() fprintf('STIMULUS STARTED.\n');  % function to call once trigger is detected
TR = 1.518;
numdummies = 5;
dummy_pause = TR*numdummies;
DisableKeysForKbCheck([]);

% tweaking
offset = [-20 -48];        % [X Y] where X and Y are the horizontal and vertical
                           % offsets to apply.  for example, [5 -10] means shift 
                           % 5 pixels to right, shift 10 pixels up.
movieflip = [0 0];         % [A B] where A==1 means to flip vertical dimension
                           % and B==1 means to flip horizontal dimension

% directories
stimulusdir = 'C:\Users\dslater\Documents\MATLAB\Software_packages\pRF_fitting\Task';         % path to directory that contains the stimulus .mat files

%%%%%%%%%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW

% set rand state
rand('state',sum(100*clock));
randn('state',sum(100*clock));

% ask the user what to run
if ~exist('subjnum','var') 
  subjnum = input('What is the subj id? ')
end
% expnum = input('What experiment (89=CCW, 90=CW, 91=expand, 92=contract, 93=multibar, 94=wedgeringmash)? ')
expnum = input('What experiment (898=LREN default, 899=LREN default + HRF mapping)? ')
disp(['Using the LREN default stimuli number ' num2str(expnum)])
runnum = input('What run number (for filename)? ')

% prepare inputs
trialparams = [];
ptonparams = {ptres,[],0};
dres = [];
frameduration = 4;
grayval = uint8(127);
iscolor = 1;
soafun = @() round(meanchange*(60/frameduration) + changeplusminus*(2*(rand-.5))*(60/frameduration));

% load specialoverlay
a1 = load(fullfile(stimulusdir,'fixationgrid.mat'));

% some prep
if ~exist('images','var')
  images = [];
  maskimages = [];
end
filename = sprintf('%s_subj%d_run%02d_exp%02d.mat',gettimestring,subjnum,runnum,expnum);

% run experiment
[images,maskimages] = ...
  showmulticlass(filename,offset,movieflip,frameduration,fixationinfo,fixationsize,tfun, ...
                 ptonparams,soafun,0,images,expnum,[],grayval,iscolor,[],[],[],dres,triggerkey, ...
                 [],trialparams,[],maskimages,a1.specialoverlay,stimulusdir,dummy_pause,TR);

%%%%%%%%%%%%%%%%%%%%%%%%%%

% KK notes:
% - remove performance check at end
% - remove resampling and viewingdistance stuff
% - remove hresgrid and vresgrid stuff
% - hardcode grid and pregenerate
% - trialparams became an internal constant
