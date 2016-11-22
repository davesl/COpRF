%% CO-pRF synthetic data example
% In this script we will look at how to generate some simple synthetic data
% and estimate pRF parameters using the CO-pRF toolbox.


%% Setup stimuli and COpRF_options

% Find example stimuli and load
ExampleStimuliPath = which('Example_COpRF_Stimuli.mat');
ExamplesDir = fileparts(ExampleStimuliPath);
load(ExampleStimuliPath)
% Load default options
COpRF_options = initCOpRF_options;
% Set the repetition time between data points to be 1 second
TR=1;


%% 1. Generate synthetic data
% We create some synthetic data using a small number of parameter
% combinations. Using the following options we create 1082 noise-free
% timeseries.

SyntheticOptions = initCOpRF_options;   % load default options - open initCOpRF_options.m for details
SyntheticOptions.Nxy = 10;              % size of x/y grid
SyntheticOptions.MaxEcc = 0.5;          % eccentricity up to 60% of stimulus coverage
SyntheticOptions.SigmaMax = 0.25;       % max rfsize is 25% of stimulus coverage
SyntheticOptions.SigmaScale = 1;        % param to control number of unique rfsizes
SyntheticOptions.Min_expts = 0.25;      % minimum compressive exponent to model
SyntheticOptions.Nexpts = 3;            % Number of exponent values to model
[SyntheticData, SyntheticParams] = Build_pRF_dictionary(stimulus,'',TR,SyntheticOptions,1);


%% 2. Create dense dictionary for model fitting
% Generating the dictionary can take approx. 5-10mins. However we only 
% need to precomute the dictionary once for a given stimuli. For a study 
% with many subjects who were each presented with the same stimuli the 
% dictionary can be saved and reused when fitting each dataset.

% % If you would like to fit the original Dumoulin and Wandell (2008) pRF
% % model set the following parameters to equal 1. By default the options
% % will model a compressive expoent and thus follow the model of Kay et
% al. (2013)
% COpRF_options.Nexpts = 1;
% COpRF_options.Min_expts = 1;

% Generate the 
[Dictionary, DictionaryParams] = Build_pRF_dictionary(stimulus,'',TR,COpRF_options,1);
COpRF_options.D = Dictionary;
COpRF_options.Dparams = DictionaryParams;

save(fullfile(ExamplesDir,'Example_COpRF_options.mat'),'COpRF_options') 
% To save time Example_COpRF_options.mat could now be loaded and used for step 3.


%% 3. CO-pRF model fitting
% To estimate model parameters from the data we need to supply stimulus and
% data variables. These should take the form of cell arrays with presented
% stimuli and recorded timeseries in a separate cell for each session.
% Stimuli should be Npixel x Npixels x Ntimepoints
% Data should be Nvoxels x Ntimepoints

data = {SyntheticData(:,3:end)}; % Create a cell array containing the single data session with time
fitHRF = 0; % Use a cannonical HRF rather than estimating from the data

% Begin the model CO-pRF fitting
[results, B, G] = Convex_pRF_fit_parallel(stimulus,data,TR,COpRF_options);


%% 4. Plot estimated parameters vs. synthetic groundtruth
% We will create some simple scatter plots to show the accuracy of the
% CO-pRF model fit.

vxs = results.vxs; % Indices for modeled timeseries in data

% Groundtruth parameter values
X_GT = SyntheticParams(vxs+2,1); % groundtruth x positions
Y_GT = SyntheticParams(vxs+2,2); % groundtruth y positions
Sigma_GT = SyntheticParams(vxs+2,3); % groundtruth rfsize positions
Exponent_GT = SyntheticParams(vxs+2,5); % groundtruth exponent positions

% Estimated parameter values
X_fit = B(vxs,1); % estimated x positions
Y_fit = B(vxs,2); % estimated y positions
Sigma_fit = B(vxs,3); % estimated rfsize positions
Exponent_fit = B(vxs,5); % estimated exponent positions

% Create scatter plots of estimated vs. groundtruth
figure;
scatter(X_fit,X_GT)
figure;
scatter(Y_fit,Y_GT)
figure;
scatter(Sigma_fit,Sigma_GT)
figure;
scatter(Exponent_fit,Exponent_GT)


%% 5. Convert from pixels to degrees
% The COpRF toolbox estimates model parameters in a stimuli defined space.
% However converting to visual angle is trivial.


%% References
% Original pRF model: Dumoulin, S.O., Wandell, B. a., 2008. Population receptive field estimates in human visual cortex. Neuroimage 39, 647–660.
% Modified CSS-pRF model: Kay, K.N., Winawer, J., Mezer, A., Wandell, B. a, 2013. Compressive spatial summation in human visual cortex. J. Neurophysiol. 110, 481–94.
% Useful code for desiging pRF fMRI experiments and nonlinear model fitting: http://cvnlab.net/analyzePRF/
