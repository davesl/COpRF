function COpRF_options = initCOpRF_options
% CO-pRF default setup options.
% Future releases will add to the number of model types and pRF modalities
% which can estimated using the toolbox.
% 
% Options are currenly set to fit CSS-pRF (Kay et al., 2013)or basic pRF 
% (Dumoulin and Wandell 2008) models.
% To fit the basic pRF model eliminate the expoent component by setting 
% options.Min_expts=1; and options.Nexpts=1;
% 
% Changing the option parameters can allow a trade-off between model
% fitting speed and estimation accuracy.

%% Basic setup
COpRF_options.vxs         = '';   % Voxel indices to fit. Regions likely to be responsive can be pre-masked to save time.
COpRF_options.maxpolydeg  = '';   % Number of polynomial degrees used for detrending
COpRF_options.hrf         = [];   % A custom hrf function can be supplied. Note: hrf points should be spaced 1-TR apart (same as data and stimulus)
COpRF_options.fitHRF      = false; % Fit HRF from data or use canonical HRF (default).
COpRF_options.Nthreads    = 1;    % Number of threads to run in parallel. Note: addtional threads can create large memory demands...

%% Parameter distributions for dictionary construction
% x-y (r-theta) setup
COpRF_options.Nxy       = 51;     % Number of points to model along x and y. e.g. options.Nxy=60 creates a 60x60 grid.
COpRF_options.MaxEcc      = 1.2;  % Max eccentricity - percentage of stimulus radius e.g. 120% by default
% rfsize (sigma) setup
COpRF_options.SigmaMax    = 0.64; % Largest receptive field to model as percentage of stimulus coverage e.g. 64% by default
COpRF_options.SigmaScale  = 0.7;  % Determines the number of receptive field sizes to model. Smaller values create larger dictionaries.
% Compressive exponent (n) setup
COpRF_options.Min_expts   = 0.1;  % Minimum exponent value to use
COpRF_options.Skew_expts  = 1;    % Exponent skew factor
COpRF_options.Nexpts      = 6;    % Number of exponent values to model

%% Fitting options
COpRF_options.Ngammas     = 10;   % Number of regularization parameters to try for the lasso
COpRF_options.Min_gamma   = 1/1000; % Minimum gamma to try... Gammas=logspace(log10(0.6),log10(options.Min_gamma),options.Ngammas);
COpRF_options.GammaStop   = true; % Will stop trying new regularization terms once the fit from the addition of atoms has consistently detriorated.
COpRF_options.SearchRadius= 1/10; % Percentage of stimulus diameter to use for second stage search radius
COpRF_options.FitModel    = true; % FitModel=true will estimate model fitness from weighted atom parameters. options. FitModel=false will estimate directly from the averaged atoms and is not recommended.
COpRF_options.doGridSearch = false; % If true the smaller sub-dictionary will be used to determine options.vxs and eliminate
COpRF_options.R2threshold = 0.1;  % Voxels with grid-search R^2 less than options.R2threshold will be masked out

%% Precomputed dictionary options
% *** Note: a precomputed dictionary must be built from the same stimuli used during the data collection ***
COpRF_options.D           = ''; % Can supply a precomputed dictionary to save time. Use Build_pRF_dictionary.m
COpRF_options.Dparams     = ''; % The model parameters used to create options.D.
