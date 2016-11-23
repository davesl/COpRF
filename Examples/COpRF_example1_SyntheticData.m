%% CO-pRF Example1: synthetic data
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

% Ignore first 2 columns - these are added to dictionaries to correct for 
% non mean-centred data
SyntheticData = SyntheticData(:,3:end);
SyntheticParams = SyntheticParams(3:end,:);

%% 2. Add Gaussian white noise
% We will add some Gaussian noise to each simulated timeseries

snr=5; % Arbitrarily choise an SNR of 5:1
for i=1:size(SyntheticData,2)
    y = SyntheticData(:,i);
    SyntheticData(:,i) = awgn(y,snr);
end

%% 3. Create dense dictionary for model fitting
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

% Create a cell array containing the single data session with time along
% second dimention
data = {SyntheticData'};

% Build the CO-pRF dictionary
[Dictionary, DictionaryParams, COpRF_options] = Build_pRF_dictionary(stimulus,data,TR,COpRF_options,1);
COpRF_options.D = Dictionary;
COpRF_options.Dparams = DictionaryParams;

% If we save COpRF_options the same options file and dictionary could be
% loaded and used for other datasets using the same stimuli.
save(fullfile(ExamplesDir,'Example_COpRF_options.mat'),'COpRF_options') 


%% 4. CO-pRF model fitting
% To estimate model parameters from the data we need to supply stimulus and
% data variables. These should take the form of cell arrays with presented
% stimuli and recorded timeseries in a separate cell for each session.
% 
% Stimuli should be Npixel x Npixels x Ntimepoints
% Data should be Nvoxels x Ntimepoints

% Begin the CO-pRF model fitting
[results, B, G] = Convex_pRF_fit_parallel(stimulus,data,TR,COpRF_options);


%% 5. Convert from pixels to degrees
% The COpRF toolbox estimates model parameters in a stimuli defined space.
% However converting to visual angle is trivial.

StimEcc = 9; % Assume stimuli extends to 9 degrees of eccentricity.
StimRes=100; % Stimulus resolution (diameter)
Cfactor = (StimEcc*2)/StimRes;

results.rfsize_degrees = results.rfsize*Cfactor;
results.ecc_degrees = results.ecc*Cfactor;
results.xpos_degrees = results.xpos*Cfactor;
results.ypos_degrees = results.ypos*Cfactor;


%% 6. Plot estimated parameters vs. synthetic groundtruth
% We will create some simple plots to show the accuracy of the CO-pRF 
% model fit.

% Define CSS-pRF model function (Kay et al., 2013)
[~,xx,yy] = makegaussian2d(StimRes,2,2,2,2);
modelfun = @(pp,ss) conv2run(posrect(pp(4)) * (ss*[vflatten(placematrix(zeros([StimRes StimRes]),makegaussian2d(StimRes,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),COpRF_options.hrf,ss(:,StimRes^2+1));
% Create matrix to regress out low frequency drift in time-series
tmatrix = Build_poly_regressors(data,TR,COpRF_options);
% prepare stimuli
SS=[];
for ii=1:length(stimulus)
    SS{ii} = squish(stimulus{ii},2)';  % frames x pixels
    SS{ii} = [SS{ii} ii*ones(size(SS{ii},1),1)];  % add a dummy column to indicate run breaks
end
SS=catcell(1,SS);

% Select 10 example timeseries
N=10;
I= [331   278    42   363   626   198   241   156   219   115];
% I=randsample(size(data{1},2),N); % Or randomize selection

% Plot data with estimated parameter predictions
figure('Position',[100 80, 800, 900])
for i=1:N
    subplot(N,1,i)
    % Plot synthetic data
    Y1=tmatrix*data{1}(I(i),:)';
    plot(Y1,'b')
    hold on
    % Plot prediction from parameter estimates
    Y2 = tmatrix*modelfun(B(I(i),:),SS);
    plot(Y2,'r','LineWidth',2)
end

% Plot the groundtruth (blue) and estimated (red) pRFs. 
% pRFs with large receptive fields, extreme eccentricity and/or small 
% compressive exponets are challenging to estimate and will match more 
% poorly to the groundtruth than pRFs which overlap entirely with the 
% stimulus area.
figure('Position',[1100 200, 600, 600]); hold on; axis square; grid on;
GT_xpos_degrees=(SyntheticParams(:,1)-(StimRes+1)/2)*Cfactor;
GT_ypos_degrees=(SyntheticParams(:,2)-(StimRes+1)/2)*Cfactor;
GT_rfsize_degrees=SyntheticParams(:,3)./sqrt(SyntheticParams(:,5))*Cfactor;
ang=0:0.01:2*pi;
xp=StimEcc*cos(ang);
yp=StimEcc*sin(ang);
plot(0+xp,0+yp,'k','LineWidth',2); % Stimulus extent
for i=1:N
    % Plot groundtruth
    x=GT_xpos_degrees(I(i));
    y=GT_ypos_degrees(I(i));
    r=GT_rfsize_degrees(I(i));
    scatter(x,y,10,'b','fill')
    xp=r*cos(ang);
    yp=r*sin(ang);
    plot(x+xp,y+yp,'b');
    % Plot estimate
    x=results.xpos_degrees(I(i));
    y=results.ypos_degrees(I(i));
    r=results.rfsize_degrees(I(i));
    scatter(x,y,10,'r','fill')
    xp=r*cos(ang);
    yp=r*sin(ang);
    plot(x+xp,y+yp,'r');
end
xlim([-StimEcc StimEcc]); ylim([-StimEcc StimEcc]);

%% References
% Original pRF model: Dumoulin, S.O., Wandell, B. a., 2008. Population receptive field estimates in human visual cortex. Neuroimage 39, 647–660.
% Modified CSS-pRF model: Kay, K.N., Winawer, J., Mezer, A., Wandell, B. a, 2013. Compressive spatial summation in human visual cortex. J. Neurophysiol. 110, 481–94.
% Useful code for desiging pRF fMRI experiments and nonlinear model fitting: http://cvnlab.net/analyzePRF/
