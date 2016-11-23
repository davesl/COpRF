function [results, B, G] = Convex_pRF_fit_parallel(stimulus,data,TR,COpRF_options)
% [results, B, G] = Convex_pRF_fit_parallel(stimulus,data,TR,COpRF_options)
% 
% Function to take pRF stimulus and data variables and fit pRF model using
% Convex optimisation algorithims.
%
% INPUT
% stimulus = A cell array of stimulus matricies. Each cell should contain
%   the stimuli for a given run number. Npixel x Npixel x Ntimepoints.
% data = A cell array of data timeseries. Run order must match with
%   stimulus. Nresponses x Ntimepoints.
% TR = the repetition time in seconds. e.g. time between time intervals in
%   stimulus and data.
% COpRF_options = COpRF_options structure. if [] then default options will
%   be used.
%
% OUTPUT
% results = structure of model fit results
% B = raw model parameter estimates (stimulus refered e.g. rfsize in pixels)
% G = the gamma sparsity regularization parameter used for the lasso fit
%
% 
%% Results description
% ~ The results output structure contains all of the pRF model parameter
%   estimates. These are: polar angle (.ang), eccentricity (.ecc), pRF size
%   (.rfsize), compressive exponent (.expt) and gain (.gain). The x- and
%   y-positions are also saved (.xpos and .ypos).
% ~ The polar angles are in degrees with 0 starting on the positive x-axis
%   and rotating anti-clockwise (e.g. 90 is along positive y-axis).
% ~ Model fitting summaries are included as: coefficient of determination
%   (.R2), sparsity regularizer (.gamma), number of non-zero atom weights
%   (.Nbeta), indices to fit voxels (.vxs), HRF function (.hrf).
% ~ The results.params variable has all of the model parameter estimates
%   and is equivalent to the B output variable. It takes the form
%   Nresponses x Nparams. The parameter ordering is Xpos, Ypos, Sigma, 
%   Gain, Exponent.
% ~ All parameters are defined in stimulus space (e.g. screen pixels). To
%   convert .ecc and .rfsize to visual angle multiply by the total visual
%   angle covered by stimulus (diameter) divided by number of pixels across
%   stimulus (diameter).
% 
%   E.g. if we have a stimulus extending to an eccentricity of 9 degrees 
%   visual angle this covers 18 degrees across the stimulus diameter. We
%   have a total stimulus pixel with of 100 and therefore multiply the
%   .rfsize and .ecc values by C=18/100.
%
%% Notes and suggestions
% ~ The fitting options can be customised before calling this function. Use
%   initCOpRF_options.m to create an options structure. See the
%   documentation of initCOpRF_options.m for details on each option.
% ~ There is no need to detrend the data as this will be performed
%   automatically using polynomial regression. The same detrending will 
%   then be applied to the dictionary.
% ~ To save time it is reccomended to provide a vector of indices 
%   specifiying which timeseries to fit. This could be done with an
%   anatomical mask (for example COpRF_options.vxs=find(MASK)).
% ~ If you set COpRF_options.doGridSearch = true; then the code will first
%   fit the model to all timeseries specified by COpRF_options.vxs using
%   the reduced size sub-dictionary (~3000 atoms). COpRF_options.vxs will
%   then be updated to include only responses where explained variance
%   R2>COpRF_options.R2threshold before running the full model fit. This
%   method can be used to estimate a mask directly from the data but will
%   add a significant amount of time to the model fitting proceedure.
% ~ Rather than generating a new dictionary for every subject in a study it
%   is possible to generate the dictionary once using
%   Build_pRF_dictionary.m and saving the output COpRF_options as a .mat
%   file. These saved options will contain the full dictionary definitons
%   and will not need to be re-computed for each subject. Note: each
%   dictionary is related to a specific stimulus configuration so if the
%   stimuli changes between subjects a new dictionary should be computed
%   each time.

%% Setup
% Transpose data if required
for i=1:length(data)
    if ~(size(stimulus{i},3)==size(data{i},2))
        data{i}=data{i}';
    end
end
% Deal with 3D data
if size(data{1},4)>1
    datasize = size(data{1});
    datasize = datasize(1:3);
    for i=1:length(data)
        data{i} = single(squish(data{i},3));
    end
else
    datasize = size(data{1},1);
end
if isempty(COpRF_options)
    COpRF_options = initCOpRF_options;
end
if isempty(COpRF_options.maxpolydeg)
    COpRF_options.maxpolydeg = cellfun(@(x) round(size(x,2)*TR/60/2),data);
end
if isempty(COpRF_options.vxs)
    COpRF_options.vxs = 1:size(data{1},1);
end

% Create matrix to regress out low frequency drift in time-series
tmatrix = Build_poly_regressors(data,TR,COpRF_options);

Gammas = logspace(log10(0.6),log10(COpRF_options.Min_gamma),COpRF_options.Ngammas); % [1/3, 1/30, 1/300, 1/3000, 1/30000];
res = sizefull(stimulus{1},2);
resmx = max(res);


%% Create dictionary
% Dparams order: X Y S G N (radius, angle, sigma, gain, exponent)
[D, Dparams, COpRF_options, model] = Build_pRF_dictionary(stimulus,data,TR,COpRF_options,0);
% Regress out polys
D(:,3:end)=tmatrix*D(:,3:end); % NB. First 2 atoms in D account for offset
D_norm = repmat( 1./sqrt( sum(D.^2) ), [size(D,1),1] );
Dn = D .* D_norm;
D_norm=D_norm(1,:);

% Free up memory
clear D
COpRF_options.D=''; COpRF_options.Dparams='';

%% Create reduced size sub-dictionary
% The sub-dictionary will have 3 values for the exponent, 21 x/y
% parameters and reduced number of rf-sizes.
COpRF_options_subdict = initCOpRF_options;
COpRF_options_subdict.Nxy = 21;
COpRF_options_subdict.Nexpts = 3; COpRF_options_subdict.Min_expts = 0.1;
COpRF_options_subdict.SigmaScale = 2;
COpRF_options_subdict.hrf = COpRF_options.hrf;
[Dsub, Dsub_params, ~] = Build_pRF_dictionary(stimulus,data,TR,COpRF_options_subdict,0);
Dsub(:,3:end)=tmatrix*Dsub(:,3:end);
Dsub_norm = repmat( 1./sqrt( sum(Dsub.^2) ), [size(Dsub,1),1] );
Dsubn = Dsub .* Dsub_norm;
Dsub_norm=Dsub_norm(1,:);

% Free up memory
clear Dsub
tmatrix = single(tmatrix);

%% Convex optimisation - initialization

if COpRF_options.Nthreads>1
    % Start parallel MATLAB to speed up execution
    partoolbox = which('matlabpool');
    if ~isempty(partoolbox)
        if matlabpool('size')==0
            matlabpool(COpRF_options.Nthreads);
        end
    end
end

% Setup variables for parfor loop
Nthread=COpRF_options.Nthreads;
vxs = COpRF_options.vxs;

% Setup data and normalise
data=catcell(2,data(:))';
data=data(:,vxs);
Data_norm = repmat( 1./sqrt( sum(data.^2) ), [size(data,1),1] );
% Remove bad data
vxs=vxs(~max(Data_norm==inf,[],1));
data=data(:,~max(Data_norm==inf,[],1));
% Regress out polys
data = tmatrix*data;
Data_norm = repmat( 1./sqrt( sum(data.^2) ), [size(data,1),1] );
% We normalise the data such that the range of gammas have a consistent
% effect of promoting atom weight sparsity for all voxels
data = data .* Data_norm;
Data_norm=Data_norm(1,:);

% Define fitting params
R = COpRF_options.SearchRadius;
FitModel = COpRF_options.FitModel;
% ** Advised to use FitModel=1 in the majority of cases **
if FitModel==0
    % FitModel==0 models fitness as the sum of weighted atoms.
    stimulus = '';
    xx = '';
    yy = '';
    % model = '';
else
    % FitModel==1 uses weighted params of atoms to cross-validate with data.
    stimulus = model.stimulus;
    xx = model.xx;
    yy = model.yy;
    res = model.res;
    resmx = model.resmx;
    hrf = model.hrf;
    modelfun = model.modelfun;
    % model = '';
end
Xpos=Dparams(:,1);
Ypos=Dparams(:,2);
GammaStop = COpRF_options.GammaStop;

% Limit vxs to timeseries above R2 threshold
if COpRF_options.doGridSearch
    
    fprintf('Create mask using sparse-grid model fit...\n');
    fprintf([num2str(size(data,2)) ' to fit...\n']);
    parfor_progress(size(data,2)); % Initialize
    tic
    R2grid=nan(size(data,2),1);
    
    % Fit timeseries using reduced size sub dictionary
    for p=1:size(data,2)
        Y = data(:,p);
        [~,~,~,~,~, R2grid(p),~,~,~] = PRF_GridSearch(Y,Dsubn,Dsub_norm,Dsub_params,Gammas,tmatrix,res,resmx,FitModel,stimulus,xx,yy,hrf,modelfun,1);
        parfor_progress; % Count
    end
    Time = toc;
    parfor_progress(0); % Clean up
    disp(['Sparse-grid fit in: ' num2str(round(Time/60)) ' minutes ' num2str(round(rem(Time,60))) ' seconds...'])
    
    % Use results to mask out unresponsive voxels
    vxs = vxs(R2grid>COpRF_options.R2threshold);
    data = data(:,R2grid>COpRF_options.R2threshold);
    Data_norm = Data_norm(R2grid>COpRF_options.R2threshold);
    disp([num2str(length(vxs)) ' passed R2 threshold masking...'])
end

% Setup chunks
ChunkSize=50;
if ChunkSize*Nthread>length(vxs); ChunkSize=round(length(vxs)/Nthread)+1; end
CInd=round(linspace(0,length(vxs),round(length(vxs)/ChunkSize)+1));

for i=1:length(CInd)-1
    vxs_chunk{i}=CInd(i)+1:CInd(i+1);
    data_chunk{i}=data(:,vxs_chunk{i});
end

clear data

% Preallocate variable memory
for i=1:length(data_chunk)
    B1{i} = nan(size(data_chunk{i},2),1);
    B2{i} = nan(size(data_chunk{i},2),1);
    B3{i} = nan(size(data_chunk{i},2),1);
    B4{i} = nan(size(data_chunk{i},2),1);
    B5{i} = nan(size(data_chunk{i},2),1);
    R_2{i} = nan(size(data_chunk{i},2),1);
    Gamma{i} = nan(size(data_chunk{i},2),1);
end

fprintf('Starting Convex model fitting...\n');
fprintf([num2str(length(vxs)) ' to fit...\n']);
parfor_progress(length(vxs)); % Initialize

tic
%% Convex optimisation - fitting (parallel and serial)
if COpRF_options.Nthreads>1
    parfor i=1:length(data_chunk)
        for p=1:size(data_chunk{i},2)
            
            % Voxel timeseries
            Y = data_chunk{i}(:,p);
            
            % Fit the model parameters
            [B1{i}(p), B2{i}(p), B3{i}(p), B4{i}(p), B5{i}(p), R_2{i}(p), Gamma{i}(p), Nbeta{i}(p), ~] ...
                = Fit_COpRF(Y,Dn,D_norm,Dparams,Dsubn,Dsub_norm,Dsub_params,Gammas,tmatrix, ...
                Xpos,Ypos,res,resmx,R,FitModel,stimulus,xx,yy,hrf,modelfun,GammaStop);
            
            parfor_progress; % Count
            
        end
    end
else
    for i=1:length(data_chunk)
        for p=1:size(data_chunk{i},2)
            
            % Voxel timeseries
            Y = data_chunk{i}(:,p);
            
            % Fit the model parameters
            [B1{i}(p), B2{i}(p), B3{i}(p), B4{i}(p), B5{i}(p), R_2{i}(p), Gamma{i}(p), Nbeta{i}(p), ~] ...
                = Fit_COpRF(Y,Dn,D_norm,Dparams,Dsubn,Dsub_norm,Dsub_params,Gammas,tmatrix, ...
                Xpos,Ypos,res,resmx,R,FitModel,stimulus,xx,yy,hrf,modelfun,GammaStop);
            
            parfor_progress; % Count
            
        end
    end
end

Time = toc;
parfor_progress(0); % Clean up
disp(['Finished Convex fit in: ' num2str(round(Time/60)) ' minutes ' num2str(round(rem(Time,60))) ' seconds...'])


%% Combine parfor chunks and prepare output
Bvx=zeros(length(vxs),5);
R2vx=zeros(length(vxs),1);
Gvx=zeros(length(vxs),1);
Nbetavx=zeros(length(vxs),1);
FSvx.R2=cell(length(vxs),1);
FSvx.Nx=cell(length(vxs),1);
FSvx.Gammas=cell(length(vxs),1);
for i=1:length(data_chunk)
    Bvx(CInd(i)+1:CInd(i+1),1)=B1{i};
    Bvx(CInd(i)+1:CInd(i+1),2)=B2{i};
    Bvx(CInd(i)+1:CInd(i+1),3)=B3{i};
    Bvx(CInd(i)+1:CInd(i+1),4)=B4{i};
    Bvx(CInd(i)+1:CInd(i+1),5)=B5{i};
    R2vx(CInd(i)+1:CInd(i+1),1)=R_2{i};
    Gvx(CInd(i)+1:CInd(i+1),1)=Gamma{i};
    Nbetavx(CInd(i)+1:CInd(i+1),1)=Nbeta{i};
end
Bvx(:,4) = Bvx(:,4)./Data_norm';

% Precompute vectors
clear results;
results.ang = nan([prod(datasize) 1]);
results.ecc = nan([prod(datasize) 1]);
results.expt = nan([prod(datasize) 1]);
results.rfsize = nan([prod(datasize) 1]);
results.R2 = nan([prod(datasize) 1]);
results.gain = nan([prod(datasize) 1]);
results.Gamma = nan([prod(datasize) 1]);
results.Nbeta = nan([prod(datasize) 1]);

% Save model fit results in results structure
results.ang(vxs,:) = mod(atan2((1+res(1))/2 - Bvx(:,1), Bvx(:,2) - (1+res(2))/2),2*pi)/pi*180;
results.ecc(vxs,:) = sqrt(((1+res(1))/2 - Bvx(:,1)).^2 + (Bvx(:,2) - (1+res(2))/2).^2);
results.expt(vxs,:) = Bvx(:,5);
results.rfsize(vxs,:) = abs(Bvx(:,3)) ./ sqrt(Bvx(:,5));
results.R2(vxs,:) = R2vx;
results.Gamma(vxs,:) = Gvx;
results.gain(vxs,:) = Bvx(:,4);
results.xpos(vxs,:) = Bvx(:,1) - (1+res(1))/2;
results.ypos(vxs,:) = Bvx(:,2) - (1+res(2))/2;
results.Nbeta(vxs,:) = Nbetavx;

results.ang = reshape(results.ang, [datasize 1]);
results.ecc = reshape(results.ecc, [datasize 1]);
results.expt = reshape(results.expt, [datasize 1]);
results.rfsize = reshape(results.rfsize, [datasize 1]);
results.R2 = reshape(results.R2, [datasize 1]);
results.gain = reshape(results.gain, [datasize 1]);
results.Gamma = reshape(results.Gamma, [datasize 1]);
results.Nbeta = reshape(results.Nbeta, [datasize 1]);

results.params = Bvx;
results.vxs = vxs;
results.hrf = COpRF_options.hrf;
results.COpRF_options = COpRF_options;

B = Bvx;
G=Gvx;

