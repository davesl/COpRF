function [D, Dparams,COpRF_options,model] = Build_pRF_dictionary(stimulus,data,TR,COpRF_options,NewD)
% Function to build appropriate dictionary for stimulus presentation.
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
% NewD = if true will compute the full dictionary regardless of whether one
%   already within the COpRF_options structure.
%
% OUTPUT
% D = the full size dictionary defined by the COpRF_options.
% Dparams = the parameter values used to build D. These are ordered: Xpos,
%   Ypos, Sigma, Gain, Exponent.
% COpRF_options = the options structure which now contains D and Dparams.
% model = definitions of the pRF model used. Note: in current release this
%   applies to the visual CSS-pRF model.

warning off

%% Setup
if isempty(COpRF_options)
    COpRF_options = initCOpRF_options;
end
res = sizefull(stimulus{1},2);
resmx = max(res);
% prepare stimuli
for ii=1:length(stimulus)
    stimulus{ii} = squish(stimulus{ii},2)';  % frames x pixels
    stimulus{ii} = [stimulus{ii} ii*ones(size(stimulus{ii},1),1)];  % add a dummy column to indicate run breaks
    stimulus{ii} = single(stimulus{ii});  % make single to save memory
end

if isempty(COpRF_options.hrf)
    if COpRF_options.fitHRF
        % Transpose data if required
        for i=1:length(data)
            if ~(size(stimulus{i},3)==size(data{i},2))
                data{i}=data{i}';
            end
        end
        [~,COpRF_options.hrf] = analyzePRFcomputeGLMdenoiseregressors(stimulus,data,TR,0,1);
    else
        COpRF_options.hrf = getcanonicalhrf(TR,TR)';
    end
end


%% PREPARE MODEL - following the CSS-pRF convention: https://github.com/kendrickkay/analyzePRF
% pre-compute some cache
[~,xx,yy] = makegaussian2d(resmx,2,2,2,2);
% define the CSS-pRF model (parameters are: Xpos Ypos Sigma Gain Exponent)
modelfun = @(pp,ss) conv2run(posrect(pp(4)) * (ss*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),COpRF_options.hrf,ss(:,prod(res)+1));

%% Setup param vecs
x_pos = linspace(-0.5,0.5,COpRF_options.Nxy);
y_pos = linspace(-0.5,0.5,COpRF_options.Nxy);
expts = linspace(1,COpRF_options.Min_expts^(1/COpRF_options.Skew_expts),COpRF_options.Nexpts).^COpRF_options.Skew_expts;

% calculate rfsizes to model
maxn = log2(resmx*COpRF_options.SigmaMax);
ss_vec = 0:COpRF_options.SigmaScale:maxn;
if ~(rem(maxn/COpRF_options.SigmaScale,1)==0); ss_vec = [ss_vec maxn]; end
ssindices = 2.^(ss_vec);

fprintf('constructing parameter grid.\n');

Count = 1;
% When using x-position and y-position
Dparams = zeros(length(x_pos)*length(y_pos)*length(ssindices)*length(expts),5);
for aa=1:length(x_pos)
    for bb=1:length(y_pos)
        for cc=1:length(ssindices)
            for dd=1:length(expts)
                
                Dparams(Count,:) = ...
                    [x_pos(aa)*res(1)*COpRF_options.MaxEcc + (res(1)+1)/2 ... % x
                    y_pos(bb)*res(2)*COpRF_options.MaxEcc + (res(2)+1)/2  ... % y
                    ssindices(cc)*sqrt(expts(dd)) 1 expts(dd)];     ... % Sigma, gain, exponent
                    
                Count = Count + 1;
            end
        end
    end
end
eccs = (((Dparams(:,1)-(res(1)+1)/2)).^2 + (Dparams(:,2)-(res(2)+1)/2).^2).^0.5;
Dparams=Dparams(eccs<=resmx/2*COpRF_options.MaxEcc,:); % Use circular x,y grid

%% Build dictionary
D = zeros(sum(cellfun(@(x) size(x,1),stimulus)),size(Dparams,1),'single');  % time x seeds
temp = catcell(1,stimulus);
fprintf('generating dictionary time-series...'); tic
if isempty(COpRF_options.D) || NewD
    
    % Only keep atoms where pRF overlaps with stimulus area
    eccs = (((Dparams(:,1)-(res(1)+1)/2)).^2 + (Dparams(:,2)-(res(2)+1)/2).^2).^0.5;
    Inds = Dparams(:,3)./sqrt(Dparams(:,5))>eccs-(resmx-1)/2;
    Dparams = Dparams(Inds,:);
    
    for ii=1:size(Dparams,1)
        D(:,ii) = modelfun(Dparams(ii,:),temp);
    end
    
    D=cat(2,ones(1,size(D,1))',D); % Add atoms to account for +ve/-ve offset
    D=cat(2,-ones(1,size(D,1))',D);
    Dparams=single(cat(1,nan(2,size(Dparams,2)),Dparams));
else
    D = COpRF_options.D;
    Dparams = COpRF_options.Dparams;
end

COpRF_options.D = D;
COpRF_options.Dparams = Dparams;

disp(['Dictionary size of ' num2str(size(Dparams,1)) ' stimuli responsive elements'])

%% Create model structure for future fitting
model.stimulus = temp;
model.xx = xx;
model.yy = yy;
model.res = res;
model.resmx = resmx;
model.hrf = COpRF_options.hrf;
% define the model (parameters are R C S G N)
model.modelfun = @(pp,dd,res,resmx,xx,yy,hrf) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));

clear temp;
fprintf('done.'); toc

warning on


