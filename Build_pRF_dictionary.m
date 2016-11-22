function [D, Dparams,options,model] = Build_pRF_dictionary(stimulus,data,TR,options,NewD)
% Function to build appropriate dictionary for stimulus presentation.
%
% ParamGridding = is cell array of parameter vectors to use in dictionary.
%           The order of the parameters is: R C S G N (radius, angle, sigma, gain, exponent)

warning off

%% Setup
if isempty(options)
    options = initOptions;
end
res = sizefull(stimulus{1},2);
resmx = max(res);
% prepare stimuli
for ii=1:length(stimulus)
    stimulus{ii} = squish(stimulus{ii},2)';  % frames x pixels
    stimulus{ii} = [stimulus{ii} ii*ones(size(stimulus{ii},1),1)];  % add a dummy column to indicate run breaks
    stimulus{ii} = single(stimulus{ii});  % make single to save memory
end

if isempty(options.hrf)
    if options.fitHRF
        [~,options.hrf] = analyzePRFcomputeGLMdenoiseregressors(stimulus,data,TR,0,1);
    else
        options.hrf = getcanonicalhrf(TR,TR)';
    end
end


%% PREPARE MODEL - following the CSS-pRF convention: https://github.com/kendrickkay/analyzePRF
% pre-compute some cache
[~,xx,yy] = makegaussian2d(resmx,2,2,2,2);
% define the CSS-pRF model (parameters are R C S G N)
modelfun = @(pp,ss) conv2run(posrect(pp(4)) * (ss*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),options.hrf,ss(:,prod(res)+1));

%% Setup param vecs
x_pos = linspace(-0.5,0.5,options.Nxy);
y_pos = linspace(-0.5,0.5,options.Nxy);
expts = linspace(1,options.Min_expts^(1/options.Skew_expts),options.Nexpts).^options.Skew_expts;

% calculate rfsizes to model
maxn = log2(resmx*options.SigmaMax);
ss_vec = 0:options.SigmaScale:maxn;
if ~(rem(maxn/options.SigmaScale,1)==0); ss_vec = [ss_vec maxn]; end
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
                    [x_pos(aa)*res(1)*options.MaxEcc + (res(1)+1)/2 ... % x
                    y_pos(bb)*res(2)*options.MaxEcc + (res(2)+1)/2  ... % y  
                    ssindices(cc)*sqrt(expts(dd)) 1 expts(dd)];     ... % Sigma, gain, exponent
                
                Count = Count + 1;
            end
        end
    end
end
eccs = (((Dparams(:,1)-(res(1)+1)/2)).^2 + (Dparams(:,2)-(res(2)+1)/2).^2).^0.5;
Dparams=Dparams(eccs<=resmx/2*options.MaxEcc,:); % Use circular x,y grid

%% Build dictionary
D = zeros(sum(cellfun(@(x) size(x,1),stimulus)),size(Dparams,1),'single');  % time x seeds
temp = catcell(1,stimulus);
fprintf('generating dictionary time-series...'); tic
if isempty(options.D) || NewD
    
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
    D = options.D;
    Dparams = options.Dparams;
end

options.D = D;
options.Dparams = Dparams;

disp(['Dictionary size of ' num2str(size(Dparams,1)) ' stimuli responsive elements'])

%% Create model structure for future fitting
model.stimulus = temp;
model.xx = xx;
model.yy = yy;
model.res = res;
model.resmx = resmx;
model.hrf = options.hrf;
% define the model (parameters are R C S G N)
model.modelfun = @(pp,dd,res,resmx,xx,yy,hrf) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));

clear temp;
fprintf('done.'); toc

warning on


