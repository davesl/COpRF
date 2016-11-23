function [B,Y2,R2,Gamma,DInds,Nx,FitSummary]=ConvexEstimate(Y,Dn,Dnorm,Dparams,Gammas,FitModel,tmatrix,stimulus,xx,yy,res,resmx,hrf,modelfun,GammaStop)
% [B,Y2,R2,Gamma,DInds,Nx,FitSummary]=ConvexEstimate(Y,Dn,Dnorm,Dparams,Gammas,FitModel,tmatrix,stimulus,xx,yy,res,resmx,hrf,modelfun,GammaStop)
% 
% Takes a set of pRF input variables and estimates a convex model fit using
% the LASSO algorithm.
% Once atoms are selected at a given Gamma value the LASSO is run again
% without regularasation over subset of selected atoms. This is a debiasing
% step to improve weight estimates (Figueiredo et al. 2007).
%
% INPUT
% Y = timeseries to fit.
% Dn = normalised dictionary
% Dnorm = a vector of weights used to normalise atoms in Dn
% Dparams = the set of parameters used to create the dictionary
% Gammas = set of gamma values to try as the sparsity regularizer. If
%   multiple values provided the estimate with greatest R2 will be chosen.
% FitModel = if true weighted average of Dparams to estimate model fitness.
%   If falue uses atom sums (not recommended)
% tmatrix = polynomial regressor matrix
% stimulus = stimulus corresponding to the timeseries Y
% xx/yy = descriptors of guassian receptive field
% res = stimulus resolution
% resmx = max(res)
% hrf = hrf function to use in pRF model
% modelfun = a handle to the pRF model function
% GammaStop = if true will stop when reducing sparsity impairs the fit
% 
% OUTPUT
% B = estimated model parameters
% Y2 = modeled timeseries
% R2 = coefficient of determination (% explained variance)
% Gamma = chosen sparsity regularizer value
% DInds = non-zero atom indices
% Nx = number of selected atoms (with non-zero weights)
% FitSummary = summary of data fit and atom selection


%% Setup
SPAMS_params.mode    = 2;
SPAMS_params.pos     = true;
SPAMS_params.lambda2 = 0;
R2_old = -inf;
R2_new  = -inf;

%% Loop over all sparsity parameters - Gammas
for j=1:length(Gammas)
    if R2_old<=R2_new
        R2_old = R2_new;
        
        % First Lasso - pick atoms
        SPAMS_params.lambda = Gammas(j);
        x = mexLasso( Y, Dn, SPAMS_params );
        DInds{j}=find(x);
        if isempty(DInds{j})
            Y2{j} = nan;
            R2(j) = -inf;
            continue;
        end
        
        % Second Lasso - no regularisation on selected atoms
        if nnz(DInds{j})>1
            SPAMS_params.lambda = 0;
            x = full(mexLasso( Y, Dn(:,DInds{j}), SPAMS_params ));
        else
            x = full(x(DInds{j}));
        end
        
        % Calculate offset and remove its influence on parameter estimation
        if nnz(DInds{j}==1 & DInds{j}==2)
            Offset(j) = -x(1)*Dnorm(1)+x(2)*Dnorm(2)*1;
            DInds{j} = DInds{j}(3:end);
            x = x(3:end);
        elseif nnz(DInds{j}==2)
            Offset(j) = x(1)*Dnorm(2);
            DInds{j} = DInds{j}([2:end]);
            x = x([2:end]);
        elseif nnz(DInds{j}==1)
            Offset(j) = -x(1)*Dnorm(1);
            DInds{j} = DInds{j}(2:end);
            x = x(2:end);
        else
            Offset(j) = 0;
        end
        
        % If no atoms other than offset continue
        if ~nnz(x)
            Y2{j} = nan;
            R2(j) = -inf;
            Nx(j) = nan;
            continue;
        end
        
        % Gain is the weighted sum of normalisation factors
        Gain = sum(Dnorm(DInds{j})'.*x);
        
        % Derive parameter estimates from weighted atoms
        x2 = x ./ ( sum(x) + eps );
        B{j}(1) = Dparams(DInds{j},1)' * x2;
        B{j}(2) = Dparams(DInds{j},2)' * x2;
        B{j}(3) = Dparams(DInds{j},3)' * x2;
        B{j}(4) = Gain;
        B{j}(5) = Dparams(DInds{j},5)' * x2;
        
        Nx(j) = nnz(x2);
        
        if ~FitModel
            % Estimate Y2 using the sum of weighted atoms.
            Y2{j} = sum( Dn(:,DInds{j}) .* repmat( x', [size(Dn,1),1] ), 2)+Offset(j);
            % Calculate best gain from signal vectors
            try        G2 = pinv(Y2{j})*Y;    catch;        G2 = 1;    end
            B{j}(4) = G2*Gain;
            % Signal vector using new gain
            Y2{j} = Y2{j}*G2;
        else
            % Estimate Y2 using the model function and estimated parameters.
            Y2{j} = tmatrix*modelfun(B{j},stimulus,res,resmx,xx,yy,hrf);
            % Calculate best gain from signal vectors
            try        G2 = pinv(Y2{j})*Y;    catch;        G2 = 1;    end
            B{j}(4) = G2*Gain;
            % Signal vector using new gain
            Y2{j} = Y2{j}*G2;
        end
        
        
        SStot = sum((Y-mean(Y)).^2);
        SSres = sum((Y-Y2{j}).^2);
        
        R2(j)  = 1 - SSres/SStot;
        
        % If GammaStop option then compare with previous R2 value
        if GammaStop
            R2_new = R2(j);
        end
        
    else
        break;
    end
end

% figure;plot(Gammas,R2)

%% Pick best estimates
[~,I]=max(R2);
FitSummary.R2 = R2;
R2=R2(I);
if R2>-inf
    FitSummary.Nx = Nx;
    FitSummary.Gammas = Gammas;
    B=B{I};
    Y2=Y2{I};
    Gamma=Gammas(I);
    DInds=DInds{I};
    Nx=Nx(I);
else
    B=nan(1,5);
    Y2=nan;
    Gamma=nan;
    DInds=nan;
    Nx=nan;
    FitSummary.Nx = Nx;
    FitSummary.Gammas = Gammas;
end
