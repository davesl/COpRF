function [B1, B2, B3, B4, B5, R_2, Gamma, Nx, FitSummary] = PRF_GridSearch(Y,Dsubn,Dsub_norm,Dsub_params,Gammas,tmatrix,res,resmx,FitModel,stimulus,xx,yy,hrf,modelfun,GammaStop)

% First use sub-dictionary to estimate x- and y-position
[BB,~,R_2,Gamma,~,Nx,FitSummary] = ConvexEstimate(Y,Dsubn,Dsub_norm,Dsub_params,Gammas,FitModel,tmatrix,stimulus,xx,yy,res,resmx,hrf,modelfun,GammaStop);

if max(isnan(BB))
    % Fitting failure (e.g. noise data) so set as nan
    R_2 = nan;
    Gamma = nan;
    Nx = nan;
    FitSummary.Nx = nan;
    FitSummary.Gammas = nan;
    FitSummary.R2 = nan;
end


B1 = BB(1);
B2 = BB(2);
B3 = BB(3);
B4 = BB(4);
B5 = BB(5);
