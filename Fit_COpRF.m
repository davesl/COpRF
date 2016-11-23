function [B1, B2, B3, B4, B5, R_2, Gamma, Nx, FitSummary] = Fit_COpRF(Y,Dn,D_norm,Dparams,Dsubn,Dsub_norm,Dsub_params,Gammas,tmatrix,Xpos,Ypos,res,resmx,R,FitModel,stimulus,xx,yy,hrf,modelfun,GammaStop)
% [B1, B2, B3, B4, B5, R_2, Gamma, Nx, FitSummary] = Fit_COpRF(Y,Dn,D_norm,Dparams,Dsubn,Dsub_norm,Dsub_params,Gammas,tmatrix,Xpos,Ypos,res,resmx,R,FitModel,stimulus,xx,yy,hrf,modelfun,GammaStop)
% 
% A wrapper for the main CO-pRF fitting proceedures.
% 1) The reduced size sub-dictionary estimates model parameters.
% 2) Estimated xpos and ypos values are used to select only atoms contained 
%    in the full ditionary within a distance R from estimated position.
% 3) The selected full size dictionary atoms are used to fit model
%    parameters.

% First use sub-dictionary to estimate x- and y-position
Bsub = ConvexEstimate(Y,Dsubn,Dsub_norm,Dsub_params,Gammas,FitModel,tmatrix,stimulus,xx,yy,res,resmx,hrf,modelfun,GammaStop);

if ~max(isnan(Bsub))
    % Find entries within radius of x/y position estimate
    I = ((Xpos-Bsub(1)).^2 + (Ypos-Bsub(2)).^2).^0.5<=resmx*R; I(1:2)=1;
    
    % Now estimate on dense dictionary centred on x/y position estimate
    [BB,~,R_2,Gamma,~,Nx,FitSummary]=ConvexEstimate(Y,Dn(:,I),D_norm(I),Dparams(I,:),Gammas,FitModel,tmatrix,stimulus,xx,yy,res,resmx,hrf,modelfun,GammaStop);
    
else
    % Fitting failure (e.g. noise data) so set as nan
    BB = Bsub;
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

