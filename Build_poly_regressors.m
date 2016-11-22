function tmatrix = Build_poly_regressors(data,tr,options)
%% Polynomial code adapted from: https://github.com/kendrickkay/analyzePRF
% Copyright (c) 2013, Kendrick Kay
% All rights reserved.
% 
% Modified by D Slater, Aug. 2016

% calc
is3d = size(data{1},4) > 1;
if is3d
  dimdata = 3;
  dimtime = 4;
  xyzsize = sizefull(data{1},3);
else
  dimdata = 1;
  dimtime = 2;
  xyzsize = size(data{1},1);
end
numtime = cellfun(@(x) size(x,2),data);
numruns = length(data);

% what polynomials should we use?
if isempty(options.maxpolydeg)
  options.maxpolydeg = cellfun(@(x) round(size(x,dimtime)*tr/60/2),data);
end
if isscalar(options.maxpolydeg)
  options.maxpolydeg = repmat(options.maxpolydeg,[1 numruns]);
end
fprintf('using the following maximum polynomial degrees: %s\n',mat2str(options.maxpolydeg));


% construct polynomial regressors matrix
polyregressors = {};
for p=1:length(data)
    if isnan(options.maxpolydeg(p))
        polyregressors{p} = zeros(numtime(p),0);
    else
        polyregressors{p} = constructpolynomialmatrix(numtime(p),0:options.maxpolydeg(p));
    end
end

% Multiplication matrix for regression
tmatrix = projectionmatrix(blkdiag(polyregressors{:}));


