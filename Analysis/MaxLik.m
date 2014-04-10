function [params,ML,vcov,g,H,exitflag,output] = MaxLik(loglikfun,params,obs,options,varargin)
% MAXLIK Maximizes a log-likelihood function
%
% PARAMS = MAXLIK(LOGLIKFUN,PARAMS,OBS)
%
% PARAMS = MAXLIK(LOGLIKFUN,PARAMS,OBS,OPTIONS)
%
% PARAMS = MAXLIK(LOGLIKFUN,PARAMS,OBS,OPTIONS,VARAGIN)
%
% [PARAMS,ML] = MAXLIK(LOGLIKFUN,PARAMS,OBS,...)
%
% [PARAMS,ML,VCOV] = MAXLIK(LOGLIKFUN,PARAMS,OBS,...)
%
% [PARAMS,ML,VCOV,G] = MAXLIK(LOGLIKFUN,PARAMS,OBS,...)
%
% [PARAMS,ML,VCOV,G,H] = MAXLIK(LOGLIKFUN,PARAMS,OBS,...)
%
% [PARAMS,ML,VCOV,G,H,EXITFLAG] = MAXLIK(LOGLIKFUN,PARAMS,OBS,...)
%
% [PARAMS,ML,VCOV,G,H,EXITFLAG,OUTPUT] = MAXLIK(LOGLIKFUN,PARAMS,OBS,...)

% Copyright (C) 2014 Christophe Gouel
% Licensed under the Expat license

%% Initialization
defaultopt = struct('ActiveParams'         , [],...
                    'cov'                  , 3,...
                    'numhessianoptions'    , struct(),...
                    'numjacoptions'        , struct(),...
                    'ParamsTransform'      , @(P) P,...
                    'ParamsTransformInv'   , @(P) P,...
                    'ParamsTransformInvDer', @(P) ones(size(P)),...
                    'solver'               , 'fminunc',...
                    'solveroptions'        , struct());
if nargin < 4 || isempty(options)
  options = defaultopt;
else
  warning('off','catstruct:DuplicatesFound')
  if isfield(options,'numhessianoptions')
    options.numhessianoptions = catstruct(defaultopt.numhessianoptions,options.numhessianoptions);
  end
  if isfield(options,'numjacoptions')
    options.numjacoptions = catstruct(defaultopt.numjacoptions,options.numjacoptions);
  end
  if isfield(options,'solveroptions')
    options.solveroptions = catstruct(defaultopt.solveroptions,options.solveroptions);
  end
  options = catstruct(defaultopt,options);
end
ActiveParams          = options.ActiveParams;
cov                   = options.cov;
ParamsTransform       = options.ParamsTransform;
ParamsTransformInv    = options.ParamsTransformInv;
ParamsTransformInvDer = options.ParamsTransformInvDer;
solver                = options.solver;

validateattributes(loglikfun,{'char','function_handle'},{},1)
validateattributes(params,{'numeric'},{'column','nonempty'},2)

if norm(ParamsTransformInv(ParamsTransform(params))-params)>=sqrt(eps)
  error('Functions to transform parameters are not inverse of each other')
end

if norm(diag(numjac(@(P) ParamsTransformInv(P),ParamsTransform(params)))-ParamsTransformInvDer(ParamsTransform(params)))>=sqrt(eps)
  error(['The function to differentiate transformed parameters does not correspond ' ...
         'to its finite difference gradient.'])
end  
  
if isa(loglikfun,'char'), loglikfun = str2func(loglikfun); end

if isempty(ActiveParams)
  ActiveParams = true(size(params));
else
  validateattributes(ActiveParams,{'logical','numeric'},{'vector','numel',numel(params)})
  ActiveParams = ActiveParams(:)~=zeros(size(params));
end

%% Functions and matrices to extract active parameters for the estimation
SelectParamsMat = zeros(length(params),sum(ActiveParams));
ind             = 1:length(params);
SelectParamsMat(sub2ind(size(SelectParamsMat),ind(ActiveParams),1:sum(ActiveParams))) = 1;
SelectParams    = @(P) ParamsTransform(params).*(~ActiveParams)+SelectParamsMat*P;

%% Maximization of the log-likelihood
Objective = @(P) -sum(loglikfun(ParamsTransformInv(SelectParams(P)),obs,varargin{:}));
problem = struct('objective', Objective,...
                 'x0'       , SelectParamsMat'*ParamsTransform(params),...
                 'solver'   , solver,...
                 'options'  , options.solveroptions);
[PARAMS,ML,exitflag,output] = feval(solver,problem);
params                      = ParamsTransformInv(SelectParams(PARAMS));
ML                          = -ML;

%% Covariance and hessian of parameters
if nargout>=4 || (nargout>=3 && any(cov==[2 3]))
  G   = numjac(@(P) loglikfun(ParamsTransformInv(SelectParams(P)),obs,varargin{:}),...
               PARAMS,options.numjacoptions);
  g   = -sum(G,1);
end

if nargout>=5 || (nargout>=3 && any(cov==[1 3]))
  H = numhessian(Objective,PARAMS,options.numhessianoptions);
  if ~all(eig(H)>=0)
    warning('Hessian is not positive definite')
  end
end

if nargout>=3
  D = diag(ParamsTransformInvDer(SelectParams(PARAMS)));
  D = D(ActiveParams,ActiveParams);
  switch cov
    case 1
      vcov = D'*(H\D);
    case 2
      vcov = D'*((G'*G)\D);
    case 3
      vcov = D'*(H\(G'*G)/H)*D;
    otherwise
      vcov = [];
  end
end
