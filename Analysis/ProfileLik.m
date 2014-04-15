function Results = ProfileLik(loglikfun,params,obs,options,varargin)
% PROFILELIK Profiles the log-likelihood 

ActiveParams0 = options.ActiveParams;
if isempty(ActiveParams0)
  ActiveParams0 = true(size(params));
else
  validateattributes(ActiveParams0,{'logical','numeric'},{'vector','numel',numel(params)})
  ActiveParams0 = ActiveParams0(:)~=zeros(size(params));
end

vec     = @(X) X(:);
N       = 10;
params0 = params;
ind     = 1:length(params);
Results = NaN(sum(ActiveParams0),N,2);

iter    = 0;
for i=ind(ActiveParams0)
  iter = iter+1;
  params               = params0;
  options.ActiveParams = ActiveParams0;
  
  grid = sort(linspace(params(i)*0.8,params(i)*1.2,N));
  options.ActiveParams(i) = 0;
  Results(iter,:,1) = grid;
  for j=1:N
    params(i) = grid(j);
    [~,Results(iter,j,2)] = MaxLik(loglikfun,params,obs,options,varargin{:});
  end
end

figure
for i=1:sum(ActiveParams0)
  subplot(ceil((sum(ActiveParams0))/ceil(sqrt(sum(ActiveParams0)))),ceil(sqrt(sum(ActiveParams0))),i)
  plot(vec(Results(i,:,1)),vec(Results(i,:,2)))
end
