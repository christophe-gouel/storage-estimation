function interp = SolveStorageEGM(model,interp,options)
% SOLVESTORAGEEGM Solve a standard storage model using endogenous grid method
  
%% Initialization
params              = num2cell(model.params);
[a, b, delta, k, r] = params{:};
beta                = (1-delta)/(1+r);

e = model.shocks.e;

demand    = @(p) (p-a)/b;
invdemand = @(d) a+b*d;

S         = interp.gridS;

Anext = gridmake(e,(1-delta)*S)*ones(2,1); % Next period availibility
W     = kron(speye(length(S)),model.shocks.w');

minA  = min(e);

Display = lower(options.Display);
MaxIter = options.MaxIter;
TolX    = options.TolX;

nA      = interp.N(1);

dis  = inf;
Iter = 0;
if isfield(interp,'s'), A = interp.s; end

% First guess
if ~isfield(interp,'cx')
  if ~exist('A','var'), A = linspace(minA,max(S),length(S)+nA); end
  PriceInterp = griddedInterpolant(A,invdemand(A),options.interp);
else
  PriceInterp = interp.cx{2};
end

if isfield(interp,'x')
  P = interp.x(:,2);
else
  P = max(PriceInterp(A),invdemand(A))';
end

%% Successive iterations
if isequal(Display,'iter'), 
  fprintf(1,'Successive approximation\n');
  fprintf(1,' Iteration  Residual\n');
end
while(dis > TolX && Iter < MaxIter)
  Iter = Iter+1;
  Pold = P;

  %% Generate endogenous grid
  EP = W*max(PriceInterp(Anext),invdemand(Anext));

  A2 = S+demand(beta*EP-k);
  A1 = linspace(minA,min(A2),nA+1)';
  A1 = A1(1:end-1);
  A  = [A1
        A2];
  
  %% New price function
  P  = invdemand([A1; A2-S]);

  dis         = norm(Pold-P);
  PriceInterp = griddedInterpolant(A,P,options.interp);

  if isequal(Display,'iter'), fprintf(1,'%9i   %7.2E\n',Iter,dis); end
end

%% Export results
StockInterp = griddedInterpolant(A,[zeros(nA,1); S],options.interp);

interp.s  = A;
interp.x  = [[zeros(size(A1)); S] P];
interp.cx = {StockInterp; PriceInterp};

if Iter==MaxIter
  warning('Failrure to solve for REE')
end