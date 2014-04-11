function [interp,exitflag] = SolveStorageDL(model,interp,options)
% SOLVESTORAGEDL Solves the storage model by fixed point as proposed by Deaton and Laroque

%% Initialization
params              = num2cell(model.params);
[a, b, delta, k, r] = params{:};
beta                = (1-delta)/(1+r);

e = model.shocks.e;

demand      = @(p) (p-a)/b;
invdemand   = @(d) a+b*d;

A           = interp.s;
W           = kron(speye(length(A)),model.shocks.w');

Display      = lower(options.Display);
InterpMethod = options.InterpMethod;
MaxIter      = options.MaxIter;
TolX         = options.TolX;

% First guess
if ~isfield(interp,'cx')
  PriceInterp = griddedInterpolant(A,max(invdemand(A),0),InterpMethod);
else
  PriceInterp = interp.cx{2};
end

if isfield(interp,'x')
  P = interp.x(:,2);
else
  P = max(PriceInterp(A),invdemand(A));
end

dis      = inf;
exitflag = 1;
Iter     = 0;

%% Successive iterations
if isequal(Display,'iter'), 
  fprintf(1,'Successive approximation\n');
  fprintf(1,' Iteration  Residual\n');
end

while(dis > TolX && Iter < MaxIter)
  Iter        = Iter + 1;
  Pold = P;
  
  %% Calculate next-period availability  
  Anext = gridmake(e,(1-delta)*max(A-demand(PriceInterp(A)),0))*ones(2,1);
  
  %% Update the price and its approximation
  P = max(invdemand(A),beta*(W*max(PriceInterp(Anext),invdemand(Anext)))-k); 
  PriceInterp = griddedInterpolant(A,P,InterpMethod);
   
  dis = norm(P-Pold);
  if isequal(Display,'iter'), fprintf(1,'%9i   %7.2E\n',Iter,dis); end
end

S                  = max(A-demand(PriceInterp(A)),0);
StockInterp        = griddedInterpolant(A,S,InterpMethod);
interp.PriceInterp = PriceInterp;
interp.x           = [S P];
interp.cx          = {StockInterp; PriceInterp};

if Iter==MaxIter
  warning('Failure to solve for REE, residuals= %4.2e',dis)
  exitflag = 0;
end
