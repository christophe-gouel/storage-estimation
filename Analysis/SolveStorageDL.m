function [interp,exitflag] = SolveStorageDL(model,interp,options)
% SOLVESTORAGEDL Solves the storage model by fixed point as proposed by Deaton and Laroque

%% Initialization
params              = num2cell(model.params);
[a, b, delta, k, r] = params{:};
beta                = (1-delta)/(1+r);

e = model.shocks.e;
w = model.shocks.w;

demand       = @(p) (p-a)/b;
invdemand    = @(d) a+b*d;

A            = interp.s;

Display      = isequal(lower(options.Display),'iter');
InterpMethod = options.InterpMethod;
MaxIter      = options.MaxIter;
TolX         = options.TolX;

n            = length(A);
kshocks      = size(model.shocks.e,1);
ind          = (1:n);
ind          = ind(ones(1,kshocks),:);
inde         = repmat(1:kshocks,1,n);

% First guess
if isfield(interp,'cx')
  PriceInterp = interp.cx{2};
else
  PriceInterp = griddedInterpolant(A,max(invdemand(A),0),InterpMethod);
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
if Display
  fprintf(1,'Successive approximation\n');
  fprintf(1,' Iteration  Residual\n');
end

while(dis > TolX && Iter < MaxIter)
  Iter = Iter + 1;
  Pold = P;
  
  %% Calculate next-period availability  
  S     = max(A-demand(PriceInterp(A)),0);
  Anext = (1-delta)*S(ind,:)+e(inde,:);

  %% Update the price and its approximation
  Pnext = max(PriceInterp(Anext),invdemand(Anext));
  P     = max(invdemand(A),beta*reshape(w'*reshape(Pnext,kshocks,n),n,1)-k); 
  PriceInterp.Values = P;
   
  dis = norm(P-Pold);
  if Display, fprintf(1,'%9i   %7.2E\n',Iter,dis); end
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
