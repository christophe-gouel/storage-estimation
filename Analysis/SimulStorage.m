function [Asim,Xsim,esim] = SimulStorage(model,interp,A0,nper)
% SIMULSTORAGE Simulate a standard storage model

params              = num2cell(model.params);
[a, b, delta, k, r] = params{:};
invdemand           = @(d) a+b*d;

nrep = length(A0);

Asim          = zeros(nrep,nper);
Asim(:,1)     = A0;

esim          = NaN(nrep,nper);
esim(:,2:end) = randn(nrep,nper-1);

P          = zeros(nrep,nper);
S          = zeros(nrep,nper);

[StockInterp,PriceInterp] = interp.cx{:};

for t=1:nper
  if t>1, Asim(:,t) = (1-delta)*S(:,t-1)+esim(:,t); end
  P(:,t) = max(PriceInterp(Asim(:,t)),invdemand(Asim(:,t)));
  S(:,t) = max(StockInterp(Asim(:,t)),0);
end

Xsim = permute(cat(3,S,P),[1 3 2]);