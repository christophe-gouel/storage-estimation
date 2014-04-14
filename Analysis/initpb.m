function [Pobs,model,interp,theta] = initpb(commodity,fgparams,r,rangeA,N,ComPrices,options)
% INITPB Initializes the problem

%% Observed prices
Pobs = table2array(ComPrices(:,commodity));

%% Definition of model structure
[model.shocks.e,model.shocks.w] = qnwnorm(20);

if isempty(fgparams), fgparams = [0.3 0.02]; end
EP     = mean(Pobs);
sigmaP = std(Pobs);
theta  = [EP -sigmaP/fgparams(1) 0 fgparams(2)*EP];
model.params = [theta r];

%% Grid on availability
interp.s = linspace(rangeA(1),rangeA(2),N)';       

%% Solve the storage model
interp = SolveStorageDL(model,interp,options);
