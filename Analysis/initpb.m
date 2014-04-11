function [Pobs,model,interp,theta] = initpb(commodity,fgparams,r,rangeA,N,ComPrices,options)
% INITPB Initializes the problem

%% Observed prices
Pobs = table2array(ComPrices(:,commodity));

%% Definition of model structure
%model.shocks.e = [-1.755 -1.045 -0.677 -0.386 -0.126 0.126 0.386 0.677 1.045 1.755]';
%model.shocks.w = ones(size(model.shocks.e))/length(model.shocks.e);
[model.shocks.e,model.shocks.w] = qnwnorm(20);

if isempty(fgparams), fgparams = [0.3 0.02]; end
EP     = mean(Pobs);
sigmaP = std(Pobs);
theta  = [EP -sigmaP/fgparams(1) 0 fgparams(2)*EP];
model.params = [theta r];

%% Grid on stocks
interp.N     = N;
interp.gridS = linspace(0,5,N(3))';       % Coarse grid far from the kink
interp.gridS = [interp.gridS; linspace(0,1,N(2))'];  % Precise grid close to the kink
interp.gridS = sort(unique(interp.gridS));

%% Solve the storage model
interp = SolveStorageEGM(model,interp,options);
