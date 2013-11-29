function [Pobs,model,interp,theta] = initpb(commodity,fgparams,r,rangeA,N,ComPrices,options)
% INITPB Initializes the problem

%% Observed prices
Pobs = table2array(ComPrices(:,commodity));

%% Definition of model object
model = recsmodel('storageexplicit.yaml');
model.shocks.e = [-1.755 -1.045 -0.677 -0.386 -0.126 0.126 0.386 0.677 1.045 1.755]';
model.shocks.w = ones(size(model.shocks.e))/length(model.shocks.e);

if isempty(fgparams), fgparams = [0.3 0.02]; end
EP     = mean(Pobs);
sigmaP = std(Pobs);
theta  = [EP -sigmaP/fgparams(1) fgparams(2)*EP];
model.params = [theta(3) 0 r theta(1) theta(2)];

[model.sss,model.xss] = recsSS(model,0,[0 theta(1)],options);

%% Definition of the interpolation structure
[interp,s] = recsinterpinit(N,rangeA(1),rangeA(2));

%% Solve for the REE
x = [zeros(size(s)) max(0,theta(1)+theta(2)*s)]; % First guess
interp =  recsSolveREE(interp,model,s,x,options);
