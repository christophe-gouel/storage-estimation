function [PL,l]= PseudoLogLik(Pobs,model,interp)
% PSEUDOLOGLIK Calculates the pseudo-likelihood for given price observations and given model solution

%% Initialization
T      = length(Pobs);
k      = size(model.shocks.e,1);
ind    = (1:T);
ind    = ind(ones(1,k),:);
s      = interp.s;
x      = interp.x;

%% Find availabilities corresponding to observed prices
Aobs    = interp1(x(:,2),s,Pobs,'spline');

%% Expectations calculation
xobs    = funeval(interp.cx,interp.fspace,Aobs);
AAobs   = Aobs(ind,:);
xxobs   = xobs(ind,:);
Anext   = model.functions.g(AAobs,xxobs,model.shocks.e(repmat(1:k,1,T),:),...
                            model.params);
Pnext   = funeval(interp.cx(:,2),interp.fspace,Anext);

EP  = reshape(model.shocks.w'*reshape(Pnext,k,T),T,1);           % E(P|Aobs)
ESP = reshape(model.shocks.w'*reshape(Pnext.^2,k,T),T,1)-EP.^2;  % Var(P|Aobs)

%% Log-pseudo-likelihood
l = -0.5*(log(2*pi)+log(ESP(1:T-1)) + ((Pobs(2:T) - EP(1:T-1)).^2)./ESP(1:T-1));
PL = sum(l);