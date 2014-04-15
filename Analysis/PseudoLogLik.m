function l = PseudoLogLik(params,Pobs,model,interp,options)
% PSEUDOLOGLIK Calculates the pseudo-likelihood for given price observations and given model solution

%% Initialization
persistent cx x
if ~isempty(cx), interp.cx = cx; end
if ~isempty(x),  interp.x  = x;  end

model.params(1:4)         = params;
par                       = num2cell(model.params);
[a, b, delta, ~, ~]       = par{:};
interp                    = SolveStorageDL(model,interp,options);
cx                        = interp.cx;
s                         = interp.s;
x                         = interp.x;
T                         = length(Pobs);
k                         = size(model.shocks.e,1);
ind                       = (1:T);
ind                       = ind(ones(1,k),:);
[~,PriceInterp]           = interp.cx{:};
demand                    = @(p) (p-a)/b;
invdemand                 = @(d) a+b*d;

%% Find availabilities corresponding to observed prices
Aobs    = max(interp1(x(:,2),s,Pobs,options.InterpMethod),demand(Pobs));

%% Expectations calculation
Sobs    = max(Aobs-demand(Pobs),0);
Sobs    = Sobs(ind,:);
Anext   = (1-delta)*Sobs+model.shocks.e(repmat(1:k,1,T),:);
Pnext   = max(PriceInterp(Anext),invdemand(Anext));

EP  = reshape(model.shocks.w'*reshape(Pnext,k,T),T,1);           % E(P|Aobs)
ESP = reshape(model.shocks.w'*reshape(Pnext.^2,k,T),T,1)-EP.^2;  % Var(P|Aobs)

%% Log-pseudo-likelihood
l = -0.5*(log(2*pi)+log(ESP(1:T-1)) + ((Pobs(2:T) - EP(1:T-1)).^2)./ESP(1:T-1));
