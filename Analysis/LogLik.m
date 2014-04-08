function l = LogLik(params,Pobs,model,interp,options)
% LOGLIK Calculates the log-likelihood for given price observations and given parameters

%% Initialization
persistent cx x
if ~isempty(cx), interp.cx = cx; end
if ~isempty(x),  interp.x  = x;  end

model.params([1 4:5]) = params([3 1:2]);
par                   = num2cell(model.params);
[k, delta, r, a, b]   = par{:}; %#ok<ASGLU>
interp                = SolveStorage(model,interp,options);
cx                    = interp.cx;
s                     = interp.s;
x                     = interp.x;
T                     = length(Pobs);
e                     = model.shocks.e;
w                     = model.shocks.w;

%% Find availabilities corresponding to observed prices
invDemandFunction = interp1(x(:,2),s,'spline','pp');
Aobs              = ppval(invDemandFunction,Pobs);

%% Residuals
Sobs       = max(funeval(cx(:,1),interp.fspace,Aobs),0);
omega      = NaN(T,1);
omega(2:T) = Aobs(2:T) - (1-delta)*Sobs(1:T-1);

%% Jacobian
pstar          = w'*funeval(cx(:,2),interp.fspace,e)*(1-delta)/(1+r)-k;
%{
There are three methods to calculate the Jacobian:
 1. Calculate the derivative of the inverse price function by finite 
    differences.
 2. Use the fact that (f^{-1})'(f(x))=1/f'(x) and calculate the inverse of the
    derivative of the price function.
 3. Calculate the derivative of the inverse price function analytically.
In all cases, we can correct for the fact that above pstar the derivative of
the inverse price function is equal to 1/b.

The methods 1 and 3 are equivalent, except that 3 is faster and more precise.
The methods 2 and 3 are equivalent for high levels of precision of the spline
approximation. We opt for method 3, which should be the fastest method.
%}
% 1. Finite-difference of the inverse price function
% J    = diag(numjac(@(P) ppval(invDemandFunction,P),Pobs));
% J(Pobs>=pstar) = 1/b;

% 2. Inverse of the derivative of the price function
% J    = 1./funeval(cx(:,2),interp.fspace,Aobs,1);
% J(Pobs>=pstar) = 1/b;

% 3. Differentiate the inverse price function
[breaks,coefs,l,order,d] = unmkpp(invDemandFunction);
dinvDemandFunction       = mkpp(breaks,repmat(order-1:-1:1,d*l,1).*coefs(:,1:order-1),d);
J                        = ones(size(Pobs))/b;
J(Pobs<pstar)            = ppval(dinvDemandFunction,Pobs(Pobs<pstar));

%% Log-pseudo-likelihood
l = -0.5*(log(2*pi)+ omega(2:T).^2)+ log(abs(J(2:T)));
