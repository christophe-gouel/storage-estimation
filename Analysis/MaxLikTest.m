function tests = MaxLikTest
tests = functiontests(localfunctions);
end

function testCompareWithOLS(testCase)

warning('off','all');
  
%% Generate data
par = [1 0.1 0.5 0.2]';
N = 500;
X = [ones(N,1) randn(N,2)];
Y = X*par(1:end-1)+par(end)*randn(N,1);

%% OLS
MLfit = fitlm(X(:,2:end),Y);

%% MLE
data = [X Y];
LogLik = @(params,obs) -0.5*(log(params(end)^2)+log(2*pi)+(obs(:,end)-obs(:,1:end-1)*params(1:end-1)).^2/params(end)^2);
par0 = par.*(1+0.2*randn(size(par)));
options = struct('cov',1,...
                 'solveroptions',struct('Display','off'));
[paropt,ML,vcov] = MaxLik(LogLik,par0,data,options);

%% Tests
verifyLessThanOrEqual(testCase,abs((MLfit.LogLikelihood-ML)/MLfit.LogLikelihood),5E-3)
verifyLessThanOrEqual(testCase,abs(MLfit.Coefficients.Estimate-paropt(1:end-1)),1E-5)
verifyLessThanOrEqual(testCase,abs((MLfit.Coefficients.SE-sqrt(diag(vcov(1:end-1,1:end-1))))./MLfit.Coefficients.SE),0.15)

options.cov = 2;
[~,~,vcov] = MaxLik(LogLik,paropt,data,options);
verifyLessThanOrEqual(testCase,abs((MLfit.Coefficients.SE-sqrt(diag(vcov(1:end-1,1:end-1))))./MLfit.Coefficients.SE),0.15)

options.cov = 3;
[~,~,vcov] = MaxLik(LogLik,paropt,data,options);
verifyLessThanOrEqual(testCase,abs((MLfit.Coefficients.SE-sqrt(diag(vcov(1:end-1,1:end-1))))./MLfit.Coefficients.SE),0.15)
end
