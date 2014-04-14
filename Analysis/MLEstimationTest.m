% MLEstimationTest

a     = 1;
b     = -2;
delta = 0;
k     = 0.02;
r     = 0.02;

model.params                    = [a b delta k r];
[model.shocks.e,model.shocks.w] = qnwnorm(20);

interp.s = linspace(min(model.shocks.e),20,1000)';

options = struct('ActiveParams' , [1 1 0 1],...
                 'Display'      , 'off',...
                 'InterpMethod' , 'linear',...
                 'MaxIter'      , 1E3  ,...
                 'TolX'         , 1E-10 ,...
                 'reesolveroptions',struct('atol',1E-10),...
                 'cov'          , 3,...
                 'solver','fminsearch',...
                 'ParamsTransformInvDer', @(P) [1; -exp(P(2)); exp(P(3)); exp(P(4))],...
                 'solveroptions',optimset('DiffMinChange', eps^(1/3),...
                                          'Display'      , 'iter',...
                                          'FinDiffType'  , 'central',...
                                          'LargeScale'   , 'off',...
                                          'MaxFunEvals'  , 2000,...
                                          'MaxIter'      , 1000,...
                                          'TolFun'       , 1e-6,...
                                          'TolX'         , 1e-7,...
                                          'UseParallel'  , 'always'),...
                 'numjacoptions',struct('FinDiffType','central'),...
                 'numhessianoptions',struct('FinDiffRelStep'  , 1E-3,...
                                            'UseParallel'     , 'always'),...
                 'T'          , 5,...
                 'UseParallel', 'never');

interp = SolveStorageDL(model,interp,options);
rng(0)
[Asim,Xsim] = SimulStorage(model,interp,0,1E5);

Pobs = squeeze(Xsim(1,2,:));

options.ParamsTransform = @(P) [P(1); log(-P(2)); log(P(3)+r); log(P(4))];
options.ParamsTransformInv = @(P) [P(1); -exp(P(2)); exp(P(3))-r; exp(P(4))];
clear LogLik
[theta,Lik,vcov,g,hess,exitflag,output] = MaxLik(@(theta,obs) LogLik(theta,obs,model,interp,options),...
                                                 [1.2 -4 0 0.05]', ...
                                                 Pobs,options);
options.solver = 'fminunc';
[theta,Lik,vcov,g,hess,exitflag,output] = MaxLik(@(theta,obs) LogLik(theta,obs,model,interp,options),...
                                                 theta, ...
                                                 Pobs,options);
