%% Initialization
ComPrices = readtable(fullfile('..','Data','ComPrices-DL1995.csv'),'ReadRowNames',true);

ComList = {'Coffee'; 'Copper'; 'Jute'; 'Maize'; 'Palmoil'; 'Sugar'; 'Tin'};

NActiveParams = 3;
Lik           = zeros(length(ComList),1);
exitflag      = zeros(length(ComList),1);
output        = cell(length(ComList),1);
V             = zeros(length(ComList),NActiveParams);
theta         = zeros(4,length(ComList));
thetainit     = zeros(4,length(ComList));
pstar         = zeros(length(ComList),1);
G             = zeros(length(ComList),1);

GridLimits = table(-5*ones(7,1),[30 40 30 40 30 20 45]','RowNames',ComList, ...
                   'VariableNames',{'Min' 'Max'});

options = struct('ActiveParams' , [1 1 0 1],...
                 'Display'      , 'off',...
                 'interp'       , 'linear',...
                 'MaxIter'      , 1E3  ,...
                 'TolX'         , 1E-10 ,...
                 'reesolveroptions',struct('atol',1E-10),...
                 'cov'          , 3,...
                 'ParamsTransformInvDer', @(P) [1; -exp(P(2)); exp(P(3)); exp(P(4))],...
                 'solveroptions',optimset('DiffMinChange', eps^(1/3),...
                                          'Display'      , 'off',...
                                          'FinDiffType'  , 'central',...
                                          'LargeScale'   , 'off',...
                                          'MaxFunEvals'  , 2000,...
                                          'MaxIter'      , 1000,...
                                          'TolFun'       , 1e-6,...
                                          'TolX'         , 1e-7,...
                                          'UseParallel'  , 'never'),...
                 'numjacoptions',struct([]),...
                 'numhessianoptions',struct('FinDiffRelStep'  , 1E-3,...
                                            'UseParallel'     , 'never'),...
                 'T'          , 5,...
                 'UseParallel', 'never');

gcp;
pctRunOnAll warning('off','backtrace');
pctRunOnAll warning('off','RECS:FailureREE');
pctRunOnAll warning('off','MATLAB:interp1:ppGriddedInterpolant');
N = [2 100 200];

%% Estimate in all situations
iter = 0;
for r=[0.02 0.05]
  options.ParamsTransform = @(P) [P(1); log(-P(2)); log(P(3)+r); log(P(4))];
  options.ParamsTransformInv = @(P) [P(1); -exp(P(2)); exp(P(3))-r; exp(P(4))];
  for solver={'fminsearch','fminunc'}
    iter = iter+1;
    options.solver = solver{:};
    parfor com=1:length(ComList)
      [Pobs,model,interp,tmp] = initpb(ComList{com},...
                                       [],...
                                       r,...
                                       GridLimits{com,:},...
                                       N, ...
                                       ComPrices,...
                                       options);
      thetainit(:,com) = tmp';
      try
        clear LogLik
        [thetatmp,Lik(com),vcov,g,~,exitflag(com),output{com}] = MaxLik(@(theta,obs) LogLik(theta,obs,model,interp,options),...
                                                          thetainit(:,com), ...
                                                          Pobs,options);
        V(com,:)          = sqrt(diag(vcov));
        theta(:,com)      = thetatmp;
        model.params(1:4) = thetatmp;
        interp            = SolveStorageEGM(model,interp,options);
        e                 = model.shocks.e;
        PriceInterp       = interp.cx{2};
        invdemand         = @(d) model.params(1)+model.params(2)*d;
        pstar(com)        = model.shocks.w'*max(PriceInterp(e),invdemand(e))*...
            (1-model.params(3))/(1+model.params(5))-model.params(4);
        G(com)            = norm(g,'inf');
      catch err
        exitflag(com) = 0;
        output{com}   = err;
        Lik(com)      = NaN;
        pstar(com)    = NaN;
        theta(:,com)  = NaN;
        V(com,:)      = NaN;
        G(com)        = NaN;
      end
    end
    Results(iter).r           = r; %#ok<*SAGROW>
    Results(iter).solver      = solver{:};
    Results(iter).exitflag    = exitflag;
    Results(iter).output      = output;
    Results(iter).Lik         = Lik;
    Results(iter).pstar       = pstar;
    Results(iter).theta       = theta';
    Results(iter).thetainit   = thetainit';
    Results(iter).V           = V;
    Results(iter).G           = G;
  end
end

%% Display estimation results
for i=1:length(Results)
  tmp = FormatResults(Results(i),options.ActiveParams);
  disp(tmp.Properties.Description)
  disp(tmp)
end
