function X = FormatResults(Results,ActiveParams)
% FORMATRESULTS Organize results in table

if nargin<2, ActiveParams = true(4,1); end

exitflag    = Results.exitflag;
solver      = Results.solver;
Lik         = Results.Lik;
r           = Results.r;
theta       = Results.theta;
V           = Results.V;
pstar       = Results.pstar;
G           = Results.G;

nparams    = sum(ActiveParams~=0);
iparams    = 1:4;
paramsname = {'a' 'b' 'delta' 'k'};
paramsname = paramsname(iparams(ActiveParams~=0));
X = cell(size(theta,1)*2,nparams+4);

i = 1:2;
for com=1:size(theta,1)
  if exitflag(com)>0
    iter = 0;
    for param=iparams(ActiveParams~=0)
      iter = iter+1;
      X(i,iter) = {num2str(theta(com,param),'%8.4f');
                    ['(' num2str(V(com,iter),'%8.4f') ')']};
    end
    X{i(1),nparams+1} = Lik(com);
    X{i(1),nparams+2} = pstar(com);
  end
  X{i(1),nparams+3} = exitflag(com);
  X{i(1),nparams+4} = G(com);
  i = i+2;
end

X = array2table(X,...
                'VariableNames',{paramsname{:} 'Lik' 'pstar' 'exitflag' 'gradient'});
X.Properties.Description = ...
    ['Optimized using ' solver ' with r=' ,num2str(r,'%4.2f')];

