function X = FormatResults(Results)
% FORMATRESULTS Organize results in table

exitflag    = Results.exitflag;
solver      = Results.solver;
PL          = Results.PL;
r           = Results.r;
theta       = Results.theta;
V           = Results.V;
pstar       = Results.pstar;
G           = Results.G;

X = cell(size(theta,1)*2,7);

i = 1:2;
for com=1:size(theta,1)
  if exitflag(com)>0
    for param=1:3
      X(i,param) = {num2str(theta(com,param),'%8.4f');
                    ['(' num2str(V(com,param),'%8.4f') ')']};
    end
    X{i(1),4} = PL(com);
    X{i(1),5} = pstar(com);
  end
  X{i(1),6} = exitflag(com);
  X{i(1),7} = G(com);
  i = i+2;
end

X = table(X(:,1),X(:,2),X(:,3),X(:,4),X(:,5),X(:,6),X(:,7),...
          'VariableNames',{'a' 'b' 'k' 'PL' 'pstar' 'exitflag' 'gradient'});
X.Properties.Description = ...
    ['Optimized using ' solver ' with r=' ,num2str(r,'%4.2f')];

