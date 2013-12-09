function [PL,theta,pstar,exitflag,output,V]= MaxPL(theta,Pobs,model,interp,options)
% MAXPL Maximizes the log pseudo-likelihood and calculates the standard deviations of the estimation for a storage model

THETA = [theta(1) log(-theta(2)) log(theta(3))];
options.optimsolveroptions.TypicalX = THETA;

[THETA,fval,exitflag,output] = feval(options.optimsolver,...
                                  @Objective,...
                                  THETA,...
                                  options.optimsolveroptions);

theta = [THETA(1) -exp(THETA(2)) exp(THETA(3))];
PL    = -fval;
pstar = (model.shocks.w'*funeval(interp.cx(:,2),interp.fspace,model.shocks.e))...
        *((1-model.params(2))/(1+model.params(3)))-theta(3);
      
if nargout==6
  G      = numjac(@LogLikComponent,THETA,options.numjacoptions);
  H      = numhessian(@Objective,THETA,options.numhessianoptions);
  Vtilde = H\(G'*G)/H;
  Dtilde = diag([1 -exp(THETA(2)) exp(THETA(3))]);
  V      = sqrt(diag(Dtilde*Vtilde*Dtilde'))';
end

function Obj = Objective(X)

  [model,interp] = SolveStorage(X,model,interp,options);
  PL  = PseudoLogLik(Pobs,model,interp);
  Obj = -PL;

end

function lnl = LogLikComponent(X)

 [model,interp] = SolveStorage(X,model,interp,options);
 [~,lnl]  = PseudoLogLik(Pobs,model,interp);

end

end