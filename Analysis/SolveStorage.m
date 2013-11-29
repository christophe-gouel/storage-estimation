function [model,interp] = SolveStorage(THETA,model,interp,options)
% SOLVESTORAGE Solves for the REE given the parameters THETA

model.params(4) = THETA(1);
model.params(5) = -exp(THETA(2));
model.params(1) = exp(THETA(3));

[interp,~,~,~,exitflag] = recsSolveREE(interp,model,[],[],options);

if exitflag~=1
  interp = rmfield(interp,{'cx','x'});
  s      = interp.s;
  theta  = [THETA(1) -exp(THETA(2)) exp(THETA(3))];
  x      = [zeros(size(s)) max(0,theta(1)+theta(2)*s)];

  [interp,~,~,~,exitflag] = recsSolveREE(interp,model,s,x,options);

  if exitflag~=1
    interp = recsFirstGuess(interp,model,[],[],[],options); % Check if this
                                                            % function works with explicit models
    interp = recsSolveREE(interp,model,s,[],options);
  end
end
