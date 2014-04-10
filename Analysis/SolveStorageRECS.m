function interp = SolveStorageRECS(model,interp,options)
% SOLVESTORAGE Solves for the REE using RECS

[interp,~,~,~,exitflag] = recsSolveREE(interp,model,[],[],options);

if exitflag~=1
  interp = rmfield(interp,{'cx','x'});
  s      = interp.s;
  x      = [zeros(size(s)) max(0,model.params(4)+model.params(5)*s)];

  [interp,~,~,~,exitflag] = recsSolveREE(interp,model,s,x,options);

  if exitflag~=1
    interp = recsFirstGuess(interp,model,[],[],[],options); % Check if this
                                                            % function works with explicit models
    interp = recsSolveREE(interp,model,s,[],options);
  end
end
