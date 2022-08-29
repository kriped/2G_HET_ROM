function [DS] = ode(t,s)


%%
load tempFile

%tc = 2*3600; % time (in hours) of control rod insertion
%CR = zeros(M).*(t<=tc) + PHID_CR_PHI.*(t>tc); %control rod step function

f = FunctionGen(M);

[DS] = eval(['[' f ']']);

end