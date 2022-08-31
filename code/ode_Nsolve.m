function [time,state_values] = ode_Nsolve()


%%
load ../data/tempFile

%tc = 2*3600; % time (in hours) of control rod insertion
%CR = zeros(M).*(t<=tc) + PHID_CR_PHI.*(t>tc); %control rod step function

f = FunctionGen(M);
ti=0; tf = 30*3600;
tspan = [ti,tf];
IC = zeros(1,M*3)+0.01;
f_handle = eval(['@(t,s)[' f ']']);

[time, state_values] = ode15s(f_handle,tspan,IC,[]);

figure()

plot(time,state_values(:,1:3:end))
%ylim([-0.1,0.1])
ylabel("Amplitude (AU)")
xlabel("Time (h)")
legend("Eq mode", "rad mode 1", "rad mode 2", "ax mode")
end