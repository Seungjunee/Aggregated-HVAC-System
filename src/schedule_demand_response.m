function [Psched,soc_sched] = schedule_demand_response(a_agg,b_agg,soc0,Paggdown,Paggup,tmax,RTP,tau,socLimit)

if nargin < 9 || isempty(socLimit)
    socLimit = 1;
end
assert(socLimit > 0, 'The symmetric SOC limit must be positive.');

Pagg = optimvar("Pagg",1,tmax);
SOC = optimvar("SOC",1,tmax+1);
prob = optimproblem;
prob.Objective = sum(Pagg*RTP) + SOC.^2*tau'; % 
prob.Constraints.cons2 = -Paggdown <= Pagg;
prob.Constraints.cons3 = Pagg <= Paggup;
prob.Constraints.cons4 = -socLimit <= SOC;
prob.Constraints.cons5 = SOC <= socLimit;
prob.Constraints.cons6 = SOC(1) == soc0;
% prob.Constraints.cons7 = SOC(end) == 0;
cons8 = optimconstr(tmax);
for t=1:tmax
    cons8(t) = SOC(t+1) == a_agg*SOC(t)+b_agg*Pagg(t); % 
end
prob.Constraints.cons8 = cons8;
[sol,~,exitflag] = solve(prob);
assert(exitflag > 0, 'Demand-response scheduling did not converge.');
Psched = sol.Pagg;
soc_sched = sol.SOC;
end
