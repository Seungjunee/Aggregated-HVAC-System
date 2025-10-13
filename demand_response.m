function [Psched,soc_sched] = demand_response(a_agg,b_agg,soc0,Paggdown,Paggup,tmax,RTP,tau)

Pagg = optimvar("Pagg",1,tmax);
SOC = optimvar("SOC",1,tmax+1);
prob = optimproblem;
prob.Objective = sum(Pagg*RTP) + SOC.^2*tau'; % 
prob.Constraints.cons2 = -Paggdown <= Pagg;
prob.Constraints.cons3 = Pagg <= Paggup;
prob.Constraints.cons4 = -5 <= SOC;
prob.Constraints.cons5 = SOC <= 5;
prob.Constraints.cons6 = SOC(1) == soc0;
% prob.Constraints.cons7 = SOC(end) == 0;
cons8 = optimconstr(tmax);
for t=1:tmax
    cons8(t) = SOC(t+1) == a_agg*SOC(t)+b_agg*Pagg(t); % 
end
prob.Constraints.cons8 = cons8;
tic
[sol,fval,exitflag,output] = solve(prob);
toc
Psched = sol.Pagg;
soc_sched = sol.SOC;
end