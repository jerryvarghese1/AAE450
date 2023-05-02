function cost = monte_carlo_fcn(file, NumBins, NumSims,table,iter,years,inf_facs)

% Read in costs and sort
[NUM,TXT,~] = xlsread(file);
item = TXT(2:end,1);
low = NUM(:,1);
avg = NUM(:,2);
high = NUM(:,3);

% Run triangular distribution fcn for each cost item
cost = zeros(NumSims,1);
for j = 1:length(item)
    cost = tri_dist_fcn(low(j),avg(j),high(j),NumSims,years,inf_facs) + cost;
end

% Plot Monte Carlo results
figure(iter)
histogram(cost/1000, NumBins)
title("Monte Carlo Cost Estimation, " + table,"FontSize",14)
xlabel("Total Cost [$M]","FontSize",12,"FontWeight","bold")
ylabel("Frequency","FontSize",12,"FontWeight","bold")
grid on

% Print data statistics
fprintf("\n%s\n",table)
fprintf("Low Cost = %.2f [$M]\n",min(cost)/1000)
fprintf("Mean Cost = %.2f [$M]\n",mean(cost)/1000)
fprintf("High Cost = %.2f [$M]\n",max(cost)/1000)