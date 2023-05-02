%% Cost Simulation %%
clear; clc; close all;

% Histogram / Simulation Properties
NumBins = 50;
NumSims = 10000;
inf_facs = readmatrix("inflation_factors.txt");

% Run Monte Carlo Simulations
cost_8 = monte_carlo_fcn("table_11_8.xlsx",NumBins,NumSims,"Table 11-8",1,8,inf_facs);
cost_9 = monte_carlo_fcn("table_11_9.xlsx",NumBins,NumSims,"Table 11-9",2,8,inf_facs);
cost_14 = monte_carlo_fcn("table_11_14.xlsx",NumBins,NumSims,"Table 11-14",3,8,inf_facs);
cost_28 = tri_dist_fcn(138341*0.8,138341,138341*1.2,NumSims,17,inf_facs(12:end));

% Table 11-28 Plot
figure(4)
histogram(cost_28/1000, NumBins)
title("Triangular Distribution Cost Estimation, Table 11-28","FontSize",14)
xlabel("Total Cost [$M]","FontSize",12,"FontWeight","bold")
ylabel("Frequency","FontSize",12,"FontWeight","bold")
grid on

% Combine Table Costs
cost = cost_8 + cost_9 + cost_14 + cost_28;

% 70th Percentile Estimate
cost_sorted = sort(cost);
fifty_perc = cost_sorted(round(length(cost_sorted)*0.5));
seventy_perc = cost_sorted(round(length(cost_sorted)*0.7));
cost_margin = seventy_perc*1.3;

leg1 = sprintf("70 %% Conf = %.2f [$M]", seventy_perc/1000);
leg2 = sprintf("50 %% Conf (avg) = %.2f [$M]", fifty_perc*1.3/1000);
leg3 = sprintf("70 %% Conf = %.2f [$M]", cost_margin/1000);

% Plot Combined Histogram
figure(5)
hold on;
xline(seventy_perc/1000,'r','LineWidth',2.5)
histogram(cost/1000, NumBins)
title("Monte Carlo Total Cost Estimation","FontSize",14)
xlabel("Total Cost [$M]","FontSize",12,"FontWeight","bold")
ylabel("Frequency","FontSize",12,"FontWeight","bold")
legend(leg1,"FontSize",10)
grid on

% Plot Combined Histogram with 30% Margin
figure(6)
hold on;
xline(cost_margin/1000,'r','LineWidth',2.5)
histogram(cost*1.3/1000, NumBins)
title("Monte Carlo Total Cost Estimation (30% Cost Margin)","FontSize",14)
xlabel("Total Cost [$M]","FontSize",12,"FontWeight","bold")
ylabel("Frequency","FontSize",12,"FontWeight","bold")
legend(leg3,"FontSize",10)
grid on

% Distribution Plot with 30% Margin
y = linspace(0,1,NumSims);
figure(7)
hold on;
xline(fifty_perc*1.3/1000,'k--','LineWidth',2.5)
xline(cost_margin/1000,'r--','LineWidth',2.5)
plot(cost_sorted*1.3/1000, y, 'b*','MarkerSize',3)
title("Monte Carlo Total Cost Distribution (30% Cost Margin)","FontSize",14)
xlabel("Total Cost [$M]","FontSize",12,"FontWeight","bold")
legend(leg2,leg3,"location","northwest","FontSize",10)
grid on

% Print Data Statistics
fprintf("\nCombined Cost Statistics\n")
fprintf("Low Cost = %.2f [$M]\n",min(cost)/1000)
fprintf("Mean Cost = %.2f [$M]\n",mean(cost)/1000)
fprintf("High Cost = %.2f [$M]\n",max(cost)/1000)
fprintf("70%% Conf Cost = %.2f [$M]\n",seventy_perc/1000)
fprintf("\nLow Cost (30%% Margin) = %.2f [$M]\n",min(cost)*1.3/1000)
fprintf("Mean Cost (30%% Margin) = %.2f [$M]\n",mean(cost)*1.3/1000)
fprintf("High Cost (30%% Margin) = %.2f [$M]\n",max(cost)*1.3/1000)
fprintf("70%% Conf Cost (30%% Margin) = %.2f [$M]\n",cost_margin/1000)