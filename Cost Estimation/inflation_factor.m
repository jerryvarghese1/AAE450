function total_cost = inflation_factor(cost,years,inf_facs)

sub_cost = cost/years;
inf_cost(1) = sub_cost*inf_facs(1);

for j = 2:years
    inf_cost(j) = inf_cost(j-1) + sub_cost*inf_facs(j);
end

total_cost = inf_cost(years);

% % Plot Running Total Cost
% plot(1:years, inf_cost/1000, 'b*-')
% title("Operations Cost Running Total (Inflation)")
% xlabel("Years of Operation")
% ylabel("Cost [$M]")
% grid on