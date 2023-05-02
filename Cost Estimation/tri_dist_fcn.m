function cost = tri_dist_fcn(low,avg,high,NumSims,years,inf_facs)

pd = makedist('Triangular','A',low,'B',avg,'C',high);

cost = random(pd,NumSims,1);
for j = 1:NumSims
    cost(j) = inflation_factor(cost(j),years,inf_facs);
end