function [LUAl, Al, L] = compute_LU(AM, P, dt, M)

L = nchoosek(M,P); % number of index subsets with nonzero probability
Sl = perms(1:M);
Sl = unique(sort(Sl(:,1:P),2), 'rows'); % all index combinations
pim = P*ones(1,M)/M;   % probability that m is an element of the selected subset
I = speye(length(AM{1}));
for l = 1:L
    Al{l} = 0*AM{1};
    for m = Sl(l,:)
        Al{l} = Al{l} + AM{m}/pim(m);
    end
    [LUAl{l}.L,LUAl{l}.U,LUAl{l}.P,LUAl{l}.Q,LUAl{l}.D] = lu(I - dt/2*Al{l});
end