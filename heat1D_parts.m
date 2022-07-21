function AM = heat1D_parts(n, dx, M)

P = (n-1)/M;

if n-1 ~= M*P
    error('n-1 should be divisable by M')
end

AM{M} = [];
for m = 1:M
    AM{m} = sparse(n,n);
    for e = ((m-1)*P+1):m*P
        ind = e:e+1;
        AM{m}(ind,ind) = AM{m}(ind,ind) + [-1 1; 1 -1]/dx^2;
    end
    AM{m}(1,:)   = 2*AM{m}(1,:);
    AM{m}(end,:) = 2*AM{m}(end,:);
end