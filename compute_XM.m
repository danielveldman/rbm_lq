function [X, duration] = compute_XM(AM, X0, B, U, tgrid, batches, pim, Mass)

tic
X = zeros(length(X0), length(tgrid));
X(:,1) = X0; N = length(X0);
dt = diff(tgrid); I = speye(length(X0)); P = size(batches,1);
for ii = 2:length(tgrid)
    batch = batches(:,ii-1);
    Aii = sparse(N,N);
    for p = batch.'
        Aii = Aii + AM{p} / pim(p);
    end
    X(:,ii) = (Mass-dt(ii-1)/2*Aii)\(Mass*X(:,ii-1)+dt(ii-1)/2*Aii*X(:,ii-1) + dt(ii-1)*B*U(:,ii-1));
end
duration = toc;