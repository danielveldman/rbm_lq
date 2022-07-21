function [phi, duration] = compute_phiM(AM, Q, X, xd, tgrid, batches, pim, Mass)

tic
dt = diff(tgrid); ndt = length(dt);
phi = zeros(size(X,1), length(dt)); N = size(X,1);
dx = X - xd(tgrid.'); I = speye(N); P = size(batches,1);

batch = batches(:,end);
Aii = sparse(N,N);
for p = batch.'
    Aii = Aii + AM{p}/pim(p);
end
phi(:,end) = (I-dt(end)/2*Aii)\((dt(end)*Q*dx(:,end)));
for ii = ndt-1:-1:1
    batch = batches(:,ii);
    Aii = sparse(N,N);
    for p = batch.'
        Aii = Aii + AM{p}/pim(p);
    end
    phi(:,ii) = (Mass-dt(ii)/2*Aii)\(Mass*phi(:,ii+1)+dt(ii)/2*Aii*phi(:,ii+1) + 2*dt(ii)*Q*dx(:,ii+1));
end
duration = toc;