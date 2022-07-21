function [phi, duration] = compute_phiMLU(A, LU, Q, X, xd, tgrid, batches)

tic
dt = diff(tgrid); ndt = length(dt);
phi = zeros(size(X,1), length(dt));
dx = X - xd(tgrid.');

m = batches(ndt);
phi(:,end) = LU{m}.Q*(LU{m}.U\(LU{m}.L\(LU{m}.P*(LU{m}.D\((dt(end)*Q*dx(:,end)))))));
for ii = ndt-1:-1:1
    m = batches(ii);
    phi(:,ii) = LU{m}.Q*(LU{m}.U\(LU{m}.L\(LU{m}.P*(LU{m}.D\(phi(:,ii+1)+dt(ii)/2*A{m}*phi(:,ii+1) + 2*dt(ii)*Q*dx(:,ii+1))))));
end
duration = toc;