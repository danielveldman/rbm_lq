function [X, duration] = compute_XMLU(AM, LU, X0, B, U, tgrid, batches)

tic
X = zeros(length(X0), length(tgrid));
X(:,1) = X0; dt = diff(tgrid);
for ii = 2:length(tgrid)
    m = batches(ii-1);
    X(:,ii) = LU{m}.Q*(LU{m}.U\(LU{m}.L\(LU{m}.P*(LU{m}.D\(X(:,ii-1)+dt(ii-1)/2*AM{m}*X(:,ii-1) + dt(ii-1)*B*U(:,ii-1))))));
end
duration = toc;