function J0 = compute_costJ(A, X0, B, u0, Q, R, xd, tgrid, Mass)

    dt = diff(tgrid);

    Aref{1} = A{1};
    for ii = 2:length(A)
        Aref{1} = Aref{1} + A{ii};
    end
    batches = ones(1,length(tgrid)-1);
    pim = 1;
    xsol0 = compute_XM(Aref, X0, B, u0, tgrid, batches, pim, Mass);
    J0 = cost_function(Q,R,xd,xsol0,u0,tgrid,dt);

end

function J = cost_function(Q,R,xd,x,u,tgrid,dt)
J = 0;
ndt = length(dt);
dx = x - xd(tgrid.');
for ii = 1:ndt
    J = J + dt(ii)*(dx(:,ii+1).'*Q*dx(:,ii+1)/2 + dx(:,ii).'*Q*dx(:,ii)/2 ...
        + u(:,ii).'*R*u(:,ii));
end
end