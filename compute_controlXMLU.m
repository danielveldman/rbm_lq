function [u0, J0, duration] = compute_controlXMLU(A, LU, X0, B, u0, Q, R, xd, tgrid, batches)

dt = diff(tgrid);
tic

max_iters = 500; tol = 1e-6;
for ii = 1:max_iters
    xsol0 = compute_XMLU(A, LU, X0, B, u0, tgrid, batches);
    J0 = cost_function(Q,R,xd,xsol0,u0,tgrid,dt);
    
    phi = compute_phiMLU(A, LU, Q, xsol0, xd, tgrid, batches);
    gradJ = B.'*phi + 2*R*u0;
    G = innerproduct(gradJ,gradJ,dt);
    
    dx = compute_XMLU(A, LU, 0*X0, B, gradJ, tgrid, batches);
    H = hessian(Q,R,dx,gradJ,dt);
    step = G/H/2;
    
%         steps = linspace(0,-step,10);
%         for kk = 1:length(steps)
%             u1 = u0 + steps(kk)*gradJ;
%             xsol1 = compute_XMLU(A, LU, X0, B, u1, tgrid, batches);
%             J(kk) = cost_function(Q,R,xd,xsol1,u1,tgrid,dt);
%             J2(kk) = cost_function(Q,R,xd,xsol0+steps(kk)*dx,u0+steps(kk)*gradJ,tgrid,dt);
%         end
%         plot(steps, J, steps, J0+G*steps+H*steps.^2, steps, J2)
    
    u1 = u0 - step*gradJ;
    xsol1 = xsol0 - step*dx;
    J1 = cost_function(Q,R,xd,xsol1,u1,tgrid,dt);
    
    if abs(J1-J0) < tol*abs(J0)
    %if norm(u1-u0) < tol*norm(u0) && norm(xsol1-xsol0) < tol*norm(xsol0)
    %norm(u1-u0)/norm(u0)
    %ii
    duration = toc;
        return;
    end
    
    u0 = u1;
    %     xsol0 = xsol1;
end

disp('not converged')
duration = toc;
u0 = [];

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

function out = innerproduct(u,v,dt)
out = 0;
for ii = 1:length(dt)
    out = out + dt(ii)*(u(:,ii).'*v(:,ii));
end
end

function H = hessian(Q,R,dx,du,dt)
H = 0;
ndt = length(dt);
for ii = 1:ndt
    H = H + dt(ii)*(dx(:,ii+1).'*Q*dx(:,ii+1)/2 + dx(:,ii).'*Q*dx(:,ii)/2 ...
        + du(:,ii).'*R*du(:,ii));
end
end