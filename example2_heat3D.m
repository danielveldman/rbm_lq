clear all
close all
clc

%% spatial and temporal grids
nx = 16;        ny = 16;        nz = 16;
Lx = 1.5;       Ly = 1.5;       Lz = 1.5;
dx = Lx/(nx-1); dy = Ly/(ny-1); dz = Lz/(nz-1);

nn = nx*ny*nz;  node_pos = reshape(1:nn, nx, ny, nz);

K = 2^6;
T = 2;
tgrid = linspace(0,T,K+1);
dt = tgrid(2)-tgrid(1);

% construct a splitting of the A matrix by randomly selecting
% interconnections
M = 3;
AM = heat3D_Arandom(nx, ny, nz, dx, dy, dz, node_pos, M);
[LUAl, Al, L] = compute_LU(AM, 1, dt, M);
                
% control operator B and control U
nodes_top = node_pos(:,:,end); nodes_bottom = node_pos(:,:,1); nodes_side = node_pos(1,:,:);
m = 1;
B = sparse(nn, m); B(nodes_top) = 1;

% initial condition
xgrid = 0:dx:Lx;  ygrid = 0:dy:Ly;  zgrid = 0:dz:Lz;
[X,Y,Z] = ndgrid(xgrid, ygrid, zgrid);
X0 = exp(-(X.^2+0*Y.^2+Z.^2)/8/Lx^2); %tanh(X/Lx-0.5)+1;% ones(nn,1);% (sin(X/Lx).*sin(Y/Ly/3) + sin(2*X/Lx).*cos(Y/2/Ly)).*cos(3*Z/Lz);
X0 = X0(:);

% reference trajectory and weighting matrices
Xd = 0*sin(3*Z/Lz); Xd = Xd(:);
xd =@(t) Xd*sin(t.'/2);
Q = sparse(nn, nn);  Q(nodes_side, nodes_side) = eye(numel(nodes_side));
R = speye(m);


%% compute the forward dynamics for zero controls
U = zeros(1, K); % zero control
n_trials = 7;   % compute n_trials random realizations 
labels = {}; 
for trial = 1:n_trials
    batches = randi([1,L], 1, K);
    X = compute_XMLU(Al, LUAl, X0, B, U, tgrid, batches);
    
    figure(1)
    hold on
    plot(tgrid, X(end,:))
    xlabel 'time'
    ylabel 'x_h(\omega,t) at x=y=z=L'
    labels = {labels{:}, ['trial ', num2str(trial)]};
end
grid on
legend(labels, 'Location','southeast','NumColumns',2)
print('example2_heat3D_states', '-djpeg', '-r300')

%% compute optimal controls

n_trials = 7;   % compute n_trials random realizations 
labels = {}; 
for trial = 1:n_trials
    batches = randi([1,L], 1, K);
    [Uopt, Jh] = compute_controlXMLU(Al, LUAl, X0, B, U, Q, R, xd, tgrid, batches);
    disp(['Jh(\omega,u^*_h(\omega)) = ', num2str(Jh)])
    
    figure(2)
    hold on
    plot(tgrid(1:end-1)+diff(tgrid)/2, Uopt)
    xlabel 'time'
    ylabel 'u^*_h(\omega,t)'
    labels = {labels{:}, ['trial ', num2str(trial)]};
end
grid on
legend(labels, 'Location','southeast','NumColumns',2)
print('example2_heat3D_control', '-djpeg', '-r300')