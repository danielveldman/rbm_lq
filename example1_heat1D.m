clear all
close all
clc

%% spatial and temporal grids
n = 61;     % number of spatial grid points
L = 1.5;
xgrid = linspace(-L,L,n);
dx = L/(n-1);

K = 2^7;   % number of temporal grid points
T = 0.5;
tgrid = linspace(0,T,K+1);
dt = tgrid(2)-tgrid(1);

% split the 1D Laplacian into three parts
M = 3;
AM = heat1D_parts(n, dx, M);  % matrix for 1D heat equation
Mass = speye(n);
[LUAl, Al, L] = compute_LU(AM, 1, dt, M);   % precompute LU factorization

% control operator B and control U
B = sparse(n, 1); B(-L/3 <= xgrid & xgrid <= 0) = 1;

% initial condition
X0 = exp(-xgrid.^2) + xgrid.^2*exp(-L^2);
% plot(xgrid, X0)

% reference trajectory and weighting matrices
xd =@(t) 0;
Q = sparse(n, n);  Q(xgrid<=0, xgrid<=0) = 100*eye(sum(xgrid<=0))*dx;
R = 1;

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
    ylabel 'x_h(\omega,t) at x=L'
    labels = {labels{:}, ['trial ', num2str(trial)]};
end
grid on
legend(labels, 'Location','southeast','NumColumns',2)
print('example1_heat1D_states', '-djpeg', '-r300')

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
print('example1_heat1D_control', '-djpeg', '-r300')