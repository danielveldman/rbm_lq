clear all
close all
clc

%% spatial grids
N = 96;
s = 0.7;
L = 5;

[A, xgrid, Mass] = create_Amatrix_fl(N, 2*L, s);
xgrid = xgrid - L;
dx = xgrid(2) - xgrid(1);

K = 2^8;
T = 0.25;
tgrid = linspace(0,T,K+1);
dt = tgrid(2)-tgrid(1);

% create a splitting of the A matrix into M parts
P = 4;
AM = fl1D_parts(A, P);
M = length(AM);
pim = ones(1,M)/M;

% control operator B and control U
B = sparse(N, 2); 
B(-L/3 <= xgrid & xgrid <= 0   , 1) = dx;
B(L/3 <= xgrid & xgrid <= 2*L/3, 2) = dx;

% initial condition
X0 = exp(-0.4^2*xgrid.^2) - exp(-0.4^2*M^2);
X0 = X0(2:end-1);

% reference trajectory and weighting matrices
xd =@(t) 0;
% Q = 100*speye(N)*dx;
Q = sparse(N, N);  Q(xgrid<=0, xgrid<=0) = 100*eye(sum(xgrid<=0))*dx;
R = 1;

%% compute the forward dynamics for zero controls
U = zeros(2, K); % zero control
n_trials = 7;   % compute n_trials random realizations 
labels = {}; 
for trial = 1:n_trials
    batches = randi([1,M], 1, K);
    X = compute_XM(AM, X0, B, U, tgrid, batches, pim, Mass);
    
    figure(1)
    hold on
    plot(tgrid, X(end,:))
    xlabel 'time'
    ylabel 'x_h(\omega,t) at x=L'
    labels = {labels{:}, ['trial ', num2str(trial)]};
end
grid on
legend(labels, 'Location','northeast','NumColumns',2)
print('example3_fractional1D_states', '-djpeg', '-r300')

%% compute optimal controls

n_trials = 7;   % compute n_trials random realizations 
labels = {}; 
for trial = 1:n_trials
    batches = randi([1,M], 1, K);
    [Uopt, Jh] = compute_controlXM(AM, X0, B, U, Q, R, xd, tgrid, batches, pim, Mass);
    disp(['Jh(\omega,u^*_h(\omega)) = ', num2str(Jh)])
    
    figure(2)
    hold on
    plot(tgrid(1:end-1)+diff(tgrid)/2, Uopt(1,:))
    xlabel 'time'
    ylabel 'u^*_{h,1}(\omega,t)'
    
    figure(3)
    hold on
    plot(tgrid(1:end-1)+diff(tgrid)/2, Uopt(2,:))
    xlabel 'time'
    ylabel 'u^*_{h,1}(\omega,t)'
    labels = {labels{:}, ['trial ', num2str(trial)]};
end
figure(2)
grid on
legend(labels, 'Location','southeast','NumColumns',2)
print('example3_fractional1D_control1', '-djpeg', '-r300')
figure(3)
grid on
legend(labels, 'Location','northeast','NumColumns',2)
print('example3_fractional1D_control2', '-djpeg', '-r300')