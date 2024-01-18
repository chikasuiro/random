%% read Excel file
T = readtable('P2.1.xlsx');
% T = readtable('P2.2.xlsx');
T = table2array(T);
nb = size(T, 1);  % n_b

%% random discretization
k = 60;
s = RandStream('dsfmt19937', 'Seed', 1234);
R = sort(randsample(s, nb, k));
nb = k;

%% set vectors of boundary values
Zeta_j   = T(:, 1) + 1i .* T(:, 2);  % vector of zeta_j
% Zeta_j   = Zeta_j(R);                % discretized
Zeta_j1  = circshift(Zeta_j, -1);    % vector of zeta_{j+1} 
Omega_j  = T(:, 3) + 1i .* T(:, 4);  % vector of Omega_j
% Omega_j  = Omega_j(R);               % discretized
Omega_j1 = circshift(Omega_j, -1);   % vector of Omega_{j+1}

%% evaluatie the stream function
z = 0.25 + 0.5i;                     % evaluation point
z = repelem(z, nb, 1);               % extend z to (nb,1)-vector
% calculate Omega-tilde by dealing BCs in vector format
val = 1 / (2 * pi * 1i) * sum( ...
    ((z - Zeta_j) .* Omega_j1 - (z - Zeta_j1) .* Omega_j) ./ (Zeta_j1 - Zeta_j) ...
    .* log ((Zeta_j1 - z) ./ (Zeta_j - z)) ...
    );
disp(imag(val));                     % output value

%% prepare internal mesh
X0 = 0.05 : 0.05 : 0.95;             % x coord of internal points
Y0 = 0.05 : 0.05 : 0.95;             % y coord of internal points
[Xint, Yint] = meshgrid(X0, Y0);     % internal points mesh
Z = Xint + 1i * Yint;                % z = x + iy
n = size(Xint, 1);                   % # of internal points along 1 axis
Z = repelem(Z, nb, 1);               % extend z to (nb,1)-vector 

%% tile boundary vectors on (n, n)
Zeta_j   = repmat(Zeta_j,   n);
Zeta_j1  = repmat(Zeta_j1,  n);
Omega_j  = repmat(Omega_j,  n);
Omega_j1 = repmat(Omega_j1, n);

%% calculate Omega values of internal points
A = ( ...
    ((Z - Zeta_j) .* Omega_j1 - (Z - Zeta_j1) .* Omega_j) ./ (Zeta_j1 - Zeta_j) ...
    .* log ((Zeta_j1 - Z) ./ (Zeta_j - Z)) ...
    );
% create group matrix
% each entry in (1:n) represents each group
G = reshape(1:1:n*n, n, n);
% extend group matrix to (nb*n, n)
% each (nb*1, 1)-vector in A is assigned to each group
G = repelem(G, nb, 1);
% calculate sum over each group, and divide by (2*pi*1i)
S = 1 / (2 * pi * 1i) * groupsummary(A, G, 'sum');

%% error estimation for P2.1
Ana = (Xint + 1i * Yint).^2;
err = sqrt(mean(abs(Ana - S).^2, 'all'));
% err = rmse(Ana, S, 'all');
disp(err);

%% prepare entire mesh
X1 = 0 : 0.05 : 1;                   % x coord of internal points
Y1 = 0 : 0.05 : 1;                   % y coord of internal points
[Xall, Yall] = meshgrid(X1, Y1);     % entire points mesh
Omega = nan(n+2) + 1i * nan(n+2);    % dummy entire Omega-tilde
Omega(2:n+1, 2:n+1) = S;             % assign internal values

%% assign boundary conditions for P2.1
% Omega_j(1:20) on y = 0
Omega(1, 1:n+1)   = Omega_j(1 : n+1);
% Omega_j(21:40) on x = 1
Omega(1:n+1, n+2) = Omega_j(n+2 : 2*(n+1));
% Omega_j(41:60) in inverse order on y = 1
Omega(n+2, 2:n+2) = Omega_j(3*(n+1) : -1 : 2*(n+1)+1);
% Omega_j(61:80) in inverse order on x = 0
Omega(2:n+2, 1)   = Omega_j(4*(n+1) : -1 : 3*(n+1)+1);

%% plot P2.1
fig_r = figure;
fig_r.PaperSize = [20 20];
contour(Xall, Yall, real(Omega), 101, 'k');
axis('square');
print(fig_r, '-dpdf', 'P2.1_r.pdf', '-fillpage');

fig_i = figure;
fig_i.PaperSize = [20 20];
contour(Xall, Yall, imag(Omega), 101, 'k');
axis('square');
print(fig_i, '-dpdf', 'P2.1_i.pdf', '-fillpage');

%% assign bounday conditions for P2.2
nlow = floor(n / 2);  % half # of points
% Omega_j(1:20) on y = 0
Omega(1, 1:n+1)           = Omega_j(1 : n+1);
% Omega_j(21:30) on x = 1 (y < 0.5)
Omega(1:nlow+1, n+2)      = Omega_j(n+2 : n+nlow+2);
% Omega_j(31:40) in reverse order on y = 0.5 (x > 0.5)
Omega(nlow+2, nlow+3:n+2) = Omega_j(2*(n+1) : -1 : n+nlow+3);
% skip Omega_j(41:42)
% Omega_j(43) at 0.5+0.5i
Omega(nlow+2, nlow+2)     = Omega_j(2*(n+1)+3);
% skip Omega_j(44:45)
% Omega_j(46:54) on x = 0.5 (y > 0.5)
Omega(nlow+3:n+1, nlow+2) = Omega_j(2*(n+1)+6 : 2*(n+1)+nlow+5);
% Omega_j(55:64) in inverse order on y = 1
Omega(n+2, 2:nlow+2)      = Omega_j(3*(n+1)+4 : -1 : 2*(n+1)+nlow+6);
% Omega_j(65:84) in inverse order on x = 0
Omega(2:n+2, 1)           = Omega_j(4*(n+1)+4 : -1 : 3*(n+1)+5);

%% plot P2.2
fig_r = figure;
fig_r.PaperSize = [20 20];
contour(Xall, Yall, real(Omega), 101, 'k');
axis('square');
hold on;
patch([0.5 1 1 0.5], [0.5 0.5 1 1], 'w', 'EdgeColor', 'None');
print(fig_r, '-dpdf', 'P2.2_r.pdf', '-fillpage');

fig_i = figure;
fig_i.PaperSize = [20 20];
contour(Xall, Yall, imag(Omega), 101, 'k');
axis('square');
hold on;
patch([0.5 1 1 0.5], [0.5 0.5 1 1], 'w', 'EdgeColor', 'None');
print(fig_i, '-dpdf', 'P2.2_i.pdf', '-fillpage');
