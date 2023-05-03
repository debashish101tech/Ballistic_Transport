L = 1; % Length of quantum well
N = 1000; % Number of points to use
x = linspace(0, L, N); % Spatial domain

% defined a square well
V0 = 1; 
V = zeros(size(x)); % Potential energy everywhere
V(x < 0.2*L | x > 0.8*L) = V0; % Define well

% Define parameters for solving Schrödinger's equation
hbar = 1; % Reduced Planck constant
m = 1; % Mass of particle
n_max = 10; % Maximum energy level to calculate

% Calculate eigen-energies and wavefunctions
E = zeros(n_max, 1);
psi = zeros(N, n_max);
for n = 1:n_max
    E(n) = (n^2*pi^2*hbar^2)/(2*m*L^2);
    k = sqrt(2*m*E(n)/hbar^2);
    psi(:, n) = sin(k.*x);
end

% Plot potential energy and wavefunctions
figure;
plot(x, V, 'k');
hold on;
plot(x, psi, 'LineWidth', 1.5);
xlabel('Position (L)');
ylabel('Wavefunction (\psi)');
title('Wavefunctions and Eigen-Energies for Quantum Well');
legend(['V(x)' repmat(['\psi_' num2str((1:n_max)')], 1, 1)], 'Location', 'NorthEast');
