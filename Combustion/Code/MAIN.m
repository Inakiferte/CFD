%=======================================
%
% Equilibrium Composirion Under 
% Incomplete Combustion 
%
% Master in Space and Aeronautical
% Engineering.
% Thermal Turbomachinery and Combustion: Assignment 2.
% By: Iñaki Fernandez.
% Last modification: 06/06/2023
%
%=======================================
clc;
clear;
close all;
header;

%=======================================
%% INPUT PARAMETERS
%=======================================
addpath('./cp/');

% Universal data
R_hat = 8.31447;          % [kJ kmol^-1 K^-1]

% Fuel Data (Methane: CH4)
a      = 1.0;
b      = 4.0;
W_fuel = 16.04;                          % [kg/kmol]
[h_fhat_fuel,s_hat_fuel,a0_fuel,a1_fuel,a2_fuel,a3_fuel,a4_fuel] = CH4_properties;

% Air Data
W_O2  = 32;                              % [kg/kmol]
W_N2  = 28;                              % [kg/kmol]
W_air = W_O2 + 3.76 * W_N2;
[h_fhat_N2,s_hat_N2,a0_N2,a1_N2,a2_N2,a3_N2,a4_N2] = N2_properties;
[h_fhat_O2,s_hat_O2,a0_O2,a1_O2,a2_O2,a3_O2,a4_O2] = O2_properties;

% Exercise Data
Tfuel_in = 300;                          % [K]
Tair_in  = 400;                          % [K]
mdot_air = 0.1;                          % [kg/s]
lamda    = 0.9;                          % [kg/s]
T0       = 25 + 273;                     % [K]

% Fuel mass flow
stoic_rel = (a + b/4.0) * (W_air / W_fuel);
mdot_fuel = mdot_air / (lamda * stoic_rel);

% Species data
W_CO2 = 44.01;                           % [kg/kmol]
W_H2O = 18.02;                           % [kg/kmol]
W_CO  = 28.01;                           % [kg/kmol]
W_H2  = 2;                               % [kg/kmol]
[h_fhat_CO,s_hat_CO,a0_CO,a1_CO,a2_CO,a3_CO,a4_CO] = CO_properties;
[h_fhat_CO2,s_hat_CO2,a0_CO2,a1_CO2,a2_CO2,a3_CO2,a4_CO2] = CO2_properties;
[h_fhat_H2,s_hat_H2,a0_H2,a1_H2,a2_H2,a3_H2,a4_H2] = H2_properties;
[h_fhat_H2O,s_hat_H2O,a0_H2O,a1_H2O,a2_H2O,a3_H2O,a4_H2O] = H2O_properties;

%=======================================
%% QUESTION 1
%=======================================
addpath('./Kp/');

% Moran & Shaphiro data
T_exp  = [298,500,1000,1200,1400,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500]; % [K]
Kp_exp = [-5.018,-2.139,-0.159,0.135,0.333,0.474,0.530,0.577,0.619,0.656,0.688,0.716,0.742,0.764,0.784,0.802,0.818,0.833,0.846,0.858,0.869,0.878,0.888,0.895,0.902]; % [#]

% Build vectors
T      = linspace(298,3500,100); % [K]
Kp_sim = zeros(1,length(T)); 

for i=1:length(T)
    % Compute gibbs free energy
    g_CO  = gibbs_free_energy(h_fhat_CO,s_hat_CO,a0_CO,a1_CO,a2_CO,a3_CO,a4_CO,T(i));
    g_H2O = gibbs_free_energy(h_fhat_H2O,s_hat_H2O,a0_H2O,a1_H2O,a2_H2O,a3_H2O,a4_H2O,T(i));
    g_CO2 = gibbs_free_energy(h_fhat_CO2,s_hat_CO2,a0_CO2,a1_CO2,a2_CO2,a3_CO2,a4_CO2,T(i));
    g_H2 = gibbs_free_energy(h_fhat_H2,s_hat_H2,a0_H2,a1_H2,a2_H2,a3_H2,a4_H2,T(i));

    % Compute Kp
    Kp_sim(i) = Kp(g_CO,1,g_H2O,1,g_CO2,1,g_H2,1,T(i));
end

% Plot
figure;
fs = 25;
plot(T, log10(Kp_sim),'Color', 'b', 'LineWidth', 2);
hold on;
plot(T_exp, Kp_exp, 's--', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r', 'Color', 'r', 'LineWidth', 2);
xlabel('T [K]', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$\log_{10}(K_{p})$', 'Interpreter', 'latex', 'FontSize', fs);
xlim([T(1) T(end)]);
set(gca, 'ticklabelinterpreter', 'latex','FontSize', fs);
legend("Computed", "Experimental", 'Location', 'east', 'Interpreter', 'latex', 'FontSize', fs, 'NumColumns', 1);
grid on;
directory = 'Plots'; % Change 'Plots' to your desired directory
if ~exist(directory, 'dir')
    mkdir(directory);
end
set(gcf, 'WindowState', 'maximized');
filename = fullfile(directory, 'Question1.png');
saveas(gcf, filename);

%=======================================
%% QUESTION 2
%=======================================
addpath('./LU_Solver/');

% Compute Kp at given temperature
Tcc    = 1400;                              %[K]
g_CO   = gibbs_free_energy(h_fhat_CO,s_hat_CO,a0_CO,a1_CO,a2_CO,a3_CO,a4_CO,Tcc);
g_H2O  = gibbs_free_energy(h_fhat_H2O,s_hat_H2O,a0_H2O,a1_H2O,a2_H2O,a3_H2O,a4_H2O,Tcc);
g_CO2  = gibbs_free_energy(h_fhat_CO2,s_hat_CO2,a0_CO2,a1_CO2,a2_CO2,a3_CO2,a4_CO2,Tcc);
g_H2   = gibbs_free_energy(h_fhat_H2,s_hat_H2,a0_H2,a1_H2,a2_H2,a3_H2,a4_H2,Tcc);
Kp_Tcc = Kp(g_CO,1,g_H2O,1,g_CO2,1,g_H2,1,Tcc);


% Initial guess for the unknowns
n1_guess = a;
n2_guess = b / 2;
n3_guess = max(0.005,(lamda - 1) * (a + b/4));
n4_guess = max(0.005,(lamda - 1) * (a + b/4));
n5_guess = 3.76 * lamda * (a + b/4);
n_guess = [n1_guess; n2_guess; n3_guess; n4_guess; n5_guess];

% LU SOLVER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
[n_new,iter,delta_k] = LU_decom(n_guess,Kp_Tcc,lamda,a,b);

% Display the results
disp('=======================================')
disp('=             QUESTION 2              =')
disp('=======================================')
disp('Converged values are:');
disp(['[CO2]: ', num2str(n_new(1)) , ' mol']);
disp(['[H2O]: ', num2str(n_new(2)) , ' mol']);
disp(['[CO]: ', num2str(n_new(3)) , ' mol']);
disp(['[H2]: ', num2str(n_new(4)) , ' mol']);
disp(['[N2]: ', num2str(n_new(5)) , ' mol']);
disp(['Number of iterations: ', num2str(iter)]);

%=======================================
%% QUESTION 3
%=======================================
addpath('./cp/');
addpath('./Energy_Equation/');

% Define constants and initial guesses
n_species = 5; 
tol = 1e-7;
max_iter = 100;
iter = 0;
dif = 1;

% Compute Reactants enthalpies sum
h_hat_fuel = h_fhat_fuel + (cp_hat_int(Tfuel_in,a0_fuel,a1_fuel,a2_fuel,a3_fuel,a4_fuel) - cp_hat_int(T0,a0_fuel,a1_fuel,a2_fuel,a3_fuel,a4_fuel));
h_hat_O2   = h_fhat_O2 + (cp_hat_int(Tair_in,a0_O2,a1_O2,a2_O2,a3_O2,a4_O2) - cp_hat_int(T0,a0_O2,a1_O2,a2_O2,a3_O2,a4_O2));
h_hat_N2   = h_fhat_N2 + (cp_hat_int(Tair_in,a0_N2,a1_N2,a2_N2,a3_N2,a4_N2) - cp_hat_int(T0,a0_N2,a1_N2,a2_N2,a3_N2,a4_N2));
sumR       = h_hat_fuel + (lamda * (a + b / 4) * h_hat_O2) + (3.76 * lamda * (a + b / 4) * h_hat_N2);

% Initial guess for the unknowns concentrations
n1_guess = a;
n2_guess = b / 2;
n3_guess = max(0.005,(lamda - 1) * (a + b/4));
n4_guess = max(0.005,(lamda - 1) * (a + b/4));
n5_guess = 3.76 * lamda * (a + b/4);
n_guess = [n1_guess; n2_guess; n3_guess; n4_guess; n5_guess];

% Initial guess for the Tcc*
Tcc_s = 1400; % [K]

while dif > tol
    iter = iter + 1;
    
    % Evaluate Kp(Tcc*) and total moles from nk*
    g_CO   = gibbs_free_energy(h_fhat_CO,s_hat_CO,a0_CO,a1_CO,a2_CO,a3_CO,a4_CO,Tcc_s);
    g_H2O  = gibbs_free_energy(h_fhat_H2O,s_hat_H2O,a0_H2O,a1_H2O,a2_H2O,a3_H2O,a4_H2O,Tcc_s);
    g_CO2  = gibbs_free_energy(h_fhat_CO2,s_hat_CO2,a0_CO2,a1_CO2,a2_CO2,a3_CO2,a4_CO2,Tcc_s);
    g_H2   = gibbs_free_energy(h_fhat_H2,s_hat_H2,a0_H2,a1_H2,a2_H2,a3_H2,a4_H2,Tcc_s);
    Kp_Tcc = Kp(g_CO,1,g_H2O,1,g_CO2,1,g_H2,1,Tcc_s);
    n_total = sum(n_guess);
    
    % Solve LU decomposition
    [n_new,iter2,delta_k] = LU_decom(n_guess,Kp_Tcc,lamda,a,b);

    % Calculate dif2
    dif2 = norm(delta_k);
    
    % Compute Tcc from energy equation
    [Tcc] = energy_equation(Tcc_s,T0,n_new,sumR);    
    
    % Update dif and guess values
    dif = dif2 + abs(Tcc - Tcc_s);
    n_guess = n_new;
    Tcc_s = Tcc;
    
    % Check for convergence or maximum iterations
    if dif < tol || iter >= max_iter
        break;
    end
end

% Display results
disp('=======================================')
disp('=             QUESTION 3              =')
disp('=======================================')
disp('Converged values are:');
disp(['[CO2]: ', num2str(n_new(1)) , ' mol']);
disp(['[H2O]: ', num2str(n_new(2)) , ' mol']);
disp(['[CO]: ', num2str(n_new(3)) , ' mol']);
disp(['[H2]: ', num2str(n_new(4)) , ' mol']);
disp(['[N2]: ', num2str(n_new(5)) , ' mol']);
disp(['Tcc: ', num2str(Tcc) , ' K']);
disp(['Number of iterations: ', num2str(iter)]);

%=======================================
%% EXTRA
%=======================================
lambda = [0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99];
Tcc_lambda = [1331.2644, 1627.3247, 1850.3644, 2025.8631, 2168.7819,2288.38,2341.3624,2380.9539];
% Plot
figure;
fs = 25;
plot(lambda,Tcc_lambda, 's--', 'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b', 'Color', 'b', 'LineWidth', 2);
ylabel('$T_{cc}$ [K]', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', fs);
xlim([lambda(1) lambda(end)]);
set(gca, 'ticklabelinterpreter', 'latex','FontSize', fs);
grid on;
directory = 'Plots'; % Change 'Plots' to your desired directory
if ~exist(directory, 'dir')
    mkdir(directory);
end
set(gcf, 'WindowState', 'maximized');
filename = fullfile(directory, 'Question3.png');
saveas(gcf, filename);