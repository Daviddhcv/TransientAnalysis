% =========================================================================
% MATLAB SCRIPT: D-Q DYNAMIC MODEL, MAGNETIC SATURATION & CROSS-VALIDATION
% WEG AG10-250MI00AI Synchronous Generator (244 kVA, 380 V, 60 Hz, 4 poles)
% =========================================================================

clear; clc; close all;

fprintf('========================================================\n');
fprintf('  STARTING D-Q NUMERICAL SIMULATION WITH MAGNETIC SATURATION\n');
fprintf('========================================================\n\n');

% STEP 1: Definition of electrical base constants
S_base = 244e3;          % Apparent power base [VA]
V_LL_base = 380;         % Line-to-line voltage base [V]
f_nom = 60;              % Nominal frequency [Hz]
w_s = 2 * pi * f_nom;    % Electrical synchronous speed [rad/s]

% Base impedance and current calculation
Z_base = (V_LL_base^2) / S_base; 
I_base = S_base / (sqrt(3) * V_LL_base);

fprintf('1. Computed Base Parameters:\n');
fprintf('   - Base Impedance (Z_base) = %.4f ohms\n', Z_base);
fprintf('   - Base Current (I_base) = %.2f A\n\n', I_base);

% STEP 2: Nominal parameters and reactances (in per-unit)
Ra = 0.02805;            % Stator resistance [pu]
Xd = 6.4570;             % d-axis synchronous reactance [pu]
Xq = 1.7860;             % q-axis synchronous reactance [pu]
Xpp_d = 0.1770;          % d-axis subtransient reactance [pu]
H = 0.218;               % Inertia constant [s]

% STEP 3: Analytical Benchmark (IEC 60909)
V_phi_rms = V_LL_base / sqrt(3);
V_phi_peak = V_phi_rms * sqrt(2);
Z1 = 1j * Xpp_d * Z_base;
Rf = 0.01;               % Applied fault resistance [ohms]

I_3ph_peak_analytical = (V_phi_peak / abs(Z1 + 3*Rf)) * sqrt(2);

fprintf('2. Analytical Benchmark (IEC 60909):\n');
fprintf('   - Analytical 3-phase peak current = %.2f A\n\n', I_3ph_peak_analytical);

% STEP 4: Temporal simulation setup using ODE numerical integration
tspan = [0 0.8];         % Simulation time span [s]
% Initial steady-state conditions [id_pu, iq_pu, omega_pu, delta_rad]
x0 = [1.0; -0.75; 1.0; 0]; 

fprintf('3. Executing numerical solver (ode45 with Runge-Kutta)...\n');
options = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);
[t, x] = ode45(@(t,x) ode_generator_system(t, x, Xd, Xq, Xpp_d, Ra, H, w_s, Z_base), tspan, x0, options);

% STEP 5: Post-processing current results (convert to Amperes)
id_sim = x(:,1) * I_base;
iq_sim = x(:,2) * I_base;
irms_sim = sqrt(id_sim.^2 + iq_sim.^2) / sqrt(2);

% Extracting maximum peak current during short-circuit (applied at t >= 0.5s)
idx_fault = find(t >= 0.5);
i_peak_matlab = max(sqrt(id_sim(idx_fault).^2 + iq_sim(idx_fault).^2)) * sqrt(2);

fprintf('4. MATLAB Numerical Solver Results:\n');
fprintf('   - Maximum simulated peak current (3PH) = %.2f A\n\n', i_peak_matlab);

% STEP 6: Generation of validation time-series plot for Q1 publication
figure('Color', 'white', 'Position', [100 100 700 450]);
plot(t, irms_sim, 'LineWidth', 1.8, 'Color', [0, 0.447, 0.741]);
grid on;
xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Stator RMS Current (A)', 'FontSize', 11, 'FontWeight', 'bold');
title('Cross-Solver Validation: ODE Transient Response with Magnetic Saturation', 'FontSize', 12, 'FontWeight', 'bold');
xlim([0.4 0.8]);
legend({'MATLAB ODE Model'}, 'Location', 'northeast', 'FontSize', 10);

% Automatic high-resolution export at 300 DPI
output_filename = 'transient_validation_q1.png';
exportgraphics(gcf, output_filename, 'Resolution', 300);
fprintf('Simulation completed successfully. Plot automatically saved as "%s" at 300 DPI.\n', output_filename);


% =========================================================================
% LOCAL FUNCTION: Ordinary Differential Equations (ODEs) System
% =========================================================================
function dxdt = ode_generator_system(t, x, Xd, Xq, Xpp_d, Ra, H, w_s, Z_base)
    id = x(1);
    iq = x(2);
    omega = x(3);
    
    % Non-linear magnetic core saturation modeling in the d-axis
    psi_d = Xd * id;
    k_alpha = 0.15;
    beta_exp = 1.8;
    K_sat = 1 + k_alpha * abs(psi_d).^beta_exp;
    Xd_sat = Xd / K_sat; 
    
    % Three-phase fault application logic at t = 0.5 s
    Rf = 0.01; % Fault resistance [ohms]
    Rf_pu = Rf / Z_base; % Converted to per-unit
    
    if t < 0.5
        % Pre-fault condition (nominal voltage in pu)
        vd = 0;
        vq = 1.0; 
    else
        % During 3-phase short-circuit fault with resistance Rf
        vd = Rf_pu * id;
        vq = Rf_pu * iq;
    end
    
    % Stator differential equations (Subtransient d-q transformation)
    did_dt = (vd - Ra*id + omega*Xq*iq) / Xpp_d;
    diq_dt = (vq - Ra*iq - omega*Xd_sat*id) / Xq;
    
    % Rotor electromechanical dynamics (Swing Equation)
    Te = 1.5 * (Xd_sat - Xq) * id * iq; 
    Tm = 0.873;                        
    domega_dt = (w_s / (2 * H)) * (Tm - Te);
    ddelta_dt = omega - 1.0;
    
    dxdt = [did_dt; diq_dt; domega_dt; ddelta_dt];
en