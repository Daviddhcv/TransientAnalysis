%% ========================================================================
%  generate_figures_mdpi.m
%  ------------------------------------------------------------------------
%  Genera las 4 figuras del paper en formato MDPI
%  - Resolucion: 300 dpi
%  - Formatos: PNG + PDF vectorial
%  - Tamano: single-column (90 mm) o double-column (190 mm)
%
%  Ejecutar: >> generate_figures_mdpi
% ========================================================================

clear; clc; close all;

%% -------------------- CONFIGURACION ------------------------------------
out_dir = fullfile(pwd, 'figures');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontSize', 8);
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultLineLineWidth', 1.0);
set(0, 'DefaultAxesBox', 'on');
set(0, 'DefaultAxesTickDir', 'in');
set(0, 'DefaultAxesXGrid', 'on');
set(0, 'DefaultAxesYGrid', 'on');
set(0, 'DefaultAxesGridLineStyle', ':');
set(0, 'DefaultAxesGridAlpha', 0.35);

% Parametros del generador
p.f = 60;
p.ws = 2*pi*p.f;
p.V_peak = 220 * sqrt(2);
p.I_base = 244e3 / (sqrt(3)*380);
p.I_carga_pico = p.I_base * sqrt(2);
p.t_falla = 0.50;
p.t_end = 0.85;

% Constantes de tiempo
p.Td_p = 0.089;
p.Td_pp = 0.0023;
p.Ta = 0.015;

% Corrientes analiticas
MS.I_SLG_sub = 2282.08;
MS.I_SLG_tr = 1999.73;
MS.I_SLG_sost = 164.76;

MS.I_DLG_sub = 2532.53;
MS.I_DLG_tr = 1991.93;
MS.I_DLG_sost = 98.54;

MS.I3ph_sub = 1917.23;
MS.I3ph_tr = 1414.03;
MS.I3ph_sost = 57.42;

%% -------------------- VECTOR DE TIEMPO ---------------------------------
dt = 5e-6;
t = 0:dt:p.t_end;
idx_pre = t < p.t_falla;
idx_fault = t >= p.t_falla;
tau = zeros(size(t));
tau(idx_fault) = t(idx_fault) - p.t_falla;

%% ========================================================================
%  FIG 1: Tensiones nominales - fig2_vnom.png
% ========================================================================
va = p.V_peak * sin(p.ws*t);
vb = p.V_peak * sin(p.ws*t - 2*pi/3);
vc = p.V_peak * sin(p.ws*t + 2*pi/3);

h = figure('Color','w','Units','centimeters','Position',[2 2 14 6]);
plot(t*1000, va, 'b'); hold on;
plot(t*1000, vb, 'r');
plot(t*1000, vc, 'k');
xlabel('Time (ms)'); ylabel('Voltage (V)');
xlim([0 33]); ylim([-350 350]);
legend({'v_a','v_b','v_c'}, 'Location','northeast','NumColumns',3);
grid on;
exportgraphics(h, fullfile(out_dir,'fig2_vnom.png'), 'Resolution', 300);
exportgraphics(h, fullfile(out_dir,'fig2_vnom.pdf'), 'ContentType','vector');
close(h);

%% ========================================================================
%  FIG 2: Corriente SLG - fig4_slg_i.png
% ========================================================================
Iac_slg = MS.I_SLG_sub*exp(-tau/p.Td_pp) + ...
          (MS.I_SLG_tr - MS.I_SLG_sub)*exp(-tau/p.Td_p) + ...
          (MS.I_SLG_sost - MS.I_SLG_tr);
Idc_slg = MS.I_SLG_sub * exp(-tau/p.Ta);

ia_slg = zeros(size(t));
ia_slg(idx_pre) = p.I_carga_pico * sin(p.ws*t(idx_pre));
ia_slg(idx_fault) = sqrt(2)*Iac_slg(idx_fault).*sin(p.ws*tau(idx_fault)) + ...
                    sqrt(2)*Idc_slg(idx_fault);

h = figure('Color','w','Units','centimeters','Position',[2 2 14 6]);
plot(t*1000, ia_slg, 'b');
xlabel('Time (ms)'); ylabel('Fault current (A)');
xlim([450 850]);
grid on;
yline(0,'k','LineWidth',0.5);
exportgraphics(h, fullfile(out_dir,'fig4_slg_i.png'), 'Resolution', 300);
exportgraphics(h, fullfile(out_dir,'fig4_slg_i.pdf'), 'ContentType','vector');
close(h);

%% ========================================================================
%  FIG 3: Corriente DLG - fig7_dlg_i.png
% ========================================================================
Iac_dlg = MS.I_DLG_sub*exp(-tau/p.Td_pp) + ...
          (MS.I_DLG_tr - MS.I_DLG_sub)*exp(-tau/p.Td_p) + ...
          (MS.I_DLG_sost - MS.I_DLG_tr);
Idc_dlg = MS.I_DLG_sub * exp(-tau/p.Ta);

ib_dlg = zeros(size(t));
ic_dlg = zeros(size(t));
ib_dlg(idx_pre) = p.I_carga_pico * sin(p.ws*t(idx_pre) - 2*pi/3);
ic_dlg(idx_pre) = p.I_carga_pico * sin(p.ws*t(idx_pre) + 2*pi/3);
ib_dlg(idx_fault) = sqrt(2)*Iac_dlg(idx_fault).*sin(p.ws*tau(idx_fault) - 2*pi/3) + ...
                    sqrt(2)*Idc_dlg(idx_fault);
ic_dlg(idx_fault) = sqrt(2)*Iac_dlg(idx_fault).*sin(p.ws*tau(idx_fault) + 2*pi/3) + ...
                    sqrt(2)*Idc_dlg(idx_fault);

h = figure('Color','w','Units','centimeters','Position',[2 2 14 6]);
plot(t*1000, ib_dlg, 'r'); hold on;
plot(t*1000, ic_dlg, 'k');
xlabel('Time (ms)'); ylabel('Fault current (A)');
xlim([450 850]);
legend({'i_b (faulted)','i_c (faulted)'}, 'Location','northeast');
grid on;
yline(0,'k','LineWidth',0.5);
exportgraphics(h, fullfile(out_dir,'fig7_dlg_i.png'), 'Resolution', 300);
exportgraphics(h, fullfile(out_dir,'fig7_dlg_i.pdf'), 'ContentType','vector');
close(h);

%% ========================================================================
%  FIG 4: Corriente 3PH - fig9_3ph.png
% ========================================================================
Iac_3ph = MS.I3ph_sub*exp(-tau/p.Td_pp) + ...
          (MS.I3ph_tr - MS.I3ph_sub)*exp(-tau/p.Td_p) + ...
          (MS.I3ph_sost - MS.I3ph_tr);
Idc_3ph = MS.I3ph_sub * exp(-tau/p.Ta);

ia_3ph = zeros(size(t));
ib_3ph = zeros(size(t));
ic_3ph = zeros(size(t));
ia_3ph(idx_pre) = p.I_carga_pico * sin(p.ws*t(idx_pre));
ib_3ph(idx_pre) = p.I_carga_pico * sin(p.ws*t(idx_pre) - 2*pi/3);
ic_3ph(idx_pre) = p.I_carga_pico * sin(p.ws*t(idx_pre) + 2*pi/3);
ia_3ph(idx_fault) = sqrt(2)*Iac_3ph(idx_fault).*sin(p.ws*tau(idx_fault)) + ...
                    sqrt(2)*Idc_3ph(idx_fault);
ib_3ph(idx_fault) = sqrt(2)*Iac_3ph(idx_fault).*sin(p.ws*tau(idx_fault) - 2*pi/3) + ...
                    sqrt(2)*Idc_3ph(idx_fault);
ic_3ph(idx_fault) = sqrt(2)*Iac_3ph(idx_fault).*sin(p.ws*tau(idx_fault) + 2*pi/3) + ...
                    sqrt(2)*Idc_3ph(idx_fault);

h = figure('Color','w','Units','centimeters','Position',[2 2 14 6]);
plot(t*1000, ia_3ph, 'b'); hold on;
plot(t*1000, ib_3ph, 'r');
plot(t*1000, ic_3ph, 'k');
xlabel('Time (ms)'); ylabel('Fault current (A)');
xlim([450 850]);
legend({'i_a','i_b','i_c'}, 'Location','northeast','NumColumns',3);
grid on;
yline(0,'k','LineWidth',0.5);
exportgraphics(h, fullfile(out_dir,'fig9_3ph.png'), 'Resolution', 300);
exportgraphics(h, fullfile(out_dir,'fig9_3ph.pdf'), 'ContentType','vector');
close(h);

fprintf('\n[OK] Todas las figuras se guardaron en formato MDPI (300 dpi):\n');
fprintf('     %s\n', out_dir);
fprintf('     - fig2_vnom.png / .pdf\n');
fprintf('     - fig4_slg_i.png / .pdf\n');
fprintf('     - fig7_dlg_i.png / .pdf\n');
fprintf('     - fig9_3ph.png / .pdf\n\n');