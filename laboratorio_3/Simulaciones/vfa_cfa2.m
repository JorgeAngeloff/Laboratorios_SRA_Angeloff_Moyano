% Parámetros 
Ad0 = 1e5;           % Ganancia DC del VFA 
f1 = 10;             % Hz (polo 1 del VFA)
f2 = 5.06e6;         % Hz (polo 2 del VFA)
Avf2_orig = 20;      % Ganancia CFA original 
Avf = 10;            % Ganancia deseada a lazo cerrado (20 dB)

% Compensador (Rx en paralelo con Cx; ese paralelo en serie con Ry)
Rx = 1e3;            % Ohm
Ry = 1e3;            % Ohm
Cx = 31e-12;         % Faradios 

% ajuste: 
Avf2_comp = 40;      % ganancia CFA ajustada en el caso compensado

% TFs elementales 
s = tf('s');

% VFA (modelo 2 polos)
w1 = 2*pi*f1;
w2 = 2*pi*f2;
Avfa = Ad0 / ((1 + s/w1)*(1 + s/w2));

% Compensador: 
Rpar = (Rx*Ry)/(Rx+Ry);       % Rx // Ry
kcomp = Ry/(Rx+Ry);           % atenuación DC del compensador
wz = 1/(Rx*Cx);               % cero angular del compensador
wp = 1/(Rpar*Cx);             % polo angular del compensador

Ac = kcomp * (1 + s/wz) / (1 + s/wp);   % TF del compensador

% Creación de bloques
Acomp_open = Avfa * Ac .* Avf2_comp;     %  (VFA * Ac * CFA)
Aopen_nocomp = Avfa * Avf2_orig;         %  sin compensador (para comparar)

% Cerrar lazo con beta = 1/Avf (topología no inversora global)
beta = 1/Avf;
T_comp = feedback(Acomp_open, beta);      % cerrado con compensador
T_nocomp = feedback(Aopen_nocomp, beta);  % cerrado sin compensador

% Barrido de frecuencia 
w = logspace(3,8,2500);   % rad/s (1 kHz a 100 MHz)
f = w/(2*pi);             % Hz

[mag_comp, phase_comp] = bode(T_comp, w);
[mag_noc,  phase_noc]  = bode(T_nocomp, w);

mag_comp = squeeze(mag_comp);
mag_noc  = squeeze(mag_noc);
phase_comp = squeeze(phase_comp);
phase_noc  = squeeze(phase_noc);

magdB_comp = 20*log10(mag_comp);
magdB_noc  = 20*log10(mag_noc);

% Determinar f_-3dB del sistema compensado 
gain0_comp_dB = 20*log10(abs(evalfr(T_comp, 0)));   % (dB)
target_dB = gain0_comp_dB - 3;
idx_3db = find(magdB_comp <= target_dB, 1, 'first');

if isempty(idx_3db)
    f_3db = NaN;
    warning('No se encontró punto -3dB dentro del rango evaluado.');
else
    f_3db = f(idx_3db);
end

% Dibujar bode 
figure('Units','normalized','Position',[0.05 0.05 0.7 0.7]);

subplot(2,1,1)
semilogx(f, magdB_noc, '--', 'LineWidth', 1.2, 'Color', [0.8 0.3 0.3]); hold on
semilogx(f, magdB_comp, '-',  'LineWidth', 1.6, 'Color', [0.1 0.6 0.1]);
grid on
xlabel('Frecuencia [Hz]')
ylabel('Magnitud [dB]')
title('Bode: Amplificador VFA + CFA (sin y con compensador) - lazo cerrado')
legend('Sin compensador (T_{noComp})','Con compensador (T_{comp})','Location','SouthWest')

% marcar ganancias DC y -3dB
yline(gain0_comp_dB, '--k', sprintf('Ganancia DC comp = %.2f dB', gain0_comp_dB), 'LabelHorizontalAlignment','left');
if ~isnan(f_3db)
    xline(f_3db, '--k', sprintf('f_{-3dB} = %.3f MHz', f_3db/1e6), 'LabelHorizontalAlignment','left');
end

% marcar fz, fp_comp y fg referencia
fz = wz/(2*pi);
fp_comp = wp/(2*pi);
fg_ref = 2e6;
xline(fz, ':m', sprintf('f_z = %.2f MHz', fz/1e6));
xline(fp_comp, ':b', sprintf('f_{p,comp} = %.2f MHz', fp_comp/1e6));
xline(fg_ref, ':c', sprintf('f_g (potencial) = %.2f MHz', fg_ref/1e6));

xlim([1e3 1e8])

subplot(2,1,2)
semilogx(f, phase_noc, '--','LineWidth',1.2, 'Color', [0.8 0.3 0.3]); hold on
semilogx(f, phase_comp, '-', 'LineWidth',1.4, 'Color', [0.1 0.6 0.1]);
grid on
xlabel('Frecuencia [Hz]')
ylabel('Fase [°]')
legend('Sin compensador','Con compensador','Location','SouthWest')
xline(fz, ':m', sprintf('f_z = %.2f MHz', fz/1e6));
xline(fp_comp, ':b', sprintf('f_{p,comp} = %.2f MHz', fp_comp/1e6));
xline(fg_ref, ':c', sprintf('f_g = %.2f MHz', fg_ref/1e6));
xlim([1e3 1e8]);


%% 
L = Acomp_open * beta;

% Vector de ventanas de tiempo (s)
t_windows = [1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3];

Mp_vals = zeros(size(t_windows));

for i = 1:length(t_windows)
    t_try = linspace(0, t_windows(i), 2000);
    y = step(T_comp, t_try);
    y = squeeze(y);
    y_final = y(end);
    y_max = max(y);
    if abs(y_final) < eps
        Mp_vals(i) = 0;
    else
        Mp_vals(i) = 100 * max(0, (y_max - y_final) / abs(y_final)); % % overshoot
    end
end

% Overshoot 
Mp_rob = max(Mp_vals);

idx_best = find(Mp_vals == Mp_rob, 1, 'first');

% Estimación de zeta y PM desde Mp (si Mp>0)
if Mp_rob > 0
    Mp_frac = Mp_rob/100;
    lnMp = log(Mp_frac);
    zeta = -lnMp / sqrt(pi^2 + lnMp^2);
    PM_rad_est = atan( 2*zeta / sqrt( sqrt(1 + 4*zeta^4) - 2*zeta^2 ) );
    PM_deg_est = PM_rad_est * 180/pi;
else
    zeta = NaN;
    PM_deg_est = NaN;
end

% Margen de fase real desde la ganancia de lazo
[GM, PM_bode, Wcg, Wcp] = margin(L);
Wcg_Hz = Wcg/(2*pi);
Wcp_Hz = Wcp/(2*pi);

% Imprimir resultados
fprintf('\n---- Resultados respuesta al escalón (robusto) ----\n');
fprintf('Mp (robusto) = %.4f %%\n', Mp_rob);
if ~isnan(zeta)
    fprintf('zeta (desde Mp) = %.4f\n', zeta);
    fprintf('PM estimado (desde Mp) = %.2f deg\n', PM_deg_est);
end
fprintf('PM (desde margin(L)) = %.2f deg\n', PM_bode);
if isfinite(Wcg_Hz)
    fprintf('Gain crossover (|L|=1) = %.3f Hz\n', Wcg_Hz);
else
    fprintf('Gain crossover (|L|=1) = Inf (no cruza en rango analizado)\n');
end

% Graficar la respuesta al escalón en la ventana con Mp máximo 
t_plot = linspace(0, t_windows(idx_best), 2000);
[y_plot, t_plot] = step(T_comp, t_plot);
y_plot = squeeze(y_plot);
figure;
plot(t_plot*1e6, y_plot, 'LineWidth', 1.4); grid on;
xlabel('Tiempo [\mus]'); ylabel('Salida (V)');
title(sprintf('Respuesta al escalón - ventana = %.3g s, Mp = %.3f %%', t_windows(idx_best), Mp_rob));
yline(y_plot(end), '--k', sprintf('Valor final = %.4g', y_plot(end)), 'LabelHorizontalAlignment','left');
text(0.05*t_plot(end)*1e6, 0.8*max(y_plot), sprintf('Mp = %.3f %%', Mp_rob));

