% Parámetros
Ad0 = 1e5;              % Ganancia DC
f1 = 10;                % Hz, primer polo
f2 = 5.06e6;            % Hz, segundo polo
Avf2 = 20;              % Ganancia fija del CFA
Avf = 10;               % Ganancia en lazo cerrado deseada

% Convertir frecuencias a rad/s
w1 = 2*pi*f1;
w2 = 2*pi*f2;

% Transferencia en lazo abierto del VFA
s = tf('s');
Avfa = Ad0 / ((1 + s/w1)*(1 + s/w2));

% Multiplicación por el CFA
Acomp = Avfa * Avf2;

% Función de transferencia en lazo cerrado
T = feedback(Acomp, 1/Avf);  % la realimentación ajustada a Avf

% Barrido de frecuencias
w = logspace(1,8,2000);   % de 10 Hz a 100 MHz
f = w/(2*pi);

[mag,phase] = bode(T, w);
mag = squeeze(mag);
phase = squeeze(phase);

% Diagrama de Bode
figure;
subplot(2,1,1)
semilogx(f,20*log10(mag),'b','LineWidth',1.5), grid on
xlabel('Frecuencia [Hz]')
ylabel('Magnitud [dB]')
title('Bode del amplificador compuesto VFA+CFA')

% Línea de referencia de ganancia
yline(20*log10(Avf),'--r','Ganancia DC','LabelHorizontalAlignment','left');

% Marcar frecuencia de -3 dB
[~,idx] = min(abs(20*log10(mag) - (20*log10(Avf)-3)));
f_3db = f(idx);
xline(f_3db,'--k',['f_{-3dB} = ' num2str(f_3db/1e6,'%.2f') ' MHz']);

subplot(2,1,2)
semilogx(f,phase,'b','LineWidth',1.5), grid on
xlabel('Frecuencia [Hz]')
ylabel('Fase [°]')

%%
clear; close all; clc;

% Parámetros
Ad0 = 1e5;
f1 = 10;              % Hz
f2 = 5.06e6;          % Hz
Avf2 = 20;            % ganancia del CFA (lineal)
Avf = 10;             % ganancia a lazo cerrado deseada (lineal)
w1 = 2*pi*f1;
w2 = 2*pi*f2;

% Modelos
s = tf('s');
Avfa = Ad0 / ((1 + s/w1)*(1 + s/w2));   % VFA (LM324) modelo 2 polos
Acomp = Avfa * Avf2;                    % VFA * CFA (en lazo abierto)
% Ganancia de lazo global 
beta = 1/Avf;
L = Acomp * beta;                       % loop gain (A*beta) - convención sin signo
Tcomp = feedback(Acomp, 1/Avf);         % closed-loop transfer (Vo/Vin)

% Respuesta al escalón
tfinal = 2e-6;     % tiempo final 
t = linspace(0, tfinal, 2000);
[y,t] = step(Tcomp, t);

% calcular sobrepico
y_final = y(end);
y_max = max(y);
Mp = (y_max - y_final) / y_final;   
Mp_pct = Mp * 100;

% Obtener zeta a partir de Mp (fórmula clásica para 2º orden)
if Mp <= 0
    zeta = NaN;
else
    lnMp = log(Mp);
    zeta = -lnMp / sqrt(pi^2 + lnMp^2);
end

% Graficar respuesta al escalón y anotar overshoot ---
figure('Units','normalized','Position',[0.2 0.2 0.5 0.5]);
plot(t*1e3, y, 'LineWidth', 1.5); grid on
xlabel('Tiempo [ms]'); ylabel('Salida (V)');
title('Respuesta al escalón del amplificador compuesto (VFA+CFA)');
% anotar
text(0.05*t(end)*1e3, 0.8*y_final, sprintf('Mp = %.2f %%', Mp_pct));
yline(y_final,'--k',sprintf('Valor final = %.4f', y_final),'LabelHorizontalAlignment','left');

