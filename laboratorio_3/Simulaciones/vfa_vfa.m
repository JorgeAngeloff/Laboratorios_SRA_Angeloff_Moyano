clc; clear; close all; 
% Bode LM324
% Parámetros LM324
Ad0 = 1e5;               % Ganancia DC
f1 = 10;                 % Hz
f2 = 5.06e6;             % Hz
w1 = 2*pi*f1;            % rad/s
w2 = 2*pi*f2;            % rad/s

% Función transferencia LM324
s = tf('s');
A_LM324 = Ad0 / ((1 + s/w1)*(1 + s/w2));

% Bode para verificar
f = logspace(0,8,1000);  % 1 Hz a 100 MHz
w = 2*pi*f;
[mag,phase] = bode(A_LM324,w);

mag = squeeze(20*log10(mag));
phase = squeeze(phase);

figure
subplot(2,1,1)
semilogx(f,mag,'b','LineWidth',1.5), grid on
xlabel('Frecuencia [Hz]'), ylabel('Magnitud [dB]')
title('Bode LM324 (Hz)')
xline(f1,'--r','f1'), xline(f2,'--r','f2')
yline(20*log10(Ad0),'--k','DC Gain')

subplot(2,1,2)
semilogx(f,phase,'m','LineWidth',1.5), grid on
xlabel('Frecuencia [Hz]'), ylabel('Fase [°]')
xline(f1,'--r','f1'), xline(f2,'--r','f2')

%% Bode del amplificador compuesto
clc; clear; close all; 

% Parámetros numéricos de la FT
num = 4.714e15;
den = [1 3.179e7 4.713e14];

% Función transferencia del amplificador compuesto
Avf_comp = tf(num, den);

% Frecuencia en Hz 
w = logspace(4, 8, 2000); % 10 kHz a 100 MHz
f = w/(2*pi);

% Respuesta en frecuencia
[mag,phase] = bode(Avf_comp, w);
mag = squeeze(mag);
phase = squeeze(phase);

% Ganancia en dB
magdB = 20*log10(mag);

% Ganancia máxima y punto -3 dB
mag_max = max(magdB);
fc_idx = find(magdB <= (mag_max-3),1,'first');
fc = f(fc_idx);

% Graficar Bode
figure;
subplot(2,1,1)
semilogx(f,magdB,'b','LineWidth',1.5), grid on
xlabel('Frecuencia [Hz]')
ylabel('Magnitud [dB]')
title('Bode del amplificador compuesto')

% Línea de referencia (20 dB)
yline(20,'--r','Avf_{ideal} = 20 dB','LabelHorizontalAlignment','left');

% Marcar frecuencia de corte
xline(fc,'--k',sprintf('f_c = %.2f MHz',fc/1e6), ...
    'LabelHorizontalAlignment','left','LabelOrientation','horizontal');

subplot(2,1,2)
semilogx(f,phase,'b','LineWidth',1.5), grid on
xlabel('Frecuencia [Hz]')
ylabel('Fase [°]')

% Marcar también la frecuencia de corte en fase
xline(fc,'--k',sprintf('f_c = %.2f MHz',fc/1e6), ...
    'LabelHorizontalAlignment','left','LabelOrientation','horizontal');
