clc;
clear;
close all;

%% Parámetros de Entrada
fp = [800 1250];     % Banda de Paso [Hz]
fs = [200 5000];     % Banda de Rechazo [Hz]
Wp = 2*pi*fp;        % Banda de Paso [rad/s]
Ws = 2*pi*fs;        % Banda de Rechazo [rad/s]
Ap = 0.25;           % Atenuación máxima en Banda de Paso [dB]
As = 30;             % Atenuación mínima en Banda de Rechazo [dB]

%% Cálculo de la Función de Transferencia
[n,Wp] = cheb1ord(Wp,Ws,Ap,As,'s');
[num,den] = cheby1(n,Ap,Wp,'s');
Filtro = tf(num,den)                     % Función de transferencia calculada
[sos,g] = tf2sos(num,den);               % Descomposición en secciones bicuadráticas

% Implementación como PasaBajo / PasaAlto
PasaBajo = tf(2*g*sos(1,1:3),sos(1,4:6))
PasaAlto = tf(1/2*sos(2,1:3),sos(2,4:6))

%% Gráficos
figure;
hold on;

% Especificaciones del Filtro
plot([fs(1)/10 fs(1) fs(1)],[-As -As -Ap],'Color','r','LineWidth',3);
plot([fs(2) fs(2) fs(2)*10],[-Ap -As -As],'Color','r','LineWidth',3);
plot([fp(1) fp(1) fp(2) fp(2)],[-As -Ap -Ap -As],'Color','g','LineWidth',3);

% Respuesta del Filtro
h = bodeplot(Filtro);
p = getoptions(h);
p.PhaseVisible = 'off';
p.FreqUnits = 'Hz';
p.Grid = 'on';
setoptions(h,p);
bode(PasaBajo);
bode(PasaAlto);