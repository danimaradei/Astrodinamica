% Graficos de los 6 elementos orbitales clasicos (COE) en funcion del
% tiempo, para verificar cuales son constantes y cual varia.
% Lee los datos calculados en Python (trayectoria_punto7.mat)

clear; clc; close all;

datos = load('trayectoria_punto7_1.mat');
t = datos.t;             % tiempo [s]
a_t = datos.a_t;         % semieje mayor [km]
e_t = datos.e_t;         % excentricidad [-]
i_t = datos.i_t;         % inclinacion [deg]
Omega_t = datos.Omega_t; % RAAN [deg]
omega_t = datos.omega_t; % argumento del perigeo [deg]
nu_t = datos.nu_t;       % anomalia verdadera [deg]

% --- 1) Semieje mayor ---
figure('Color', 'w');
plot(t, a_t, 'b-', 'LineWidth', 1.5);
xlabel('Tiempo [s]'); ylabel('a [km]');
title('Semieje mayor vs. tiempo');
grid on;

% --- 2) Excentricidad ---
figure('Color', 'w');
plot(t, e_t, 'b-', 'LineWidth', 1.5);
xlabel('Tiempo [s]'); ylabel('e [-]');
title('Excentricidad vs. tiempo');
grid on;

% --- 3) Inclinacion ---
figure('Color', 'w');
plot(t, i_t, 'b-', 'LineWidth', 1.5);
xlabel('Tiempo [s]'); ylabel('i [deg]');
title('Inclinacion vs. tiempo');
grid on;

% --- 4) RAAN ---
figure('Color', 'w');
plot(t, Omega_t, 'b-', 'LineWidth', 1.5);
xlabel('Tiempo [s]'); ylabel('\Omega [deg]');
title('RAAN vs. tiempo');
grid on;

% --- 5) Argumento del perigeo ---
figure('Color', 'w');
plot(t, omega_t, 'b-', 'LineWidth', 1.5);
xlabel('Tiempo [s]'); ylabel('\omega [deg]');
title('Argumento del perigeo vs. tiempo');
grid on;

% --- 6) Anomalia verdadera ---
figure('Color', 'w');
plot(t, nu_t, 'b-', 'LineWidth', 1.5);
xlabel('Tiempo [s]'); ylabel('\nu [deg]');
title('Anomalia verdadera vs. tiempo');
grid on;
