% Grafico 3D de la trayectoria orbital del punto 7
% Lee los datos calculados en Python (trayectoria_punto7.mat)

clear; clc; close all;

datos = load('trayectoria_punto7.mat');
X = datos.X;      % matriz 6xN: filas 1-3 = posicion, filas 4-6 = velocidad
Rt = datos.Rt;     % radio terrestre [km]

Xp = X(1, :);
Yp = X(2, :);
Zp = X(3, :);

figure('Color', 'w');
hold on;

% --- Tierra: esfera a escala ---
[xe, ye, ze] = sphere(40);
surf(Rt*xe, Rt*ye, Rt*ze, ...
    'FaceColor', [0.35 0.55 0.85], 'FaceAlpha', 0.35, 'EdgeColor', 'none');

% --- Ecuador, como referencia (z=0) ---
theta = linspace(0, 2*pi, 200);
plot3(Rt*cos(theta), Rt*sin(theta), zeros(size(theta)), ...
    'b--', 'LineWidth', 1);

% --- Trayectoria del satelite ---
plot3(Xp, Yp, Zp, 'r-', 'LineWidth', 1.8);

% --- Puntos notables ---
plot3(Xp(1), Yp(1), Zp(1), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
text(Xp(1), Yp(1), Zp(1), '  Perigeo (t=0)');

r_norm = sqrt(Xp.^2 + Yp.^2 + Zp.^2);
[~, idx_apogeo] = max(r_norm);
plot3(Xp(idx_apogeo), Yp(idx_apogeo), Zp(idx_apogeo), ...
    'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
text(Xp(idx_apogeo), Yp(idx_apogeo), Zp(idx_apogeo), '  Apogeo');

xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
title('Trayectoria orbital 3D');
axis equal;
grid on;
view(-50, 20);   % azimut, elevacion (mismo angulo que la version en Python)
legend({'Tierra', 'Ecuador', 'Trayectoria (1 periodo)'}, 'Location', 'northoutside');

hold off;
