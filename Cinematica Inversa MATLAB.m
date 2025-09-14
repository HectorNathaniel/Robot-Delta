% Cinemática Inversa para Robot Paralelo Tipo Delta
% Autor: Héctor Nathaniel López Ménez
% Fecha: 12/07/2025
% Descripción: El código calcula y verifica la posición del efector 
% final de un robot Delta mediante cinemática inversa y directa a partir
% de sus parámetros geométricos.

%% Limpieza de entorno
% Se eliminan variables, cierran ventanas y limpian la consola para 
% asegurar que la ejecución sea desde cero sin interferencias previas.

dbclear all; close all; clear; clc;

%% Definición de constantes geométricas
% Se declaran variables globales que representan las dimensiones físicas 
% del robot paralelo tipo Delta.

global L1 L2 ra rb r theta;

L1 = 0.17;        % Longitud del eslabón superior [m]
L2 = 0.32513;     % Longitud del eslabón inferior [m]
ra = 0.116;       % Radio de la plataforma fija [m]
rb = 0.05713;     % Radio del efector final [m]

r  = ra - rb;     % Diferencia de radios (distancia base a efector)
theta = [0; 2*pi/3; -2*pi/3];    % Ángulos de montaje de los motores en la base [rad]

%% Cargar archivos de datos experimentales
% El usuario selecciona archivos .txt con las posiciones deseadas del efector final (X, Y, Z)
% y ángulos phi de motores. 

[fileNamePos, pathPos] = uigetfile('*.txt','Seleccionar archivo de posiciones deseadas');
if ~isequal(fileNamePos,0)
    opts = detectImportOptions(fullfile(pathPos,fileNamePos),'FileType','text');
    opts.DataLines = [2 Inf]; % Ignora encabezado
    Tpos = readtable(fullfile(pathPos,fileNamePos), opts);
    Pd_list = table2array(Tpos(:,1:3));
    fprintf('Posiciones cargadas: %d muestras\n', size(Pd_list,1));
else
    fprintf('Omitida carga de posiciones.\n');
end

[fileNamePhi, pathPhi] = uigetfile('*.txt','Seleccionar archivo de ángulos phi');
if ~isequal(fileNamePhi,0)
    optsPhi = detectImportOptions(fullfile(pathPhi,fileNamePhi),'FileType','text');
    optsPhi.DataLines = [2 Inf];
    Tphi = readtable(fullfile(pathPhi,fileNamePhi), optsPhi);
    phi_list = table2array(Tphi(:,1:3));
    fprintf('Ángulos phi cargados: %d muestras\n', size(phi_list,1));
else
    fprintf('Omitida carga de ángulos phi.\n');
end

%% Determinar rango de trabajo teórico
% El rango de trabajo en X,Y está limitado por:
%   R_max = r + L1
% El rango en Z está limitado por:
%   0 <= Z <= L1

R_max = r + L1;
fprintf('Rango teórico: X,Y en [-%.3f, %.3f]  Z en [0, %.3f]\n',R_max,R_max,L1);

%% Solicitar posición final al usuario
% Se piden coordenadas finales del efector en (X,Y,Z).

prompt = {'X_final [m]:','Y_final [m]:','Z_final [m]:'};
defaul = {'0','0','0.23'};
dlg = inputdlg(prompt,'Posición final',1,defaul);

if isempty(dlg)
    error('No se introdujo posición.');
end

PosF = str2double(dlg);

% Validación de que la posición esté dentro del rango teórico
if any(PosF(1:2)<-R_max) || any(PosF(1:2)>R_max) || PosF(3)<0 || PosF(3)>L1
    error('Posición fuera de rango teórico.');
end

%% Cinemática Inversa
% Se calcula el ángulo de cada motor (phi) necesario para que el efector 
% alcance la posición deseada final PosF.

phiF = compute_phi(PosF);

fprintf('\nCinemática Inversa:');
for i=1:3
    fprintf(' φ%d = %.4f rad (%.2f°)\n',i,phiF(i),rad2deg(phiF(i)));
end

%% Cinemática Directa
% A partir de los ángulos calculados (phiF), se calcula la posición (X,Y,Z) 
% con cinemática directa y se verifica contra la posición deseada.

P_check = CD(phiF);
fprintf('\nCinemática Directa (verificación): X=%.4f, Y=%.4f, Z=%.4f\n',P_check);

fprintf('\nProceso completado.\n');

%%FUNCIONES AUXILIARES

%% Cinemática Inversa
% Entrada: posición deseada Pd = [Xp,Yp,Zp]
% Salida: ángulos phi de cada motor

function phi = compute_phi(Pd)
    global L1 L2 r theta;
    Xp=Pd(1); Yp=Pd(2); Zp=Pd(3);

    % Ecuaciones derivadas de la geometría del robot paralelo tipo Delta
    E = 2*L1*(r - Xp*cos(theta) - Yp*sin(theta));
    F = -2*L1*Zp*ones(3,1);
    G = Xp.^2 + Yp.^2 + Zp.^2 + r^2 + L1^2 - L2^2 - 2*r*(Xp*cos(theta) + Yp*sin(theta));

    % Discriminante
    D = E.^2 + F.^2 - G.^2;
    if any(D<0)
        error('Posición fuera de alcance. D negativo.');
    end

    % Fórmula de la cinemática inversa
    % phi = 2*atan( (-F - sqrt(D)) / (G - E) )
    phi = 2*atan((-F - sqrt(D)) ./ (G - E));
end

%% Cinemática Directa
% Entrada: ángulos de cada motor (phi)
% Salida: posición (X,Y,Z) del efector

function P = CD(phi)
    global L1 L2 r theta;

    % Cálculo de posiciones intermedias (eslabones superiores)
    X = (r + L1*cos(phi)) .* cos(theta);
    Y = (r + L1*cos(phi)) .* sin(theta);
    Z = L1*sin(phi);

    % Existen varios casos numéricos para resolver la intersección de esferas.
    % Caso 1: Z casi iguales (plano paralelo)
    if abs(Z(3)-Z(1))<1e-4 && abs(Z(2)-Z(1))<1e-4
        M = [2*(X(3)-X(1)),2*(Y(3)-Y(1));
             2*(X(2)-X(1)),2*(Y(2)-Y(1))];
        b = [X(3)^2+Y(3)^2 - (X(1)^2+Y(1)^2);
             X(2)^2+Y(2)^2 - (X(1)^2+Y(1)^2)];
        w = M\b;
        Xc=w(1); Yc=w(2);
        Zc=Z(1)+sqrt(L2^2-(Xc-X(1))^2-(Yc-Y(1))^2);
        P=[Xc;Yc;Zc];

    % Caso 2: condición de simetría especial
    elseif abs(X(2)-X(3))<1e-5 && abs(Z(3)-Z(2))<1e-5 && abs(X(3)-X(1))>1e-5
        a=-(Z(3)-Z(1))/(X(3)-X(1));
        b=(X(3)^2-X(1)^2+Y(3)^2-Y(1)^2+Z(3)^2-Z(1)^2)/(2*(X(3)-X(1)));
        A=a^2+1; B=2*a*(b-X(3))-2*Z(3); 
        C=(b-X(3))^2+Y(3)^2-L2^2+Z(3)^2;
        Zm=(-B+sqrt(B^2-4*A*C))/(2*A);
        Xm=a*Zm+b; Ym=0;
        P=[Xm;Ym;Zm];

    % Caso 3: caso general
    else
        a11=2*(X(3)-X(1)); a21=2*(X(3)-X(2));
        a12=2*(Y(3)-Y(1)); a22=2*(Y(3)-Y(2));
        a13=2*(Z(3)-Z(1)); a23=2*(Z(3)-Z(2));
        b1=X(3)^2+Y(3)^2+Z(3)^2-(X(1)^2+Y(1)^2+Z(1)^2);
        b2=X(3)^2+Y(3)^2+Z(3)^2-(X(2)^2+Y(2)^2+Z(2)^2);

        % Sistema de ecuaciones lineales reducido
        a1=a11/a13 - a21/a23;
        a2=a12/a13 - a22/a23;
        a3=b2/a23 - b1/a13;

        a4=-a2/a1; a5=-a3/a1;
        a6=-(a21*a4 + a22)/a23; a7=(b2 - a21*a5)/a23;

        % Resolución de intersección
        A=a4^2+1+a6^2;
        B=2*a4*(a5-X(1))-2*Y(1)+2*a6*(a7-Z(1));
        C=a5*(a5-2*X(1))+a7*(a7-2*Z(1))+X(1)^2+Y(1)^2+Z(1)^2-L2^2;

        Yp=(-B+sqrt(B^2-4*A*C))/(2*A);
        Xp=a4*Yp+a5; Zp=a6*Yp+a7;

        Ym=(-B-sqrt(B^2-4*A*C))/(2*A);
        Xm=a4*Ym+a5; Zm=a6*Ym+a7;

        % Selección de la solución válida (Z>0)
        if Zp>0
            P=[Xp;Yp;Zp];
        else
            P=[Xm;Ym;Zm];
        end
    end
end
