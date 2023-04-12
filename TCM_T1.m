clear all
close all
% clc

%% Configuracion de las graficas
linewidth = 3;
fontsize = 20;
line = ["-", "--", "-.", ":", "-", "--", "-.", ":"];

%% Config. Solver
options = optimoptions('fsolve', 'Display','off', 'FunctionTolerance',1e-3);


%% Data
Ne = 50;     % Nº elementos
N = Ne+1;   % Nº nodos

% Ctes
sigma = 5.670373e-8;    % [W/m2K4]


% Condiciones contorno
T0 = 273.15;   % [K]
Tf = T0;       % [K]

% Heaters
Q_heater = 50;  % [W]
N_heater = 2;   % Nº heaters

Lx1_ht = 0.05;      % [m]   Separacion1 heater
Lx2_ht = 0.2;       % [m]   Separacion2 heater
Lx_ht = 0.1;        % [m]   Tamaño heater
Lz_ht = 0.1;        % [m]   Tamaño heater

% Propiedades materiales
Lx = 0.5;   % [m]
Lz = 0.3;   % [m]
Ly_CFRP = 2.5e-3;   % [m]
Ly_hon = 15e-3;     % [m]

rho_CFRP = 1650;    % [kg/m3]
rho_hon = 40;       % [kg/m3]

ce_CFRP = 930;      % [J/kg K]
ce_hon = 890;       % [J/kg K]

Kx_CFRP = 130;  % [W/(m K)]
Ky_CFRP = 1;    % [W/(m K)]
Kz_CFRP = 130;  % [W/(m K)]
Kx_hon = 0.7;   % [W/(m K)]
Ky_hon = 1.5;   % [W/(m K)]
Kz_hon = 0.7;   % [W/(m K)]

eps_CFRP = 0.8; % Emisividad


%% Apartado 1

% Conductividad efectiva
Keff_lon = (Kz_CFRP*Ly_CFRP*2+Kx_hon*Ly_hon) / (Ly_CFRP*2+Ly_hon);
Keff_tran = (Ly_CFRP*2+Ly_hon) / (Ly_CFRP/Ky_CFRP*2 + Ly_hon/Ky_hon);

V = Lx*Lz*(Ly_CFRP*2+Ly_hon);

q1 = - Keff_lon * T0;
q2 = - (Keff_lon * T0 + Q_heater*N_heater / V * Lx^2/2 + q1) / Lx;

X = linspace(0,Lx, N);
dx = X(2)-X(1); 

% Analitico
T_1_ana = - (Q_heater*N_heater/V*X.^2/2 + q2*X + q1) / Keff_lon;
% graphT_x(1, X, T_1_ana-T0, 'T [ºC]', 'X [m]', 'Caso 1. Analítico', linewidth, fontsize)

% Numerico
[Q_ht_1] = Heaters_1(N, dx, Q_heater*N_heater, T0, Tf, V, Keff_lon);
[MT_1] = TempSystem(N);

 T_1_num = MT_1\Q_ht_1;
% T_1_num = inv(MT_1)*Q_ht_1;

% graphT_x(2, X, T_1_num-T0, 'T [ºC]', 'X [m]', 'Caso 1. Numérico', linewidth, fontsize)

graphT_x_3(3, X, T_1_ana-T0, X, T_1_num-T0, 'T [ºC]', 'X [m]', 'Caso 1. Analítico y Numérico', ...
    ["Caso 1. Analítico", "Caso 1. Numérico"], linewidth, fontsize)

%% Apartado 2
N2 = N*10;
V_ht_2 = Lx_ht*Lz*(Ly_CFRP*2+Ly_hon);
X2 = linspace(0,Lx, N2);
dx2 = X2(2)-X2(1); 

% Zonas heaters
Lx12 = Lx1_ht;
Lx23 = (Lx1_ht+Lx_ht);
Lx34 = (Lx1_ht+Lx_ht+Lx2_ht);
Lx45 = (Lx1_ht+Lx_ht+Lx2_ht+Lx_ht);

% dx = 0.01;
X_Zone1 = [0:dx:Lx12];
X_Zone2 = [Lx12:dx:Lx23];
X_Zone3 = [Lx23:dx:Lx34];
X_Zone4 = [Lx34:dx:Lx45];
X_Zone5 = [Lx45:dx:Lx];

Q_heat_2 = Q_heater;
q_ht_2 = Q_heat_2/V_ht_2;
A = (Ly_CFRP*2+Ly_hon)*Lz; %Area plano YZ

% Analitico

% Zona 1
C1 = Q_heat_2/(Keff_lon*A);
C2 = T0;

T_2_ana_Zone1 = C1 * X_Zone1 + C2;

% Zona 2
C3 = Q_heat_2/(Keff_lon*A) + q_ht_2/(Keff_lon)* Lx12;
C4 = T_2_ana_Zone1(end) + q_ht_2*Lx12^2/(Keff_lon*2) - C3 * Lx12;

T_2_ana_Zone2 = -q_ht_2*X_Zone2.^2/(Keff_lon*2) + C3*X_Zone2 + C4;

%Zona 3 T=cte
T_2_ana_Zone3 = T_2_ana_Zone2(end)*ones(1,length(X_Zone3));

% Zona 4
C5 = -Q_heat_2/(Keff_lon*A) + q_ht_2/(Keff_lon)* Lx45;
C6 = T_2_ana_Zone3(end) + q_ht_2*Lx34^2/(Keff_lon*2) - C5 * Lx34;

T_2_ana_Zone4 = -q_ht_2*X_Zone4.^2/(Keff_lon*2) + C5*X_Zone4 + C6;

% Zone 5
C7 = -Q_heat_2/(Keff_lon*A);
C8 = T0+Q_heat_2/(Keff_lon*A)*Lx;

T_2_ana_Zone5 = C7 * X_Zone5 + C8;

T_2_ana = CatTemp(T_2_ana_Zone1, T_2_ana_Zone2, T_2_ana_Zone3, T_2_ana_Zone4, T_2_ana_Zone5);

% graphT_x(4, X, T_2_ana-T0, 'T [ºC]', 'X [m]', 'Caso 2. Analítico', linewidth, fontsize)
% graphT_x_2(5, X, T_1_ana-T0, X, T_2_ana-T0, 'T [ºC]', 'X [m]', 'Casos 1 y 2. Analítico', ...
%     ["Caso 1. Analítico", "Caso 2. Analítico"], linewidth, fontsize)

% Numerico
[Q_ht_2] = Heaters_2(N2, Q_heat_2, Lx, Lx1_ht, Lx2_ht, Lx_ht, dx2, V_ht_2, Keff_lon, T0);
[MT_2] = TempSystem(N2);

T_2_num = MT_2\Q_ht_2;
%  T_2_num = (MT_2)^(-1)*Q_ht_2;

% graphT_x(6, X2, T_2_num-T0, 'T [ºC]', 'X [m]', 'Caso 2. Numérico', linewidth, fontsize)
% graphT_x_2(7, X, T_1_num-T0, X2, T_2_num-T0, 'T [ºC]', 'X [m]', 'Casos 1 y 2. Numérico', ...
%     ["Caso 1. Numérico", "Caso 2. Numérico"], linewidth, fontsize)

graphT_x_3(8, X, T_2_ana-T0, X2, T_2_num-T0, 'T [ºC]', 'X [m]', 'Caso 2. Analítico y Numérico', ...
    ["Caso 2. Analítico", "Caso 2. Numérico"], linewidth, fontsize)


%% Apartado 3. Radiacion
Nx3=100; 
X = linspace(0,Lx,Nx3);
dx = X(2)-X(1); 

Nx1_3 = int32( Lx1_ht / Lx *Nx3 )+1;
Nx_ht_3 = int32( (Lx1_ht+Lx_ht/2) / Lx *Nx3 );
Nx2_3 = int32( (Lx1_ht+Lx_ht) / Lx *Nx3 )+1;
Nx_sens_3 = int32(Nx3/2);

Q_ht3 = Q_heater;
V_ht_3 = Lx_ht*Lz*(Ly_CFRP*2+Ly_hon);
V3 = Lx*Lz*(Ly_CFRP*2+Ly_hon);

% Q_ht = ones(length(X),1)*-0.2019;
[Q_ht_3] = Heaters_2(length(X), Q_ht3, Lx, Lx1_ht, Lx2_ht, Lx_ht, dx, V_ht_3, Keff_lon, T0);
[Q_rad_3] = Rad_Case3(length(X), sigma, eps_CFRP, Lz, Lx, Lx1_ht, Lx2_ht, Lx_ht, dx, Keff_lon, Ly_CFRP, Ly_hon, T0);
[T3_num] = CalculateSystem(options, X, Q_ht_3, Q_rad_3, T0);
[T3_num_lin] = CalculateSystemLin(options, X, Q_ht_3, Q_rad_3, T0);

Q_rad_3_000 = zeros(Nx3,1);
[T3_num_SinRad] = CalculateSystem(options, X, Q_ht_3, Q_rad_3_000, T0);

% graphT_x(9, X, T3_num-T0, 'T [ºC]', 'X [m]', 'Caso 3. Numérico', linewidth, fontsize)   
% graphT_x(10, X, T3_num_lin-T0, 'T [ºC]', 'X [m]', 'Caso 3. Numérico linealizado', linewidth, fontsize)   
graphT_x_3(11, X, T3_num-T0, X, T3_num_lin-T0, 'T [ºC]', 'X [m]', 'Caso 3. Numérico lineal. y sin lineal.', ... 
    ["Caso 3. Num. sin lineal.", "Caso 3. Num. lineal."], linewidth, fontsize)

graphT_x_3(12, X, T3_num_SinRad-T0, X, T3_num-T0, 'T [ºC]', 'X [m]', 'Caso 3. Con/Sin Radiación', ... 
    ["Sin radiación ", "Radiación"], linewidth, fontsize)

%% Apartado 4. Transitorio
dt = 60;
Time_max = 1200;
Time = [0:dt:Time_max];

V_CFRP = 2*Lx*Lz*Ly_CFRP;
m_CFRP = V_CFRP*rho_CFRP;

V_hon = Lx*Lz*Ly_hon;
m_hon = V_hon*rho_hon;

fv_CFRP = V_CFRP/(V_CFRP+V_hon);
fv_hon = V_hon/(V_CFRP+V_hon);

fm_CFRP = m_CFRP/(m_CFRP+m_hon);
fm_hon = m_hon/(m_CFRP+m_hon);

% Capacidad termica especifica equivalente [J/kg K]
Ce_e = ce_CFRP*fm_CFRP + ce_hon*fm_hon;
rho_e = rho_CFRP*fv_CFRP + rho_hon*fv_hon;
% rho_e = rho_CFRP*fm_CFRP + rho_hon*fm_hon;

% Solve Estacionario
Q_ht_4 = Q_ht_3;
Q_rad_4 = Q_rad_3;
[T4_num_est] = CalculateSystem(options, X, Q_ht_4, Q_rad_4, T0);

% Solve Transitorio
T4_num(1,:) = ones(1,length(X))*T0;
for i=2:length(Time)
    T4_num(i,:) = CalculateSystem_4(options, X, Q_ht_4, Q_rad_4, T0, rho_e, Ce_e, dt, Keff_lon, dx, A, T4_num(i-1,:));
end

X_s = int32(length(X)/2); % Position sensor
X_heat1 = int32((Lx1_ht+Lx_ht/2)/Lx *length(X)); % Pos. middle heater 1
X_heat2 = int32((Lx1_ht+Lx_ht+Lx2_ht+Lx_ht/2)/Lx *length(X)); % Pos. middle heater 2
T_max4_trans = max(T4_num)-T0; % Max temperature in every step

graphT_t(15, X, T4_num_est-T0, X, T4_num-T0, 'T [ºC]', 'X [m]', 'Caso 4. Transitorio', Time, linewidth, fontsize)
graphTpoint_t(16, Time, [T4_num(:,1), T4_num(:,Nx1_3), T4_num(:,Nx_ht_3), T4_num(:,Nx2_3), T4_num(:,Nx_sens_3)]-T0, ...
    'T [ºC]', 'Time [s]', 'Caso4. Evolución temporal temperatura', ... 
    ["X = 0 [m]", "X = 0.05 [m]", "X = 0.1 [m]", "X = 0.15 [m]", "X = 0.25 [m]"], linewidth, fontsize)

%% Apartado 5
Idx_sensor = int32( length(X)/2 );

Q_ht_5 = Q_ht_4;
Q_rad_5 = Q_rad_4;
T5_num(1,:) = ones(1,length(X))*T0;

E = 0;  % Energia Heaters
count=0;
SensON = false;

for i=2:length(Time)
    T_sensor = T5_num(i-1,Idx_sensor);
    
    if T_sensor > T0+12
        SensON = false;
        E = E + count;
        count = 0;
        Q_ht_5real = zeros(length(X),1);
    end
    
    if T_sensor < T0+10
        SensON = true;
        Q_ht_5real = Q_ht_5;
    end
    
    if SensON
        count = count + dt*Q_heater;
    end
    
    T5_num(i,:) = CalculateSystem_4(options, X, Q_ht_5real, Q_rad_5, T0, rho_e, Ce_e, dt, Keff_lon, dx, A, T5_num(i-1,:));
end

%Potencia media Heaters
Pmed5 = E/Time_max;

graphT_t(17, X, T4_num_est-T0, X, T5_num-T0, 'T [ºC]', 'X [m]', 'Caso 5. Transitorio con termostato', Time, linewidth, fontsize)
graphTpoint_t(18, Time, [T5_num(:,1), T5_num(:,Nx1_3), T5_num(:,Nx_ht_3), T5_num(:,Nx2_3), T5_num(:,Nx_sens_3)]-T0, ...
    'T [ºC]', 'Time [s]', 'Caso5. Evolución temporal temperatura', ... 
    ["X = 0 [m]", "X = 0.05 [m]", "X = 0.1 [m]", "X = 0.15 [m]", "X = 0.25 [m]"], linewidth, fontsize)

%% Apartado 6 Modelo 2D

% Definicion de malla
Nx6=30;
Ny6=14;
X6 = linspace(0,Lx,Nx6);
Y6 = linspace(0,(Ly_CFRP*2+Ly_hon),Ny6);

% Calculos
Q_ht_6 = Heaters2D(X6, Y6, Q_heater, Lx1_ht, Lx2_ht, Lx, Lx_ht, Ly_CFRP, Ly_hon, Lz, Kx_CFRP, Kx_hon, Ky_CFRP, Ky_hon, T0);
Q_rad_6 = Rad2D(X6, Y6, sigma, eps_CFRP, Lx1_ht, Lx2_ht, Lx, Lx_ht, Ly_CFRP, Ly_hon, Lz, Kx_CFRP, Kx_hon, Ky_CFRP, Ky_hon, T0); 
[Kx_6, Ky_6] = K2D(X6, Y6, Ly_CFRP, Ly_hon, Kx_CFRP, Kx_hon, Ky_CFRP, Ky_hon);

[T6] = Solve2D(options, X6, Y6, Q_ht_6, Q_rad_6, Kx_6, Ky_6, T0);

% Plot
Plot2D(20, T6-T0, "Modelo 2D Plano XY. T [ºC]")

Plot2D(21, Q_ht_6, "Modelo 2D Plano XY. Q_{ht}")
Plot2D(22, Q_rad_6, "Modelo 2D Plano XY. Q_{rad}")

graphT_x_4(23, X, T3_num-T0, X6, T6(1,:)-T0, X6, T6(int32(Ny6/2),:)-T0, X6, T6(end,:)-T0, ...
    'T [ºC]', 'X [m]', 'Temperatura caras modelo 2D', ...
    ["1D", "2D Y = 0 [m]", "2D Y = 0.01 [m]", "2D Y = 0.02 [m]"], linewidth, fontsize)



%% Apartado 7 Modelo 3D

% Definicion de malla
Nx7=10;
Ny7=6;
Nz7=int32(Lz/Lx*Nx7);
X7 = linspace(0,Lx,Nx7);
Y7 = linspace(0,(Ly_CFRP*2+Ly_hon),Ny7);
Z7 = linspace(0,Lz,Nz7);

Ly = Ly_CFRP*2+Ly_hon;
Ny1 = int32( Ly_CFRP / Ly *Ny7 )+1;
Ny2 = int32( (Ly_CFRP+Ly_hon) / Ly *Ny7 );
Ny3 = int32(Ny7/2);

Nz1 = int32( (Lz-Lz_ht)/2/Lz *Nz7 )+1;
Nz2 = int32( ((Lz-Lz_ht)/2+Lz_ht)/Lz *Nz7 );
Nz3 = int32( Nz7/2 );

% Calculos
Q_ht_7 = Heaters3D(X7, Y7, Z7, Q_heater, Lx1_ht, Lx2_ht, Lx, Lz, Lx_ht, Lz_ht, Ly_CFRP, Ly_hon, Kx_CFRP, Kx_hon, Ky_CFRP, Ky_hon, T0);
Q_rad_7 = Rad3D(X7, Y7, Z7, sigma, eps_CFRP, Lx1_ht, Lx2_ht, Lx, Lz, Lx_ht, Lz_ht, Ly_CFRP, Ly_hon, Kx_CFRP, Kx_hon, Ky_CFRP, Ky_hon, T0);
[Kx_7, Ky_7, Kz_7] = K3D(X7, Y7, Z7, Ly_CFRP, Ly_hon, Kx_CFRP, Kx_hon, Ky_CFRP, Ky_hon, Kz_CFRP, Kz_hon);
[T7] = Solve3D(options, X7, Y7, Z7, Q_ht_7, Q_rad_7, Kx_7, Ky_7, Kz_7, T0);

%Plot
Plot3D(30, T7-T0, 'Modelo 3D. T [ºC]')

Plot3D(31, Q_ht_7, 'Modelo 3D. Q_{ht7}')
Plot3D(32, Q_rad_7, 'Modelo 3D. Q_{rad7}')

Plot3DSectionsXZ(33, T7, T0, Nz7, Ny2, Ny3)
Plot3DSectionsXY(34, T7, T0, Nz1, Nz2, Nz3)

%% Save data

% save('Model_1D_Est3_50W.mat', 'X', 'T3_num', 'T0')
% save('Model_1D_Tran4_50W.mat', 'X', 'Time', 'T4_num', 'T4_num_est', 'T0')
% save('Model_1D_TranTermo5_50W.mat', 'X', 'Time', 'T5_num', 'T0', 'Pmed5')
% save('Model_2D_100x60.mat','Q_ht_6', 'Q_rad_6', 'T6')
% save('Model_3D_40x15x24.mat','Q_ht_7', 'Q_rad_7', 'T7')



%% Functions 3D

% Solve 3D
function [Temp] = Solve3D(options, X, Y, Z, Q_ht, Q_rad, Kx, Ky, Kz, T0)
    Nx = length(X);
    Ny = length(Y);
    Nz = length(Z);
    dx = X(2)-X(1);
    dy = Y(2)-Y(1);
    dz = Z(2)-Z(1);
    
    for m=1:Nx
        for n=1:Ny
            for c=1:Nz
                % Calculate Conditions
                Nodos_X0 = m==1;
                Nodos_XLx = m==Nx;
                Nodos_Y0 = n==1 && (m~=1 || m~=Nx);
                Nodos_YLy = n==Ny && (m~=1 || m~=Nx);
                Nodos_Z0 = c==1 && (m~=1 || m~=Nx);
                Nodos_ZLz = c==Nz && (m~=1 || m~=Nx);
                
                % Defnine equations
                
                % Cara X=0
                if Nodos_X0
                    a{n,m,c} = @(T) T(n,m,c)-T0;
                    
                % Cara X=Lx
                elseif Nodos_XLx 
                    a{n,m,c} = @(T) T(n,m,c)-T0;
                    
                % Cara Inferior Y=0
                elseif Nodos_Y0 
                   
                    % Arista Inferior en Z=0
                    if Nodos_Z0 
                        a{n,m,c} = @(T) Kx(n,m,c)*( T(n,m-1,c)+T(n,m+1,c)-2*T(n,m,c) )/dx^2 + Ky(n,m,c)*( T(n+1,m,c)-T(n,m,c) )/dy^2 + Kz(n,m,c)*( T(n,m,c+1)-T(n,m,c) )/dz^2 + Q_ht(n,m,c) + Q_rad(n,m,c)*(T(n,m,c)^4-T0^4);
                        
                    % Arista Inferior en Z=Lz
                    elseif Nodos_ZLz 
                        a{n,m,c} = @(T) Kx(n,m,c)*( T(n,m-1,c)+T(n,m+1,c)-2*T(n,m,c) )/dx^2 + Ky(n,m,c)*( T(n+1,m,c)-T(n,m,c) )/dy^2 + Kz(n,m,c)*( T(n,m,c-1)-T(n,m,c) )/dz^2 + Q_ht(n,m,c) + Q_rad(n,m,c)*(T(n,m,c)^4-T0^4);
                       
                    % Resto de Cara Inferior Y=0
                    else 
                        a{n,m,c} = @(T) Kx(n,m,c)*( T(n,m-1,c)+T(n,m+1,c)-2*T(n,m,c) )/dx^2 + Ky(n,m,c)*( T(n+1,m,c)-T(n,m,c) )/dy^2 + Kz(n,m,c)*( T(n,m,c-1)+T(n,m,c+1)-2*T(n,m,c) )/dz^2 + Q_ht(n,m,c) + Q_rad(n,m,c)*(T(n,m,c)^4-T0^4);
                        
                    end
                        
                % Cara Superior Y=Ly
                elseif Nodos_YLy 
                    
                    % Arista Superior en Z=0
                    if Nodos_Z0 
                        a{n,m,c} = @(T) Kx(n,m,c)*( T(n,m-1,c)+T(n,m+1,c)-2*T(n,m,c) )/dx^2 + Ky(n,m,c)*( T(n-1,m,c)-T(n,m,c) )/dy^2 + Kz(n,m,c)*( T(n,m,c+1)-T(n,m,c) )/dz^2 + Q_ht(n,m,c) + Q_rad(n,m,c)*(T(n,m,c)^4-T0^4);
                        
                    % Arista Superior en Z=Lz
                    elseif Nodos_ZLz 
                        a{n,m,c} = @(T) Kx(n,m,c)*( T(n,m-1,c)+T(n,m+1,c)-2*T(n,m,c) )/dx^2 + Ky(n,m,c)*( T(n-1,m,c)-T(n,m,c) )/dy^2 + Kz(n,m,c)*( T(n,m,c-1)-T(n,m,c) )/dz^2 + Q_ht(n,m,c) + Q_rad(n,m,c)*(T(n,m,c)^4-T0^4);
                        
                    % Resto de Cara Superior Y=Ly
                    else 
                        a{n,m,c} = @(T) Kx(n,m,c)*( T(n,m-1,c)+T(n,m+1,c)-2*T(n,m,c) )/dx^2 + Ky(n,m,c)*( T(n-1,m,c)-T(n,m,c) )/dy^2 + Kz(n,m,c)*( T(n,m,c-1)+T(n,m,c+1)-2*T(n,m,c) )/dz^2 + Q_ht(n,m,c) + Q_rad(n,m,c)*(T(n,m,c)^4-T0^4);
                        
                    end
                    
                % Cara Z=0 
                elseif Nodos_Z0
                    a{n,m,c} = @(T) Kx(n,m,c)*( T(n,m-1,c)+T(n,m+1,c)-2*T(n,m,c) )/dx^2 + Ky(n,m,c)*( T(n-1,m,c)+T(n+1,m,c)-2*T(n,m,c) )/dy^2 + Kz(n,m,c)*( T(n,m,c+1)-T(n,m,c) )/dz^2 + Q_ht(n,m,c) + Q_rad(n,m,c)*(T(n,m,c)^4-T0^4);
                 
                % Cara Z=Lz
                elseif Nodos_ZLz
                    a{n,m,c} = @(T) Kx(n,m,c)*( T(n,m-1,c)+T(n,m+1,c)-2*T(n,m,c) )/dx^2 + Ky(n,m,c)*( T(n-1,m,c)+T(n+1,m,c)-2*T(n,m,c) )/dy^2 + Kz(n,m,c)*( T(n,m,c-1)-T(n,m,c) )/dz^2 + Q_ht(n,m,c) + Q_rad(n,m,c)*(T(n,m,c)^4-T0^4);
                    
                % Nodos interiores
                else 
                    a{n,m,c} = @(T) Kx(n,m,c)*( T(n,m-1,c)+T(n,m+1,c)-2*T(n,m,c) )/dx^2 + Ky(n,m,c)*( T(n-1,m,c)+T(n+1,m,c)-2*T(n,m,c) )/dy^2 + Kz(n,m,c)*( T(n,m,c-1)+T(n,m,c+1)-2*T(n,m,c) )/dz^2;
                end
                
            end
        end
    end
    
    T0_ = ones(Ny,Nx,Nz)*T0;
    Temp = fsolve(@CalcTemps, T0_, options);
    
    function fun = CalcTemps(T)
        for mm=1:Nx
            for nn=1:Ny
                for cc=1:Nz
                    fun(nn,mm,cc) = a{nn,mm,cc}(T);
                end
            end
        end
    end
end

% Matrix Thermal Conductivity 3D (K)
function [Kx, Ky, Kz] = K3D(X, Y, Z, Ly_CFRP, Ly_hon, Kx_CFRP, Kx_hon, Ky_CFRP, Ky_hon, Kz_CFRP, Kz_hon)
    dx = X(2)-X(1);
    dy = Y(2)-Y(1);    
    dz = Z(2)-Z(1);
    
    Nx = length(X);
    Ny = length(Y);
    Nz = length(Z);
    
    Ly = Ly_CFRP*2 + Ly_hon;
    
    Ny1 = int32( Ly_CFRP / Ly *Ny )+1;
    Ny2 = int32( (Ly_CFRP+Ly_hon) / Ly *Ny );
    
    % Matrix K 
    Kx = zeros(Ny,Nx,Nz); 
    Ky = zeros(Ny,Nx,Nz); 
    Kz = zeros(Ny,Nx,Nz); 
    
    % Values
    Kx(1:Ny1,      :,:) = Kx_CFRP;
    Kx(Ny1+1:Ny2-1,:,:) = Kx_hon;
    Kx(Ny2:end,    :,:) = Kx_CFRP;

    Ky(1:Ny1,      :,:) = Ky_CFRP;
    Ky(Ny1+1:Ny2-1,:,:) = Ky_hon;
    Ky(Ny2:end,    :,:) = Ky_CFRP;
    
    Kz(1:Ny1,      :,:) = Kz_CFRP;
    Kz(Ny1+1:Ny2-1,:,:) = Kz_hon;
    Kz(Ny2:end,    :,:) = Kz_CFRP;
end

% Matrix Radiation 3D
function [Q_rad] = Rad3D(X, Y, Z, sigma, epsilon, Lx1_ht, Lx2_ht, Lx, Lz, Lx_ht, Lz_ht, Ly_CFRP, Ly_hon, Kx_CFRP, Kx_hon, Ky_CFRP, Ky_hon, T0)
    dx = X(2)-X(1);    
    dy = Y(2)-Y(1);
    dz = Z(2)-Z(1);

    Nx = length(X);
    Ny = length(Y);
    Nz = length(Z);
    
    Ly = Ly_CFRP*2 + Ly_hon;
    
    Nx1 = int32( Lx1_ht / Lx *Nx )+1;
    Nx2 = int32( (Lx1_ht+Lx_ht) / Lx *Nx )+1;
    Nx3 = int32( (Lx1_ht+Lx_ht+Lx2_ht) / Lx *Nx );
    Nx4 = int32( (Lx1_ht+Lx_ht+Lx2_ht+Lx_ht) / Lx *Nx );
    
    Nz1 = int32( (Lz-Lz_ht)/2/Lz *Nz )+1;
    Nz2 = int32( ((Lz-Lz_ht)/2+Lz_ht)/Lz *Nz );
    
    % Matrix Rad 
    Q_rad = zeros(Ny,Nx,Nz); 

    % Calor Rad. Cara Inf
    Q_rad(1,      1:end,      1:end) = 1;

    % Calor Rad. Cara Sup
    Q_rad(end,    1:Nx1-1,    1:end) = 1;
    Q_rad(end,    Nx2+1:Nx3-1,1:end) = 1;
    Q_rad(end,    Nx4+1:end,  1:end) = 1;
    Q_rad(end,    Nx1:Nx2,    1:Nz1-1) = 1;
    Q_rad(end,    Nx1:Nx2,    Nz2+1:end) = 1;
    Q_rad(end,    Nx3:Nx4,    1:Nz1-1) = 1;
    Q_rad(end,    Nx3:Nx4,    Nz2+1:end) = 1;

    q_rad = -epsilon*sigma/(dy);
    Q_rad = q_rad * Q_rad;
end

% Matrix Q_heaters 3D
function [Q_ht] = Heaters3D(X, Y, Z, Q_heater, Lx1_ht, Lx2_ht, Lx, Lz, Lx_ht, Lz_ht, Ly_CFRP, Ly_hon, Kx_CFRP, Kx_hon, Ky_CFRP, Ky_hon, T0)
    dx = X(2)-X(1);
    dy = Y(2)-Y(1);    
    dz = Z(2)-Z(1);
    
    Nx = length(X);
    Ny = length(Y);
    Nz = length(Z);
    
    Ly = Ly_CFRP*2 + Ly_hon;
    
    Nx1 = int32( Lx1_ht / Lx *Nx )+1;
    Nx2 = int32( (Lx1_ht+Lx_ht) / Lx *Nx )+1;
    Nx3 = int32( (Lx1_ht+Lx_ht+Lx2_ht) / Lx *Nx );
    Nx4 = int32( (Lx1_ht+Lx_ht+Lx2_ht+Lx_ht) / Lx *Nx );
    
    Nz1 = int32( (Lz-Lz_ht)/2/Lz *Nz )+1;
    Nz2 = int32( ((Lz-Lz_ht)/2+Lz_ht)/Lz *Nz );
    
    % Matrix Heater 
    Q_ht = zeros(Ny,Nx,Nz); 
    
    % Calor heaters en cara superior
    Q_ht(end,Nx1:Nx2,Nz1:Nz2) = 1;
    Q_ht(end,Nx3:Nx4,Nz1:Nz2) = 1;
    
    V = Lx_ht*Lz_ht*dy;
    q = Q_heater/V;
    Q_ht = q*Q_ht;
   
    % Contorno
    Q_ht(:,1,:) = T0;
    Q_ht(:,end,:) = T0;
end

% Plot points 3D
function Plot3D(f, M, tit)
    px = size(M,1);
    py = size(M,2);
    pz = size(M,3);

    [Y,X,Z] = ndgrid(1:px, 1:py, 1:pz);
    
    figure(f)
    scatter3(X(:), Z(:), Y(:), 500, M(:), 'filled')
    colorbar;
    xlabel('X');
    ylabel('Z');
    zlabel('Y');
    title(tit);
end 

function Plot3DSectionsXZ(f, T7, T0, Nz7, Ny2, Ny3)
    figure(f)
    subplot(2,2,1)
    imagesc(reshape(T7(end,:,:)-T0,[],Nz7)')
    xlabel('Nodos X');
    ylabel('Nodos Z');
    title('Y = 0.02 [m]')
    colorbar;

    subplot(2,2,2)
    imagesc(reshape(T7(Ny2,:,:)-T0,[],Nz7)')
    xlabel('Nodos X');
    ylabel('Nodos Z');
    title('Y = 0.0175 [m]')
    colorbar;

    subplot(2,2,3)
    imagesc(reshape(T7(Ny3,:,:)-T0,[],Nz7)')
    xlabel('Nodos X');
    ylabel('Nodos Z');
    title('Y = 0.01 [m]')
    colorbar;

    subplot(2,2,4)
    imagesc(reshape(T7(1,:,:)-T0,[],Nz7)')
    xlabel('Nodos X');
    ylabel('Nodos Z');
    title('Y = 0 [m]')
    colorbar;
end

function Plot3DSectionsXY(f, T7, T0, Nz1, Nz2, Nz3)
    figure(f)
    subplot(2,2,1)
    imagesc(T7(:,:,1)-T0)
    set(gca,'YDir','normal')
    xlabel('Nodos X');
    ylabel('Nodos Y');
    title('Z = 0 [m]')
    colorbar;

    subplot(2,2,2)
    imagesc(T7(:,:,Nz1)-T0)
    set(gca,'YDir','normal')
    xlabel('Nodos X');
    ylabel('Nodos Y');
    title('Z = 0.1 [m]')
    colorbar;

    subplot(2,2,3)
    imagesc(T7(:,:,Nz3)-T0)
    set(gca,'YDir','normal')
    xlabel('Nodos X');
    ylabel('Nodos Y');
    title('Z = 0.15 [m]')
    colorbar;

    subplot(2,2,4)
    imagesc(T7(:,:,Nz2)-T0)
    set(gca,'YDir','normal')
    xlabel('Nodos X');
    ylabel('Nodos Y');
    title('Z = 0.2 [m]')
    colorbar;
end

%% Funtions 2D
% Solve 2D
function [Temp] = Solve2D(options, X, Y, Q_ht, Q_rad, Kx, Ky, T0)
    Nx = length(X);
    Ny = length(Y);
    dx = X(2)-X(1);
    dy = Y(2)-Y(1);
    
    for m=1:Nx
        for n=1:Ny
            % Calculate Conditions
            Nodos_X0 = m==1;
            Nodos_XLx = m==Nx;
            Nodos_Y0 = n==1 && (m~=1 || m~=Nx);
            Nodos_YLy = n==Ny && (m~=1 || m~=Nx);
            
            % Defnine equations
            if Nodos_X0
                a{n,m} = @(T) T(n,m)-T0;
            elseif Nodos_XLx 
                a{n,m} = @(T) T(n,m)-T0;
            elseif Nodos_Y0 % Cara Inferior: Radiacion
                a{n,m} = @(T) Kx(n,m)*(T(n,m-1)+T(n,m+1)-2*T(n,m))/dx^2 + Ky(n,m)*(T(n+1,m)-T(n,m))/dy^2 + Q_rad(n,m)*(T(n,m)^4-T0^4);
            elseif Nodos_YLy % Cara Superior: Radiacion + Heaters
                a{n,m} = @(T) Kx(n,m)*(T(n,m-1)+T(n,m+1)-2*T(n,m))/dx^2 + Ky(n,m)*(T(n-1,m)-T(n,m))/dy^2 + Q_ht(n,m) + Q_rad(n,m)*(T(n,m)^4-T0^4);
            else % Nodos interior
                a{n,m} = @(T) Kx(n,m)*( T(n,m-1)+T(n,m+1)-2*T(n,m) )/dx^2 + Ky(n,m)*( T(n-1,m)+T(n+1,m)-2*T(n,m) )/dy^2;
%                 a{n,m} = @(T) ( T(n,m-1)+T(n,m+1)-2*T(n,m) )/dx^2 + ( T(n-1,m)+T(n+1,m)-2*T(n,m) )/dy^2;
            end
        end
    end
    
    T0_ = ones(Ny,Nx)*T0;
    Temp = fsolve(@CalcTemps, T0_, options);
    
    function fun = CalcTemps(T)
        for mm=1:Nx
            for nn=1:Ny
                fun(nn,mm) = a{nn,mm}(T);
            end
        end
    end
end

% Matrix Thermal Conductivity 2D (K)
function [Kx, Ky] = K2D(X, Y, Ly_CFRP, Ly_hon, Kx_CFRP, Kx_hon, Ky_CFRP, Ky_hon)
    Nx = length(X);
    Ny = length(Y);
    
    Ly = Ly_CFRP*2 + Ly_hon;
    
    Ny1 = int32( Ly_CFRP / Ly *Ny )+1;
    Ny2 = int32( (Ly_CFRP+Ly_hon) / Ly *Ny );
    
    
     % Matrix K
     Kx = zeros(Ny,Nx); 
     Ky = zeros(Ny,Nx); 
     
     % Values
     Kx(1:Ny1,      :) = Kx_CFRP;
     Kx(Ny1+1:Ny2-1,:) = Kx_hon;
     Kx(Ny2:end,    :) = Kx_CFRP;
     
     Ky(1:Ny1,      :) = Ky_CFRP;
     Ky(Ny1+1:Ny2-1,:) = Ky_hon;
     Ky(Ny2:end,    :) = Ky_CFRP;
end

% Matrix Radiation 2D
function [Q_rad] = Rad2D(X, Y, sigma, epsilon, Lx1_ht, Lx2_ht, Lx, Lx_ht, Ly_CFRP, Ly_hon, Lz, Kx_CFRP, Kx_hon, Ky_CFRP, Ky_hon, T0) 
    dx = X(2)-X(1);
    dy = Y(2)-Y(1);
    
    Nx = length(X);
    Ny = length(Y);
    
    Ly = Ly_CFRP*2 + Ly_hon;
    
    Nx1 = int32( Lx1_ht / Lx *Nx )+1;
    Nx2 = int32( (Lx1_ht+Lx_ht) / Lx *Nx )+1;
    Nx3 = int32( (Lx1_ht+Lx_ht+Lx2_ht) / Lx *Nx );
    Nx4 = int32( (Lx1_ht+Lx_ht+Lx2_ht+Lx_ht) / Lx *Nx );
    
    Ny1 = int32( Ly_CFRP / Ly *Ny )+1;
    Ny2 = int32( (Ly_CFRP+Ly_hon) / Ly *Ny );
    
    % Matrix Rad 
    Q_rad = zeros(Ny,Nx); 
    
    % Calor Rad. Cara Inf
    Q_rad(1,    1+1:end-1) = 1;
    
    % Calor Rad. Cara Sup
    Q_rad(end,    1+1:Nx1-1) = 1;
    Q_rad(end,    Nx2+1:Nx3-1) = 1;
    Q_rad(end,    Nx4+1:end-1) = 1;
    
    q_rad = -epsilon*sigma / dy;
    
    Q_rad = q_rad * Q_rad;
end

% Matrix Q_heaters 2D
function [Q_ht] = Heaters2D(X, Y, Q_heater, Lx1_ht, Lx2_ht, Lx, Lx_ht, Ly_CFRP, Ly_hon, Lz, Kx_CFRP, Kx_hon, Ky_CFRP, Ky_hon, T0)
    dx = X(2)-X(1);
    dy = Y(2)-Y(1);
    
    Nx = length(X);
    Ny = length(Y);
    
    Ly = Ly_CFRP*2 + Ly_hon;
    
    Nx1 = int32( Lx1_ht / Lx *Nx )+1;
    Nx2 = int32( (Lx1_ht+Lx_ht) / Lx *Nx )+1;
    Nx3 = int32( (Lx1_ht+Lx_ht+Lx2_ht) / Lx *Nx );
    Nx4 = int32( (Lx1_ht+Lx_ht+Lx2_ht+Lx_ht) / Lx *Nx );
    
    Ny1 = int32( Ly_CFRP / Ly *Ny )+1;
    Ny2 = int32( (Ly_CFRP+Ly_hon) / Ly *Ny );
    
   
    % Matrix Heater 
    Q_ht = zeros(Ny,Nx); 
    
    % Calor heaters
    Q_ht(end,Nx1:Nx2) = 1;
    Q_ht(end,Nx3:Nx4) = 1;
    
    V = Lz*dy*dx;
    q = Q_heater/(Lx_ht/dx) /V;
    Q_ht = q*Q_ht;
   
    % Contorno
    Q_ht(:,1) = T0;
    Q_ht(:,end) = T0;
end

% Plot 2D
function Plot2D(f, M, tit)
    figure(f);
    imagesc(M);
    set(gca,'YDir','normal')
    colorbar;
    % shading flat
    title(tit);
    xlabel('Nodos X');
    ylabel('Nodos Y');

end

%% Functions
% Matrix Temperatures, Contour
function [MT] = TempSystem(N)

    % Matrix system
    MT = -2*eye(N);
    MT(1,1) = 1;
    MT(end,end) = 1;
    
    for i=1:N-2
        MT(i+1,i) = 1;
        MT(i+1,i+2) = 1;
    end
end

% Vector Heaters Apartado 1
function [Q_ht] = Heaters_1(N, dx, q_ht, T0, Tf, V, Keff_lon)
    % Vector Heater 
    Q_ht = zeros(N,1);
    
    Q_ht(1,1) = T0;
    Q_ht(end,1) = Tf;
    
    Q_ht(2:end-1,1) = -dx^2*q_ht / (V*Keff_lon);
end


%  Vector Heaters Apartado 2
function [Q_ht] = Heaters_2(N, q_ht, Lx, Lx1_ht, Lx2_ht, Lx_ht, dx, V, Keff_lon, T0)
    N1 = round( Lx1_ht / Lx *N );
    N2 = round( (Lx1_ht+Lx_ht) / Lx *N );
    N3 = round( (Lx1_ht+Lx_ht+Lx2_ht) / Lx *N );
    N4 = round( (Lx1_ht+Lx_ht+Lx2_ht+Lx_ht) / Lx *N );

    % Vector Heater 
    Q_ht = zeros(N,1);
    
    Q_ht(N1:N2,1) = 1;
    Q_ht(N3:N4,1) = 1;
    Q_ht = -dx^2*q_ht / (V*Keff_lon) * Q_ht;
    
    Q_ht(1,1) = T0;
    Q_ht(end,1) = T0;
end

% Function enlazar tempertaturas analitico apartado 2
function [T] = CatTemp(T1, T2, T3, T4, T5)
    T1 = T1(1:end-1);
    T2 = T2(1:end-1);
    T3 = T3(1:end-1);
    T4 = T4(1:end-1);
    
    T = horzcat(T1, T2, T3, T4, T5);
end

% Radiation Vertex
function [Q_rad] = Rad_Case3(N, sigma, epsilon, Lz, Lx, Lx1_ht, Lx2_ht, Lx_ht, dx, Keff_lon, Ly_CFRP, Ly_hon, T0)
    N1 = round( Lx1_ht / Lx *N );
    N2 = round( (Lx1_ht+Lx_ht) / Lx *N );
    N3 = round( (Lx1_ht+Lx_ht+Lx2_ht) / Lx *N );
    N4 = round( (Lx1_ht+Lx_ht+Lx2_ht+Lx_ht) / Lx *N );
    Ly = Ly_CFRP*2+Ly_hon;

    % Vector Heater 
    Q_rad = zeros(N,1);
    
    Q_rad(1+1:N1-1,1) = 1;
    Q_rad(N2+1:N3-1,1) = 1;
    Q_rad(N4+1:end-1,1) = 1;
%     Q_rad = -epsilon*sigma*dx*2/Keff_lon * Q_rad;
%     Q_rad = -epsilon*sigma*dx^2/(Keff_lon*0.02) * Q_rad;
    Q_rad = -epsilon*sigma*dx^2*2/(Keff_lon*Ly) * Q_rad;
end


% Solve nodal system
function [T] = CalculateSystem_OLD(X, Q_in_out, Q_rad, T0)
    T_ = sym('a', [1, length(X)]);

    fun = zeros(1,length(X));
    fun(1,1) = T0;
    fun(1,end) = T0;

    eq(1,1) = fun(1,1) == T_(1);

    for i=2:length(X)-1
        eq(1,i) = fun(1,i) == T_(i-1)-2*T_(i)+T_(i+1) - Q_in_out(i) + Q_rad(i)*(T_(i)^4 - T0^4);
    end

    eq(1,end+1) = fun(1,end) == T_(end);

%     T0_ = ones(1,length(X))*T0;
%     Tsym_ = vpasolve(eq, T_, T0);
    Tsym_ = fsolve(eq, T0);

    for i=1:length(fun)
        aux = struct2cell(Tsym_);
        T(i) = vpa(aux(i));
    end
end

function [Temp] = CalculateSystem(options, X, Q_in_out, Q_rad, T0)

    len = length(X);
    a{1} = @(T) T(1)-T0;
    for i=2:len-1
        a{i} = @(T) T(i-1)-2*T(i)+T(i+1) - Q_in_out(i) + Q_rad(i)*(T(i)^4 - T0^4);
%         a{i} = @(T) T(i-1)-2*T(i)+T(i+1) - Q_in_out(i) + Q_rad(i)*( 4*T0^3*(T(i)-T0) ); % linealizado
    end
    a{end+1} = @(T) T(i+1)-T0;

    T0_ = ones(1,len)*T0;
    
    Temp = fsolve(@CalcTemps, T0_, options);
    
    
    function fun = CalcTemps(T)
        for j=1:len
            fun(j) = a{j}(T);
        end
    end
end

function [Temp] = CalculateSystemLin(options, X, Q_in_out, Q_rad, T0)

    len = length(X);
    a{1} = @(T) T(1)-T0;
    for i=2:len-1
        a{i} = @(T) T(i-1)-2*T(i)+T(i+1) - Q_in_out(i) + Q_rad(i)*( 4*T0^3*(T(i)-T0) ); % linealizado
    end
    a{end+1} = @(T) T(i+1)-T0;

    T0_ = ones(1,len)*T0;
    
    Temp = fsolve(@CalcTemps, T0_, options);
    
    
    function fun = CalcTemps(T)
        for j=1:len
            fun(j) = a{j}(T);
        end
    end
end

% Solve Transient Response
function [Temp] = CalculateSystem_4(options, X, Q_in_out, Q_rad, T0, rho, ce, dt, k, dx, A, T_old)
   
    len = length(X);
    a{1} = @(T) T(1)-T0;
    for i=2:len-1
        a{i} = @(T) ( T(i-1)-2*T(i)+T(i+1) )/dx^2 - Q_in_out(i)/dx^2 + Q_rad(i)/dx^2*(T(i)^4 - T0^4) - rho*ce/k *(T(i)-T_old(i))/dt;
%         a{i} = @(T) ( T(i-1)-2*T(i)+T(i+1) )/dx^2 - Q_in_out(i)/dx^2 + Q_rad(i)/dx^2*( 4*T0^3*(T(i)-T0) ) - rho*ce/k *(T(i)-T_old(i))/dt; % linealizado
    
    end
    
    a{end+1} = @(T) T(i+1)-T0;

    T0_ = ones(1,len)*500;
    Temp = fsolve(@CalcTemps, T0_, options);
    
    
    function fun = CalcTemps(T)
        for j=1:len
            fun(j) = a{j}(T);
        end
    end
end



%% Functions graphics
% T(x)
function graphT_x(f, X, T, ylab, xlab, tit, linewidth, fontsize)

    figure(f);
    hold on;
    grid on;
    plot(X, T, 'LineWidth', linewidth);
    hold off;

%     xlim( rad2deg([min(alpha), max(alpha)]) )
%     ylim([0, max(m)*1.3])
    ylabel(ylab);
    xlabel(xlab);
    title(tit);
    set(gca, 'FontSize', fontsize);
end

function graphT_x_2(f, X1, T1, X2, T2, ylab, xlab, tit, leg, linewidth, fontsize)

    figure(f);
    grid on;
    plot(X1, T1, 'LineWidth', linewidth);
    hold on;
    plot(X2, T2, 'LineWidth', linewidth);
    hold off;
    grid on;

%     xlim( rad2deg([min(alpha), max(alpha)]) )
%     ylim([0, max(m)*1.3])
    ylabel(ylab);
    xlabel(xlab);
    title(tit);
    legend(leg)
    set(gca, 'FontSize', fontsize);
end

function graphT_x_3(f, X, T1, X2, T2, ylab, xlab, tit, leg, linewidth, fontsize)

    figure(f);
    grid on;
    plot(X, T1, 'LineWidth', linewidth);
    hold on;
    plot(X2, T2, "-square", 'markersize', 5);
    hold off;
    grid on;

%     xlim( rad2deg([min(alpha), max(alpha)]) )
%     ylim([0, max(m)*1.3])
    ylabel(ylab);
    xlabel(xlab);
    title(tit);
    legend(leg)
    set(gca, 'FontSize', fontsize);
end

function graphT_x_4(f, X1, T1, X2, T2, X3, T3, X4, T4, ylab, xlab, tit, leg, linewidth, fontsize)

    figure(f)
    plot(X1, T1, 'LineWidth', linewidth);
    hold on;
    plot(X2, T2, "--", 'LineWidth', linewidth);
    hold on;
    plot(X3, T3, ":", 'LineWidth', linewidth);
    hold on;
    plot(X4, T4, "-.", 'LineWidth', linewidth);
    hold off;
    grid on;

%     xlim( rad2deg([min(alpha), max(alpha)]) )
%     ylim([0, max(m)*1.3])
    ylabel(ylab);
    xlabel(xlab);
    title(tit);
    legend(leg)
    set(gca, 'FontSize', fontsize);
end


% Temperature evolution function of Time
function graphT_t(f, X0, T0, X, T, ylab, xlab, tit, Time, linewidth, fontsize)
    line = ["-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":", "-", "--", "-.", ":"];

    figure(f);
    for i=1:length(T(:,1))
        plot(X, T(i,:), line(i), 'linewidth',linewidth/2);
        hold on;
    end
    plot(X0, T0, '--', 'color',"k", 'linewidth',linewidth)
    hold off;
    grid on;

%     xlim( rad2deg([min(alpha), max(alpha)]) )
%     ylim([0, max(m)*1.3])
    ylabel(ylab);
    xlabel(xlab);
    title(tit);
    labels = [compose('Time =  %d', Time'); 'Estacionario'];
    legend(labels,'location','best')
    set(gca, 'FontSize', fontsize);
end

% PointTemperature vs Time
function graphTpoint_t(f, Time, T, ylab, xlab, tit, leg, linewidth, fontsize)
    line = ["-", "--", "-.", ":", "-", "--", "-.", ":"];

    figure(f);
    for i=1:length(T(1,:)) % Filas: Evol. Temp. Columnas: Nodos
        plot(Time, T(:,i), line(i), 'LineWidth', linewidth);
        hold on;
    end
    hold off;
    grid on;

%     xlim( rad2deg([min(alpha), max(alpha)]) )
%     ylim([0, max(m)*1.3])
    ylabel(ylab);
    xlabel(xlab);
    title(tit);
    legend(leg)
    set(gca, 'FontSize', fontsize);
end
