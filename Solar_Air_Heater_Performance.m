Radi = input('Enter the solar radiation (W/m2) = ');
Ac   = input('Enter the collector area (m2) = ');
T1   = input('Enter the initial (inlet) temperature (C) = ');
Cp   = input('Enter the specific heat of air (J/kg.K) = ');
effi = input('Enter the thermal efficiency (fraction, 0-1) = ');
rho  = input('Enter the air density (kg/m3) = ');
A    = input('Enter the air passage area (m2) = ');
V0   = input('Enter the maximum air velocity (m/s) = ');
Dh   = input('Enter the hydraulic diameter (m) = ');
mu   = input('Enter the dynamic viscosity (Pa.s) = ');
K    = input('Enter the thermal conductivity (W/m.K) = ');
L    = input('Enter the air passage length (m) = ');
% Mass flow rate range
m_max = rho*A*V0;
if m_max > 0.01
    m_dot = 0.01:0.005:m_max;
else
    m_dot = linspace(0.2*m_max, m_max, 8);   % fallback range
end
% Thermal performance
T2 = T1 + effi*Radi*Ac ./ (m_dot*Cp);        % outlet temperature
Qu = m_dot .* Cp .* (T2 - T1);               % useful heat gain (W)
% Flow properties (vectors)
V  = m_dot ./ (rho*A);                       % velocity for each m_dot
Pr = (mu*Cp)/K;
Re = (rho*V*Dh)/mu;
DeltaT = T2 - T1;
PumpingPower = del_P .* (m_dot ./ rho);

Nu    = zeros(size(Re));
h     = zeros(size(Re));
f     = zeros(size(Re));
del_P = zeros(size(Re));
% Values at the edges of the transition region (for interpolation)
Nu_lam  = 3.66;        f_lam  = 64/2300;
Nu_turb = 0.023*4000^0.8*Pr^0.4;
f_turb  = 0.3164*4000^(-0.25);
for i = 1:length(Re)
    if Re(i) < 2300                          % Laminar
        Nu(i) = 3.66;
        f(i)  = 64/Re(i);
    elseif Re(i) > 4000                      % Turbulent
        Nu(i) = 0.023*Re(i)^0.8*Pr^0.4;
        f(i)  = 0.3164*Re(i)^(-0.25);
    else                                     % Transition (interpolated)
        w = (Re(i)-2300)/(4000-2300);
        Nu(i) = (1-w)*Nu_lam + w*Nu_turb;
        f(i)  = (1-w)*f_lam  + w*f_turb;
    end
    h(i)     = Nu(i)*K/Dh;
    del_P(i) = f(i)*(L/Dh)*(rho*V(i)^2/2);
end
% Results table
disp('---------------------------------------------------------------');
disp('m_dot      V        Re        Nu        h          DeltaP');
disp('(kg/s)    (m/s)     (-)       (-)     (W/m2K)      (Pa)');
disp('---------------------------------------------------------------');
for i = 1:length(m_dot)
    fprintf('%.3f     %.2f    %.0f     %.2f     %.2f     %.2f\n', ...
        m_dot(i), V(i), Re(i), Nu(i), h(i), del_P(i));
end
% Plots
figure(1);
plot(m_dot, T2, '-o');
xlabel('Mass Flow Rate (kg/s)');
ylabel('Outlet Temperature (^\circC)');
title('Outlet Temperature vs Mass Flow Rate');
grid on;

figure(2);
plot(m_dot, Qu, '-o');
xlabel('Mass Flow Rate (kg/s)');
ylabel('Useful Heat Gain (W)');
title('Useful Heat Gain vs Mass Flow Rate');
grid on;

figure(3);
plot(m_dot, Re, '-o');
xlabel('Mass Flow Rate (kg/s)');
ylabel('Reynolds Number');
title('Reynolds Number vs Mass Flow Rate');
grid on;

figure(4);
plot(m_dot, h, '-o');
xlabel('Mass Flow Rate (kg/s)');
ylabel('Heat Transfer Coefficient (W/m^2 K)');
title('Heat Transfer Coefficient vs Mass Flow Rate');
grid on;

figure(5);
plot(m_dot, del_P, '-o');
xlabel('Mass Flow Rate (kg/s)');
ylabel('Pressure Drop (Pa)');
title('Pressure Drop vs Mass Flow Rate');
grid on;

figure(6);
plot(m_dot, eta, '-o');
xlabel('Mass Flow Rate (kg/s)');
ylabel('Thermal Efficiency');
title('Thermal Efficiency vs Mass Flow Rate');
grid on;

figure(7);
plot(m_dot, DeltaT, '-o');
xlabel('Mass Flow Rate (kg/s)');
ylabel('Temperature Rise (°C)');
title('Temperature Rise vs Mass Flow Rate');
grid on;

figure(8);
plot(m_dot, PumpingPower, '-o');
xlabel('Mass Flow Rate (kg/s)');
ylabel('Pumping Power (W)');
title('Pumping Power vs Mass Flow Rate');
grid on;