%This code is designed to study Acuistic Leaky-Wave Antennas (ALWA), based on the theory of composite right left handed materials (CRLH)
%The model proposed is of cilindrical waveguide with axisymmetric open channels
%Suplemental material for:
%"Educational Open Source Kit for the Evaluation of Acoustic Leaky Wave Antennas with Metamaterials"
%Eduardo Romero-Vivas, Javier Romero-Vivas, Omar A. Bustamante, Braulio Leon-Lopez
%JASA Eduaction in Acoustics
%%Version 1.2, May 2021, Octave/Matlab
%
%
%
%ALWA parameters
%a - waveguide radius
%b - shunt width
%l - shunt length
%h - membrane thickness
%N - unit cell number
%d - unit cell length
%L - total ALWA length

close all
clear all

%%%%%%%%%%%%%%%%%%%%%%% Frequencies of study
f_beg = 1700;       %first frequency of study
f_end = 5000;       %last frequency of study
res   = 1000;       %number of samples between frequencies interval

f_th_list = linspace(f_beg,f_end,res);

%%%%%%%%%%%%%%%%%%%%%%% Characteristics of the medium

rho    = 0.9402;                   % air density STP reference - La Paz BCS, Mexico
c      = 342.4;                    % free-space sound velocity

%%%%%%%%%%%%%%%%%%%%%%% ALWA parameters

h        = 0.000067;        % membrane thickness
E        = 3.6e9;             % membrane Young's modulus
rho_mem  = 1370;              % membrane density
v_mem    = 0.33;            % membrane Poisson's ratio

a      = 0.0039;          % waveguide radius
b      = 0.0004;          % shunt width
l      = 0.0198;          % shunt length

S_a    = pi * (a^2);      % waveguide transversal area

N = 20;           % unit cell num
d = 0.0124;       % unit cell length
L = N * d;        % total ALWA length


%%%%%%%%%%%%%%%%%%%%%%% Elements of the impedance Z_se

M_wg = (rho/S_a) * (d-h);                                              % mass of the waveguide section
M_mem  = 1.8830 * ((rho_mem*h)/(pi * a^2));                            % mass of the membrane
C_mem     = (pi * a^6 ) / (196.51 * ( (E * h^3)/(12*(1-(v_mem^2)) ) ) );  % compliance of the membrane


%%%%%%%%%%%%%%%%%%%%%%% Elements of the admittance Y_sh

C_wg       = (S_a/(rho * c^2)) * (d-h);       % compliance of the waveguide section

f_zero          = 1 /(2*pi * sqrt( (M_wg + M_mem) * C_mem ) );   % resonance frequency of Z_se

M_shunt  = (rho/(2*pi*b)) * log(1 + l/a);            % shunt mass of the shunt
C_shunt     = ( 1 / ( 4 * pi^2 * (f_zero)^2 * M_shunt) ) - C_wg;  % shunt compliance of the shunt


%%%%%%%%%%%%%%%%%%%%%%% Transfer function

for num=1:length(f_th_list)
  for f=f_th_list(num)
    omega = 2 * pi * f;
    k     = omega/343.4;


    %%%%%%%%%%%%%%%%%%%%%%% Series impedance Z_se
    Z_se = 1i * ( (omega*(M_wg + M_mem)) - (1/(omega*C_mem)));

    %%%%%%%%%%%%%%%%%%%%%%% Shunt admittance Y_sh
    Y_sh = 1i * ( (omega*(C_shunt + C_wg)) - (1/(omega*M_shunt))  );


    %%%%%%%%%%%%%%%%%%%%%%% ABCD matrix  - Equation 15 from supplementary 1
    %   Y/2 --- Z --- Y/2  cell configuration

    A = 1 + ((Z_se * Y_sh) / 2);
    B = Z_se;
    C = Y_sh * (1 + ((Z_se * Y_sh) / 4));
    D = 1 + ((Z_se * Y_sh) / 2);

            %ABCD = [A, B, C, D]

    %%%%%%%%%%%%%%%%%%%%%%% Bloch parameters

    gamma_bloch = (acosh(A)) / d;          %Equation 1 from supplementary 1
    Z_bloch     = sqrt(B/C);               %Equation 2 from supplementary 1
    gamma_bloch_list(num)=gamma_bloch;
    Z_bloch_list(num)=Z_bloch;

  end
end


%%%%%%%%%%%%%%%%%%%%%%% Offset freq

aRef = f_zero;                        % Reference to set the transition frequency f_zero within f_th_list
aDiff = abs(f_th_list - aRef);        % Calculate the difference
[minVal, offset] = min(aDiff);        % Find closest value to aRef


%%%%%%%%%%%%%%%%%%%%%%% Alpha Bloch - attenuator factor
alpha_list = real(gamma_bloch_list);

%%%%%%%%%%%%%%%%%%%%%%% Beta Bloch - phase constant
beta_list        = imag(gamma_bloch_list);
beta_LH          = -(beta_list(1:offset));    % LH region
beta_RH          = beta_list(offset+1:end);   % RH region
beta_ABCD_CRLH   = [beta_LH, beta_RH];        % full CRLH

for i=1:length(beta_list)
  beta_d_ABCD(i) = d*beta_ABCD_CRLH(i);
end


%%%%%%%%%%%%%%%%%%%%%%% Wavenumber k = omega/c
for num=1:length(f_th_list)
  k_list(num) = ( (2*pi*f_th_list(num))/c );
end

k_d_list = k_list * d;
k_L =-1 * k_list(1:offset);
k_R = k_list(offset+1:end);


%%%%%%%%%%%%%%%%%%%%%%% Cutoff frequencies

f_L  = 1 / (2*pi*sqrt( M_shunt * C_mem ));                       % LH resonance frequency
f_R  = 1 /(2*pi*sqrt( (M_wg + M_mem) * (C_shunt + C_wg) ));      % RH resonance frequency

fcL  = f_R * abs( 1 - (sqrt( 1 + ( f_L / f_R ) )) );     % LH cut off frequency
fcR  = f_R   *  ( 1 + (sqrt( 1 + ( f_L / f_R ) )) );     % RH cut off frequency

f1_list      = abs(beta_LH-k_L) < 1;      % finding f1 
f1_idx_list  = find (f1_list);
f1_idx       = int16(mean(f1_idx_list));
f1           = f_th_list(f1_idx(1));

f2_list      = abs(beta_RH-k_R) < 1;      % finding f2 
f2_idx_list  = find (f2_list);
f2_idx       = int16(mean(f2_idx_list));
f2           = f_th_list(offset + f2_idx(1));

%%%%%%%%%% vertical bars for cutoff frequencies 

fl_x = [fcL fcL];    
fl_y = [-4 4];

fr_x = [fcR fcR];
fr_y = [-4 4];

f1_x = [f1 f1];
f1_y = [-4 4];

f2_x = [f2 f2];
f2_y = [-4 4];

f_zero_x = [f_zero f_zero];
f_zero_y = [-4 4];


%%%%%%%%%%%%%%%%%%%%%%% Dispersion diagram

figure
plot(f_th_list, beta_d_ABCD, 'color','blue', 'linewidth', 2)
xlabel ('Frequency [Hz]');
ylabel ('\beta * d [rad]');
xlim([1800 5000]);
ylim([-3.2 3.2])
legend('location', 'eastoutside');
set( legend, 'fontsize', 13 );
hold on
title ('Dispersion diagram');

plot(f_th_list, k_d_list, 'color','k', 'linewidth', 0.5);
plot(f_th_list, -k_d_list, 'color','k', 'linewidth', 0.5);

plot(fl_x, fl_y, 'color','k', 'Marker','.', 'linewidth', 0.1)
text(fcL + 25, 2.5, 'f_{cL}', 'fontsize', 15)

plot(fr_x, fr_y, 'color','k', 'Marker','.', 'linewidth', 0.1)
text(fcR + 25, 2.5, 'f_{cR}', 'fontsize', 15)

plot(f1_x, f1_y, 'color','k', 'Marker','.', 'linewidth', 0.1)
text(f1 + 25, 2.5, 'f_1', 'fontsize', 15)

plot(f2_x, f2_y, 'color','k', 'Marker','.', 'linewidth', 0.1)
text(f2 + 25, 2.5, 'f_2', 'fontsize', 15)

plot(f_zero_x, f_zero_y, 'color','k','Marker','.', 'linewidth', 0.1)
text(f_zero + 25, 2.5, 'f_0', 'fontsize', 15)
set( gca, 'fontsize', 13 );
legend('CRLH', '\pm k*d')

hold off

%%%%%%%%%%%%%%%%%%%%%% Angular response

k_full     = [-k_list(1:offset), k_list(offset + 1:end)];
for it=1:length(k_full)
  theta_list(it) = real(asin(beta_ABCD_CRLH(it)/k_full(it)) * 180/pi);
end


figure
plot(f_th_list(1:offset), -theta_list(1:offset), 'color','red', 'linewidth', 3);
hold on
xlabel ('Frequency [Hz]');
ylabel ('\theta [deg]');
xlim([f1 f2]);
ylim([-90 90])
hold on
title ('Angular response');
plot(f_th_list(offset:end),theta_list(offset:end), 'color','blue','linewidth', 3);

%%%%%%%%%%%% Experimental data

plot(2859, -60, '-s', 'MarkerSize',12);
plot(3046, -30, '^', 'MarkerSize',12);
plot(3234, 0, 'x', 'MarkerSize',12);
plot(3421, 30, 'o', 'MarkerSize',12);
plot(3656, 60, '+', 'MarkerSize',12);


%%%%%%%%%%%%%%%%%%%%%%% HPBW

for k=1:length(k_list)
  lambda_list(k) = (2*pi)/k_list(k);
end

for delta = 1:length(theta_list)
  delta_theta(delta) = [1/( ( L / lambda_list(delta) )  * cos(theta_list(delta)* pi/180))];
end

delta_theta_deg = (delta_theta  * 180/pi);
hpbw = delta_theta_deg/2;

for i = 1:length(theta_list)
  dr(i) = [(1.02*lambda_list(i))/ L ];
end

dr_deg = (dr  * 180/pi);


plot(f_th_list(1:offset), -theta_list(1:offset) + hpbw(1:offset), 'color','black');
plot(f_th_list(1:offset), -theta_list(1:offset) - hpbw(1:offset), 'color','black');
plot(f_th_list(offset:end), theta_list(offset:end) + hpbw(offset:end), 'color','black');
plot(f_th_list(offset:end), theta_list(offset:end) - hpbw(offset:end), 'color','black');


legend('location', 'northwest');
set( legend, 'fontsize', 13 );
grid on
set( gca, 'fontsize', 13 );
legStr = { '\theta LH', '\theta RH', 'Exp:-60^{\circ}', 'Exp:-30^{\circ}', 'Exp:0^{\circ}',  'Exp:30^{\circ}', 'Exp:60^{\circ}', 'HPBW' };
legend( legStr );
legend('location', 'eastoutside');


%%%%%%%%%%%%%%%%%%%% Print data

fprintf('LH cut off frequency fcL = %d Hz\n', int16(fcL))
fprintf('Backfire cut off frequency f1 = %d Hz\n', int16(f1))
fprintf('Transition frequency f_zero = %d Hz\n', int16(f_zero))
fprintf('Endfire cut off frequency f2 = %d Hz\n', int16(f2))
fprintf('RH cut off frequency fcR = %d Hz\n', int16(fcR))
