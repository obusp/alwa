%This code is designed to study Acuistic Leaky-Wave Antennas (ALWA), based on the theory of composite right left handed materials (CRLH)
%The model proposed is of cilindrical waveguide with axisymmetric open channels
%Suplemental material for:
%"Educational Open Source Kit for the Evaluation of Acoustic Leaky Wave Antennas with Metamaterials"
%Eduardo Romero-Vivas, Javier Romero-Vivas, Omar Bustamante, Braulio Leon-Lopez
%JASA Eduaction in Acoustics
%%Version 1.1, December 2021, Octave/Matlab
%
%f_beg - first frequency of study
%f_end - last frequency of study
%res - number of samples between frequencies interval
%
%rho - air density
%c - speed of sound
%
%ALWA parameters
%a - waveguide radius
%b - shunt width
%L - shunt length
%h - membrane thickness
%N - unit cell number
%d - unit cell length
%l - total ALWA length

close all
clear all

%%%%%%%%%%%%%%%%%%%%%%% Frequencies of study
f_beg = 1700;
f_end = 5000;
res   = 1000;

f_th_list = linspace(f_beg,f_end,res);

%%%%%%%%%%%%%%%%%%%%%%% Characteristics of the medium

rho    = 0.9402;                   %% air density STP reference - La Paz BCS, Mexico
c      = 342.4;                    %% speed of sound

%%%%%%%%%%%%%%%%%%%%%%% LWA parameters

a      = 0.0039;          %% waveguide radius
b      = 0.0004;          %% shunt width
L      = 0.0198;          %% shunt length
h      = 0.000067;        %% membrane thickness

S_a    = pi * (a^2);         %% waveguide transversal area

N = 20;           %% unit cell num
d = 0.0124;       %% unit cell length
l = N * d;        %% total ALWA length


%%%%%%%%%%%%%%%%%%%%%%% Elements of the impedance Z_th

mass_tl = (rho/S_a) * (d-h);          % mass of the waveguide section
mass_m  = 1.8830 * ((1370*h)/(pi * a^2));              % mass of the membrane
c_m     = (pi * a^6 ) / (196.51 * ( (3e9 * h^3)/(12*(1-0.33^2) ) ) );       % compliance of the membrane


%%%%%%%%%%%%%%%%%%%%%%% Elements of the admittance Y_th

c_tl       = (S_a/(rho * c^2)) * (d-h);       % compliance of the waveguide section

f_zero         = 1 /(2*pi * sqrt( (mass_tl + mass_m) * c_m ) )   % resonance frequency of Z_th

mass_slit  = (rho/(2*pi*b)) * log(1 + L/a);            % shunt mass of the slit
c_slit     = ( 1 / ( 4 * pi^2 * (f_zero)^2 * mass_slit) ) - c_tl;  % shunt compliance of the slit


%%%%%%%%%%%%%%%%%%%%%%% Cutoff frequencies

w_R        = 1 /sqrt( (mass_tl + mass_m) * (c_slit + c_tl) );
w_cut_R    = 2 * w_R;
f_cut_R    = w_cut_R / (2*pi);

w_L        = 1 / sqrt( mass_slit * c_m );
w_cut_L    = w_L / 2;
f_cut_L    = w_cut_L / (2*pi);

w_shunt    = 1 /sqrt( mass_slit * (c_slit + c_tl) );  % resonance frequency of Y_th
f_shunt    = w_shunt / (2*pi);

w_series   = 1 /sqrt( (mass_tl + mass_m) * c_m );
f_series   = w_series / (2*pi);

f_L        = ( w_L/(2*pi) );
f_R        = ( w_R/(2*pi) );

fcL        = f_R * abs( 1 - (sqrt( 1 + ( f_L / f_R ) )) );
fcR        = f_R   *  ( 1 + (sqrt( 1 + ( f_L / f_R ) )) );


%%%%%%%%%%%%%%%%%%%%%%% Transfer function

for num=1:length(f_th_list)
  for f=f_th_list(num)
    omega = 2 * pi * f;
    k     = omega/343.4;


    %%%%%%%%%%%%%%%%%%%%%%% Series impedance Z_th
    Z_th = 1i * ( (omega*(mass_tl + mass_m)) - (1/(omega*c_m)));

    %%%%%%%%%%%%%%%%%%%%%%% Shunt admittance Y_th
    Y_th = 1i * ( (omega*(c_slit + c_tl)) - (1/(omega*mass_slit))  );


    %%%%%%%%%%%%%%%%%%%%%%% ABCD matrix
    %   Y/2 --- Z --- Y/2 cell configuration

    A = 1 + ((Z_th * Y_th) / 2);
    B = Z_th;
    C = Y_th * (1 + ((Z_th * Y_th) / 4));
    D = 1 + ((Z_th * Y_th) / 2);

    ABCD = [A, B; C, D];
    ABCD_N = ABCD^N;

    A_N = ABCD_N(1);
    B_N = ABCD_N(3);
    C_N = ABCD_N(2);
    D_N = ABCD_N(4);

    %%%%%%%%%%%%%%%%%%%%%%% Bloch constant (Bongard, 2010)
    gamma_bloch = (acosh(A)) / d;
    Z_bloch     = sqrt(B/C);
    gamma_bloch_list(num)=gamma_bloch;
    Z_bloch_list(num)=Z_bloch;


    %%%%%%%%%%%%%%%%%%%%%%% Scattering parameters
    bulk_m = rho * c^2;        % bulk modulus
    Zc     = sqrt(mass_tl/(c_slit + c_tl));     % impedance

    S_11   = (A_N + (B_N/Zc) - (C_N*Zc) - D_N ) / (A_N + (B_N/Zc) + (C_N*Zc) + D_N );
    S_21   = 2 / (A_N + (B_N/Zc) + (C_N*Zc) + D_N );

    S_11_list(num)=S_11;
    S_21_list(num)=S_21;

  end
end


%%%%%%%%%%%%%%%%%%%%%%% Offset freq

aRef = f_zero;                 % Reference to set the phase origin whether its zero or the closest value
aDiff = abs(f_th_list - aRef); % Calculate the difference
[minVal, offset] = min(aDiff); % Find closest value to aRef


%%%%%%%%%%%%%%%%%%%%%%% Alpha - attenuator factor
alpha_list = real(gamma_bloch_list);

%%%%%%%%%%%%%%%%%%%%%%% Beta - bloch constant through ABCD matrix
beta_list        = imag(gamma_bloch_list);
beta_LH          = -(beta_list(1:offset));
beta_RH          = beta_list(offset+1:end);
beta_ABCD_CRLH   = [beta_LH, beta_RH];

for i=1:length(beta_list)
  beta_d_ABCD(i) = d*beta_ABCD_CRLH(i);
end




%%%%%%%%%%%%%%%%%%%%%%%
for ps21=1:length(S_21_list)
  phi_S21_list(ps21)=angle(S_21_list(ps21));
end

%%%%%%%%%%%%%%%%%%%%%%% Unwrap

phi_S21_lh_list = flip(phi_S21_list(1:offset));

phi_S21_lh=unwrap(phi_S21_lh_list);
phi_S21_lh=flip(phi_S21_lh);
phi_S21_rh=unwrap(phi_S21_list(offset:end));


%%%%%%%%%%%%%%%%%%%%%%% Beta through scattering parameters
phi_S21_LH = -phi_S21_lh/l;
phi_S21_RH = -phi_S21_rh/l;

%%%%%%%%%%%%%%%%%%%%%%% Wavenumber k = omega/c
for num=1:length(f_th_list)
  k_list(num) = ( (2*pi*f_th_list(num))/c );
end

k_d_list = k_list * d;
k_L =-1 * k_list(1:offset);
k_R = k_list(offset+1:end);


%%%%%%%%%%%%%%%%%%%%%%% Cutoff frequencies
fl_x = [fcL fcL];
fl_y = [-4 4];

fr_x = [fcR fcR];
fr_y = [-4 4];

f1_list = abs(phi_S21_LH-k_L) < 1;
f1_idx = find (f1_list);

f2_list = abs(phi_S21_RH(2:end)-k_R) < 1;
f2_idx = find (f2_list);

f1   = f_th_list(f1_idx(1))
f2   = f_th_list(offset + f2_idx(1))

f1_x = [f1 f1];
f1_y = [-4 4];

f2_x = [f2 f2];
f2_y = [-4 4];

f_zero_x = [f_zero f_zero];
f_zero_y = [-4 4];

%%test

figure
%plot(f_th_list, beta_d_ABCD, 'color','blue', 'linewidth', 2, ';CRLH;')
%%%%AQUI
plot(f_th_list, beta_d_ABCD, 'color','blue', 'linewidth', 2)
xlabel ('Frequency [Hz]');
ylabel ('\Delta \phi_{S21} = \beta * d [rad]');
xlim([1800 5000]);
ylim([-3.2 3.2])
legend('location', 'eastoutside');
set( legend, 'fontsize', 13 );
hold on
title ('Dispersion diagram');
%%plot(f_th_list, k_d_list, 'color','k', 'linewidth', 0.5, ';k*d;');
%%%octave

plot(f_th_list, k_d_list, 'color','k', 'linewidth', 0.5);
plot(f_th_list, -k_d_list, 'color','k', 'linewidth', 0.5);

%plot(fl_x, fl_y, 'color','k', '--', 'linewidth', 0.1)
plot(fl_x, fl_y, 'color','k', 'Marker','.', 'linewidth', 0.1)
text(fcL + 25, 2.5, 'fcL', 'fontsize', 15)

%plot(fr_x, fr_y, 'color','k', '--', 'linewidth', 0.1)
plot(fr_x, fr_y, 'color','k', 'Marker','.', 'linewidth', 0.1)
text(fcR + 25, 2.5, 'fcR', 'fontsize', 15)

%plot(f1_x, f1_y, 'color','k', '--', 'linewidth', 0.1)
plot(f1_x, f1_y, 'color','k', 'Marker','.', 'linewidth', 0.1)
text(f1 + 25, 2.5, 'f_1', 'fontsize', 15)

%plot(f2_x, f2_y, 'color','k', '--', 'linewidth', 0.1)
plot(f2_x, f2_y, 'color','k', 'Marker','.', 'linewidth', 0.1)
text(f2 + 25, 2.5, 'f_2', 'fontsize', 15)

%plot(f_zero_x, f_zero_y, 'color','k', '--', 'linewidth', 0.1)
plot(f_zero_x, f_zero_y, 'color','k','Marker','.', 'linewidth', 0.1)
text(f_zero + 25, 2.5, 'f_0', 'fontsize', 15)
set( gca, 'fontsize', 13 );
legend('CRLH', 'k*d')

hold off

%%%%%%%%%%%%%%%%%%%%%%% Theta

k_full     = [-k_list(1:offset), k_list(offset + 1:end)];
for it=1:length(k_full)
  theta_list(it) = (asin(beta_ABCD_CRLH(it)/k_full(it)) * 180/pi );
end


figure
%plot(f_th_list(1:offset), -theta_list(1:offset), 'color','red', 'linewidth', 3, ';\theta CRLH;');
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

%legStr = { '\theta LH', '\theta RH', 'Exp:-60ª', 'Exp:-30ª', 'Exp:0ª',  'Exp:30ª', 'Exp:60ª' };
%legend( legStr );

%%%%%%%%%%%%%%%%%%%%%%% HPBW

for k=1:length(k_list)
  lambda_list(k) = (2*pi)/k_list(k);
end

for delta = 1:length(theta_list)
  delta_theta(delta) = [1/( ( l / lambda_list(delta) )  * cos(theta_list(delta)*pi/180) )];
end

delta_theta_deg = delta_theta * 180/pi;
hpbw = delta_theta_deg/2;

%plot(f_th_list(1:offset), -theta_list(1:offset) + hpbw(1:offset), 'color','black', ';HPBW;');
%AQUI
plot(f_th_list(1:offset), -theta_list(1:offset) + hpbw(1:offset), 'color','black');
plot(f_th_list(1:offset), -theta_list(1:offset) - hpbw(1:offset), 'color','black');
plot(f_th_list(offset:end), theta_list(offset:end) + hpbw(offset:end), 'color','black');
plot(f_th_list(offset:end), theta_list(offset:end) - hpbw(offset:end), 'color','black');


legend('location', 'northwest');
set( legend, 'fontsize', 13 );
grid on
set( gca, 'fontsize', 13 );
legStr = { '\theta LH', '\theta RH', 'Exp:-60ª', 'Exp:-30ª', 'Exp:0ª',  'Exp:30ª', 'Exp:60ª', 'HPBW' };
legend( legStr );