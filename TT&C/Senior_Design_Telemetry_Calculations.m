%Senior Design Telemetry Calculations
clc
clear

%Designing for Table 16-16 Small User
Ground_diam  = 10;  %m   -Design Choice
Uplink_freq = 14;   %Ghz   -Design Choice
Uplink_wavelength = 2.998*10^3 / (Uplink_freq * 10^9);

Uplink_beamwidth = 21 / Uplink_freq / Ground_diam;   %deg   -Calculation

Ground_antenna_eff = 0.55;  %Efficiency    -Assumption
Uplink_gain = 20.4 + 20*log10(Uplink_freq) + 20*log10(Ground_diam) + 10*log10(Ground_antenna_eff);  %Uplink gain dBi  -calculated

Uplink_transmit_power = 200;    %W  -Input
Loss_Backoff_and_line = 4; %dB(positive value, subtracted later)   -Assume

Uplink_EIRP = 10*log10(Uplink_transmit_power) + Uplink_gain - Loss_Backoff_and_line;    %dB, EIRP of uplink(DSN)  -calculated

Uplink_Range = 38000;   %km     -Design Choice/Mission vaule
Uplink_space_loss = 92.45 + 20*log10(Uplink_Range) + 20*log10(Uplink_freq);
Loss_atmosphere_Up = 10;    %dB (positive)   -unsure of how to find
Net_path_loss_up = Uplink_space_loss + Loss_atmosphere_Up;

SC_antenna_efficiency = 0.7; %efficiency  -assumption
Sat_coverage_area = 13.30;      %deg^2  -unsure of calculation
G_ideal = 10*log10(41253/Sat_coverage_area) + 10*log10(SC_antenna_efficiency);   %db  -calculation
Sat_line_loss = 2;  %db  -assumption
Sat_receive_power = Uplink_EIRP + G_ideal  - Net_path_loss_up - Sat_line_loss;



T_a_est = 290;  %K  -estimate
T_L_est = 290;  %K  -estimate
T_o = 290;  %K  -reference temp, given
L =  0.5;   %Ratio of output power to input power, estimate
F_db = 2.5;   %db   -estimate
F = 10^(F_db/10);   %calculation

T_s_nonDB = T_a_est * L + (1-L)*T_L_est + (F-1)*T_o;
T_s = 10*log10(T_s_nonDB);  %
Noise_bandwidth_B = 3;  %db  -assumption

G_T = G_ideal - T_s;    %db   -calculation



%%
Downlink_freq = 11.7;   %GHz
SC_diam = 1.2;  %m

















