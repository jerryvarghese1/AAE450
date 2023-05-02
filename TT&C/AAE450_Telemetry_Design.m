%% Sizing Mass versus Antenna Diameter
clc
clear
close all

antenna_diam = [2.5, 4, 3.66, 4.8, 3.1];
antenna_mass = [21.3, 100.6, 53, 255.55, 36.24];
mass_log = log10(antenna_mass);

p = polyfit(antenna_diam, mass_log, 1)
x_range = 0:0.1:6;
y_range = p(1) * x_range + p(2);
max_weight = 106.913;   %Max weight of all telecommunications systems
assumption_component_weight = 80;
antenna_mass = max_weight - assumption_component_weight
dish_percentage = 100/118; %Based on previous mission architecture for percentage of dish to overall telecommunicaiton mass
weight_limit = log10(300) * ones(1, length(x_range))%log10(max_weight*dish_percentage) * ones(1, length(x_range));

figure(1)
hold on
plot(antenna_diam, mass_log, 'kx')
plot(x_range, y_range, 'b')
plot(x_range, weight_limit, 'r')
axis([1, 6, 0, 3.5])
grid on
title('Relating mass of antenna vs S/C Antenna Diameter')
xlabel('Antenna Diameter (m)')
ylabel('Antenna Mass (log(kg))')
legend('Individual Data Points', 'Sizing Estimate', 'Antenna Weight Limit')

%% Antenna plots
% clc
% clear
% close all
%Equation: [LM_uplink,  LM_downlink] = LinkBudgetCalculation(sc_diam, sc_power, sc_freq, sc_Data_rate, DSN_diam, DSN_power, DSN_freq, DSN_Data_rate)
%LM_uplink: Link Margin (dB) of Telemetry Uplink(DSN-->S/C) 
%LM_downlink: Link Margin (dB) of Telemetry Downlink(S/C-->DSN)
%sc_diam: Diameter of S/C Antenna, m
%sc_power: Power of S/C Antenna, W
%sc_freq: Sending Frequency of S/C Antenna during Downlink, GHz
%sc_Data_rate: Data rate of S/C Antenna during Downlink, bits/sec
%DSN_diam: Diameter of DSN, m --> Either 34m(X, Ka) or 70m (X)
%DSN_power: Power of DSN Antenna, kW
%DSN_freq: Sending Frequency of DSN during Uplink, GHz
%DSN_Data_rate: Data rate of DSN Antenna during Uplink, bits/sec

%Both Link Margins need to be greater than zero for link to be closed but
%need to account for margin of error
%All spacecraft values are chosen
%DSN Values have certain combinations depending on which antenna is used
%34m X Band: 8.45-8.5 GHz, 20 kW, 500-32000 bit/sec data rate
%34m Ka Band: 31.8–32.3 GHz, 0.8 kW, 500-32000 bit/sec data rate
%70m X Band: 8.45-8.5 GHz, 20 kW, 500-32000 bit/sec data rate
%70m X Band High Power:  8.45-8.5 GHz, 400 kW, 500-32000 bit/sec data rate


%% Plotting SC Transmit Power vs SC Data Rate at 3dB Margin and unspecified Error
clc
clear all
close all
SC_Antenna_1 = 3.5; %m
SC_Antenna_2 = 3.6; %m
SC_Antenna_3 = 3.7; %m
SC_Antenna_4 = 3.8; %m
SC_Antenna_5 = 3.9; %m

SC_power = 100:175;
Margin = 3; %Margin of Link Margin
ML_min = 0;  %Minimum acceptable Link Margin (dB)
[C_N0_uplink_1,  C_N0_downlink_1] = LinkBudget_Varying_Datarate(SC_Antenna_1, SC_power, 34.7, 34, 20, 8.5);
[C_N0_uplink_2,  C_N0_downlink_2] = LinkBudget_Varying_Datarate(SC_Antenna_2, SC_power, 34.7, 34, 20, 8.5);
[C_N0_uplink_3,  C_N0_downlink_3] = LinkBudget_Varying_Datarate(SC_Antenna_3, SC_power, 34.7, 34, 20, 8.5);
[C_N0_uplink_4,  C_N0_downlink_4] = LinkBudget_Varying_Datarate(SC_Antenna_4, SC_power, 34.7, 34, 20, 8.5);
[C_N0_uplink_5,  C_N0_downlink_5] = LinkBudget_Varying_Datarate(SC_Antenna_5, SC_power, 34.7, 34, 20, 8.5);

%Because we vary only power of spacecraft antenna, only the downlink data
%rate will change and uplink data rate is the maximum value 
R_uplink_1 = 10.^(1/10 * (C_N0_uplink_1 - Margin - ML_min));   %DSN Bit rate bits/sec
R_downlink_1 = 10.^(1/10 * (C_N0_downlink_1 - Margin - ML_min));    %DSN Bit rate bits/sec
R_uplink_2 = 10.^(1/10 * (C_N0_uplink_2 - Margin - ML_min));   %DSN Bit rate bits/sec
R_downlink_2 = 10.^(1/10 * (C_N0_downlink_2 - Margin - ML_min));    %DSN Bit rate bits/sec
R_uplink_3 = 10.^(1/10 * (C_N0_uplink_3 - Margin - ML_min));   %DSN Bit rate bits/sec
R_downlink_3 = 10.^(1/10 * (C_N0_downlink_3 - Margin - ML_min));    %DSN Bit rate bits/sec
R_uplink_4 = 10.^(1/10 * (C_N0_uplink_4 - Margin - ML_min));   %DSN Bit rate bits/sec
R_downlink_4 = 10.^(1/10 * (C_N0_downlink_4 - Margin - ML_min));    %DSN Bit rate bits/sec
R_uplink_5 = 10.^(1/10 * (C_N0_uplink_5 - Margin - ML_min));   %DSN Bit rate bits/sec
R_downlink_5 = 10.^(1/10 * (C_N0_downlink_5 - Margin - ML_min));    %DSN Bit rate bits/sec

figure(1)
hold on
% plot(SC_power, R_uplink_1*ones(1, length(SC_power)), 'r--')
plot(SC_power, R_downlink_1, 'r')
% plot(SC_power, R_uplink_2*ones(1, length(SC_power)), 'b--')
plot(SC_power, R_downlink_2, 'b')
% plot(SC_power, R_uplink_3*ones(1, length(SC_power)), 'g--')
plot(SC_power, R_downlink_3, 'g')
% plot(SC_power, R_uplink_4*ones(1, length(SC_power)), 'k--')
plot(SC_power, R_downlink_4, 'k')
% plot(SC_power, R_uplink_5*ones(1, length(SC_power)), 'c--')
plot(SC_power, R_downlink_5, 'c')
grid on
title({'Determining Data Rate for Telemetry while varying S/C Transmit Power', '3 dB Margin, Link Margin = 1 dB, DSN 34m, 34.7GHz downlink, 8.5 GHz Uplink'})
xlabel('SC Antenna Power (W)')
ylabel('Data Rate (bits/sec)')
% legend({'DSN Max Uplink Data Rate (bits/sec), SC diam = 2m','Spacecraft Max Downlink Data Rate (bits/sec), SC diam = 2m', 'DSN Max Uplink Data Rate (bits/sec), SC diam = 2.4m','Spacecraft Max Downlink Data Rate (bits/sec), SC diam = 2.4m','DSN Max Uplink Data Rate (bits/sec), SC diam = 2.8m','Spacecraft Max Downlink Data Rate (bits/sec), SC diam = 2.8m','DSN Max Uplink Data Rate (bits/sec), SC diam = 3.2m','Spacecraft Max Downlink Data Rate (bits/sec), SC diam = 3.2m','DSN Max Uplink Data Rate (bits/sec), SC diam = 3.6m','Spacecraft Max Downlink Data Rate (bits/sec), SC diam = 3.6m'}, 'location', 'northwest')
legend({'Spacecraft Max Downlink Data Rate (bits/sec), SC diam = 2m','Spacecraft Max Downlink Data Rate (bits/sec), SC diam = 2.4m','Spacecraft Max Downlink Data Rate (bits/sec), SC diam = 2.8m','Spacecraft Max Downlink Data Rate (bits/sec), SC diam = 3.2m','Spacecraft Max Downlink Data Rate (bits/sec), SC diam = 3.6m'}, 'location', 'northwest')

%% 
clc
clear
close all
sc_diam = 3.5:0.01:4;    %m
sc_power_1 = 100; %W
sc_power_2 = 110; %W
sc_power_3 = 120; %W
sc_power_4 = 130; %W
sc_power_5 = 140; %W
sc_power_6 = 150; %W
sc_power_7 = 160; %W

sc_freq = 34.7; %Ghz
sc_Data_rate = 5000; %bits/sec
DSN_diam = 34;
DSN_power = 20; %KW
DSN_freq = 8.5;
DSN_Data_rate = 10000;

[LM_uplink_1,  LM_downlink_1] = LinkBudgetCalculation(sc_diam, sc_power_1, sc_freq, sc_Data_rate, DSN_diam, DSN_power, DSN_freq, DSN_Data_rate);
[LM_uplink_2,  LM_downlink_2] = LinkBudgetCalculation(sc_diam, sc_power_2, sc_freq, sc_Data_rate, DSN_diam, DSN_power, DSN_freq, DSN_Data_rate);
[LM_uplink_3,  LM_downlink_3] = LinkBudgetCalculation(sc_diam, sc_power_3, sc_freq, sc_Data_rate, DSN_diam, DSN_power, DSN_freq, DSN_Data_rate);
[LM_uplink_4,  LM_downlink_4] = LinkBudgetCalculation(sc_diam, sc_power_4, sc_freq, sc_Data_rate, DSN_diam, DSN_power, DSN_freq, DSN_Data_rate);
[LM_uplink_5,  LM_downlink_5] = LinkBudgetCalculation(sc_diam, sc_power_5, sc_freq, sc_Data_rate, DSN_diam, DSN_power, DSN_freq, DSN_Data_rate);
[LM_uplink_6,  LM_downlink_6] = LinkBudgetCalculation(sc_diam, sc_power_6, sc_freq, sc_Data_rate, DSN_diam, DSN_power, DSN_freq, DSN_Data_rate);
[LM_uplink_7,  LM_downlink_7] = LinkBudgetCalculation(sc_diam, sc_power_7, sc_freq, sc_Data_rate, DSN_diam, DSN_power, DSN_freq, DSN_Data_rate);

figure(2)
hold on
grid on
plot(sc_diam, LM_uplink_1)
plot(sc_diam, LM_uplink_2)
plot(sc_diam, LM_uplink_3)
plot(sc_diam, LM_uplink_4)
plot(sc_diam, LM_uplink_5)
plot(sc_diam, LM_uplink_6)
plot(sc_diam, LM_uplink_7)
title({'Telemecommunications Uplink Link Margin varying SC Diameter and SC Power', 'SC Freq: 8.5 GHz, SC R: 325 b/s, DSN 34m 20kW, , DSN Freq: 8.5 GHz, DSN R: 10000 b/s'})
xlabel('SC Antenna Diameter (m)')
ylabel('Uplink Link Margin (dB)')
legend('SC Power = 100 W', 'SC Power = 110 W', 'SC Power = 120 W', 'SC Power = 130 W', 'SC Power = 140 W', 'SC Power = 150 W','SC Power = 160 W')  

figure(3)
hold on
grid on
plot(sc_diam, LM_downlink_1)
plot(sc_diam, LM_downlink_2)
plot(sc_diam, LM_downlink_3)
plot(sc_diam, LM_downlink_4)
plot(sc_diam, LM_downlink_5)
plot(sc_diam, LM_downlink_6)
plot(sc_diam, LM_downlink_7)
title({'Telemecommunications Downlink Link Margin varying SC Diameter and SC Power', 'SC Freq: 8.5 GHz, SC R: 50 b/s, DSN 34m 20kW, , DSN Freq: 8.5 GHz, DSN R: 10000 b/s'})
xlabel('SC Antenna Diameter (m)')
ylabel('Downlink Link Margin (dB)')
legend('SC Power = 100 W', 'SC Power = 110 W', 'SC Power = 120 W', 'SC Power = 130 W', 'SC Power = 140 W', 'SC Power = 150 W','SC Power = 160 W')  

%% Bit Error rate for individual uplink and Donwnlink
clc
clear
close all

sc_diam = 3.8;  %kg
sc_power = 110;  %W
sc_freq = 34.7; %GHz
sc_Data_rate = 4000;    %bits/sec'
[EB_NO_uplink,  EB_NO_downkink] = LinkBudget_BER(sc_diam, sc_power, sc_freq, sc_Data_rate, 34, 20, 8.5, 10000);



%Data Modulation Options:
%Type: QPSK  Bits per Symbol: 2  Code Rate: 1/2  Eb/No_req: 1.05
%Type: QPSK  Bits per Symbol: 2  Code Rate: 3/4  Eb/No_req: 2.31
%Type: QPSK  Bits per Symbol: 2  Code Rate: 5/9  Eb/No_req: 2.99
%Type: QPSK  Bits per Symbol: 2  Code Rate: 9/10  Eb/No_req: 3.89
%Type: 8PSK  Bits per Symbol: 3  Code Rate: 5/6  Eb/No_req: 5.41
%Type: 8PSK  Bits per Symbol: 3  Code Rate: 9/10  Eb/No_req: 6.7
%Type: 16PSK  Bits per Symbol: 4  Code Rate: 5/6  Eb/No_req: 6.42
%Type: 16PSK  Bits per Symbol: 4  Code Rate: 9/10  Eb/No_req: 7.61

%Modulation Options for Data. Type_Code Rate = [Bits per symbol, Eb/No_req]
QPSK_1_2 = [2, 1.05];
QPSK_3_4 = [2, 2.31];
QPSK_5_9 = [2, 2.99];
QPSK_9_10 = [2, 3.89];
PSK8_5_6 = [3, 5.41];
PSK8_9_10 = [3, 6.7];
PSK16_5_6 = [4, 6.42];
PSK16_9_10 = [4, 7.61];

%Can be any one of the options above
Modulation_choice = QPSK_5_9;

m_1 = Modulation_choice(1);
Eb_No_req_1 = Modulation_choice(2);

BER_up_1 = 1/m_1 * erfc(sqrt(m_1 * EB_NO_uplink)*sin(pi/(2^m_1)));
BER_down_1 = 1/m_1 * erfc(sqrt(m_1 * (EB_NO_downkink))*sin(pi/(2^m_1)));
LM_up_1 = EB_NO_uplink - Eb_No_req_1;
LM_down_1 = EB_NO_downkink - Eb_No_req_1;


%BER for Uplink and Downlink
figure(1)
hold on
plot(sc_power, ones(1, length(sc_power))*(BER_up_1))
plot(sc_power, (BER_down_1))
xlabel('SC Transmit Power (W)')
ylabel('Bit Error Rate (BER)')
title('Comparing Bit Error Rate to SC Transmit Power for various Modulations')
legend('Uplink', 'Downlink')
grid on
set(gca, 'YScale', 'log')

% %Link Margins
figure(2)
hold on
plot(sc_power, ones(1, length(sc_power))*LM_up_1)
plot(sc_power, LM_down_1)
grid on
xlabel('SC Transmit Power (W)')
ylabel('Link Margin (db)')
title('Comparing Link Margin to SC Transmit Power for various Modulations')
legend('Uplink', 'Downlink')

%% Comparing Bit Error Rates and ML for Uplink only and Downlink Only
clc
clear
close all

sc_diam = 3.8;  %m
sc_power = 100:160;  %W
sc_freq = 34.7; %GHz
sc_Data_rate = 4000;    %bits/sec'

[EB_NO_uplink,  EB_NO_downkink] = LinkBudget_BER(sc_diam, sc_power, sc_freq, sc_Data_rate, 34, 20, 8.5, 10000);

%Data Modulation Options:
%Type: QPSK  Bits per Symbol: 2  Code Rate: 1/2  Eb/No_req: 1.05
%Type: QPSK  Bits per Symbol: 2  Code Rate: 3/4  Eb/No_req: 2.31
%Type: QPSK  Bits per Symbol: 2  Code Rate: 5/9  Eb/No_req: 2.99
%Type: QPSK  Bits per Symbol: 2  Code Rate: 9/10  Eb/No_req: 3.89
%Type: 8PSK  Bits per Symbol: 3  Code Rate: 5/6  Eb/No_req: 5.41
%Type: 8PSK  Bits per Symbol: 3  Code Rate: 9/10  Eb/No_req: 6.7
%Type: 16PSK  Bits per Symbol: 4  Code Rate: 5/6  Eb/No_req: 6.42
%Type: 16PSK  Bits per Symbol: 4  Code Rate: 9/10  Eb/No_req: 7.61

%Modulation Options for Data. Type_Code Rate = [Bits per symbol, Eb/No_req]
QPSK_1_2 = [2, 1.05];
QPSK_3_4 = [2, 2.31];
QPSK_5_9 = [2, 2.99];
QPSK_9_10 = [2, 3.89];
PSK8_5_6 = [3, 5.41];
PSK8_9_10 = [3, 6.7];
PSK16_5_6 = [4, 6.42];
PSK16_9_10 = [4, 7.61];

%Can be any one of the options above
Modulation_choice_1 = QPSK_1_2;
Modulation_choice_2 = QPSK_3_4;
Modulation_choice_3 = QPSK_5_9;
Modulation_choice_4 = QPSK_9_10;
Modulation_choice_5 = PSK8_5_6;
Modulation_choice_6 = PSK8_9_10;
Modulation_choice_7 = PSK16_5_6;
Modulation_choice_8 = PSK16_9_10;

% EB_NO_downkink = 10.^(EB_NO_downkink/10);
m_1 = Modulation_choice_1(1);
Eb_No_req_1 = Modulation_choice_1(2);
BER_up_1 = 1/m_1 * erfc(sqrt(m_1 * EB_NO_uplink)*sin(pi/(2^m_1)));
BER_down_1 = 1/m_1 * erfc(sqrt(m_1 * (EB_NO_downkink))*sin(pi/(2^m_1)));
LM_up_1 = EB_NO_uplink - Eb_No_req_1;
LM_down_1 = EB_NO_downkink - Eb_No_req_1;

m_2 = Modulation_choice_2(1);
Eb_No_req_2 = Modulation_choice_2(2);
BER_up_2 = 1/m_2 * erfc(sqrt(m_2 * EB_NO_uplink)*sin(pi/(2^m_2)));
BER_down_2 = 1/m_2 * erfc(sqrt(m_2 * EB_NO_downkink)*sin(pi/(2^m_2)));
LM_up_2 = EB_NO_uplink - Eb_No_req_2;
LM_down_2 = EB_NO_downkink - Eb_No_req_2;

m_3 = Modulation_choice_3(1);
Eb_No_req_3 = Modulation_choice_3(2);
BER_up_3 = 1/m_3 * erfc(sqrt(m_3 * EB_NO_uplink)*sin(pi/(2^m_3)));
BER_down_3 = 1/m_3 * erfc(sqrt(m_3 * EB_NO_downkink)*sin(pi/(2^m_3)));
LM_up_3 = EB_NO_uplink - Eb_No_req_3;
LM_down_3 = EB_NO_downkink - Eb_No_req_3;

m_4 = Modulation_choice_4(1);
Eb_No_req_4 = Modulation_choice_4(2);
BER_up_4 = 1/m_4 * erfc(sqrt(m_4 * EB_NO_uplink)*sin(pi/(2^m_4)));
BER_down_4 = 1/m_4 * erfc(sqrt(m_4 * EB_NO_downkink)*sin(pi/(2^m_4)));
LM_up_4 = EB_NO_uplink - Eb_No_req_4;
LM_down_4 = EB_NO_downkink - Eb_No_req_4;

m_5 = Modulation_choice_5(1);
Eb_No_req_5 = Modulation_choice_5(2);
BER_up_5 = 1/m_5 * erfc(sqrt(m_5 * EB_NO_uplink)*sin(pi/(2^m_5)));
BER_down_5 = 1/m_5 * erfc(sqrt(m_5 * EB_NO_downkink)*sin(pi/(2^m_5)));
LM_up_5 = EB_NO_uplink - Eb_No_req_5;
LM_down_5 = EB_NO_downkink - Eb_No_req_5;

m_6 = Modulation_choice_6(1);
Eb_No_req_6 = Modulation_choice_6(2);
BER_up_6 = 1/m_6 * erfc(sqrt(m_6 * EB_NO_uplink)*sin(pi/(2^m_6)));
BER_down_6 = 1/m_6 * erfc(sqrt(m_6 * EB_NO_downkink)*sin(pi/(2^m_6)));
LM_up_6 = EB_NO_uplink - Eb_No_req_6;
LM_down_6 = EB_NO_downkink - Eb_No_req_6;

m_7 = Modulation_choice_7(1);
Eb_No_req_7 = Modulation_choice_7(2);
BER_up_7 = 1/m_7 * erfc(sqrt(m_7 * EB_NO_uplink)*sin(pi/(2^m_7)));
BER_down_7 = 1/m_7 * erfc(sqrt(m_7 * EB_NO_downkink)*sin(pi/(2^m_7)));
LM_up_7 = EB_NO_uplink - Eb_No_req_7;
LM_down_7 = EB_NO_downkink - Eb_No_req_7;

m_8 = Modulation_choice_8(1);
Eb_No_req_8 = Modulation_choice_8(2);
BER_up_8 = 1/m_8 * erfc(sqrt(m_8 * EB_NO_uplink)*sin(pi/(2^m_8)));
BER_down_8 = 1/m_8 * erfc(sqrt(m_8 * EB_NO_downkink)*sin(pi/(2^m_8)));
LM_up_8 = EB_NO_uplink - Eb_No_req_8;
LM_down_8 = EB_NO_downkink - Eb_No_req_8;

%BER for Various Downlinks
figure(1)
loglog(sc_power, (BER_down_1))
hold on
loglog(sc_power, (BER_down_5))
loglog(sc_power, (BER_down_7))
xlabel('SC Transmit Power (W)')
ylabel('Bit Error Rate (BER)')
title({'Comparing Downlink Bit Error Rate to SC Transmit Power for various Modulations', '3.2m, 34.7 GHz, 2000 b/s to 34m DSN'})
legend('QPSK', '8-PSK', '16-PSK', 'location', 'best')
grid on

% %Link Margins for Downlinks
figure(2)
hold on
% plot(sc_power, ones(1, length(sc_power))*LM_up_1)
plot(sc_power, LM_down_1, 'r:')
plot(sc_power, LM_down_2, 'r-.')
plot(sc_power, LM_down_3, 'r--')
plot(sc_power, LM_down_4, 'r')
plot(sc_power, LM_down_5, 'b--')
plot(sc_power, LM_down_6, 'b')
plot(sc_power, LM_down_7, 'g--')
plot(sc_power, LM_down_8, 'g')
grid on
xlabel('SC Transmit Power (W)')
ylabel('Link Margin (db)')
title({'Comparing Downlink Link Margin to SC Transmit Power for various Modulations', '3.8m, 34.7 GHz, 4000 b/s to 34m DSN'})
legend('QPSK 1/2', 'QPSK 3/4', 'QPSK 5/9', 'QPSK 9/10', '8-PSK 5/6', '8-PSK 9/10', '16-PSK 5/6', '16-PSK 9/10', 'location', 'best')

%BER for Various Uplinks
figure(3)
hold on
plot(sc_power, ones(1, length(sc_power))*(BER_up_1))
plot(sc_power, ones(1, length(sc_power))*(BER_up_5))
plot(sc_power, ones(1, length(sc_power))*(BER_up_7))
xlabel('SC Transmit Power (W)')
ylabel('Bit Error Rate (BER)')
title({'Comparing Uplink Bit Error Rate to SC Transmit Power for various Modulations', '34m DSN, 8.5 GHz, 10000b/s to 3.2m, 34.7 GHz, 2000 b/s to 34m DSN'})
legend('QPSK', '8-PSK', '16-PSK', 'location', 'best')
grid on
set(gca, 'YScale', 'log')

% %Link Margins for Uplink
figure(4)
hold on
plot(sc_power, ones(1, length(sc_power))*LM_up_1)
plot(sc_power, ones(1, length(sc_power))*LM_up_1, 'r:')
plot(sc_power, ones(1, length(sc_power))*LM_up_2, 'r-.')
plot(sc_power, ones(1, length(sc_power))*LM_up_3, 'r--')
plot(sc_power, ones(1, length(sc_power))*LM_up_4, 'r')
plot(sc_power, ones(1, length(sc_power))*LM_up_5, 'b--')
plot(sc_power, ones(1, length(sc_power))*LM_up_6, 'b')
plot(sc_power, ones(1, length(sc_power))*LM_up_7, 'g--')
plot(sc_power, ones(1, length(sc_power))*LM_up_8, 'g')
grid on
xlabel('SC Transmit Power (W)')
ylabel('Link Margin (db)')
title({'Comparing Downlink Link Margin to SC Transmit Power for various Modulations', '3.8m, 34.7 GHz, 4000 b/s to 34m DSN'})
legend('QPSK 1/2', 'QPSK 3/4', 'QPSK 5/9', 'QPSK 9/10', '8-PSK 5/6', '8-PSK 9/10', '16-PSK 5/6', '16-PSK 9/10', 'location', 'best')

%% Link Margin vs location
clc
clear
close all
range = 0.1:0.1:21;

HGA_diam = 3.8;
N_MGA = 70;
LGA_max_power = 150;
sc_power = 110;
Ka = 34.7;
X = 8.45;
LGA_R = 1;
MGA_R = 6;
HGA_R = 4000;

[LGA_up_34,  LGA_down_34] = LGA_range(LGA_max_power, X, LGA_R, 34, 20, 8.5, 500);
[LGA_up_70,  LGA_down_70] = LGA_range(LGA_max_power, X, LGA_R, 70, 20, 8.5, 500);
% [LGA_up_70_high,  LGA_down_70_high] = LGA_range(LGA_max_power, X, LGA_R, 70, 400, 8.5, 500);
% [LGA_up_70_34,  LGA_down_70_34] = LGA_range(LGA_max_power, X, LGA_R, 70+34, 20+20, 8.5, 500);
% [LGA_up_70_34_34,  LGA_down_70_34_34] = LGA_range(LGA_max_power, X, LGA_R, 70+2*34, 20+2*20, 8.5, 500);
% [LGA_up_70_high_34_34,  LGA_down_70_high_34_34] = LGA_range(LGA_max_power, X, LGA_R, 70+2*34, 400+2*20, 8.5, 500);
% [LGA_up_Worst_case,  LGA_down_Worst_case] = LGA_range(LGA_max_power, X, LGA_R, 2*(70+2*34), 2*(400+2*20), 8.5, 500);
[MGA_up,  MGA_down] = MGA_range(N_MGA, sc_power, X, MGA_R, 70, 20, 8.5, 200);
[HGA_up,  HGA_down] = HGA_range(HGA_diam, sc_power, Ka, HGA_R, 34, 20, 8.5, 10000);

figure(1)
hold on
plot(range, LGA_up_70, 'g')
plot(range, LGA_down_70, 'k')
plot(range, MGA_up, 'b')
plot(range, MGA_down, 'r')
xlabel('Spacecraft Distance from Earth (Au)')
ylabel('Uplink Link Margin (dB)')
title('MGA and LGA Link Margin enroute to Uranus')
grid on
legend('LGA Uplink', 'LGA Downlink', 'MGA Uplink', 'MGA Downlink')
axis([0, 21, -5, 55])

figure(2)
hold on
plot(range, HGA_up, 'b')
plot(range, HGA_down, 'r')
xlabel('Spacecraft Distance from Earth (Au)')
ylabel('Uplink Link Margin (dB)')
title('HGA Link Margin enroute to Uranus')
grid on
legend('HGA Uplink', 'HGA Downlink')
axis([0, 21, -5, 55])

% %uplink
% figure(1)
% hold on
% plot(range, HGA_up, 'k--')
% plot(range, MGA_up, 'k:')
% plot(range, LGA_up_34)
% plot(range, LGA_up_70)
% plot(range, LGA_up_70_high)
% plot(range, LGA_up_70_34)
% plot(range, LGA_up_70_34_34)
% plot(range, LGA_up_70_high_34_34)
% plot(range, LGA_up_Worst_case)
% xlabel('SC Distance from Earth (Au)')
% ylabel('Uplink Link Margin (dB)')
% title('Plotting Link Margin for Uplink over the course of the journey')
% grid on
% axis([0, 21, -5, 65])
% grid on
% legend('HGA(110W, Ka, 4000b/s) to DSN 34m,', 'MGA(110W, X, 6b/s) to DSN 70m', ...
%     'LGA (150W, X, 1b/s) to 34m', 'LGA to 70m', 'LGA to 70m(400kW)', ...
%     'LGA to 70m + 34m', 'LGA to 70m + 2*34m', 'LGA to 70m(400kW) + 2*34m',...
%     'Worst Case: LGA to 2 DSN locations with 70m(400kW) + 2*34m', 'location', 'best')
% 
% 
% %downlink
% figure(2)
% hold on
% plot(range, HGA_down, 'k--')
% plot(range, MGA_down, 'k:')
% plot(range, LGA_down_34, 'm')
% plot(range, LGA_down_70, 'r')
% plot(range, LGA_down_70_high, 'g')
% plot(range, LGA_down_70_34, 'c')
% plot(range, LGA_down_70_34_34, 'b')
% plot(range, LGA_down_70_high_34_34, 'm--')
% plot(range, LGA_down_Worst_case, 'r--')
% xlabel('SC Distance from Earth (Au)')
% ylabel('Downlink Link Margin (dB)')
% title('Plotting Link Margin for Downlink over the course of the journey')
% grid on
% axis([0, 21, -5, 65])
% grid on
% legend('HGA(110W, Ka, 4000b/s) to DSN 34m,', 'MGA(110W, X, 6b/s) to DSN 70m', ...
%     'LGA (150W, X, 1b/s) to 34m', 'LGA to 70m', 'LGA to 70m(400kW)', ...
%     'LGA to 70m + 34m', 'LGA to 70m + 2*34m', 'LGA to 70m(400kW) + 2*34m',...
%     'Worst Case: LGA to 2 DSN locations with 70m(400kW) + 2*34m', 'location', 'best')
% 



%% Data rate calclation
clc

Data_amount = 2;   %GB
Data_bits = Data_amount*10^9;   %bits
Transmit_days = 13;
Transmit_hrs_per_day = 14;
transmit_sec = Transmit_days * Transmit_hrs_per_day * 3600;
Max_data_rate = Data_bits/transmit_sec;




