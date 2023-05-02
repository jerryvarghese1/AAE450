%% Plotting SC Transmit Power vs SC Data Rate at 3dB Margin and unspecified Error
clc
clear all
close all
SC_Antenna_1 = 2; %m
SC_Antenna_2 = 2.4; %m
SC_Antenna_3 = 2.8; %m
SC_Antenna_4 = 3.2; %m
SC_Antenna_5 = 3.6; %m

SC_power = 0:1:200;
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
