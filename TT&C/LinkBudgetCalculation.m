function [LM_uplink,  LM_downlink] = LinkBudgetCalculation(sc_diam, sc_power, sc_freq, sc_Data_rate, DSN_diam, DSN_power, DSN_freq, DSN_Data_rate)
    %Inital Variables
    T0 = 290; %"Noise figure are always quoted at this reference temp"-SMAD
    eta_ka = 0.55; %overall efficiency for TWTA (SMAD page 637)
    R = 21;  %Range to Uranus in AU
    R = R * 1.49*10^8; %range to Uranus in km 
%     Latm = -10; %atmospheric loss  (clear day) (unit: dB)
 
    Latm_uplink = -1;
    %atmospheric loss depends on what band we are transmitting in, with ka
    %band being affected more than x band
    if sc_freq > 10
        Latm_downlink = -5;
    else
        Latm_downlink = -1;
    end
    Cov_DSN = 120; %Coverage area for DSN (I believe it should be 180?)

    Rb_SC = sc_Data_rate;%linspace(100, 10000, 100); %bits per seconds rate for S/C (500<->32000)-estimated from Juno
    Rb_SC_dB = 10*log10(Rb_SC); %Data rate in dB-Hz
    Rb_DSN = DSN_Data_rate;%linspace(100, 32000, 100); %bits per seconds rate for DSN (500<->32000)
    Rb_DSN_db = 10*log10(Rb_DSN); %Data rate in dB-Hz

    % Know inputs Ground Station
    dia_g = DSN_diam; %DSN Antennas either 34m or 70m
    Ptx_g = DSN_power*10^3; %transmit power of ground (current estimate from SMAD) (W)
    Lin_g = -1; %line loss for the ground station (Known), dB
    Lbackoff_g = -1; %backoff loss (assumed from SMAD table 16-13), dB
    Louttx_g = Lin_g + Lbackoff_g; %backoff + Lineloss ground station (calc), dB
    eta_g = .7; % Antenna Efficiency ground station (Known)
    Ltxrx_g = -0.1; %Pointing loss ground station (assume), ddB
    Ta_g = 160; %Antenna Noise Temp (assume) Kelvin
    Tl_g = 10; %Noise temperture  ground station (assume)
    T_rec_g = 250;
    F_g = 1 + T_rec_g / T0; %GS Receiver Noise (assume)(1 + T/290), assuming T = 290

    % Know inputs Spacecraft
    dia_sc = sc_diam; %Antennas on Spacecraft Bus (m)
    Ptx_sc = sc_power; %51:150; %transmit power of SC (current estimate from SMAD) (W)
    Lin_sc = -1; %Line loss spacecraft (estimated)
    Lbackoff_sc = -1; %backoff loss (assumed from SMAD table 16-13)
    Louttx_sc = Lin_sc + Lbackoff_sc; %backoff + line loss spacecraft (calc)

    etLin_sca_sc = .7; % Antenna Efficiency spacecraft (Assume)
    Ltxrx_sc = -1; %Pointing loss spacecraft (assume)
    Ta_s = 3000; %Antenna Noise Temp (assume) Kelvin
    Tl_sc = 200; %Noise temperture (input from thermals)
    T_sc = 70;  %Spacecraft avg. Temp (K)
    F_sc = 1 + T_sc/T0 ; %SC Receiver Noise (assume)

    % Uplink(DSN to Spacecraft)

    f_up = DSN_freq; %uplink for X band

    %%%Gateway Termnial Type-DSN
    beamwidth_up = 21./(f_up*dia_g); %half-power beamwidth (deg)

    %gain of DSN
    Gain_gs = 20.4 + 20*log10(f_up) + 20*log10(dia_g) + 10*log10(eta_g);   %dBi

    %Gain of transmit antenna (dBi)
    Gr_sc_up = 20.4 + 20*log10(f_up) + 20*log10(dia_sc) + 10*log10(eta_ka); 

    %Decibels Conversion
    P_ref = 1; %in SMAD values are consistent to 1W for all units
    PtxdB_g = 10*log10(Ptx_g./P_ref);
    EIRP_g = PtxdB_g + Gain_gs + Louttx_g; %equation 16-20 SMAD

    %%%Space loss Calcuations
    Ls_up = -(92.45 + 20*log10(R) +  20*log10(f_up)); %Space loss GS to SC (unit: dB)
    Net_loss_up = Ls_up + Latm_uplink; %Atm loss + Space Loss

    %Total Power at receiver
    C_up = EIRP_g + Gr_sc_up + Net_loss_up + Lin_sc + Ltxrx_sc;

    %Effective system noise temperture
    L_up = 0.5; %Assumed loss ratio due to components upstream of receiver and downstream of antenna
    Ts_up = Ta_s * L_up + (1-L_up) * Tl_sc + (F_sc-1) * T0;    %K
    Ts_up_dB = 10*log10(Ts_up); %dbK

    %Receiver G/T
    G_T_sc = Gr_sc_up - Ts_up_dB; %(unit: dB/k)

    Lcomb_up = Ls_up + Latm_uplink + Louttx_g + Ltxrx_sc + Ltxrx_g;

    %Receiver C/N0 (k = 228.6 is Boltzman constant in dBJ/K)
    C_N0_up = EIRP_g + G_T_sc + Lcomb_up + 228.6;

    %Available Eb/N0
    Eb_N0_up = C_N0_up - Rb_DSN_db;

    % Downlink, S/C to DSN
    f_down = sc_freq; %downlink for KA band

    %User Termnial Type
    beamwidth_down = 21./(f_down*dia_sc); %half-power beamwidth (deg)

    %Satellite antenna type
    Gr_sc_down = 20.4 + 20*log10(f_down) + 20*log10(dia_sc) + 10*log10(eta_ka);

    %Decibels Conversion
    PtxdB_sc = 10*log10(Ptx_sc/P_ref); %decibal value for Sc power

    EIRP_sc = PtxdB_sc + Gr_sc_down + Louttx_sc; %%%%% EIRP of our spacecraft

    %Propagation Range and Total Net Loss
    Ls_down = -(92.45 + 20*log10(R) + 20*log10(f_down)); %space loss
    Net_Loss_down = Latm_downlink + Ls_down; %Total Loss through space
    
    %Gain of receive antenna-Ground Station(dBi)
    Gr_recieve_down = 20.4 + 20*log10(f_down) + 20*log10(dia_g) + 10*log10(eta_g); 

    %Total Power at receiver
    C_down = EIRP_sc + Gr_recieve_down + Net_Loss_down + Lin_g + Ltxrx_g+Ltxrx_sc; 

    %Effective system noise temperture
    L_down = 0.5; %Assumed loss ratio due to components upstream of receiver and downstream of antenna
    Ts_down = Ta_g*L_down + (1-L_down)*Tl_g + (F_g-1)*T0; %K
    Ts_down = 10*log10(Ts_down);    %dB-K
    
    %Receiver G/T
    G_T_g = Gr_recieve_down - Ts_down;
    Lcomb_down = Net_Loss_down + Louttx_sc + Ltxrx_g+ Ltxrx_sc

    %Receiver C/N0 (k = 228.6 is Boltzman constant in dBJ/K)
    C_N0_down = EIRP_sc + G_T_g + Lcomb_down + 228.6;  

    %Available Eb/N0
    Eb_N0_down = C_N0_down - Rb_SC_dB;

    % Link Margin
    %end to end Eb/N0 (eq. 16-32)
%     end_Eb_N0 = 10 * log10((10.^-(Eb_N0_up/10)) + (10.^-(Eb_N0_down/10)))
    % modem_loss = -1.2; %modem implementation loss SMAD table 16-11
    req_Eb_No = -2.99; %Required Eb/N0 SMAD table 16-12
    LM_uplink = Eb_N0_up  + req_Eb_No; %Link Margin Value
    LM_downlink = Eb_N0_down + req_Eb_No;
%     PtxdB_sc
%     Gr_sc_down
%     Gain_gs
%     Lcomb_down
%     Rb_SC_dB
%     Ts_down
   
%     PtxdB_sc+Gr_sc_down+Gain_gs+Lcomb_down-Rb_SC_dB-Ts_down+228.6
end