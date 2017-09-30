%% This Matlab Script simulates the BER for FBMC and OFDM, including channel estimation.
%{
% Ronald Nissel, rnissel@nt.tuwien.ac.at
% (c) 2016 by Institute of Telecommunications, TU Wien
% www.tc.tuwien.ac.at 

% Pilot symbol aided channel estimation is based on
% R. Nissel, M. Rupp, "On Pilot-Symbol Aided Channel Estimation in FBMC-OQAM", IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP 2016)

% BER and BEP in FBMC and OFDM: 
% R. Nissel, M. Rupp, "Bit Error Probability for Pilot-Symbol Aided Channel Estimation in FBMC-OQAM", IEEE International Conference on Communications (ICC 2016)
% Theoretical BEP calculation for higher modulation order than 4QAM, see
% R. Nissel, M. Lerch, M. Simko, M. Rupp, "Bit Error Probability for Pilot-Symbol-Aided OFDM Channel Estimation in Doubly-Selective Channels", International ITG Workshop on Smart Antennas (WSA 2014)
% and 
% R. Nissel, S. Caban, M. Rupp, "Closed-Form Capacity Expression for Low Complexity BICM with Uniform Inputs",IEEE International Symposium on Personal, Indoor, and Mobile Radio Communications (PIMRC2015)

% For perfect channel knowledge, FBMC and OFDM have the same BER (same number of subcarriers). However, FBMC requires less guard bands (and also does not include a CP) so that the spectral efficiency is higher!
% We assume the same number of subcarriers (OFDM, FBMC) and no CP to make the BER comparision fair (same number of transmitted bits).
%}
clear;clc;
close all;
addpath('./Theory');

M_SNR_OFDM_dB =[0:5:30];

NrRepetitions = 1000;
NrTime=50;
QAM_ModulationOrder = 16; % 4 16 64 128 ...

%% FBMC Object
FBMC = Modulation.FBMC(...
    16,...                          % Number subcarriers
    8,...                          % Number FBMC symbols
    15e3,...                        % Subcarrier spacing (Hz)
    15e3*14*12,...                  % Sampling rate (Samples/s)
    15e3*20,...                     % Intermediate frequency first subcarrier (Hz)
    false,...                       % Transmit real valued signal 
    'Hermite-OQAM',...              % Prototype filter (Hermite, PHYDYAS, RRC) and OQAM or QAM, 
    8, ...                          % Overlapping factor (corresponding to the prototype filter length)
    0, ...                          % Initial phase shift
    true ...                        % Polyphase implementation
    );

%% OFDM Object
OFDM = Modulation.OFDM(...
    16,...                          % Number subcarriers
    4,...                          % Number OFDM Symbols
    15e3,...                        % Subcarrier spacing (Hz)
    15e3*14*12,...                  % Sampling rate (Samples/s)
    15e3*20,...                     % Intermediate frequency first subcarrier (Hz)
    false,...                       % Transmit real valued signal
    0, ...                          % Cyclic prefix length (s), LTE: 1/15e3/14
    (8-1/2)*1/15e3*1/2 ...          % Zero guard length (s)
    );

%% PAM and QAM Object
PAM = Modulation.SignalConstellation(sqrt(QAM_ModulationOrder),'PAM');
QAM = Modulation.SignalConstellation(QAM_ModulationOrder,'QAM');

%% Channel Estimation Objects
ChannelEstimation_FBMC = ChannelEstimation.PilotSymbolAidedChannelEstimation(...
    'Diamond',...                           % Pilot pattern
    [...                                    % Matrix that represents the pilot pattern parameters
    FBMC.Nr.Subcarriers,...                 % Number of subcarriers
    4; ...                                  % Pilot spacing in the frequency domain
    FBMC.Nr.MCSymbols,...                   % Number of FBMC/OFDM Symbols
    5 ...                                   % Pilot spacing in the time domain
    ],...                                   
    'linear'...                             % Interpolation(Extrapolation) method 'linear','spline','FullAverage,'MovingBlockAverage',...
    );

%% Imaginary Interference Cancellation Objects
AuxiliaryMethod = ChannelEstimation.ImaginaryInterferenceCancellationAtPilotPosition(...
    'Auxiliary', ...                                    % Cancellation method
    ChannelEstimation_FBMC.GetAuxiliaryMatrix(2), ...   % PilotMatrix
    FBMC.GetFBMCMatrix, ...                             % Imaginary interference matrix
    16, ...                                             % Cancel 28 closest interferers
    2 ...                                               % Pilot to data power offset
    );

BER_FBMC_Aux = nan(length(M_SNR_OFDM_dB),NrRepetitions);
MSE_FBMC_Aux = nan(length(M_SNR_OFDM_dB),NrRepetitions);
Time_FBMC_Aux = nan(length(M_SNR_OFDM_dB),NrRepetitions);
for i_rep = 1:NrRepetitions
for i_SNR = 1:length(M_SNR_OFDM_dB)

        SNR_OFDM_dB = M_SNR_OFDM_dB(i_SNR);
        Pn_time = FBMC.PHY.SamplingRate/(FBMC.PHY.SubcarrierSpacing*FBMC.Nr.Subcarriers)*10^(-SNR_OFDM_dB/10);
        Tp=0;
        Tpd=0;
        %Generate data(time domain)
    for t = 1:NrTime
        if mod(t,2)==1
            Tp=Tp+1;
            [BinaryDataStream_FBMC_Aux_signal(:,Tp),xP_FBMC(:,Tp),x_FBMC_Aux(:,:,t),s_FBMC_Aux(:,t)]= FBMC_signal(AuxiliaryMethod,FBMC,PAM,ChannelEstimation_FBMC);
            index(Tp)=t;
            %have pilots
        else
            Tpd=Tpd+1;
            [BinaryDataStream_FBMC_Aux_data(:,Tpd),x_FBMC_Aux(:,:,t),s_FBMC_Aux(:,t)]= FBMC_data(AuxiliaryMethod,FBMC,PAM);
        %pure data
        end
    end
        %% Channel (doubly flat fading and AWGN in accordance with our testbed measurements!)     
        [h,~] = Jakes_Flat(FBMC,NrTime);
        h=abs(h);
        h=h./norm(h,'inf');
        % h = sqrt(1/2)*(randn(NrTime,1)+1j*randn(NrTime,1));
        n_FBMC = sqrt(1/2)*sqrt(Pn_time/2)*(randn(size(s_FBMC_Aux))+1j*randn(size(s_FBMC_Aux)));
        
        Tp = 0;
    for t = 1:NrTime
        r_FBMC_Aux(:,t) = h(t).*s_FBMC_Aux(:,t)+n_FBMC(:,t);
        %Demodulate FBMC signal
        y_FBMC_Aux(:,:,t) = FBMC.Demodulation(r_FBMC_Aux(:,t));
        %LS channel estimates at pilot positions
        if(mod(t,2)==1)
            y_FBMC_Aux_temp = y_FBMC_Aux(:,:,t);
            Tp = Tp+1;
            hP_LS_FBMC_Aux(Tp) =mean(y_FBMC_Aux_temp(ChannelEstimation_FBMC.PilotMatrix==1)./xP_FBMC(:,Tp)/...
                sqrt(AuxiliaryMethod.PilotToDataPowerOffset*AuxiliaryMethod.DataPowerReduction));
        end
    end
        %% Channel Estimation using Interpolation and calculate MSE
        tic;
        h_FBMC_Aux = interp1(index,hP_LS_FBMC_Aux,(1:NrTime),'linear','extrap');
        Time_FBMC_Aux(i_SNR,i_rep)=toc;
        MSE_FBMC_Aux(i_SNR,i_rep)=var(h-(h_FBMC_Aux).');
            %h_FBMC_Aux = svminterp(index,hP_LS_FBMC_Aux,h,NrTime);
        %tic;
        %h_FBMC_Aux = svminterp_real(index,hP_LS_FBMC_Aux,h,NrTime);
        %Time_FBMC_Aux(i_SNR,i_rep)=toc;
        %MSE_FBMC_Aux(i_SNR,i_rep)=var(h-(h_FBMC_Aux));
        %Equalized received symbols at data position
        Tp=0;Tpd=0;
        for t = 1:NrTime
            if mod(t,2)==1
                y_FBMC_Aux_temp=y_FBMC_Aux(:,:,t);
                y_EQ_FBMC_Aux_temp = real(y_FBMC_Aux_temp(AuxiliaryMethod.PilotMatrix==0)./h_FBMC_Aux(t)...
                    /sqrt(AuxiliaryMethod.DataPowerReduction));
                %Detect BitStream
                Tp=Tp+1;
                DetectedBitStream_FBMC_Aux_signal(:,Tp) = PAM.Symbol2Bit(real(y_EQ_FBMC_Aux_temp(:)));
                BER_FBMC_Aux_signal(Tp)=mean(BinaryDataStream_FBMC_Aux_signal(:,Tp)~=DetectedBitStream_FBMC_Aux_signal(:,Tp));
            else
                y_FBMC_Aux_temp=y_FBMC_Aux(:,:,t);
                y_EQ_FBMC_Aux_temp = real(y_FBMC_Aux_temp./h_FBMC_Aux(t));
                %Detect BitStream
                Tpd=Tpd+1;
                DetectedBitStream_FBMC_Aux_data(:,Tpd) = PAM.Symbol2Bit(real(y_EQ_FBMC_Aux_temp(:)));
                BER_FBMC_Aux_data(Tp)=mean(BinaryDataStream_FBMC_Aux_data(:,Tp)~=DetectedBitStream_FBMC_Aux_data(:,Tp));
             end
        end
        %% Calculate BER
        BER_FBMC_Aux(i_SNR,i_rep)=mean([BER_FBMC_Aux_data,BER_FBMC_Aux_signal]);
        
end
    if mod(i_rep,100)==0
       disp([int2str(i_rep/NrRepetitions*100) '%']);
    end
end

%% Theoretical BEP for perfect channel knowledge 
% BEP_4QAM = 1/2-1./(2*sqrt(2*(1+10.^(-M_SNR_OFDM_dB/10))-1));
M_SNR_OFDM_dB_morePoints = min(M_SNR_OFDM_dB):0.5:max(M_SNR_OFDM_dB);
BEP_perfect = BitErrorProbabilityRayleighDoublyFlat(M_SNR_OFDM_dB_morePoints,QAM.SymbolMapping,QAM.BitMapping);

%% Plot MSE
figure();
semilogy(M_SNR_OFDM_dB,trimmean(MSE_FBMC_Aux',2),'red -o');
hold on;
%semilogy(M_SNR_OFDM_dB_morePoints,BEP_perfect','black');
xlabel('SNR for OFDM (dB)'); 
ylabel('MSE');
legend('Simulation: FBMC Auxiliary','Location','SouthWest');
grid on;
%% Plot BER and BEP
figure();
semilogy(M_SNR_OFDM_dB,trimmean(BER_FBMC_Aux',2),'red -o');
hold on;
%semilogy(M_SNR_OFDM_dB_morePoints,BEP_perfect','black');
xlabel('SNR for OFDM (dB)'); 
ylabel('BER');
legend('Simulation: FBMC Auxiliary','Location','SouthWest');
grid on;

%% Plot H
figure();
%plot(abs(h),'red -');
plot(h,'red -');
hold on;
%plot(abs(h_FBMC_Aux),'green --');
%plot(angle(h),'blue -.');
%plot(angle(h_FBMC_Aux),'black :');
plot(real(h_FBMC_Aux),'black :');
xlabel('Time'); 
%ylabel('Magnitude and Phase');
ylabel('Magnitude');
%legend('h-Mag','h-FBMC-Aux-Mag',...
%    'h-Pha','h-FBMC-Aux-Pha','Location','SouthWest');
legend('h','h-FBMC-Aux');
grid on;

%% Plot Pilot Pattern
figure();
ChannelEstimation_FBMC.PlotPilotPattern(AuxiliaryMethod.PilotMatrix)
title('FBMC Auxiliary');

%% Calculate and Plot Expected Transmit Power Over Time
[Power_FBMC_Aux,t_FBMC] = FBMC.PlotTransmitPower(AuxiliaryMethod.PrecodingMatrix*AuxiliaryMethod.PrecodingMatrix');
[Power_OFDM,t_OFDM] = OFDM.PlotTransmitPower;
figure();
plot(t_FBMC,Power_FBMC_Aux,'red');
hold on;
plot(t_OFDM,Power_OFDM,'black ');
legend({'FBMC Auxiliary','OFDM'});
ylabel('Transmit Power');
xlabel('Time(s)');

%% Calculate Power Spectral Density
[PSD_FBMC_Aux,t_FBMC] = FBMC.PlotPowerSpectralDensity(AuxiliaryMethod.PrecodingMatrix*AuxiliaryMethod.PrecodingMatrix');
[PSD_OFDM,t_OFDM] = OFDM.PlotPowerSpectralDensity;
figure();
plot(t_FBMC,10*log10(PSD_FBMC_Aux),'red');
hold on;
plot(t_OFDM,10*log10(PSD_OFDM),'black ');
legend({'FBMC Auxiliary','OFDM'});
ylabel('Power Spectral Density (dB)');
xlabel('Frequency (Hz)');

save('temp_Cubic_real_50.mat');