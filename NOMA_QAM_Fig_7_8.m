%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB script is used to generate Figure 7 and 8 in the paper:
%
% E. M. Almohimmah and M. T. Alresheedi, "Error Analysis of NOMA-Based VLC
% Systems With Higher Order Modulation Schemes," in IEEE Access, vol. 8, 
% pp. 2792-2803, 2020, doi: 10.1109/ACCESS.2019.2962331.
%
% Download article: https://ieeexplore.ieee.org/document/8943113
%
% This is version 1.0 (Last edited: 2021-01-04)
%
% License: This code is licensed under the GPLv2 license. If you in any 
% way use this code for research that results in publications, please cite
% our paper as described above.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SER of the first and second user vs the PA coefficient with various SNR values.

clear; clc;
set(groot, 'defaultLegendInterpreter','latex');
style = ['o', '*', 'x', 's', '^'];

% colors
red = [1 0 0];
green = [0 1 0];
blue = [0 0 1];
dred = [217 83 25]/255;
dgreen = [0 127 0]/255;
dblue = [0 114 189]/255;
lred = [255 153 200]/255;
lgreen = [119 172 48]/255;
lblue = [77 190 238]/255;
orange = [237 177 32]/255;
purple = [126 47 142]/255;
color1 = [dblue; dred; lgreen; orange; purple];
color2 = [blue; red; dgreen; orange; purple];


SNR = [22 24 26 28]; 
a1 = 0.5:0.001:1; % for theoretical SER
pa1 = 0.5:0.03:1; % for simulation

M1 = 4;
M2 = 16;

for n = 1:length(SNR)
    disp(['snrdB = ' num2str(SNR(n))]); % display current SNR
    
    % Theory SER
    for i = 1:length(a1)
        [ser1t(i), ser2t(i)] = Theoretical_SER(M1,M2,a1(i),SNR(n));
    end
    
    % Simulation SER
    for i = 1:length(pa1)
        disp(['a1 = ' num2str(pa1(i))]); % display current a1
        [ser1s(i), ser2s(i)] = Simulated_SER(M1,M2,pa1(i),SNR(n));
    end

    figure(1)
    semilogy(a1, ser1t, '-', 'color', color1(n,:), 'linewidth', 1.5);
    hold on;
    semilogy(pa1, ser1s, style(n), 'color', color2(n,:), 'linewidth',1);
    
    figure(2)
    semilogy(a1,ser2t,'-','color', color1(n,:),'linewidth',1); hold on;
    semilogy(pa1, ser2s, style(n), 'color', color2(n,:), 'linewidth',1);
end

figure(1)
grid on; xlabel('$a_1$', 'Interpreter', 'latex');
ylabel('$SER_1$', 'Interpreter', 'latex');
axis([0.5 1 1e-7 1]);
legend(['SNR = ' num2str(SNR(1)) ' (Theo.)'], ...
       ['SNR = ' num2str(SNR(1)) ' (Simu.)'], ...
       ['SNR = ' num2str(SNR(2)) ' (Theo.)'], ...
       ['SNR = ' num2str(SNR(2)) ' (Simu.)'], ...
       ['SNR = ' num2str(SNR(3)) ' (Theo.)'], ...
       ['SNR = ' num2str(SNR(3)) ' (Simu.)'], ...
       ['SNR = ' num2str(SNR(4)) ' (Theo.)'], ...
       ['SNR = ' num2str(SNR(4)) ' (Simu.)'], ...
        'location','ne');
   
figure(2)
grid on; xlabel('$a_1$', 'Interpreter', 'latex'); 
ylabel('$SER_2$', 'Interpreter', 'latex');
axis([0.5 1 1e-8 1]);

legend(['SNR = ' num2str(SNR(1)) ' (Theo.)'], ...
       ['SNR = ' num2str(SNR(1)) ' (Simu.)'], ...
       ['SNR = ' num2str(SNR(2)) ' (Theo.)'], ...
       ['SNR = ' num2str(SNR(2)) ' (Simu.)'], ...
       ['SNR = ' num2str(SNR(3)) ' (Theo.)'], ...
       ['SNR = ' num2str(SNR(3)) ' (Simu.)'], ...
       ['SNR = ' num2str(SNR(4)) ' (Theo.)'], ...
       ['SNR = ' num2str(SNR(4)) ' (Simu.)'], ...
        'location','sw');

%% Theoretical SER1 and SER2
function [ser1, ser2] = Theoretical_SER(M1,M2,a1,SNR_dB)
SNR = db2pow(SNR_dB);
a2 = 1 - a1;
k1 = 1/sqrt(2/3*(M1-1)); % scaling factor to normalize the power of s1
k2 = 1/sqrt(2/3*(M2-1)); % scaling factor to normalize the power of s2
yy = sqrt(2*SNR);

m = 1:sqrt(M2);
A = 2*m-1-sqrt(M2);    % the alphabet (integers) for s2

% calculate SER1
Pe = zeros(1,length(SNR));
for n = 1:length(SNR)
    y = yy(n);
    for i = 1:sqrt(M2)
        for j = 1:sqrt(M2)
            Pe(n) = Pe(n) + (M1-sqrt(M1))*qfunc(y*(k1*sqrt(a1)+A(i)*k2*sqrt(a2))) - (M1-2*sqrt(M1)+1)*qfunc(y*(k1*sqrt(a1)+A(i)*k2*sqrt(a2)))*qfunc(y*(k1*sqrt(a1)+A(j)*k2*sqrt(a2)));
        end
    end
end

ser1 = 4*Pe/M1/M2;

% Calculate SER2
% First calculate ser2 if s1 was detected correctly
ser2cS1 = 4*(1 - 1/sqrt(M2))*qfunc(yy*k2*sqrt(a2)) - 4*(1 - 2/sqrt(M2)+1/M2)*qfunc(yy*k2*sqrt(a2)).^2;
% Now calculate ser2
ser2 = (1 - ser1).*ser2cS1 + ser1;
end

%% Approximation SER1 and SER2
function [ser1, ser2] = Approximated_SER(M1,M2,a1,SNR_dB)
SNR = db2pow(SNR_dB);
a2 = 1 - a1;
k1 = 1/sqrt(2/3*(M1-1)); % scaling factor to normalize the power of s1
k2 = 1/sqrt(2/3*(M2-1)); % scaling factor to normalize the power of s2
yy = sqrt(2*SNR);

m = 1:sqrt(M2);
A = 2*m-1-sqrt(M2);    % the alphabet (integers) for s2

% calculate ser1
Pe = zeros(1,length(SNR));
for n = 1:length(SNR)
    y = yy(n);
    Pe(n) = sum(qfunc(y*(k1*sqrt(a1)+A*k2*sqrt(a2))));
end

ser1 = 4*Pe*(sqrt(M1) - 1)/sqrt(M1*M2);

% Calculate ser2
% First calculate ser2 if s1 was detected correctly
ser2cS1 = 4*(1 - 1/sqrt(M2))*qfunc(yy*k2*sqrt(a2));
% Now calculate ser2
ser2 = (1 - ser1).*ser2cS1 + ser1;
end

%% Simulated SER1 and SER2
function [ser1, ser2] = Simulated_SER(M1,M2,a1,SNR_dB)
a2 = 1 - a1;

% ACO-OFDM params
fftSize     = 512;               % fft size
cpSize      = 0;                 % cyle prefix size
nSubcar     = fftSize/4;         % number of subcarrier
modParams = [nSubcar,cpSize];

frmLen = 100;   % number of birs per frame
maxErrs = 1e3;  % target number of errors at each Es/No
maxBits = 1e7;  % maximum number of Bits at each Es/No

% QAM modulator and demodulator
QAMMod1 = comm.RectangularQAMModulator(M1,'NormalizationMethod', ...
    'Average Power','AveragePower',a1);
QAMMod2 = comm.RectangularQAMModulator(M2,'NormalizationMethod', ...
    'Average Power','AveragePower',a2);
QAMDemod1 = comm.RectangularQAMDemodulator(M1,'NormalizationMethod', ...
    'Average Power','AveragePower',a1);
QAMDemod2 = comm.RectangularQAMDemodulator(M2,'NormalizationMethod', ...
    'Average Power','AveragePower',a2);

% pre-allocate memory for ser1 and ser2
ser1 = zeros(length(SNR_dB),1);
ser2 = zeros(length(SNR_dB),1);
errorRate1 = comm.ErrorRate;
errorRate2 = comm.ErrorRate;

for n = 1:length(SNR_dB)
    snr = SNR_dB(n);
%     disp(['snrdB = ' num2str(snr)]); % display current SNR
    
    errorStats1 = zeros(1,3);
    errorStats2 = zeros(1,3);
    
    while errorStats1(2) < maxErrs && errorStats1(3) < maxBits
        
        % Generate random data
        dataIn_1 = randi([0 M1-1],nSubcar,1);
        dataIn_2 = randi([0 M2-1],nSubcar,1);
        
        % Modulate signals using M-QAM with power allocation
        dataMod_1 = QAMMod1(dataIn_1);
        dataMod_2 = QAMMod2(dataIn_2);
        
        % Total electrical signal
        x = dataMod_1 + dataMod_2;
        
        in = aco_ofdm_modulator(x, modParams);
        
        % The equalized received signal
        rxSig = awgn(in,snr,'measured');
        
        out = aco_ofdm_demodulator(rxSig, modParams);
        
        dataOut1 = QAMDemod1(out);
        errorStats1 = errorRate1(dataIn_1,dataOut1);
        
        y2 = out - QAMMod1(dataOut1); % Interference cancelation
        dataOut2 = QAMDemod2(y2);
        errorStats2 = errorRate2(dataIn_2,dataOut2);
    end
    
    ser1(n) = errorStats1(1);
    ser2(n) = errorStats2(1);
    reset(errorRate1);
    reset(errorRate2);
end
end

%% ACO-OFDM Modulation
function acoMod = aco_ofdm_modulator(in, modulator_param)
% get parameters
nSubcar   = modulator_param(1);
cpSize    = modulator_param(2);
fftSize   = nSubcar*4; % it's in ACO-OFDM

if mod(length(in),nSubcar) == 0 
    
    Nsym = length(in)/nSubcar;
    blkSize = fftSize+cpSize;
    acoMod = zeros(Nsym*blkSize,1);  %pre-allocate memeory
    
    for i = 1:Nsym
  
        in_temps = in(1+(i-1)*nSubcar:i*nSubcar);
        
        % pilot insertion and hermetian sysmetry to ensure a real data after FFT
        pilot_ins_data = zeros(fftSize,1);
        for j = 1:nSubcar
            pilot_ins_data(j*2) = in_temps(j);
            pilot_ins_data(fftSize - j*2+2) = conj(in_temps(j));
        end

        % fourier transform time doamain data
        IFFT_data =sqrt(fftSize)*ifft(pilot_ins_data);
       
        % remove negative part
        IFFT_data(IFFT_data<0) = 0;
        
        % add cycle prefix 
        acoMod((i-1)*blkSize+1:i*blkSize) = [IFFT_data(end-cpSize+1:end); IFFT_data];
        
    end
else
    error('ACO-OFDM wrong input size for modulation');
end
end

%% ACO-OFDM Demodulation 
function out = aco_ofdm_demodulator(in, demodParams)
% in: input time-domain data, 1-D column array
% out: subcarrier data, dimension (nSubcar,nOfdmSymbol)

nSubcar = demodParams(1);
cpSize  = demodParams(2);

fftSize = nSubcar*4; % FFT size (this relation is only true in ACO-OFDM)
blkSize = fftSize+cpSize;  % OFDM symbol length

if mod(length(in),blkSize) == 0
    nOfdmSymbol = length(in)/blkSize;    
    in          = reshape(in,blkSize,nOfdmSymbol);
    inFFT       = 1/sqrt(fftSize)*fft(in(cpSize+1:end,:));     % cycle prefix remove and conversion to frequency domain
    out         = 2*inFFT(2:2:2*nSubcar,:);    % Pilot remove and composenate the 1/2 attenuation in ACO-OFDM
    out         = reshape(out,nSubcar*nOfdmSymbol,1);
else
    error('ACO-OFDM wrong input size for demodulation');
end
end
