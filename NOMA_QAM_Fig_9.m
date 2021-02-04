%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB script is used to generate Figure 9 in the paper:
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

%% SER of the second user vs. the PA coefficient with variousmodulation schemes

clear; clc;
set(groot, 'defaultLegendInterpreter','latex');
style = ["-.", "--", ".-", ":", "-"];

SNR_dB = [20 26 31 26 31];
a1 = 0.5:0.001:1;
M1 = [4  4  4 16 16];
M2 = [4 16 64  4 16];

for k = 1:length(SNR_dB)
    ser2 = zeros(1,length(a1));
    for i = 1:length(a1)
        [~, ser2(i)] = Theoretical_SER(M1(k),M2(k),a1(i),SNR_dB(k));        
    end
    
    semilogy(a1,ser2,style(k),'linewidth',1.5,'markersize',10, ...
        'markerindices',1:10:length(a1)); hold on;
end

axis([0.5 1 1e-6 1]); grid on;
xlabel('$a_1$', 'Interpreter', 'latex');
ylabel('$SER_2$', 'Interpreter', 'latex');
legend('$M_1=4$, $M_2=4$, $SNR=20$', '$M_1=4$, $M_2=16$, $SNR=26$', ...
    '$M_1=4$, $M_2=64$, $SNR=31$', '$M_1=16$, $M_2=4$, $SNR=26$', ...
    '$M_1=16$, $M_2=16$, $SNR=31$','location','sw')


%% Theoretical SER1 and SER2
function [ser1, ser2] = Theoretical_SER(M1,M2,a1,SNR_dB)
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
    for i = 1:sqrt(M2)
        for j = 1:sqrt(M2)
            Pe(n) = Pe(n) ...
                + (M1-sqrt(M1))*qfunc(y*(k1*sqrt(a1)+A(i)*k2*sqrt(a2))) ...
                - (M1-2*sqrt(M1)+1) ...
                * qfunc(y*(k1*sqrt(a1)+A(i)*k2*sqrt(a2))) ...
                * qfunc(y*(k1*sqrt(a1)+A(j)*k2*sqrt(a2)));
        end
    end
end

ser1 = 4*Pe/M1/M2;

% Calculate ser2
% First calculate ser2 if s1 was detected correctly
ser2cS1 = 4*(1 - 1/sqrt(M2))*qfunc(yy*k2*sqrt(a2)) - 4*(1 - 2/sqrt(M2)+1/M2)*qfunc(yy*k2*sqrt(a2)).^2;
% Now calculate ser2
ser2 = (1 - ser1).*ser2cS1 + ser1;
end
