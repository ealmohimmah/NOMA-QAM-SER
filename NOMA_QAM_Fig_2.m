%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB script is used to generate Figure 2 in the paper:
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
clear;
clc;

p2 = 0.12;
p1 = 1 - p2;
M1 = 4;
M2 = 16;
M = M1*M2;

x1 = sqrt(p1)*qammod((0:M1-1)',M1,'UnitAveragePower',true);
x2 = sqrt(p2)*qammod((0:M2-1)',M2,'UnitAveragePower',true);

% All possible superimposed symbols
xt(:,1) = x1(1) + x2;
xt(:,2) = x1(2) + x2;
xt(:,3) = x1(3) + x2;
xt(:,4) = x1(4) + x2;
xt = xt(:);

% User 1 constellation
figure;
axis(gca,'equal');
axis([-1.2 1.2 -1.2 1.2]);
grid on;
box on;
hold on;
plot(-1.5:1.5,zeros(1,4),'k','linewidth',1);
plot(zeros(1,4),-1.5:1.5,'k','linewidth',1);
viscircles([0 0],sqrt(p1),'LineStyle','--','linewidth',0.5,'color','b');
plot(x1,'bo','linewidth',0.5,'markersize',6,'MarkerFaceColor','blue','MarkerEdgeColor','black');

% User 2 constellation
figure;
axis(gca,'equal');
axis([-1.2 1.2 -1.2 1.2]);
grid on;
box on;
hold on;
plot(-1.5:1.5,zeros(1,4),'k','linewidth',1);
plot(zeros(1,4),-1.5:1.5,'k','linewidth',1);
viscircles([0 0],sqrt(p2),'LineStyle','--','linewidth',0.5,'color','r');
plot(x2,'ro','linewidth',0.5,'markersize',6,'MarkerFaceColor','red','MarkerEdgeColor','black');

% Superimposed signal constellation
figure;
hold on;
axis(gca,'equal');
axis([-1.2 1.2 -1.2 1.2]);
grid on;
box on;
plot(-1.5:1.5,zeros(1,4),'k','linewidth',1);
plot(zeros(1,4),-1.5:1.5,'k','linewidth',1);
viscircles([0 0],sqrt(p1),'LineStyle','--','linewidth',0.5,'color','b');
viscircles([sqrt(p1)*cosd(45) sqrt(p1)*cosd(45)],sqrt(p2),'LineStyle','--','linewidth',0.5);
plot(xt,'go','linewidth',0.5,'markersize',6,'MarkerFaceColor','green','MarkerEdgeColor','black');
