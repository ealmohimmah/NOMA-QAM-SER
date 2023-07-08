%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB script is used to generate Figure 3 in the paper:
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
clear
clc

M1 = 16;
P = 2/3*(M1-1);
g = sqrt(M1)-1;
refpts = qammod(0:M1-1,M1);
s1 = refpts(real(refpts) > 0 & imag(refpts) > 0);

% line([-sqrt(M) sqrt(M)], [0 0], 'color', 'k', 'linewidth', 1);hold on;
% line([0 0], [-sqrt(M) sqrt(M)], 'color', 'k', 'linewidth', 1);

% Plot Decision Region for
% interior symbol
rectangle('Position',[0 0 2 2], 'linestyle','-.', 'FaceColor',[252 71 201 75]/255,'EdgeColor','k','LineWidth',0.5); hold on;
% edge symbol
rectangle('Position',[2 0 4 2], 'linestyle','none', 'FaceColor',[71 201 252 75]/255,'EdgeColor','k','LineWidth',0.5); hold on;
% corner symbol
rectangle('Position',[2 2 4 4], 'linestyle','-.', 'FaceColor',[201 252 71 75]/255,'EdgeColor','k','LineWidth',0.5); hold on;

plot(s1, 'kx', 'linewidth', 1, 'MarkerSize', 12); hold on;

xlabel('In-Phase');
ylabel('Quadrature');
grid off;
% box on;
axis equal
axis([0 sqrt(M1) 0 sqrt(M1)]);

ax = gca;
ax.XTick = [0 1 2 3];
ax.YTick = [0 1 2 3];
ax.TickLabelInterpreter='latex';
ax.XTickLabel = {'$0$', '$\xi_1\sqrt{a_{1}E_{s}}$', '$2\xi_1\sqrt{a_{1}E_{s}}$', '$3\xi_1\sqrt{a_{1}E_{s}}$'};
ax.YTickLabel = {'$0$', '$\xi_1\sqrt{a_{1}E_{s}}$', '$2\xi_1\sqrt{a_{1}E_{s}}$', '$3\xi_1\sqrt{a_{1}E_{s}}$'};

M2 = 16;
P2 = 0.05;
g = -1;
s2 = sqrt(P2)*qammod(0:M2-1,M2);

x1 = s1(4) + s2;
x2 = s1(2) + s2;
x3 = s1(1) + s2;

plot(x1, 'o', 'linewidth', 0.5, 'MarkerSize', 6, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
plot(x2, 'o', 'linewidth', 0.5, 'MarkerSize', 6, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k');
plot(x3, 'o', 'linewidth', 0.5, 'MarkerSize', 6, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');


A2 = [-3 -1 1 3];
idx = 1;
for i = 1:4
    yp = 1+sqrt(P2)*A2(i)-0.15;
    for j = 1:4
        xp = 1+sqrt(P2)*A2(j)-0.1;
        str = ['$x_{' num2str(idx) '}$']; 
        text(xp, yp, str, 'interpreter', 'latex');
        idx = idx + 1;
    end
end

idx = 1;
for i = 1:4
    yp = 1+sqrt(P2)*A2(i)-0.15;
    for j = 1:4
        xp = 3+sqrt(P2)*A2(j)-0.1;
        str = ['$t_{' num2str(idx) '}$']; 
        text(xp, yp, str, 'interpreter', 'latex');
        idx = idx + 1;
    end
end

idx = 1;
for i = 1:4
    yp = 3+sqrt(P2)*A2(i)-0.15;
    for j = 1:4
        xp = 3+sqrt(P2)*A2(j)-0.1;
        str = ['$z_{' num2str(idx) '}$']; 
        text(xp, yp, str, 'interpreter', 'latex');
        idx = idx + 1;
    end
end

text(0.9, 0.85, '$s_{1}$', 'interpreter', 'latex', 'color', 'k', 'fontsize',12);
text(2.9, 0.85, '$s_{2}$', 'interpreter', 'latex', 'color', 'k', 'fontsize',12);
text(2.9, 2.85, '$s_{3}$', 'interpreter', 'latex', 'color', 'k', 'fontsize',12);
text(0.9, 2.85, '$s_{4}$', 'interpreter', 'latex', 'color', 'k', 'fontsize',12);


