%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB script is used to generate Figure 4 in the paper:
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

sty = ['x' 'o' 's' '+'];

p1 = [85/100 16/21 9/14];
p2 = 1 - p1;
M1 = 4;
M2 = 16;
M  = M1*M2;

% Axes limits
minX = -1.5; maxX = 1.5; minY = -1.5; maxY = 1.5;

for i = 1:length(p1)
    figure(i);
    hold on; axis(gca,'equal'); axis([minX maxX minY maxY]);
    set(gca,'xticklabel',{[]},'yticklabel',{[]});
    xlabel('In-Phase'); ylabel('Quadrature'); grid on; box on;
    % plot the axes
    line([minX maxX],[0 0],'Color','k','LineStyle','-','linewidth',1);
    line([0 0],[minY maxY],'Color','k','LineStyle','-','linewidth',1);
end

for i = 1:length(p1)
    x1 = sqrt(p1(i)).*qammod((0:M1-1)',M1,'UnitAveragePower',true);
    x2 = sqrt(p2(i)).*qammod((0:M2-1)',M2,'UnitAveragePower',true);
    
    for j = 1:M1
        for k = 1:M2
            x(j,k) = x1(j) + x2(k);
        end
        figure(i)
        c = plot(x(j,:),sty(j),'linewidth',1.5,'markersize',8);
%         scatter(real(x(j,:)),imag(x(j,:)),'filled','MarkerFaceAlpha',4/8,...
%             'MarkerEdgeColor','k');
    end
    str1 = ['$p_{1}$' ' = ' num2str(p1(i))];
    str2 = ['$p_{2}$' ' = ' num2str(p2(i))];
    txt = text(-1,1.6,str1,'interpreter','latex','fontsize',14);
    txt = text(0.20,1.6,str2,'interpreter','latex','fontsize',14);
    
end



