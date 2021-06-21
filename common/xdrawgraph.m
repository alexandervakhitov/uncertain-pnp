%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This toolbox illustrates how to use the Uncertainty-aware PnPL
% algorithms described in:
%
%       A. Vakhitov, L. Ferraz, A. Agudo, F. Moreno-Noguer
%       Uncertainty-Aware Camera Pose Estimation from Points and Lines 
%       CVPR 2021
%
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the version 3 of the GNU General Public License
% as published by the Free Software Foundation.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.       
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xdrawgraph(xs,yrange,method_list,field,ti,lx,ly,save_pref,do_legend, w, h)
%the legend is at upper right in default
if (nargin < 9)
   do_legend = 0; 
   save_pref = '';
end

box('on');
hold('all');

p= zeros(size(method_list));
mnames = {};
for i= 1:length(method_list)
    mnames{i} = method_list(i).name;
    p(i)= plot(xs,method_list(i).(field),method_list(i).marker,...
        'color',method_list(i).color,...
        'markerfacecolor',method_list(i).markerfacecolor,...
        'displayname',method_list(i).name, ...
        'LineWidth',2,'MarkerSize',8);
%        );
end
if length(xs)~=1
    ylim(yrange);
    xlim(xs([1 end]));
end

if length(xs)>5       
    set(gca,'xtick',xs(1:round(length(xs)/5):length(xs)));
else
    set(gca,'xtick',xs);
end
title(ti,'FontSize',16, 'FontWeight', 'normal');
xlabel(lx,'FontSize',16);
ylabel(ly,'FontSize',16);

tix=get(gca,'xtick')';
a = get(gca,'XTickLabel');
set(gca, 'XTickLabel', a, 'fontsize', 16);

if isa(xs, 'int32')
    set(gca,'xticklabel',num2str(tix,'%d'))
else
    set(gca,'xticklabel',num2str(tix,'%.1f'))
end
tiy=get(gca,'ytick')';
set(gca,'yticklabel',num2str(tiy,'%.1f'))



set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf, 'PaperPositionMode', 'auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
if ~exist('./figures/', 'dir')
   mkdir('./figures')
end
print(gcf, ['figures/' save_pref field '.pdf'], '-dpdf','-r0','-fillpage');

if (do_legend == 2)
   legend(p); 
end
if (do_legend == 1)
    hold off;    
    fig2 = figure('position',[100,200,w,h]);
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf, 'PaperPositionMode', 'auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);    
    legendflex(p, mnames, 'ref', fig2, 'box', 'off', 'padding', [0 0 0], 'FontSize', 16, 'nrow', 2, 'anchor', [2,2]);
    set(gcf, 'PaperPositionMode', 'auto');
    print(gcf, ['figures\' save_pref 'legend.pdf'], '-dpdf', '-r0', '-fillpage');
end


return
