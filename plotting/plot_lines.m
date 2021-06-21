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
%
% This script contains the plotting for the 3D noise for points and lines experiments 
% from the main paper (Figure 3, bottom row)
% To reproduce, run the script after main_ordinary_3d_lines_noise_3d.m

function plot_lines(matpath, save_pref)
    if nargin == 0
        matpath = 'results/ordinary3Dlines';
        save_pref = 'ord3d_lines';
    end    

    load(matpath);

    npts = int32(npts);

    close all;

    i= 0; w= 375; h= 330; wnarr = 340; h_norm=300;

    wl = 600;
    hl = 110;
    dol = 1;


    yrange = [0 3];
    figure('color','w','position',[w*i,100,w,h]);%i=i+1;
    xdrawgraph(npts,yrange,method_list,'deleted_mean_r','Mean Rot.',...
        'Number of Points','Rotation Err. (deg.)', save_pref, dol, wl, hl);
    figure('color','w','position',[w*i,100+h,wnarr,h]);i=i+1;
    xdrawgraph(npts,yrange,method_list,'deleted_med_r','Median Rot.',...
        'Number of Points','', save_pref, dol, wl, hl);
    set(gcf,'Units','Inches');    
        pos = get(gcf,'Position');
        set(gcf, 'PaperPositionMode', 'auto', 'PaperUnits','Inches','PaperSize',1.25*[pos(3), pos(4)]);    
        print(gcf, ['figures/med_rot_lines_2d.pdf'], '-dpdf', '-r0');

    figure('color','w','position',[w*i,100,w,h]);%i=i+1;
    yrange = [0 1.5];
    xdrawgraph(npts,yrange,method_list,'deleted_mean_t','Mean Translation',...
        'Number of Points','Translation Err. (%)', save_pref, dol, wl, hl);
    figure('color','w','position',[w*i,100+h,wnarr,h]);i=i+1;
    xdrawgraph(npts,yrange,method_list,'deleted_med_t','Median Translation',...
        'Number of Points','', save_pref, dol, wl, hl);

    set(gcf,'Units','Inches');    
        pos = get(gcf,'Position');
        set(gcf, 'PaperPositionMode', 'auto', 'PaperUnits','Inches','PaperSize',1.25*[pos(3), pos(4)]);    
        print(gcf, ['figures/med_trans_lines_2d.pdf'], '-dpdf', '-r0');
    

    
    yrange = [0 15];
    figure('color','w','position',[w*i,100+h,w,h_norm]);i=i+1;
    xdrawgraph(npts,yrange,method_list,'deleted_med_c','Points+Lines',...
        'Number of Points','Med. runtime (ms)', save_pref, dol, wl, hl);
    
%     title_figure('Points+Lines, 3D+2D Noise', 'title_lines', w, h, npts, method_list);
%     title_figure('Points, 3D+2D Noise', 'title_pt3d', w, h, npts, method_list);
%     title_figure('Points, 2D Noise', 'title_pt2d', w, h, npts, method_list);
%     title_figure('KITTI', 'title_kitti', w, h, npts, method_list);
%     title_figure('TUM', 'title_tum', w, h, npts, method_list);
end

function title_figure(title, save_pref, w, h, npts, method_list)
    i = 1;    
    figure('color','w','position',[w*i,100+h,w,h]);i=i+1;
    
    plot(npts,method_list(i).deleted_med_c,method_list(i).marker,...
    'color',method_list(i).color,...
    'markerfacecolor',method_list(i).markerfacecolor,...
    'displayname',method_list(i).name, ...
    'LineWidth',2,'MarkerSize',8);
    ylabel(title,'FontSize',16,'FontWeight', 'bold');
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf, 'PaperPositionMode', 'auto', 'PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(gcf, ['figures/' save_pref '.pdf'], '-dpdf','-r0','-fillpage');
end