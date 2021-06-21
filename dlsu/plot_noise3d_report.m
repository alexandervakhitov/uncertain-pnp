function plot_noise3d_report
    load('noise3d_report.mat');
    metrics_labels = {'rot', 'trans'};
    linestyles = {'-', '-', '--', '-', '--', '-'};
    colors = {'r', 'g', 'g', 'b', 'b', 'k'};    
    y_low_high = [[0,5]; [0,80]; ];  
    close all;
    lw = 2;
    for rep_id = 1:2
        figure(rep_id);
        hold on;
        for meth_id = 1:length(method_labels)
%             plot(noises_3d, report(:, meth_id, rep_id), [linestyles{meth_id} colors{meth_id}], 'Linewidth', lw, 'DisplayName', method_labels{meth_id}); 
            neg = -0.1*report(:, meth_id, 3*rep_id);
            pos = +0.1*report(:, meth_id, 3*rep_id);
            errorbar(noises_3d, report(:, meth_id, 3*rep_id-2), neg, pos, [linestyles{meth_id} colors{meth_id}], 'Linewidth', lw, 'DisplayName', method_labels{meth_id});             
        end
        ylabel(metrics_labels{rep_id});
        xlabel('3D Noise level., m.');
        ylim(y_low_high(rep_id,:)); 
        legend('Location', 'best');
        
        set(gcf,'Units','Inches');    
        pos = get(gcf,'Position');
        set(gcf, 'PaperPositionMode', 'auto', 'PaperUnits', 'Inches','PaperSize',1.25*[pos(3), pos(4)]);    
        print(gcf, ['figs/noise3d/' metrics_labels{rep_id} '.pdf'], '-dpdf', '-r0');        
    end

end