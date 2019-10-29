clear all;
close all;
clc;

fs = 2;
fsu = 4;
fss = 2;
bw = fsu;
subjects = [1, 5, 8, 9, 12, 16];
fsy = 2; fsu = 4;
count = 0;
dddata = load('..\UT_Dallas_data\s.mat');
mini_emo_start_i = 4;
cog_stress_start_i = 5;
relax_start_i = 6;
emo_stress_start_i = 7;
final_relax_start_i = 8;

addpath('..\Brewer_Map_Color_library');
for sub = subjects
    
    count = count+1;
    load(['result_stress_',num2str(sub),'.mat']);
    u = data.uj_whole(1:end-end_);
    u = u(:)';
    u(u<=0) = NaN;

    fs = 8;
    s = dddata.s;
    xp = s(sub).y(mini_emo_start_i:final_relax_start_i);

    xp(end) = xp(end) - 1;
    xp = xp - xp(1);
    xp = [xp(1:2) 0 xp(3:end)];
    xp(3) = xp(2) + 3 * 60 * fs;
    xp = xp / bw * 2;


    signal = s(sub).x(s(sub).y(mini_emo_start_i):(s(sub).y(final_relax_start_i) - 1));
    signal = downsample(signal,4);
    ty = (0:(length(signal) - 1)) / fsy;
    tu = (0:(length(u) - 1)) / fsu;
    
    

    fig = figure('units','normalized','outerposition',[0 0 0.49 0.850*0.85]);
    left_color = [.5 .5 0]*0;
    right_color = [0 .5 .5]*0;
    set(fig,'defaultAxesColorOrder',[left_color; right_color]);
    offset = 0.1;
    
    % title('Participant 1')
    
    %% plot stress state
    subplot(513);
    plot(xK, 'linewidth', 1.5,'Color',[54 84 33]/255); hold on;  
    ylabel({'State','z_{j|J}'});
    
    lcl_x = norminv(0.025, xK, sqrt(vK));
    ucl_x = norminv(0.975, xK, sqrt(vK));
    
    plot(lcl_x, 'c', 'linewidth', 0.8); hold on;
    plot(ucl_x, 'c', 'linewidth', 0.8); hold on;

    fill([(1:k) (k:(-1):1)], [lcl_x fliplr(ucl_x)], 'c', 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    
    
    %xlim([0 (N + 1)]); 

    ylim([min(lcl_x)*1.1 (min(lcl_x)+(max(ucl_x)-min(lcl_x))*1.2)]);
    yticks([-6 -4 0 +4]);
    yl = ylim;
    patch([xp(1), xp(2), xp(2), xp(1)]/fsy/2  , [yl(1) yl(1) yl(2) yl(2)], cr, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch([xp(2), xp(3), xp(3), xp(2)]/fsy/2  , [yl(1) yl(1) yl(2) yl(2)], cg, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    patch([xp(3), xp(4), xp(4), xp(3)]/fsy/2  , [yl(1) yl(1) yl(2) yl(2)], cc, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch([xp(4), xp(5), xp(5), xp(4)]/fsy/2  , [yl(1) yl(1) yl(2) yl(2)], cb, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch([xp(5), xp(6), xp(6), xp(5)]/fsy/2 , [yl(1) yl(1) yl(2) yl(2)], cv, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(gca,'xticklabel',[]);
    % title('Participant 6', 'FontWeight', 'Normal');
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    xlim([0 (N + 1)]);
    plot(find(n == 0), yl(2) * ones(length(find(n == 0))), 's', 'MarkerFaceColor', cv,'Color',cv); hold on; grid;
    plot(find(n == 1), yl(2) * ones(length(find(n == 1))), 's','MarkerFaceColor', cg, 'Color',cg);

    
    

    
     %% plot probability of a peak       
    subplot(514);

    plot(lclK, 'g', 'linewidth', 1.5); hold on;
    plot(uclK, 'g', 'linewidth', 1.5); hold on;
    plot(pK, 'k', 'linewidth', 1.5); hold on;
    plot([0, (N + 1)], [chance_prob, chance_prob], 'k--', 'linewidth', 1);

    fill([(1:k) (k:(-1):1)], [uclK(1:k) lclK(k:(-1):1)], 'g', ...
        'EdgeColor', 'none', 'FaceAlpha', 0.15);

    ylabel({'Probability','p_{j|J}'});
    xlim([0 (N + 1)]); 
    % xlim([0 600]);
    % ylim([0 1.1]);
    yl = 1.1 * ylim;
    ylim(yl);
    patch([xp(1), xp(2), xp(2), xp(1)]  / fsy /2, [0 0 yl(2) yl(2)], cr, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch([xp(2), xp(3), xp(3), xp(2)]  / fsy /2, [0 0 yl(2) yl(2)], cg, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    patch([xp(3), xp(4), xp(4), xp(3)]  / fsy /2, [0 0 yl(2) yl(2)], cc, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch([xp(4), xp(5), xp(5), xp(4)]  / fsy /2, [0 0 yl(2) yl(2)], cb, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch([xp(5), xp(6), xp(6), xp(5)]  / fsy /2, [0 0 yl(2) yl(2)], cv, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    set(gca,'xticklabel',[]);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    plot(find(n == 0), yl(2) * ones(length(find(n == 0))), 's', 'MarkerFaceColor', cv,'Color',cv); hold on; grid;
    plot(find(n == 1), yl(2) * ones(length(find(n == 1))), 's','MarkerFaceColor', cg, 'Color',cg);
    hold on; 
    
    breaks = ceil(xp / fsu);
    pos = xK((breaks(1) + 1):breaks(4));
    neg = xK((breaks(4) + 1):breaks(5));
    [fpr, tpr, T, auc, opt] = perfcurve([ones(1, length(pos)) zeros(1, length(neg))], [pos neg], 1);
    
   
    % for median based threshold high arousal index (HAI)
    io_certainty = 1 - normcdf(prctile(xK, 50) * ones(1, length(xK)), xK, sqrt(vK));
    

    %% plot high arousal index (HAI)
    subplot(515);
    v = [0 0.9; ((N + 1) * bw / fsy) 0.9; ((N + 1) * bw / fsy) 1; 0 1];
    c = [1 (220 / 255) (220 / 255); 1 (220 / 255) (220 / 255); 1 0 0; 1 0 0];
    f = [1 2 3 4];
     patch('Faces',f,'Vertices',v, 'FaceVertexCData', c, 'FaceColor','interp', ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5); hold on;
    ylim([0 1]); grid; xlim([0 (N + 1)]);
    ylabel({'High','Arousal','Index'});
    xlabel('Time (seconds)');
    plot(io_certainty, 'linewidth', 1.5, 'Color', cv);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);

    u = data.uj_whole(1:end-end_);
    u = u(:)';
    u(u<=0) = NaN;
    
    
    %% plot raw skin conductance data
    subplot(511);
    
    title(['Participant ',num2str(count)]);
    hold on;
    tu = (0:length(u)-1)/fsu;
    %yyaxis right,
    c = brewermap(9,'Accent');hold on;
   % stem(tu, u, '.', 'linewidth', 1.5, 'Markersize', 15, 'Color', c(5,:)); hold on; ylabel('Amplitude');


    hold on; plot(ty, signal, 'b', 'linewidth', 1.5,'Color', c(1,:)*0.6); hold on; ylabel({'Skin','Conductance','(\muS)'});
    yl = 1.1 * ylim;
    ylim(yl);
    set(gca,'xticklabel',[]);

    patch([xp(1), xp(2), xp(2), xp(1)]  / fsy / 2, [0 0 yl(2) yl(2)], cr, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch([xp(2), xp(3), xp(3), xp(2)]  / fsy / 2, [0 0 yl(2) yl(2)], cg, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    patch([xp(3), xp(4), xp(4), xp(3)]  / fsy / 2, [0 0 yl(2) yl(2)], cc, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch([xp(4), xp(5), xp(5), xp(4)]  / fsy / 2, [0 0 yl(2) yl(2)], cb, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch([xp(5), xp(6), xp(6), xp(5)]  / fsy / 2, [0 0 yl(2) yl(2)], cv, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    % xt = get(gca, 'XTick');
    xlim([ty(1) ty(end)])
    set(gca, 'FontSize', 16);
    
    %% plot neural stimuli
    subplot(512);
    
    %title(['Participant ',num2str(count)]);

    tu = (0:length(u)-1)/fsu;
    %yyaxis right,
    c = brewermap(9,'Accent');hold on;
    stem(tu, u, '.', 'linewidth', 1.5, 'Markersize', 15, 'Color', [0 0 0]); hold on; ylabel('Amplitude');
    set(gca,'xticklabel',[]);
    yl = 1.1 * ylim;
    ylim(yl);
    patch([xp(1), xp(2), xp(2), xp(1)]  / fsy / 2, [0 0 yl(2) yl(2)], cr, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch([xp(2), xp(3), xp(3), xp(2)]  / fsy / 2, [0 0 yl(2) yl(2)], cg, 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    patch([xp(3), xp(4), xp(4), xp(3)]  / fsy / 2, [0 0 yl(2) yl(2)], cc, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch([xp(4), xp(5), xp(5), xp(4)]  / fsy / 2, [0 0 yl(2) yl(2)], cb, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    patch([xp(5), xp(6), xp(6), xp(5)]  / fsy / 2, [0 0 yl(2) yl(2)], cv, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
    %yyaxis left, hold on; plot(ty, signal, 'b', 'linewidth', 1.5,'Color', c(1,:)*0.6); hold on; ylabel('SC (\muS)');
    % xt = get(gca, 'XTick');
    xlim([ty(1) ty(end)])
    set(gca, 'FontSize', 16);
    
    
    saveas(gcf,['s',num2str(sub)],'epsc');

end
