% Confirming visual lateralization via model selection
% _
% Bayesian Model Selection for orientation pop-out processing
%
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% Date  : 12/05/2017, 05:50


%%% SPECIFY DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% BMS maps
BMS_dir = 'C:\Joram\projects\MACS\DataSets\Bogler_et_al_2013\analyses\model_selection\MS_SAL_01_GLM_III-l_vs_GLM_III-r\';
BMS_str = {strcat(BMS_dir,'GLM_III-l_model_LFM.nii');
           strcat(BMS_dir,'GLM_III-r_model_LFM.nii')};

% localizer masks
mask_dir = 'C:\Joram\projects\MACS\DataSets\Bogler_et_al_2013\analyses\model_selection\_masks_\';
mask_str = {strcat(mask_dir,'V4_left.nii');
            strcat(mask_dir,'V4_right.nii')};

% set bin centers
n_mask = length(mask_str);
n_vols = length(BMS_str);
x_vals = [0.05:0.05:0.95];


%%% PLOT HISTOGRAMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open figure
hf = figure('Name', 'Likeliest Frequencies', 'Color', [1 1 1], 'Position', [50 50 500 1000]);
hl = 'ABC';

% show histograms
for i = 1:length(mask_str)
    
    % load volumes
    [out, V] = MF_analyze('sum(V,2)', [BMS_str; mask_str(i)]);
    mods1 = {'GLM with left OC only', 'GLM with right OC only'};
    mods2 = {'left OC only', 'right OC only'};
    
    % calculate bins
    n = zeros(length(BMS_str),length(x_vals));
    for j = 1:n_vols
        n(j,:) = hist(V(j,~isnan(V(end,:))),x_vals);
    end;
    
    % calculate wins
    WM(:,i) = zeros(length(BMS_str),1);
    [C,I] = max(V(1:n_vols,~isnan(V(end,:))));
    for j = 1:n_vols
        WM(j,i) = sum(I==j)/numel(C);
    end;
    
    % plot histogram
    subplot(3,1,i);
    hsp = bar(x_vals, n', 1.7, 'grouped');
    for j = 1:n_vols
        set(hsp(j),'FaceColor',((3-j)/3)*[0 1 0]);
    end;
    grid on;
    set(gca, 'box', 'on');
    axis([0 1 0 (4.1/3)*max(n(:))]);
    legend(mods1{1},mods1{2}, 'Location', 'NorthEast');
    set(gca,'XTick',[0:0.1:1],'XTickLabel',cellstr(num2str([0:0.1:1]'))');
    xlabel('likeliest frequency', 'FontSize', 12);
    ha = ylabel('number of voxels', 'FontSize', 12);
    if i == 1, hb = title('Voxels from left V4 only', 'FontSize', 24); end;
    if i == 2, hb = title('Voxels from right V4 only', 'FontSize', 24); end;
    ca = get(ha,'Position');
    cb = get(hb,'Position');
    text(ca(1),cb(2),hl(i), 'FontSize', 24, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Bottom');
    
end;

% plot bar plot
subplot(3,1,3);
hold on;
for j = 1:n_vols
    bar([0 2]+j,WM(j,:)*100,0.3,'FaceColor',((3-j)/3)*[0 1 0]);
end;
grid on;
set(gca, 'box', 'on');
axis([(1-0.5) (2*n_vols+0.5) 0 100]);
set(gca,'XTick',[1:2*n_vols],'XTickLabel',[mods2 mods2]);
xlabel('left V4 only                            right V4 only', 'FontSize',  12);
h3a = ylabel('proportion of voxels', 'FontSize', 12);
h3b = title('V4 model preferences', 'FontSize', 24);
c3a = (1-0.5) - (0-ca(1))/(1-0)*((2*n_vols+0.5)-(1-0.5));
c3b = get(h3b,'Position');
text(c3a,c3b(2),hl(3), 'FontSize', 24, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Bottom');


% save figure
% saveas(hf,'LF_hists.fig');
% save('LF_hists.mat','x_vals');