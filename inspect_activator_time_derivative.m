clear; clc;

inp_fold = 'E:\Asef_Cdc42_Rac1_model\ruffling_differentiator\2D_dynamic_cell\cell_R_90_K1_edge_1.4_alpha_A_30_A_act_0.7_gamma_A_0.1';
inp_fold = 'E:\Asef_Cdc42_Rac1_model\ruffling_differentiator\2D_dynamic_cell\cell_R_90_K1_edge_1.4_alpha_A_50_A_act_0.1_gamma_A_0.2';
load(fullfile(inp_fold, 'parameters.mat'));
load(fullfile(inp_fold, 'data_inspection.mat'));

d = 10;

k4 = 0.01;


F1_max = -Inf;
F1_min = Inf;

mkdir(fullfile(inp_fold, 'F1'));

fig = figure();
for i = 1:3450
    disp(i);
    load(fullfile(inp_fold, 'mat', strcat(num2str(i), '.mat')));
    
    [Im, i_start, i_end, j_start, j_end] = crop_frame(Im_L, crop_d);
    Im_outline = outline_8p(~Im_L);
    K1_L = K1_basal*ones(size(Im_L));
    K1_L(Im_outline==1) = K1_edge;
    
    As = As_L(j1-d:j2+d, i1-d:i2+d);
    A = A_L(j1-d:j2+d, i1-d:i2+d);
    Cs = Cs_L(j1-d:j2+d, i1-d:i2+d);
    C = C_L(j1-d:j2+d, i1-d:i2+d);
    K1 = K1_L(j1-d:j2+d, i1-d:i2+d);
    s = size(As);
    
    F1 = (k1 + gamma1*As.^n1./(K1.^n1+beta1*Cs.^n1+As.^n1)).*A - k2.*As;
    F2 = (k3 + gamma3*As.^n3).*C-k4.*Cs;
    
    F1_min_i = min(F1(Im==1));
    F1_max_i = max(F1(Im==1));
    
    if F1_min_i < F1_min
        F1_min = F1_min_i;
    end
    if F1_max_i > F1_max
        F1_max = F1_max_i;
    end
    
    clf;
    set(fig, 'Position', [50 50 2*s(2) 2*s(1)]);
    hold on;
    axis ij;
    axis off;
    axis image;
    imagesc(F1);
    colormap(jet);
    colorbar;
    caxis([-0.04, 0.25]);
    drawnow;
    %saveas(fig, fullfile(inp_fold, 'F1', strcat(num2str(i), '.png')));
end
