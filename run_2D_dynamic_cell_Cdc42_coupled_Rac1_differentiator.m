 clear; clc;

rng('shuffle');

numt = 5000000; %number of time steps

dt = 0.001; %time step
dx = 0.02; %spatial step

S = [1000, 1000]; %size of the large canvas for cell simulation

mask_file = 'E:\Asef_Cdc42_Rac1_model\Rac1_regulator_coupled_Cdc42\2D_dynamic_cell_differentiator\Cell 3 Rac1 ratio_mask.mat';
load(mask_file);

Im_L = zeros(S);
Im_L(S(1)/2 - size(mask, 1)/2:S(1)/2 + size(mask, 1)/2 - 1, S(2)/2 - size(mask, 2)/2:S(2)/2 + size(mask, 2)/2 - 1) = mask;
Im_outline_L = outline_8p(~Im_L);
Im_stimulus_L = zeros(S);
Im_stimulus_L(end/2-120-10:end/2-120+10,end/2+30-10:end/2+30+10) = 1;
Im_stimulus_L(end/2-20-10:end/2-20+10,end/2+65-10:end/2+65+10) = 1;
Im_stimulus_L(end/2+120-10:end/2+120+10,end/2-15-10:end/2-15+10) = 1;
Im_stimulus_L(Im_L==0) = 0;

Im_nucl_L = Im_L - Im_outline_L; %nucleous
Im_pnucl_L = outline_8p(Im_nucl_L); %perinuclear region
Im_rb_L = Im_pnucl_L; %barier for retraction

Rs_init = 0.067885;
R_init = 0.93211;
Is_init = 0.25568;
I_init = 2.7443;
Cs_init = 0.058767;
C_init = 0.94123;
Cs_stimulus = 1;

Rs_L = Rs_init*ones(S).*Im_L;
R_L = R_init*ones(S).*Im_L;
Is_L = Is_init*ones(S).*Im_L;
I_L = I_init*ones(S).*Im_L;
Cs_L = Cs_init*ones(S).*Im_L;
C_L = C_init*ones(S).*Im_L;
Cs_L(Im_stimulus_L==1) = Cs_stimulus;
Cs_overhead = (Cs_stimulus-Cs_init)*sum(Im_stimulus_L);
C_L(C_L>0) = C_init - Cs_overhead/sum(Im_L);


k1 = 0.005;
gamma1 = 1.5;
alpha1 = 0.85;
K1_basal = 2.1;
K1_edge = 1.4;
K1_L = K1_basal*ones(size(Im_L));
K1_L(Im_outline_L==1) = K1_edge;
beta1 = 0.5;
n1 = 2;

k2 = 0.1;

k3 = 0.00001;
gamma3 = 0.2;
n3 = 2;

k4 = 0.01;

k5 = 0.05;
gamma5 = 1.5;
alpha5 = 0.1;
K5_basal = 0.5;
K5_edge = 0.45;
K5_L = K5_basal*ones(size(Im_L));
K5_L(Im_outline_L==1) = K5_edge;
n5 = 2;
k6 = 1.15;

DRs = 0.001;
DR = 0.1/3;
DIs = 0.003;
DI = 0.1/3;
DCs = 0.001;
DC = 0.1/3;

alpha_V = 0.001; %parameter for the volume factor in the model of moving cell (sharpnes of the sigmoid function, that describes deviation from the initial value of volume)
V0 = sum(Im_L(:))*1; %initial volume (area) of the cell
beta_V = 0.5; %parameter for the volume factor in the model of moving cell (scales pobability of protrusion and retraction from beta_V to 1-beta_V, introduced to make the volume factor weaker then the activ factor)
gamma_V = 0.5;

alpha_A = 100; %parameter for the actin factor in the model of moving cell (sharpnes of the sigmoid function, that describes deviation from the "default" value of A)
A_act = 0.1; %parameter for the actin factor in the model of moving cell ("default" value of A starting from which cell begins to protrude)
%parameters for the actin factor in the model of moving cell (increase/decrease of protrusion/retraction probability with the respect of basal level when A increases)
%see figure for explanation
beta_A = 0;
gamma_A = 0.3;

g = 2; %g parameter for the geometric factor in the model of moving cell (sensitivity to the curvature)
k = 3; %k parameter for the geometric factor in the model of moving cell (smoothnes of boundary)

mov_T = 30; %period of cell outline update
save_T = 1000; %period of saving cell state

crop_d = 2; %width of the outline while cropping the frame adound cell

n_alf_R = 5; %noise value
n_alf_I = 5;
n_alf_C = 5;

[Im, i_start, i_end, j_start, j_end] = crop_frame(Im_L, crop_d);
Im_outline = Im_outline_L(i_start:i_end, j_start:j_end);
Im_nucl = Im_nucl_L(i_start:i_end, j_start:j_end);
Im_rb = Im_rb_L(i_start:i_end, j_start:j_end);
Rs = Rs_L(i_start:i_end, j_start:j_end);
R = R_L(i_start:i_end, j_start:j_end);
Is = Is_L(i_start:i_end, j_start:j_end);
I = I_L(i_start:i_end, j_start:j_end);
Cs = Cs_L(i_start:i_end, j_start:j_end);
C = C_L(i_start:i_end, j_start:j_end);
K1 = K1_L(i_start:i_end, j_start:j_end);
K5 = K5_L(i_start:i_end, j_start:j_end);

U = U_matrix(Im); %supplementary matrix for laplacian computation

s = size(Rs);

save_fold = 'E:\Asef_Cdc42_Rac1_model\Rac1_regulator_coupled_Cdc42\2D_dynamic_cell_differentiator';
subfold = strcat('gamma1_', num2str(gamma1), '_alpha1_', num2str(alpha1), ...
    '_gamma5_', num2str(gamma5), '_alpha5_', num2str(alpha5), ...
    '_gamma_A_', num2str(gamma_A), '_A_act_', num2str(A_act), ...
    '_beta_A_', num2str(beta_A));

save_fold = fullfile(save_fold, subfold);
mkdir(save_fold);
mkdir(fullfile(save_fold, 'Rs_L'));
mkdir(fullfile(save_fold, 'Cs_L'));
mkdir(fullfile(save_fold, 'mat'));
mkdir(fullfile(save_fold, 'colored'));

save(fullfile(save_fold, 'parameters.mat'), ...
    'numt', 'dt', 'dx', 'S', 'Im_L', 'Im_outline_L', 'Im_stimulus_L', ...
    'Im_nucl_L', 'Im_rb_L', 'Rs_init', 'R_init', 'Is_init', 'I_init', 'Cs_init', ...
    'C_init', 'Cs_stimulus', 'Rs_L', 'R_L', 'Is_L', 'I_L', 'Cs_L', 'C_L', ...
    'Cs_overhead', 'k1', 'gamma1', 'alpha1', 'K1_basal', 'K1_edge', 'K1_L', ...
    'beta1', 'n1', 'k2', 'k3', 'gamma3', 'n3', 'k4', 'k5', 'gamma5', 'alpha5', ...
    'K5_basal', 'K5_edge', 'K5_L', 'n5', 'k6', 'DRs', 'DR', 'DIs', 'DI', ...
    'DCs', 'DC', 'alpha_V', 'V0', 'beta_V', 'gamma_V', 'alpha_A', 'A_act', ...
    'beta_A', 'gamma_A', 'g', 'k', 'mov_T', 'save_T', 'crop_d', 'n_alf_R', ...
    'n_alf_I', 'n_alf_C', 'Im', 'i_start', 'i_end', 'j_start', 'j_end', ...
    'Im', 'Im_outline', 'Im_nucl', 'Im_rb', 'Rs', 'R', 'Is', 'I', 'Cs', ...
    'C', 'K1', 'K5', 'U', 's');

V_vals = linspace(0.5*V0, 2*V0, 1000);
A_vals = 0:0.01:0.25;

%protrusion
vol_p_prot = 0.5 + beta_V - (beta_V+gamma_V)*(1-1./(1+exp(alpha_V*(V_vals-V0))));
act_p_prot = 0.5 - beta_A + (beta_A+gamma_A)*(1-1./(1+exp(alpha_A*(A_vals-A_act))));

%retraction
vol_p_ret = 0.5 - beta_V + (beta_V+gamma_V)*(1-1./(1+exp(alpha_V*(V_vals-V0))));
act_p_ret = 0.5 + beta_A - (beta_A+gamma_A)*(1-1./(1+exp(alpha_A*(A_vals-A_act))));

fig = figure('Position', [50 50 1000 400]);
subplot(1,2,1);
hold on;
grid off;
box on;
ylim([0, 1]);
plot(A_vals, act_p_prot, 'Color', [1 0 0], 'LineWidth', 3);
plot(A_vals, act_p_ret, 'Color', [0 1 0], 'LineWidth', 3);
xlim([min(A_vals), max(A_vals)]);
legend({'protrusion','retraction'});
xlabel('A');
set(gca, 'LineWidth', 2);
title('actin factor');

subplot(1,2,2);
hold on;
grid off;
box on;
plot(V_vals, vol_p_prot, 'Color', [1 0 0], 'LineWidth', 3);
plot(V_vals, vol_p_ret, 'Color', [0 1 0], 'LineWidth', 3);
xlim([min(V_vals), max(V_vals)]);
ylim([0, 1]);
legend({'protrusion','retraction'});
xlabel('V');
set(gca, 'LineWidth', 2);
title('volume factor');

saveas(fig, fullfile(save_fold, 'weigthts_plot.png'));
% close(fig);


fig = figure('Position', [50 50 1200 800]);
for j = 1:numt
    
    %protrude-retract
    if mod(j, mov_T) == 0
        [Im, Rs, R, Is, I, Cs, C] = protrude_extrapolate_diff(Im, Rs, R, Is, I, Cs, C, F1, g, k, V0, alpha_V, beta_V, gamma_V, A_act, alpha_A, beta_A, gamma_A);
        [Im, Rs, R, Is, I, Cs, C] = retract_reduce_diff(Im_rb, Im, Rs, R, Is, I, Cs, C, F1, g, k, V0, alpha_V, beta_V, gamma_V, A_act, alpha_A, beta_A, gamma_A);
        Im_L = return_frame_to_canvas(Im, S, i_start, i_end, j_start, j_end);
        Rs_L = return_frame_to_canvas(Rs, S, i_start, i_end, j_start, j_end);
        R_L = return_frame_to_canvas(R, S, i_start, i_end, j_start, j_end);
        Is_L = return_frame_to_canvas(Is, S, i_start, i_end, j_start, j_end);
        I_L = return_frame_to_canvas(I, S, i_start, i_end, j_start, j_end);
        Cs_L = return_frame_to_canvas(Cs, S, i_start, i_end, j_start, j_end);
        C_L = return_frame_to_canvas(C, S, i_start, i_end, j_start, j_end);
        [Im, i_start, i_end, j_start, j_end] = crop_frame(Im_L, crop_d);
        Im_outline_L = outline_8p(~Im_L);
        Im_nucl = Im_nucl_L(i_start:i_end, j_start:j_end);
        Im_rb = Im_rb_L(i_start:i_end, j_start:j_end);
        Rs = Rs_L(i_start:i_end, j_start:j_end);
        R = R_L(i_start:i_end, j_start:j_end);
        Is = Is_L(i_start:i_end, j_start:j_end);
        I = I_L(i_start:i_end, j_start:j_end);
        Cs = Cs_L(i_start:i_end, j_start:j_end);
        C = C_L(i_start:i_end, j_start:j_end);
        U = U_matrix(Im);
        K1_L = K1_basal*ones(size(Im_L));
        K1_L(Im_outline_L==1) = K1_edge;
        K1 = K1_L(i_start:i_end, j_start:j_end);
        K5_L = K5_basal*ones(size(Im_L));
        K5_L(Im_outline_L==1) = K5_edge;
        K5 = K5_L(i_start:i_end, j_start:j_end);
        
        im1 = Rs;
        im1(Im==1) = (im1(Im==1)-min(im1(Im==1)))/(max(im1(Im==1))-min(im1(Im==1)));
        
        im2 = Is;
        im2(Im==1) = (im2(Im==1)-min(im2(Im==1)))/(max(im2(Im==1))-min(im2(Im==1)));
        
        im3 = Cs;
        im3(Im==1) = (im3(Im==1)-min(im3(Im==1)))/(max(im3(Im==1))-min(im3(Im==1)));
        
        clf;
        subplot(2,3,1);
        hold on;
        colormap hot;
        axis ij;
        axis off;
        axis equal;
        imagesc(Rs);
        colorbar;
        title('Rs');
        
        subplot(2,3,4);
        hold on;
        colormap hot;
        axis ij;
        axis off;
        axis equal;
        imagesc(im1);
        colorbar;
        title('Rs scaled');
        
        subplot(2,3,2);
        hold on;
        colormap hot;
        axis ij;
        axis off;
        axis equal;
        imagesc(Is);
        colorbar;
        title('Is');
        
        subplot(2,3,5);
        hold on;
        colormap hot;
        axis ij;
        axis off;
        axis equal;
        imagesc(im2);
        colorbar;
        title('Is scaled');
        
        subplot(2,3,3);
        hold on;
        colormap hot;
        axis ij;
        axis off;
        axis equal;
        imagesc(Cs);
        colorbar;
        title('Cs');
        
        subplot(2,3,6);
        hold on;
        colormap hot;
        axis ij;
        axis off;
        axis equal;
        imagesc(im3);
        colorbar;
        title('Cs scaled');
        drawnow;
    end
    
    noise_R = randn(size(Rs))*n_alf_R;
    noise_I = randn(size(Is))*n_alf_I;
    noise_C = randn(size(Cs))*n_alf_C;
    
    Rs_new = Rs + (DRs*laplacian_DT(Rs,dx,U) + noise_R)*dt;
    Rs_new(Im==0)=0;
    Rs_neg = Rs_new;
    Rs_neg(Rs_neg > 0) = 0;
    Rs_new(Rs_new < 0) = 0;
    
    R_new = R + (DR*laplacian_DT(R,dx,U) - noise_R)*dt;
    R_new(Im==0)=0;
    R_neg = R_new;
    R_neg(R_neg > 0) = 0;
    R_new(R_new < 0) = 0;
    
    Rs = Rs_new + R_neg;
    R = R_new + Rs_neg;
    
    Is_new = Is + (DIs*laplacian_DT(Is,dx,U) + noise_I)*dt;
    Is_new(Im==0)=0;
    Is_neg = Is_new;
    Is_neg(Is_neg > 0) = 0;
    Is_new(Is_new < 0) = 0;
    
    I_new = I + (DI*laplacian_DT(I,dx,U) - noise_I)*dt;
    I_new(Im==0)=0;
    I_neg = I_new;
    I_neg(I_neg > 0) = 0;
    I_new(I_new < 0) = 0;
    
    Is = Is_new + I_neg;
    I = I_new + Is_neg;
    
    Cs_new = Cs + (DCs*laplacian_DT(Cs,dx,U) + noise_C)*dt;
    Cs_new(Im==0)=0;
    Cs_neg = Cs_new;
    Cs_neg(Cs_neg > 0) = 0;
    Cs_new(Cs_new < 0) = 0;
    
    C_new = C + (DC*laplacian_DT(C,dx,U) - noise_C)*dt;
    C_new(Im==0)=0;
    C_neg = C_new;
    C_neg(C_neg > 0) = 0;
    C_new(C_new < 0) = 0;
    
    Cs = Cs_new + C_neg;
    C = C_new + Cs_neg;
    
    %react
    F1 = (k1 + (gamma1+alpha1*Cs.^n1).*Rs.^n1./(K1.^n1+beta1*Is.^n1+Rs.^n1)).*R - k2.*Rs;
    F2 = (k3 + gamma3*Rs.^n3).*I-k4.*Is;
    F3 = (k5 + (gamma5+alpha5*Rs.^n5).*Cs.^n5./(K5.^n5+Cs.^n5)).*C - k6.*Cs;
    Rs = Rs + F1*dt;
    R = R - F1*dt;
    Is = Is + F2*dt;
    I = I - F2*dt;
    Cs = Cs + F3*dt;
    C = C - F3*dt;
    
    %save
    if mod(j, save_T) == 0
        disp(j/save_T);
        
        imwrite(Rs_L, fullfile(save_fold, 'Rs_L', strcat(num2str(j/save_T), '.tif')));
        imwrite(Cs_L, fullfile(save_fold, 'Cs_L', strcat(num2str(j/save_T), '.tif')));
        save(fullfile(save_fold, 'mat', strcat(num2str(j/save_T), '.mat')), ...
            'Rs_L', 'R_L', 'Is_L', 'I_L', 'Cs_L', 'C_L', 'Im_L');
        saveas(fig, fullfile(save_fold, 'colored', strcat(num2str(j/save_T), '.png')));
    end
    
end
