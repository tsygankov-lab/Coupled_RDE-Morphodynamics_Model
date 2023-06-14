clear; clc;

rng('shuffle');

numt = 5000000; %number of time steps

dt = 0.001; %time step
dx = 0.02; %spatial step

P_total_conc = 1;
I_total_conc = 3;
C_total_conc = 1;
R_total_conc = 1;

S = [1000, 1000]; %size of the large canvas for cell simulation
cell_R = 80; %radius of the cell

Im_L = generate_regular_polyhedron(S, S/2, cell_R, 500);
Im_outline = outline_8p(~Im_L); %cell_outline

Ps_init = 0;
Is_init = 1;
Cs_init = 0;
Rs_init = 0;

Ps_L = Ps_init*ones(S).*Im_L;
P_L = (P_total_conc-Ps_init)*ones(S).*Im_L;
Is_L = Is_init*ones(S).*Im_L;
I_L = (I_total_conc-Is_init)*ones(S).*Im_L;
Cs_L = Cs_init*ones(S).*Im_L;
C_L = (C_total_conc-Cs_init)*ones(S).*Im_L;
Rs_L = Rs_init*ones(S).*Im_L;
R_L = (R_total_conc-Rs_init)*ones(S).*Im_L;

k1 = 0.005;
gamma1 = 2;
K1_basal = 2.1;
K1_edge = 1.4;
K1_L = K1_basal*ones(size(Im_L));
K1_L(Im_outline==1) = K1_edge;
beta1 = 0.5;
n1 = 2;
k2 = 0.1;

k3 = 0.00001;
gamma3 = 0.2;
n3 = 2;
k4 = 0.01;

k5 = 0.005;
gamma5 = 0.1;
alpha5 = 2;
K5_basal = 2.1;
K5_edge = 1.4;
K5_L = K5_basal*ones(size(Im_L));
K5_L(Im_outline==1) = K5_edge;
n5 = 2;
k6 = 0.1;

k7 = 0.005;
gamma7 = 0.1;
alpha7 = 2;
K7_basal = 2.1;
K7_edge = 1.4;
K7_L = K7_basal*ones(size(Im_L));
K7_L(Im_outline==1) = K7_edge;
n7 = 2;
k8 = 0.1;

DPs = 0.001;
DP = 0.1/3;
DIs = 0.003;
DI = 0.1/3;
DCs = 0.001;
DC = 0.1/3;
DRs = 0.001;
DR = 0.1/3;

alpha_V = 0.001; %parameter for the volume factor in the model of moving cell (sharpnes of the sigmoid function, that describes deviation from the initial value of volume)
V0 = sum(Im_L(:))*1; %initial volume (area) of the cell
beta_V = 0.5; %parameter for the volume factor in the model of moving cell (scales pobability of protrusion and retraction from beta_V to 1-beta_V, introduced to make the volume factor weaker then the activ factor)
gamma_V = 0.5;

alpha_A = 100; %parameter for the actin factor in the model of moving cell (sharpnes of the sigmoid function, that describes deviation from the "default" value of A)
A_act = 0.1; %parameter for the actin factor in the model of moving cell ("default" value of A starting from which cell begins to protrude)
%parameters for the actin factor in the model of moving cell (increase/decrease of protrusion/retraction probability with the respect of basal level when A increases)
%see figure for explanation
beta_A = 0;
gamma_A = 0.1;
%gamma_A = 0.2;
%gamma_A = 0.15;

g = 2; %g parameter for the geometric factor in the model of moving cell (sensitivity to the curvature)
k = 3; %k parameter for the geometric factor in the model of moving cell (smoothnes of boundary)

sub_iter_N = 1; %number of iterations for sepparate reaction and diffusion steps (due to operator splitting)
mov_T = 50; %period of cell outline update
save_T = 1000; %period of saving cell state

crop_d = 2; %width of the outline while cropping the frame adound cell

n_alf = 5; %noise value

[Im, i_start, i_end, j_start, j_end] = crop_frame(Im_L, crop_d);
Ps = Ps_L(i_start:i_end, j_start:j_end);
P = P_L(i_start:i_end, j_start:j_end);
Is = Is_L(i_start:i_end, j_start:j_end);
I = I_L(i_start:i_end, j_start:j_end);
Cs = Cs_L(i_start:i_end, j_start:j_end);
C = C_L(i_start:i_end, j_start:j_end);
Rs = Rs_L(i_start:i_end, j_start:j_end);
R = R_L(i_start:i_end, j_start:j_end);
K1 = K1_L(i_start:i_end, j_start:j_end);
K5 = K5_L(i_start:i_end, j_start:j_end);
K7 = K7_L(i_start:i_end, j_start:j_end);

U = U_matrix(Im); %supplementary matrix for laplacian computation

s = size(Rs);

save_fold = 'E:\Asef_Cdc42_Rac1_model\PIP_regulator_Cdc42_Rac1\2D_dynamic_model_differentiator';
subfold = strcat('cell_R_', num2str(cell_R), '_k5_', num2str(k5), ...
    '_gamma5_', num2str(gamma5), '_alpha5_', num2str(alpha5), ...
    '_K5_basal_', num2str(K5_basal), '_K5_edge_', num2str(K5_edge), ...
    '_k6_', num2str(k6), '_gamma7_', num2str(gamma7), ...
    '_alpha7_', num2str(alpha7), '_K7_basal_', num2str(K7_basal), ...
    '_K7_edge_', num2str(K7_edge), '_k8_', num2str(k8), ...
    '_alpha_A_', num2str(alpha_A), '_A_act_', num2str(A_act), ...
    '_gamma_A_', num2str(gamma_A));
subfold = 'simultaneous';

save_fold = fullfile(save_fold, subfold);
mkdir(save_fold);
mkdir(fullfile(save_fold, 'Ps_L'));
mkdir(fullfile(save_fold, 'Cs_L'));
mkdir(fullfile(save_fold, 'Rs_L'));
mkdir(fullfile(save_fold, 'mat'));
mkdir(fullfile(save_fold, 'colored'));

save(fullfile(save_fold, 'parameters.mat'), ...
    'numt', 'dt', 'dx', 'P_total_conc', 'I_total_conc', 'C_total_conc', ...
    'R_total_conc', 'S', 'cell_R', 'Im_L', 'Im_outline', 'Ps_init', ...
    'Is_init', 'Cs_init', 'Rs_init', 'Ps_L', 'P_L', 'Is_L', 'I_L', ...
    'Cs_L', 'C_L', 'Rs_L', 'R_L', 'k1', 'gamma1', 'K1_basal', 'K1_edge', ...
    'K1_L', 'beta1', 'n1', 'k2', 'k3', 'gamma3', 'n3', 'k4', 'k5', 'gamma5', ...
    'alpha5', 'K5_basal', 'K5_edge', 'K5_L', 'n5', 'k6', 'k7', 'gamma7', ...
    'alpha7', 'K7_basal', 'K7_edge', 'K7_L', 'n7', 'k8', 'DPs', 'DP', 'DIs', ...
    'DI', 'DCs', 'DC', 'DRs', 'DR', 'alpha_V', 'V0', 'beta_V', 'gamma_V', ...
    'alpha_A', 'A_act', 'beta_A', 'gamma_A', 'g', 'k', 'sub_iter_N', 'mov_T', ...
    'save_T', 'crop_d', 'n_alf', 'Im', 'i_start', 'i_end', 'j_start', 'j_end', ...
    'Ps', 'P', 'Is', 'I', 'Cs', 'C', 'Rs', 'R', 'K1', 'K5', 'K7', 'U', 's');

V_vals = linspace(0.5*V0, 2*V0, 1000);
%R_vals = -0.01:0.001:0.04;
R_vals = -0.01:0.001:0.5;

%protrusion
vol_p_prot = 0.5 + beta_V - (beta_V+gamma_V)*(1-1./(1+exp(alpha_V*(V_vals-V0))));
act_p_prot = 0.5 - beta_A + (beta_A+gamma_A)*(1-1./(1+exp(alpha_A*(R_vals-A_act))));

%retraction
vol_p_ret = 0.5 - beta_V + (beta_V+gamma_V)*(1-1./(1+exp(alpha_V*(V_vals-V0))));
act_p_ret = 0.5 + beta_A - (beta_A+gamma_A)*(1-1./(1+exp(alpha_A*(R_vals-A_act))));

fig = figure('Position', [50 50 1000 400]);
subplot(1,2,1);
hold on;
grid off;
box on;
ylim([0, 1]);
plot(R_vals, act_p_prot, 'Color', [1 0 0], 'LineWidth', 3);
plot(R_vals, act_p_ret, 'Color', [0 1 0], 'LineWidth', 3);
xlim([min(R_vals), max(R_vals)]);
legend({'protrusion','retraction'});
xlabel('R');
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
        [Im, Ps, P, Is, I, Cs, C, Rs, R] = protrude_extrapolate_diff(Im, Ps, P, Is, I, Cs, C, Rs, R, F4, g, k, V0, alpha_V, beta_V, gamma_V, A_act, alpha_A, beta_A, gamma_A);
        [Im, Ps, P, Is, I, Cs, C, Rs, R] = retract_reduce_diff(Im, Ps, P, Is, I, Cs, C, Rs, R, F4, g, k, V0, alpha_V, beta_V, gamma_V, A_act, alpha_A, beta_A, gamma_A);
        Im_L = return_frame_to_canvas(Im, S, i_start, i_end, j_start, j_end);
        Ps_L = return_frame_to_canvas(Ps, S, i_start, i_end, j_start, j_end);
        P_L = return_frame_to_canvas(P, S, i_start, i_end, j_start, j_end);
        Is_L = return_frame_to_canvas(Is, S, i_start, i_end, j_start, j_end);
        I_L = return_frame_to_canvas(I, S, i_start, i_end, j_start, j_end);
        Cs_L = return_frame_to_canvas(Cs, S, i_start, i_end, j_start, j_end);
        C_L = return_frame_to_canvas(C, S, i_start, i_end, j_start, j_end);
        Rs_L = return_frame_to_canvas(Rs, S, i_start, i_end, j_start, j_end);
        R_L = return_frame_to_canvas(R, S, i_start, i_end, j_start, j_end);
        [Im, i_start, i_end, j_start, j_end] = crop_frame(Im_L, crop_d);
        Ps = Ps_L(i_start:i_end, j_start:j_end);
        P = P_L(i_start:i_end, j_start:j_end);
        Is = Is_L(i_start:i_end, j_start:j_end);
        I = I_L(i_start:i_end, j_start:j_end);
        Cs = Cs_L(i_start:i_end, j_start:j_end);
        C = C_L(i_start:i_end, j_start:j_end);
        Rs = Rs_L(i_start:i_end, j_start:j_end);
        R = R_L(i_start:i_end, j_start:j_end);
        U = U_matrix(Im);
        Im_outline = outline_8p(~Im_L);
        K1_L = K1_basal*ones(size(Im_L));
        K1_L(Im_outline==1) = K1_edge;
        K1 = K1_L(i_start:i_end, j_start:j_end);
        
        K5_L = K5_basal*ones(size(Im_L));
        K5_L(Im_outline==1) = K5_edge;
        K5 = K5_L(i_start:i_end, j_start:j_end);
        
        K7_L = K7_basal*ones(size(Im_L));
        K7_L(Im_outline==1) = K7_edge;
        K7 = K7_L(i_start:i_end, j_start:j_end);
        
        im1 = Ps;
        im1(Im==1) = (im1(Im==1)-min(im1(Im==1)))/(max(im1(Im==1))-min(im1(Im==1)));
        
        im2 = Cs;
        im2(Im==1) = (im2(Im==1)-min(im2(Im==1)))/(max(im2(Im==1))-min(im2(Im==1)));
        
        im3 = Rs;
        im3(Im==1) = (im3(Im==1)-min(im3(Im==1)))/(max(im3(Im==1))-min(im3(Im==1)));
        
        clf;
        subplot(2,3,1);
        hold on;
        colormap hot;
        axis ij;
        axis off;
        axis equal;
        imagesc(Ps);
        colorbar;
        title('Ps');
        
        subplot(2,3,4);
        hold on;
        colormap hot;
        axis ij;
        axis off;
        axis equal;
        imagesc(im1);
        colorbar;
        title('Ps scaled');
        
        subplot(2,3,2);
        hold on;
        colormap hot;
        axis ij;
        axis off;
        axis equal;
        imagesc(Cs);
        colorbar;
        title('Cs');
        
        subplot(2,3,5);
        hold on;
        colormap hot;
        axis ij;
        axis off;
        axis equal;
        imagesc(im2);
        colorbar;
        title('Cs scaled');
        
        subplot(2,3,3);
        hold on;
        colormap hot;
        axis ij;
        axis off;
        axis equal;
        imagesc(Rs);
        colorbar;
        title('Rs');
        
        subplot(2,3,6);
        hold on;
        colormap hot;
        axis ij;
        axis off;
        axis equal;
        imagesc(im3);
        colorbar;
        title('Rs scaled');
        drawnow;
    end
    
    %diffuse with noise
    for i = 1:sub_iter_N
        noise_P = randn(size(Ps))*n_alf;
        noise_I = randn(size(Is))*n_alf;
        noise_C = randn(size(Cs))*n_alf;
        noise_R = randn(size(Rs))*n_alf;
        
        Ps_new = Ps + (DPs*laplacian_DT(Ps,dx,U) + noise_P)*dt;
        Ps_new(Im==0)=0;
        Ps_neg = Ps_new;
        Ps_neg(Ps_neg > 0) = 0;
        Ps_new(Ps_new < 0) = 0;
        
        P_new = P + (DP*laplacian_DT(P,dx,U) - noise_P)*dt;
        P_new(Im==0)=0;
        P_neg = P_new;
        P_neg(P_neg > 0) = 0;
        P_new(P_new < 0) = 0;
        
        Ps = Ps_new + P_neg;
        P = P_new + Ps_neg;
        
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
    end
    
    %react
    for i = 1:sub_iter_N
        F1 = (k1 + gamma1*Ps.^n1./(K1.^n1 + beta1*Is.^n1 + Ps.^n1)).*P - k2.*Ps;
        F2 = (k3 + gamma3*Ps.^n3).*I - k4*Is;
        F3 = (k5 + (gamma5 + alpha5*Ps.^n5).*Cs.^n5./(K5.^n5 + Cs.^n5)).*C - k6*Cs;
        F4 = (k7 + (gamma7 + alpha7*Ps.^n7).*Rs.^n7./(K7.^n7 + Rs.^n7)).*R - k8*Rs;
        
        Ps = Ps + F1*dt;
        P = P - F1*dt;
        Is = Is + F2*dt;
        I = I - F2*dt;
        Cs = Cs + F3*dt;
        C = C - F3*dt;
        Rs = Rs + F4*dt;
        R = R - F4*dt;
    end
   
    %save
    if mod(j, save_T) == 0
        disp(j/save_T);
        imwrite(Ps, fullfile(save_fold, 'Ps_L', strcat(num2str(j/save_T), '.tif')));
        imwrite(Cs, fullfile(save_fold, 'Cs_L', strcat(num2str(j/save_T), '.tif')));
        imwrite(Rs, fullfile(save_fold, 'Rs_L', strcat(num2str(j/save_T), '.tif')));
        save(fullfile(save_fold, 'mat', strcat(num2str(j/save_T), '.mat')), ...
            'Ps_L', 'P_L', 'Is_L', 'I_L', 'Cs_L', 'C_L', 'Rs_L', 'R_L', 'Im_L');
        saveas(fig, fullfile(save_fold, 'colored', strcat(num2str(j/save_T), '.png')));
    end
end

