clear; clc;

rng('shuffle');

numt = 5000000; %number of time steps

dt = 0.001; %time step
dx = 0.02; %spatial step

A_total_conc = 1;
C_total_conc = 3;

S = [1000, 1000]; %size of the large canvas for cell simulation
cell_R = 90; %radius of the cell

Im_L = generate_regular_polyhedron(S, S/2, cell_R, 500);
Im_outline = outline_8p(~Im_L); %cell_outline

As_init = 0;
Cs_init = 1;

As_L = As_init*ones(S).*Im_L;
A_L = (A_total_conc-As_init)*ones(S).*Im_L;
Cs_L = Cs_init*ones(S).*Im_L;
C_L = (C_total_conc-Cs_init)*ones(S).*Im_L;

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

DAs = 0.001;
DA = 0.1/3;
DCs = 0.003;
DC = 0.1/3;

alpha_V = 0.001; %parameter for the volume factor in the model of moving cell (sharpnes of the sigmoid function, that describes deviation from the initial value of volume)
V0 = sum(Im_L(:))*1; %initial volume (area) of the cell
beta_V = 0.5; %parameter for the volume factor in the model of moving cell (scales pobability of protrusion and retraction from beta_V to 1-beta_V, introduced to make the volume factor weaker then the activ factor)
gamma_V = 0.5;

alpha_A = 50; %parameter for the actin factor in the model of moving cell (sharpnes of the sigmoid function, that describes deviation from the "default" value of A)
A_act = 0.1; %parameter for the actin factor in the model of moving cell ("default" value of A starting from which cell begins to protrude)
%parameters for the actin factor in the model of moving cell (increase/decrease of protrusion/retraction probability with the respect of basal level when A increases)
%see figure for explanation
beta_A = 0;
gamma_A = 0.2;

% alpha_A = 30;
% A_act = 0.7;
% gamma_A = 0.1;

g = 2; %g parameter for the geometric factor in the model of moving cell (sensitivity to the curvature)
k = 3; %k parameter for the geometric factor in the model of moving cell (smoothnes of boundary)

sub_iter_N = 1; %number of iterations for sepparate reaction and diffusion steps (due to operator splitting)
mov_T = 50; %period of cell outline update
save_T = 1000; %period of saving cell state

crop_d = 2; %width of the outline while cropping the frame adound cell

n_alf = 5; %noise value

[Im, i_start, i_end, j_start, j_end] = crop_frame(Im_L, crop_d);
As = As_L(i_start:i_end, j_start:j_end);
A = A_L(i_start:i_end, j_start:j_end);
Cs = Cs_L(i_start:i_end, j_start:j_end);
C = C_L(i_start:i_end, j_start:j_end);
K1 = K1_L(i_start:i_end, j_start:j_end);

U = U_matrix(Im); %supplementary matrix for laplacian computation

s = size(As);

save_fold = 'E:\Asef_Cdc42_Rac1_model\ruffling_differentiator\2D_dynamic_cell';
subfold = strcat('cell_R_', num2str(cell_R), '_K1_edge_', num2str(K1_edge), ...
    '_alpha_A_', num2str(alpha_A), '_A_act_', num2str(A_act), ...
    '_gamma_A_', num2str(gamma_A));

save_fold = fullfile(save_fold, subfold);
mkdir(save_fold);
mkdir(fullfile(save_fold, 'As_L'));
mkdir(fullfile(save_fold, 'mat'));
mkdir(fullfile(save_fold, 'colored'));

save(fullfile(save_fold, 'parameters.mat'), ...
    'numt', 'dt', 'dx', 'A_total_conc', 'C_total_conc', 'S', 'cell_R', ...
    'Im_L', 'Im_outline', ...
    'As_L', 'A_L', 'Cs_L', 'C_L', 'k1', 'gamma1', ...
    'K1_edge', 'K1_basal', 'K1_L', 'K1', 'beta1', ...
    'n1', 'k2', 'k3', 'gamma3', 'n3', 'k4', 'DAs', 'DA', 'DCs', 'DC', 'alpha_V', ...
    'V0', 'beta_V', 'gamma_V', 'alpha_A', 'A_act', 'beta_A', 'gamma_A', ...
    'g', 'k', 'sub_iter_N', 'mov_T', 'save_T', 'crop_d', 'n_alf', 'Im', ...
    'i_start', 'i_end', 'j_start', 'j_end', ...
    'As', 'A', 'Cs', 'C', 'U', 's');

V_vals = linspace(0.5*V0, 2*V0, 1000);
A_vals = -0.1:0.01:0.5;
%A_vals = 0:0.01:2;

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


fig = figure('Position', [50 50 800 800]);
for j = 1:numt
    
    %protrude-retract
    if mod(j, mov_T) == 0
%         %[Im, As, A, Cs, C] = protrude_extrapolate(Im, As, A, Cs, C, g, k, V0, alpha_V, beta_V, gamma_V, A_act, alpha_A, beta_A, gamma_A);
%         %[Im, As, A, Cs, C] = retract_reduce(Im, As, A, Cs, C, g, k, V0, alpha_V, beta_V, gamma_V, A_act, alpha_A, beta_A, gamma_A);
%         [Im, As, A, Cs, C] = protrude_extrapolate_diff(Im, As, A, Cs, C, F1, g, k, V0, alpha_V, beta_V, gamma_V, A_act, alpha_A, beta_A, gamma_A);
%         [Im, As, A, Cs, C] = retract_reduce_diff(Im, As, A, Cs, C, F1, g, k, V0, alpha_V, beta_V, gamma_V, A_act, alpha_A, beta_A, gamma_A);
%         Im_L = return_frame_to_canvas(Im, S, i_start, i_end, j_start, j_end);
%         As_L = return_frame_to_canvas(As, S, i_start, i_end, j_start, j_end);
%         A_L = return_frame_to_canvas(A, S, i_start, i_end, j_start, j_end);
%         Cs_L = return_frame_to_canvas(Cs, S, i_start, i_end, j_start, j_end);
%         C_L = return_frame_to_canvas(C, S, i_start, i_end, j_start, j_end);
%         [Im, i_start, i_end, j_start, j_end] = crop_frame(Im_L, crop_d);
%         As = As_L(i_start:i_end, j_start:j_end);
%         A = A_L(i_start:i_end, j_start:j_end);
%         Cs = Cs_L(i_start:i_end, j_start:j_end);
%         C = C_L(i_start:i_end, j_start:j_end);
%         U = U_matrix(Im);
%         Im_outline = outline_8p(~Im_L);
%         K1_L = K1_basal*ones(size(Im_L));
%         K1_L(Im_outline==1) = K1_edge;
%         K1 = K1_L(i_start:i_end, j_start:j_end);
        
        im1 = As;
        im1(Im==1) = (im1(Im==1)-min(im1(Im==1)))/(max(im1(Im==1))-min(im1(Im==1)));
        
        im2 = Cs;
        im2(Im==1) = (im2(Im==1)-min(im2(Im==1)))/(max(im2(Im==1))-min(im2(Im==1)));
        
        clf;
        subplot(2,2,1);
        hold on;
        colormap hot;
        axis ij;
        axis off;
        axis equal;
        imagesc(As);
        colorbar;
        title('As');
        
        subplot(2,2,2);
        hold on;
        colormap hot;
        axis ij;
        axis off;
        axis equal;
        imagesc(im1);
        colorbar;
        title('As scaled');
        
        subplot(2,2,3);
        hold on;
        colormap hot;
        axis ij;
        axis off;
        axis equal;
        imagesc(Cs);
        colorbar;
        title('Cs');
        
        subplot(2,2,4);
        hold on;
        colormap hot;
        axis ij;
        axis off;
        axis equal;
        imagesc(im2);
        colorbar;
        title('Cs scaled');
        drawnow;
    end
    
    %diffuse with noise
    for i = 1:sub_iter_N
        noise_A = randn(size(As))*n_alf;
        noise_C = randn(size(Cs))*n_alf;
        
        As_new = As + (DAs*laplacian_DT(As,dx,U) + noise_A)*dt;
        As_new(Im==0)=0;
        As_neg = As_new;
        As_neg(As_neg > 0) = 0;
        As_new(As_new < 0) = 0;
        
        A_new = A + (DA*laplacian_DT(A,dx,U) - noise_A)*dt;
        A_new(Im==0)=0;
        A_neg = A_new;
        A_neg(A_neg > 0) = 0;
        A_new(A_new < 0) = 0;
        
        As = As_new + A_neg;
        A = A_new + As_neg;
        
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
    end
    
    %react
    for i = 1:sub_iter_N
        F1 = (k1 + gamma1*As.^n1./(K1.^n1+beta1*Cs.^n1+As.^n1)).*A - k2.*As;
        F2 = (k3 + gamma3*As.^n3).*C-k4.*Cs;
        As = As + F1*dt;
        A = A - F1*dt;
        Cs = Cs + F2*dt;
        C = C - F2*dt;
    end
   
    %save
    if mod(j, save_T) == 0
        disp(j/save_T);
        imwrite(As, fullfile(save_fold, 'As_L', strcat(num2str(j/save_T), '.tif')));
        save(fullfile(save_fold, 'mat', strcat(num2str(j/save_T), '.mat')), ...
            'As_L', 'A_L', 'Cs_L', 'C_L', 'Im_L');
        saveas(fig, fullfile(save_fold, 'colored', strcat(num2str(j/save_T), '.png')));
    end
    
end
