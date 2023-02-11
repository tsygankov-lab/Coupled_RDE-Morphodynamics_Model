clear; clc;

root_fold = 'E:\Asef_Cdc42_Rac1_model\ruffling_differentiator\2D_dynamic_cell';

subfolds = dir(root_fold);
subfolds_ids = [subfolds.isdir];
subfolds = {subfolds.name}';
subfolds = subfolds(subfolds_ids);
subfolds = subfolds(3:end);

subfolds = {'cell_R_90_K1_edge_1.4_alpha_A_50_A_act_0.1_gamma_A_0.2'};

for sf_id = 1:length(subfolds)
    subfold = subfolds{sf_id};
    disp(subfold);
    files = dir(fullfile(root_fold, subfold, 'mat', '*.mat'));
    files = {files.name};
    N = length(files);
    disp(N);
    
    V_vals = zeros(1, N);
    
    A_vals = zeros(1, N);
    C_vals = zeros(1, N);
    
    As_min_vals = zeros(1, N);
    As_max_vals = zeros(1, N);
    As_ampl_vals = zeros(1, N);
    
    A_min_vals = zeros(1, N);
    A_max_vals = zeros(1, N);
    A_ampl_vals = zeros(1, N);
    
    Cs_min_vals = zeros(1, N);
    Cs_max_vals = zeros(1, N);
    Cs_ampl_vals = zeros(1, N);
    
    C_min_vals = zeros(1, N);
    C_max_vals = zeros(1, N);
    C_ampl_vals = zeros(1, N);
    
    load(fullfile(root_fold, subfold, 'mat', strcat(num2str(1), '.mat')));
    S = size(Im_L);
    
    Im_L_all = zeros(S(1), S(2));
    
    for i = 1:N
        load(fullfile(root_fold, subfold, 'mat', strcat(num2str(i), '.mat')));
        Im = Im_L;
        As = As_L;
        A = A_L;
        Cs = Cs_L;
        C = C_L;
        V_vals(i) = sum(Im(:));
        A_vals(i) = sum(As(:)+A(:));
        C_vals(i) = sum(Cs(:)+C(:));
        As_min_vals(i) = min(As(Im==1));
        As_max_vals(i) = max(As(Im==1));
        As_ampl_vals(i) = max(As(Im==1))-min(As(Im==1));
        A_min_vals(i) = min(A(Im==1));
        A_max_vals(i) = max(A(Im==1));
        A_ampl_vals(i) = max(A(Im==1))-min(A(Im==1));
        Cs_min_vals(i) = min(Cs(Im==1));
        Cs_max_vals(i) = max(Cs(Im==1));
        Cs_ampl_vals(i) = max(Cs(Im==1))-min(Cs(Im==1));
        C_min_vals(i) = min(C(Im==1));
        C_max_vals(i) = max(C(Im==1));
        C_ampl_vals(i) = max(C(Im==1))-min(C(Im==1));
        
        Im_L_all = Im_L_all + Im_L;
    end
    
    fig = figure('Position', [50 50 700 700]);
    hold on;
    grid on;
    box on;
    plot(V_vals, 'LineWidth', 2, 'Color', [1 0 0]);
    title('volume');
    xlabel('frame number');
    saveas(fig, fullfile(root_fold, subfold, 'volume.png'));
    close(fig);
    
    fig = figure('Position', [50 50 700 700]);
    hold on;
    grid on;
    box on;
    plot(A_vals, 'LineWidth', 2, 'Color', [1 0 0]);
    title('A total mass');
    xlabel('frame number');
    saveas(fig, fullfile(root_fold, subfold, 'A_total_mass.png'));
    close(fig);
    
    fig = figure('Position', [50 50 700 700]);
    hold on;
    grid on;
    box on;
    plot(A_vals./V_vals, 'LineWidth', 2, 'Color', [1 0 0]);
    title('A total concentration');
    xlabel('frame number');
    saveas(fig, fullfile(root_fold, subfold, 'A_total_concentration.png'));
    close(fig);
    
    fig = figure('Position', [50 50 700 700]);
    hold on;
    grid on;
    box on;
    plot(C_vals, 'LineWidth', 2, 'Color', [1 0 0]);
    title('C total mass');
    xlabel('frame number');
    saveas(fig, fullfile(root_fold, subfold, 'C_total_mass.png'));
    close(fig);
    
    fig = figure('Position', [50 50 700 700]);
    hold on;
    grid on;
    box on;
    plot(C_vals./V_vals, 'LineWidth', 2, 'Color', [1 0 0]);
    title('C total concentration');
    xlabel('frame number');
    saveas(fig, fullfile(root_fold, subfold, 'C_total_concentration.png'));
    close(fig);
    
    fig = figure('Position', [50 50 700 700]);
    hold on;
    grid on;
    box on;
    plot(As_min_vals, 'LineWidth', 2, 'Color', [1 0 0]);
    title('As min');
    xlabel('frame number');
    saveas(fig, fullfile(root_fold, subfold, 'As_min.png'));
    close(fig);
    
    fig = figure('Position', [50 50 700 700]);
    hold on;
    grid on;
    box on;
    plot(As_max_vals, 'LineWidth', 2, 'Color', [1 0 0]);
    title('As max');
    xlabel('frame number');
    saveas(fig, fullfile(root_fold, subfold, 'As_max.png'));
    close(fig);
    
    fig = figure('Position', [50 50 700 700]);
    hold on;
    grid on;
    box on;
    plot(As_ampl_vals, 'LineWidth', 2, 'Color', [1 0 0]);
    title('As amplitude');
    xlabel('frame number');
    saveas(fig, fullfile(root_fold, subfold, 'As_amplitude.png'));
    close(fig);
    
    fig = figure('Position', [50 50 700 700]);
    hold on;
    grid on;
    box on;
    plot(log(As_ampl_vals), 'LineWidth', 2, 'Color', [1 0 0]);
    title('log(As amplitude)');
    xlabel('frame number');
    saveas(fig, fullfile(root_fold, subfold, 'As_log_amplitude.png'));
    close(fig);

    fig = figure('Position', [50 50 700 700]);
    hold on;
    grid on;
    box on;
    plot(Cs_min_vals, 'LineWidth', 2, 'Color', [1 0 0]);
    title('Cs min');
    xlabel('frame number');
    saveas(fig, fullfile(root_fold, subfold, 'Cs_min.png'));
    close(fig);
    
    fig = figure('Position', [50 50 700 700]);
    hold on;
    grid on;
    box on;
    plot(Cs_max_vals, 'LineWidth', 2, 'Color', [1 0 0]);
    title('Cs max');
    xlabel('frame number');
    saveas(fig, fullfile(root_fold, subfold, 'Cs_max.png'));
    close(fig);
    
    fig = figure('Position', [50 50 700 700]);
    hold on;
    grid on;
    box on;
    plot(Cs_ampl_vals, 'LineWidth', 2, 'Color', [1 0 0]);
    title('Cs amplitude');
    xlabel('frame number');
    saveas(fig, fullfile(root_fold, subfold, 'Cs_amplitude.png'));
    close(fig);
    
    fig = figure('Position', [50 50 700 700]);
    hold on;
    grid on;
    box on;
    plot(log(Cs_ampl_vals), 'LineWidth', 2, 'Color', [1 0 0]);
    title('log(Cs amplitude)');
    xlabel('frame number');
    saveas(fig, fullfile(root_fold, subfold, 'Cs_log_amplitude.png'));
    close(fig);
    
    fig = figure('Position', [50 50 500 500]);
    hold on;
    axis ij;
    axis off;
    colormap(hot);
    imagesc(Im_L_all);
    set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
    set(gca,'position',[0 0 1 1]);
    set(gca, 'XLim', [0 S(2)], 'YLim', [0 S(1)]);
    set(fig,'PaperUnits','points');
    set(fig,'PaperSize',[S(2) S(1)]);
    set(fig,'PaperPosition',[0 0 S(2) S(1)]);
    drawnow;
    saveas(fig, fullfile(root_fold, subfold, 'Im_L_all.png'));
    close(fig);
    
    ids_i = find(sum(Im_L_all,1));
    ids_j = find(sum(Im_L_all,2));
    
    i1 = ids_i(1);
    i2 = ids_i(end);
    j1 = ids_j(1);
    j2 = ids_j(end);
    
    save(fullfile(root_fold, subfold, 'data_inspection.mat'), ...
        'V_vals', 'A_vals', 'C_vals', 'As_min_vals', 'As_max_vals', ...
        'As_ampl_vals', 'A_min_vals', 'A_max_vals', 'A_ampl_vals', ...
        'Cs_min_vals', 'Cs_max_vals', 'Cs_ampl_vals', 'C_min_vals', ...
        'C_max_vals', 'C_ampl_vals', 'i1', 'i2', 'j1', 'j2');

end
