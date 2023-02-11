clear; clc;

root_fold = 'E:\Asef_Cdc42_Rac1_model\ruffling_differentiator\2D_dynamic_cell';

subfolds = dir(root_fold);
subfolds_ids = [subfolds.isdir];
subfolds = {subfolds.name}';
subfolds = subfolds(subfolds_ids);
subfolds = subfolds(3:end);

subfolds = {'cell_R_90_K1_edge_1.4_alpha_A_50_A_act_0.1_gamma_A_0.2'};

T = 1;
d = 30;

for sf_id = 1:length(subfolds)
    subfold = subfolds{sf_id};
    disp(subfold);
    files = dir(fullfile(root_fold, subfold, 'mat', '*.mat'));
    files = {files.name};
    N = length(files);
    
    load(fullfile(root_fold, subfold, 'parameters.mat'));
    load(fullfile(root_fold, subfold, 'data_inspection.mat'));
    
    s = [j2-j1+2*d+1, i2-i1+2*d+1];
    
    mkdir(fullfile(root_fold, subfold, 'As_jet_crop'));
    As_min = min(As_min_vals);
    As_max = max(As_max_vals);
    
    files = dir(fullfile(root_fold, subfold, 'As_jet_crop', '*.png'));
    file_1 = 1;
%     if ~isempty(files)
%         file_1 = 1 + fix(length(files)/T)*T;
%     end
    
    
    fig = figure('Position', [50 50 2*s(2) 2*s(1)]);
    for i = file_1:T:N
        load(fullfile(root_fold, subfold, 'mat', strcat(num2str(i), '.mat')));
        
%         Im = Im_L(i1-d:i2+d, j1-d:j2+d);
%         Im_nucl = Im_nucl_L(i1-d:i2+d, j1-d:j2+d);
%         As = As_L(i1-d:i2+d, j1-d:j2+d);
        
        Im = Im_L(j1-d:j2+d, i1-d:i2+d);
        As = As_L(j1-d:j2+d, i1-d:i2+d);
        
        Im_all_cell = Im;
        Im_outline = outline_8p(~Im_all_cell);
        
        im = As;
        im(im>0) = (im(im>0)-As_min)/(As_max-As_min);
        im = round(im*255);
        im_rgb = ind2rgb(im, jet);
        im_r = im_rgb(:,:,1);
        im_g = im_rgb(:,:,2);
        im_b = im_rgb(:,:,3);
        im_r(Im==0) = 0; im_g(Im==0) = 0; im_b(Im==0) = 0;
        im_rgb = cat(3, im_r, im_g, im_b);
        
        clf;
        hold on;
        axis ij;
        axis off;
        image(im_rgb);
        set(gca, 'box','off','XTickLabel',[],'XTick',[],'YTickLabel',[],'YTick',[])
        set(gca,'position',[0 0 1 1]);
        set(gca, 'XLim', [0 s(2)], 'YLim', [0 s(1)]);
        set(fig,'PaperUnits','points');
        set(fig,'PaperSize',[s(2) s(1)]);
        set(fig,'PaperPosition',[0 0 s(2) s(1)]);
        set(gca,'Color','k');
        fig.Color = 'black';
        fig.InvertHardcopy = 'off';
        drawnow;
        saveas(fig, fullfile(root_fold, subfold, 'As_jet_crop', strcat(num2str(i), '.png')));
       
    end
    close(fig);
end
