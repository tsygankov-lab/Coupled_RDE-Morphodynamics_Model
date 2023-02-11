clear; clc;

root_fold = 'E:\Asef_Cdc42_Rac1_model\ruffling_differentiator\2D_dynamic_cell';

subfolds = dir(root_fold);
subfolds_ids = [subfolds.isdir];
subfolds = {subfolds.name}';
subfolds = subfolds(subfolds_ids);
subfolds = subfolds(3:end);

subfolds = {'cell_R_90_K1_edge_1.4_alpha_A_50_A_act_0.1_gamma_A_0.2'};

plot_subfolds = {'As_jet_crop'};

T = 1;
T = 5;

for sf_id = 1:length(subfolds)
    subfold = subfolds{sf_id};
    disp(subfold);
    
    for p_sf_id = 1:length(plot_subfolds)
        plot_subfold = plot_subfolds{p_sf_id};
        disp(plot_subfold);

        imageNames = dir(fullfile(root_fold, subfold, plot_subfold, '*.png'));
        outputVideo = VideoWriter(fullfile(root_fold, subfold, strcat(plot_subfold, '_T_', num2str(T), '.avi')));
        outputVideo.FrameRate = 30;
        outputVideo.Quality = 100;
        open(outputVideo);

        start_f = 1;
        files = dir(fullfile(root_fold, subfold, plot_subfold, '*.png'));
        files = {files.name}';
        %end_f = length(files);
        %end_f = 4991;
        %end_f = str2num(strrep(files{end}, '.png', ''));
        end_f = start_f;
        for f_id = 1:length(files)
            f = files{f_id};
            f = str2num(strrep(f, '.png', ''));
            if end_f < f
                end_f = f;
            end
        end
        
        for j = start_f:T:end_f
            img = imread(fullfile(root_fold, subfold, plot_subfold, strcat(num2str(j),'.png')));
            img = imresize(img, 1);
            writeVideo(outputVideo,img);
        end
        close(outputVideo);

    end
end
