function data_exploration(path_img, path_gt, save_dir)
% File: data_exploration.m
% Esegue EDA e salva figure di raw images+istogrammi e 3D volume
%
% INPUT:
% path_img - percorso file NIfTI immagini MRI
% path_gt  - percorso file NIfTI ground truth
% save_dir - cartella dove salvare le figure

    vol_4d = double(niftiread(path_img));
    vol_gt = double(niftiread(path_gt));

    slice = round(size(vol_4d,3)/2); % slice centrale per visualizzazione
    seq_names = {'FLAIR','T1','T1c','T2'};

    % RAW images + istogrammi
    fig = figure('Visible','off');
    for i = 1:4
        img = vol_4d(:,:,slice,i);
        brain_mask = img>0;

        subplot(2,4,i); imshow(mat2gray(img)); title(seq_names{i});
        subplot(2,4,i+4); histogram(img(brain_mask)); title(seq_names{i} + " Histogram");
    end
    sgtitle("Raw Images + Histograms - Slice "+slice);
    saveas(fig, save_dir + "raw_histogram.png");
    close(fig);

    % 3D volume FLAIR
    fig3 = figure('Visible','off');
    volshow(vol_4d(:,:,:,1)); % FLAIR
    saveas(fig3, save_dir + "3d_volume.png");
    close(fig3);

end