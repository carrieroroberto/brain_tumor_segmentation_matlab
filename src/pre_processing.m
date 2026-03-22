function [imgs_proc, mask_gt, seed_map] = pre_processing(path_img, path_gt, filename)

    pre_processing_dir = "results/plots/" + filename + "/pre_processing/";
    mkdir(pre_processing_dir);

    vol_4d = double(niftiread(path_img));
    vol_gt = double(niftiread(path_gt));
    
    slice_scores = squeeze(sum(sum(vol_gt > 0, 1), 2));
    [~, z_best] = max(slice_scores);
    
    raw_fl = vol_4d(:,:,z_best,1); raw_t1 = vol_4d(:,:,z_best,2);
    raw_t1c = vol_4d(:,:,z_best,3); raw_t2 = vol_4d(:,:,z_best,4);
    mask_gt = vol_gt(:,:,z_best) > 0;
    
    mask_brain = imfill(raw_fl > 0, "holes");
    
    img_fl = raw_fl; img_t1 = raw_t1; img_t1c = raw_t1c; img_t2 = raw_t2;
    img_fl(mask_brain) = (img_fl(mask_brain) - mean(img_fl(mask_brain))) / std(img_fl(mask_brain));
    img_t1c(mask_brain) = (img_t1c(mask_brain) - mean(img_t1c(mask_brain))) / std(img_t1c(mask_brain));
    img_t2(mask_brain) = (img_t2(mask_brain) - mean(img_t2(mask_brain))) / std(img_t2(mask_brain));
    
    img_fl(~mask_brain)=0; img_t1c(~mask_brain)=0; img_t2(~mask_brain)=0; img_t1(~mask_brain)=0;
    zsc_fl = mat2gray(img_fl); zsc_t1c = mat2gray(img_t1c); zsc_t2 = mat2gray(img_t2); zsc_t1 = mat2gray(img_t1);
    
    imgs_proc.flair = adapthisteq(medfilt2(zsc_fl, [3 3]), "ClipLimit", 0.015);
    imgs_proc.t1    = adapthisteq(medfilt2(zsc_t1, [3 3]), "ClipLimit", 0.015);
    imgs_proc.t1c   = adapthisteq(medfilt2(zsc_t1c, [3 3]), "ClipLimit", 0.015);
    imgs_proc.t2    = adapthisteq(medfilt2(zsc_t2, [3 3]), "ClipLimit", 0.015);
    
    imgs_proc.fus2 = 0.6*imgs_proc.flair + 0.4*imgs_proc.t1c;
    imgs_proc.fus3 = 0.54*imgs_proc.flair + 0.16*imgs_proc.t1c + 0.30*imgs_proc.t2;
    
    seed_map = imgaussfilt(imgs_proc.fus3, 3);

    fig_eda = figure("Visible", "off", "Position", [100, 100, 1600, 800], "Color", "w");
    raws = {raw_fl, raw_t1, raw_t1c, raw_t2}; titles_raw = ["FLAIR Raw", "T1 Raw", "T1c Raw", "T2 Raw"];
    for k = 1:4
        subplot(2, 4, k); imshow(mat2gray(raws{k})); title(titles_raw(k), "FontWeight", "bold");
        subplot(2, 4, k+4); histogram(raws{k}(mask_brain), 50); xlabel("Intensità"); ylabel("Freq");
    end
    sgtitle("Esplorazione Dati Grezzi (Raw) - " + filename, "FontSize", 16, "FontWeight", "bold");
    saveas(fig_eda, pre_processing_dir + "raw_histogram.png");
    close(fig_eda);

    plot_pre_step(raw_fl, zsc_fl, imgs_proc.flair, mask_brain, "FLAIR", pre_processing_dir + "flair_preprocessing.png");
    plot_pre_step(raw_t1c, zsc_t1c, imgs_proc.t1c, mask_brain, "T1c", pre_processing_dir + "t1c_preprocessing.png");
    plot_pre_step(raw_t2, zsc_t2, imgs_proc.t2, mask_brain, "T2", pre_processing_dir + "t2_preprocessing.png");
    plot_fus2_step(imgs_proc.flair, imgs_proc.t1c, imgs_proc.fus2, mask_brain, pre_processing_dir + "fus2_preprocessing.png");
    plot_fus3_step(imgs_proc.flair, imgs_proc.t1c, imgs_proc.t2, imgs_proc.fus3, mask_brain, pre_processing_dir + "fus3_preprocessing.png");

    fig_seed = figure("Visible", "off", "Color", "w");
    imshow(seed_map, []); colormap(jet); colorbar; title("Fus3 Seed Map (Gauss Smoothed)");
    saveas(fig_seed, pre_processing_dir + "seed_map.png"); close(fig_seed);
end

function plot_pre_step(raw, zsc, clahe, mask, name, save_path)
    fig = figure("Visible", "off", "Position", [100, 100, 1200, 800], "Color", "w");
    imgs = {mat2gray(raw), zsc, clahe}; titles = ["1. Raw", "2. Z-Score", "3. Median + CLAHE"];
    for k = 1:3
        subplot(2, 3, k); imshow(imgs{k}); title(titles(k), "FontWeight", "bold");
        subplot(2, 3, k+3); histogram(imgs{k}(mask), 50); xlabel("Intensità");
    end
    sgtitle("Preprocessing Step-by-Step: " + name, "FontSize", 16, "FontWeight", "bold");
    saveas(fig, save_path); close(fig);
end

function plot_fus2_step(fl, t1c, fus, mask, save_path)
    fig = figure("Visible", "off", "Position", [100, 100, 1200, 800], "Color", "w");
    imgs = {fl, t1c, fus}; titles = ["1. FLAIR (CLAHE)", "2. T1c (CLAHE)", "3. Fus2 (0.6*FL + 0.4*T1c)"];
    for k = 1:3
        subplot(2, 3, k); imshow(imgs{k}); title(titles(k), "FontWeight", "bold");
        subplot(2, 3, k+3); histogram(imgs{k}(mask), 50); xlabel("Intensità");
    end
    sgtitle("Generazione Fusione Bimodale (Fus2)", "FontSize", 16, "FontWeight", "bold");
    saveas(fig, save_path); close(fig);
end

function plot_fus3_step(fl, t1c, t2, fus, mask, save_path)
    fig = figure("Visible", "off", "Position", [100, 100, 1600, 800], "Color", "w");
    imgs = {fl, t1c, t2, fus}; titles = ["1. FLAIR", "2. T1c", "3. T2", "4. Fus3 (Multimodale)"];
    for k = 1:4
        subplot(2, 4, k); imshow(imgs{k}); title(titles(k), "FontWeight", "bold");
        subplot(2, 4, k+4); histogram(imgs{k}(mask), 50); xlabel("Intensità");
    end
    sgtitle("Generazione Fusione Multimodale (Fus3)", "FontSize", 16, "FontWeight", "bold");
    saveas(fig, save_path); close(fig);
end