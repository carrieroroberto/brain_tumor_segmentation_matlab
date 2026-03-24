function [imgs_proc, mask_gt, seed_map] = pre_processing(path_img, path_gt, filename)
% File: pre_processing.m
% Esegue la standardizzazione dei dati e il miglioramento del segnale.
% Applica normalizzazione Z-Score, filtraggio spaziale.
% Genera inoltre combinazioni multimodali (Early Fusion) e si auto-documenta.
%
% INPUT:
% path_img - percorso al volume NIfTI (4D)
% path_gt  - percorso alla Ground Truth NIfTI (3D)
% filename - identificativo del paziente (usato per l'albero di salvataggio)
%
% OUTPUT:
% imgs_proc - struct con le configurazioni preprocessate (flair, t1, t1c, t2, fus2, fus3)
% mask_gt   - maschera GT corrispondente alla slice estratta
% seed_map  - mappa di probabilità spaziale sfocata per il Region Growing

    % inizializzazione della directory per l'auto-documentazione grafica
    pre_processing_dir = "results/plots/" + filename + "/pre_processing/";
    mkdir(pre_processing_dir);
    
    % lettura dei dati volumetrici
    vol_4d = double(niftiread(path_img));
    vol_gt = double(niftiread(path_gt));
    
    % estrazione euristica della slice con maggiore area tumorale
    slice_scores = squeeze(sum(sum(vol_gt > 0, 1), 2));
    [~, z_best] = max(slice_scores);
    
    img_fl = vol_4d(:,:, z_best, 1);
    img_t1 = vol_4d(:,:, z_best, 2);
    img_t1c = vol_4d(:,:, z_best, 3);
    img_t2 = vol_4d(:,:, z_best, 4);
    mask_gt = vol_gt(:,:, z_best) > 0;
    
    % segmentazione elementare del background
    mask_brain = imfill(img_fl > 0, "holes");
    
    % normalizzazione Z-Score per canale (allinea i domini di intensità tra scanner)
    img_fl(mask_brain) = (img_fl(mask_brain) - mean(img_fl(mask_brain))) / std(img_fl(mask_brain));
    img_t1c(mask_brain) = (img_t1c(mask_brain) - mean(img_t1c(mask_brain))) / std(img_t1c(mask_brain));
    img_t2(mask_brain) = (img_t2(mask_brain) - mean(img_t2(mask_brain))) / std(img_t2(mask_brain));
    
    % soppressione formale del background post-normalizzazione
    img_fl(~mask_brain) = 0;
    img_t1c(~mask_brain) = 0;
    img_t2(~mask_brain) = 0;
    img_t1(~mask_brain) = 0;
    
    % scaling lineare nel range [0, 1] per le operazioni morfologiche
    img_fl_norm = mat2gray(img_fl);
    img_t1c_norm = mat2gray(img_t1c);
    img_t2_norm = mat2gray(img_t2);
    img_t1_norm = mat2gray(img_t1);
    
    % riduzione rumore (Filtro Mediano)
    img_fl_filt = medfilt2(img_fl_norm, [3 3]);
    img_t1_filt = medfilt2(img_t1_norm, [3 3]);
    img_t1c_filt = medfilt2(img_t1c_norm, [3 3]);
    img_t2_filt = medfilt2(img_t2_norm, [3 3]);
    
    % salvataggio output strutturato
    imgs_proc.flair = img_fl_filt;
    imgs_proc.t1 = img_t1_filt;
    imgs_proc.t1c = img_t1c_filt;
    imgs_proc.t2 = img_t2_filt;
    
    % fusione multimodale (media pesata basata sulle caratteristiche cliniche)
    imgs_proc.fus2 = 0.6*img_fl_filt + 0.4*img_t1c_filt;
    imgs_proc.fus3 = 0.54*img_fl_filt + 0.16*img_t1c_filt + 0.30*img_t2_filt;
    
    % generazione della mappa termica mediante smoothing Gaussiano
    seed_map = imgaussfilt(imgs_proc.fus3, 3);
    
    % ================= LOGICA DI PLOTTING E AUTO-DOCUMENTAZIONE =================
    
    fig_eda = figure("Visible", "off");
    imgs = {img_fl, img_t1, img_t1c, img_t2};
    titles = ["FLAIR Originale", "T1 Originale", "T1c Originale", "T2 Originale"];
    
    for k = 1:4
        subplot(2, 4, k);
        imshow(imgs{k});
        title(titles(k));
        subplot(2, 4, k+4);
        histogram(imgs{k}(mask_brain));
        xlabel("Intensità");
        ylabel("Conteggio Pixel");
    end
    sgtitle("Analisi Sequenze Grezze - Paziente: " + filename);
    saveas(fig_eda, pre_processing_dir + "raw_histogram.png");
    close(fig_eda);
    
    plot_pre_step(img_fl, img_fl_norm, img_fl_filt, mask_brain, "FLAIR", filename, pre_processing_dir);
    plot_pre_step(img_t1c, img_t1c_norm, img_t1c_filt, mask_brain, "T1c", filename, pre_processing_dir);
    plot_pre_step(img_t2, img_t2_norm, img_t2_filt, mask_brain, "T2", filename, pre_processing_dir);
    plot_fus2_step(imgs_proc.flair, imgs_proc.t1c, imgs_proc.fus2, mask_brain, filename, pre_processing_dir);
    plot_fus3_step(imgs_proc.flair, imgs_proc.t1c, imgs_proc.t2, imgs_proc.fus3, mask_brain, filename, pre_processing_dir);
    
    fig_seed = figure("Visible", "off");
    imshow(seed_map, []);
    colormap(jet);
    colorbar;
    title("Fusion3 Seed Map (Sfocatura Gaussiana) - Paziente: " + filename);
    saveas(fig_seed, pre_processing_dir + "seed_map.png");
    close(fig_seed);
end

% --- Funzioni accessorie di plotting ---

function plot_pre_step(raw, norm, filt, mask, seq, filename, save_path)
    fig = figure("Visible", "off");
    imgs = {raw, norm, filt};
    titles = ["Sequenza Grezza", "Normalizzazione Z-Score", "Filtro Rumore Mediano"];
    
    for k = 1:3
        subplot(2, 3, k);
        imshow(imgs{k});
        title(titles(k));
        subplot(2, 3, k+3);
        histogram(imgs{k}(mask));
        xlabel("Intensità");
    end
    sgtitle("Fasi Pre-Processing " + seq + " - Paziente: " + filename);
    saveas(fig, save_path + seq + "_pre_processing.png");
    close(fig);
end

function plot_fus2_step(fl, t1c, fus, mask, filename, save_path)
    fig = figure("Visible", "off");
    imgs = {fl, t1c, fus};
    titles = ["FLAIR Pre-Processata", "T1c Pre-Processata", "Fusion2 (FLAIR+T1c)"];
    
    for k = 1:3
        subplot(2, 3, k);
        imshow(imgs{k});
        title(titles(k));
        subplot(2, 3, k+3);
        histogram(imgs{k}(mask));
        xlabel("Intensità");
    end
    sgtitle("Generazione Fusione Bimodale - Paziente: " + filename);
    saveas(fig, save_path + "Fusion2_pre_processing.png");
    close(fig);
end

function plot_fus3_step(fl, t1c, t2, fus, mask, filename, save_path)
    fig = figure("Visible", "off");
    imgs = {fl, t1c, t2, fus};
    titles = ["FLAIR Pre-Processata", "T1c Pre-Processata", "T2 Pre-Processata", "Fusion3 (FLAIR+T1c+T2)"];
    
    for k = 1:4
        subplot(2, 4, k);
        imshow(imgs{k});
        title(titles(k));
        subplot(2, 4, k+4);
        histogram(imgs{k}(mask));
        xlabel("Intensità");
    end
    sgtitle("Generazione Fusione Multimodale - Paziente: " + filename);
    saveas(fig, save_path + "Fusion3_pre_processing.png");
    close(fig);
end