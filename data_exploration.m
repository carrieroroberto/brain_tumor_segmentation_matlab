% File: data_exploration.m
% Esegue l'Exploratory Data Analysis (EDA) su un singolo paziente rappresentativo.
% Fornisce visualizzazioni comparative delle distribuzioni di intensità (istogrammi)
% e un'analisi qualitativa 3D per comprendere il contesto spaziale del tumore.

clear; clc; close all; 

%% ---------------------------- PARAMETRI ---------------------------------
id = "006"; % identificativo del paziente campione
filename = "BRATS_" + id + ".nii.gz"; 
path_img = "dataset/Task01_BrainTumour/imagesTr/" + filename; 
path_gt  = "dataset/Task01_BrainTumour/labelsTr/" + filename; 
slice = 78; % slice assiale ottimale per la visualizzazione bidimensionale
disp("Inizio analisi esplorativa sul paziente: " + filename);

%% -------------------------- CARICAMENTO DATI ----------------------------
% lettura dei volumi NIfTI in formato double precision
vol_4d = double(niftiread(path_img)); % dimensione: [H x W x slices x 4 canali]
vol_gt = double(niftiread(path_gt)); 

%% ----------------- ESTRAZIONE SLICE E SEQUENZA MRI ---------------------
img_fl  = vol_4d(:,:,slice,1); % FLAIR: evidenzia edema
img_t1  = vol_4d(:,:,slice,2); % T1: dettaglio anatomico base
img_t1c = vol_4d(:,:,slice,3); % T1c: evidenzia barriera emato-encefalica (core)
img_t2  = vol_4d(:,:,slice,4); % T2: dettaglio fluidi e confini

% creazione della maschera del tessuto per l'isolamento del segnale utile
brain_mask = img_fl > 0; 

%% ---------------------- PLOT 1: SEQUENZA + ISTOGRAMMI ----------------
fig1 = figure("Name", "Analisi Sequenze MRI"); 
imgs = {img_fl, img_t1, img_t1c, img_t2}; 
titles = ["FLAIR", "T1", "T1c", "T2"]; 

for i = 1:4
    % riga superiore: visualizzazione spaziale della risonanza
    subplot(2, 4, i); 
    imshow(mat2gray(imgs{i})); 
    title(titles(i)); 
    
    % riga inferiore: distribuzione statistica delle intensità
    subplot(2, 4, i + 4);
    data = imgs{i}(brain_mask); % isolamento dei pixel appartenenti al cervello
    
    histogram(data); 
    xlabel("Intensità"); 
    ylabel("Numero Pixel"); 
end
sgtitle("Analisi Sequenze MRI - Paziente: " + filename);

%% -------------------------- PLOT 2: VOLUME 3D ---------------------------
% attivazione del viewer OpenGL per il rendering volumetrico
volshow(vol_4d(:,:,:,1));
disp("Analisi esplorativa completata.");