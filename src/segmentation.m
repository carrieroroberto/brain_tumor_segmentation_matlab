function mask_clean = segmentation(I_proc, seed_map, seq_name, filename)
% File: segmentation.m
% Implementa l'architettura algoritmica ibrida per l'estrazione della ROI tumorale.
% Coordina la sequenza: Otsu -> Region Growing -> Trasformata Watershed.
%
% INPUT:
% I_proc   - immagine MRI bidimensionale post-equalizzazione
% seed_map - mappa di intensità sfocata per la localizzazione automatica
% seq_name - stringa indicante la configurazione elaborata
% filename - nome del file NIfTI analizzato
%
% OUTPUT:
% mask_clean - maschera finale pulita, pronta per la valutazione quantitativa

    segmentation_dir = "results/plots/" + filename + "/segmentation/";
    mkdir(segmentation_dir);
    
    % ricerca del picco di intensità sulla seed map per le coordinate iniziali
    [~, maxIdx] = max(seed_map(:));
    [seedY, seedX] = ind2sub(size(I_proc), maxIdx);
    
    % segmentazione multilivello (Otsu) per isolare l'area di confidenza tumorale
    levels = multithresh(I_proc, 2);
    mask_otsu = I_proc > levels(2);
    
    % calcolo della soglia dinamica (6° percentile) e avvio della Region Growing
    thresh_rg = prctile(I_proc(mask_otsu), 6);
    mask_rg = region_growing(I_proc, seedY, seedX, thresh_rg);
    
    % costruzione del bacino topografico tramite trasformata della distanza inversa
    D = -bwdist(~mask_rg);
    % imposizione dei minimi locali per forzare il Watershed sui marker di Otsu
    D_mod = imimposemin(D, mask_otsu);
    % restrizione del dominio per prevenire l'allagamento dello sfondo (leakage)
    D_mod(~mask_rg) = -Inf;
    
    % esecuzione dell'algoritmo spartiacque
    L = watershed(D_mod);
    
    % rimozione dei confini di separazione
    mask_raw = mask_rg;
    mask_raw(L == 0) = 0;
    
    % isolamento esclusivo della componente morfologica contenente il seed iniziale
    mask_raw = bwselect(mask_raw, seedX, seedY, 8);
    
    % ================= LOGICA DI PLOTTING E AUTO-DOCUMENTAZIONE =================
    
    fig = figure("Visible", "off");
    
    subplot(2, 2, 1);
    imshow(I_proc, []);
    hold on;
    visboundaries(mask_otsu, "Color", "y");
    title("Otsu");
    
    subplot(2, 2, 2);
    imshow(I_proc, []);
    hold on;
    visboundaries(mask_rg, "Color", "c");
    title("Region Growing");
    
    subplot(2, 2, 3);
    imshow(I_proc, []);
    hold on;
    visboundaries(mask_raw, "Color", "m");
    title("Watershed Marcato");
    
    subplot(2, 2, 4);
    imshow(I_proc, []);
    hold on;
    visboundaries(mask_raw, "Color", "r");
    title("Maschera Finale");
    
    sgtitle("Segmentazione " + seq_name + " - Paziente: " + filename);
    saveas(fig, segmentation_dir + seq_name + "_segmentation.png");
    close(fig);
    
    % delega del passo finale di regolarizzazione alla funzione dedicata
    mask_clean = post_processing(mask_raw, I_proc, seq_name, filename);
end