function mask_clean = post_processing(mask_raw, I_proc, seq_name, filename)
% File: post_processing.m
% Applica operazioni di algebra morfologica per la regolarizzazione dei contorni,
% il riempimento delle cavità e la rimozione del rumore. Genera il grafico "Prima/Dopo".
%
% INPUT:
% mask_raw - maschera binaria grezza derivante dall'algoritmo spartiacque
% I_proc   - immagine strutturale usata come livello di background per i plot
% seq_name - configurazione in elaborazione
% filename - nome del file in elaborazione
%
% OUTPUT:
% mask_clean - maschera binaria clinicamente robusta post-operazioni

    post_processing_dir = "results/plots/" + filename + "/post_processing/";
    mkdir(post_processing_dir);
    
    % definizione dell'elemento strutturante (disco simmetrico)
    se = strel("disk", 3);
    
    % operazione di chiusura per saldare frammenti al confine
    mask_clean = imclose(mask_raw, se);
    % riempimento del core necrotico centrale (spesso ipointenso e ignorato)
    mask_clean = imfill(mask_clean, "holes");
    % filtraggio conservativo per preservare solo la massa connessa dominante
    mask_clean = bwareafilt(mask_clean, 1);
    
    % ================= LOGICA DI PLOTTING E AUTO-DOCUMENTAZIONE =================
    
    fig = figure("Visible", "off");
        
    subplot(1, 2, 1);
    imshow(I_proc, []);
    hold on;
    visboundaries(mask_raw, "Color", "r");
    title("Prima - Maschera Watershed Grezza");
        
    subplot(1, 2, 2);
    imshow(I_proc, []);
    hold on;
    visboundaries(mask_clean, "Color", "g");
    title("Dopo - Maschera con Operazioni Morfologiche");
        
    sgtitle("Post-Processing " + seq_name + " - Paziente: " + filename);
    saveas(fig, post_processing_dir + seq_name + "_post_processing.png");
    close(fig);
    
end