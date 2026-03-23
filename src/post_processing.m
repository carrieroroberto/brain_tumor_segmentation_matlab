function mask_clean = post_processing(mask_raw, I_proc, seq_name, filename)
    
    post_processing_dir = "results/plots/" + filename + "/post_processing/";
    mkdir(post_processing_dir);

    se = strel("disk", 3);
    mask_clean = imclose(mask_raw, se);
    mask_clean = imfill(mask_clean, "holes");
    mask_clean = bwareafilt(mask_clean, 1);
    
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