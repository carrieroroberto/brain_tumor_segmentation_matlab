function mask_clean = segmentation(I_proc, seed_map, seq_name, filename)

    segmentation_dir = "results/plots/" + filename + "/segmentation/";
    mkdir(segmentation_dir);

    [~, maxIdx] = max(seed_map(:));
    [seedY, seedX] = ind2sub(size(I_proc), maxIdx);
    
    levels = multithresh(I_proc, 2);
    mask_otsu = I_proc > levels(2);
    
    thresh_rg = prctile(I_proc(mask_otsu), 6);
    mask_rg = region_growing(I_proc, seedY, seedX, thresh_rg);
    
    D = -bwdist(~mask_rg);
    D_mod = imimposemin(D, mask_otsu);
    D_mod(~mask_rg) = -Inf;
    L = watershed(D_mod);
    
    mask_raw = mask_rg;
    mask_raw(L == 0) = 0;
    mask_raw = bwselect(mask_raw, seedX, seedY, 8);
    
    fig = figure("Visible", "off");
    
    subplot(1, 4, 1);
    imshow(I_proc, []);
    hold on;
    visboundaries(mask_otsu, "Color", "y");
    title("Otsu");
    
    subplot(1, 4, 2);
    imshow(I_proc, []);
    hold on;
    visboundaries(mask_rg, "Color", "c");
    title("Region Growing");
    
    subplot(1, 4, 3);
    imshow(I_proc, []);
    hold on;
    visboundaries(mask_raw, "Color", "m");
    title("Watershed Marcato");
    
    subplot(1, 4, 4);
    imshow(I_proc, []);
    hold on;
    visboundaries(mask_raw, "Color", "r");
    title("Risultato Finale");
    
    sgtitle("Segmentazione " + seq_name + " - Paziente: " + filename);
    saveas(fig, segmentation_dir + seq_name + "_segmentation.png");
    close(fig);

    mask_clean = post_processing(mask_raw, I_proc, seq_name, filename);

end