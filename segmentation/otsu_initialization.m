function [seed_point, marker_int, marker_ext, img_roi] = otsu_initialization(img)
    thresholds = multithresh(img, 4);
    quantized = imquantize(img, thresholds);
    all_bright_spots = (quantized >= 4);
    
    se_break = strel('disk', 4);
    broken_spots = imopen(all_bright_spots, se_break);
    stats = regionprops(broken_spots, img, 'Area', 'MeanIntensity', 'PixelIdxList', 'Solidity');
    
    if isempty(stats)
        se_break = strel('disk', 2);
        broken_spots = imopen(all_bright_spots, se_break);
        stats = regionprops(broken_spots, img, 'Area', 'MeanIntensity', 'PixelIdxList', 'Solidity');
    end
    
    if isempty(stats)
        se_break = strel('disk', 1);
        broken_spots = all_bright_spots;
        stats = regionprops(broken_spots, img, 'Area', 'MeanIntensity', 'PixelIdxList', 'Solidity');
    end
    
    scores = [stats.Area] .* [stats.MeanIntensity] .* ([stats.Solidity].^2);
    [~, target_idx] = max(scores);
    
    temp_mask = false(size(img));
    temp_mask(stats(target_idx).PixelIdxList) = true;
    temp_mask = imfill(temp_mask, 'holes');
    
    core_mask = imerode(temp_mask, strel('disk', 3));
    if sum(core_mask(:)) == 0
        core_mask = temp_mask;
    end
    
    idx = find(core_mask);
    [r, c] = ind2sub(size(img), idx(round(end/2)));
    seed_point = [r, c];
    
    tumor_mask = bwselect(broken_spots, c, r, 4);
    tumor_mask = imdilate(tumor_mask, se_break);
    tumor_mask = imfill(tumor_mask, 'holes');
    tumor_mask = tumor_mask & all_bright_spots;
    
    img_roi = img;
    img_roi(~tumor_mask) = 0;
    
    marker_int = imerode(tumor_mask, strel('disk', 3));
    marker_ext = ~imdilate(tumor_mask, strel('disk', 15)) & (img > 0.05);
end