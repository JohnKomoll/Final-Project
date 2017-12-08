function [num_cells] = cell_counter( video )

num_cells = [];
counter = 0;
while hasFrame(video)
    counter = counter + 1;
    
    vidframe = im2double(readFrame(video));
    vidframe = imadjust(vidframe(:,:,1));
    
    % Smooth image
    rad = 5;
    sigma = 5;
    fgauss = fspecial('gaussian', rad, sigma);
    vidframe = imfilter(vidframe, fgauss);
    
    % Subtract background
    img_noise = imopen(vidframe, strel('disk', 100));
    vidframe_sb = imsubtract(vidframe, img_noise);
    
    % Take a mask to find nuclei
    mask = vidframe_sb > 0.08;
    cell_data = regionprops(mask, vidframe_sb, 'Area');
    areas = [cell_data.Area];
    
    % Filter out bad areas
    cells = length(areas);
    num_cells(counter) = cells;
    
end

end