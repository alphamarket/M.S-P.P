function Image=im2gray(Image)
    % converts an image to its gray scale image
    try
        % try rgb2gray.
        Image = rgb2gray(Image);
    catch E
        % if already in grayscale?
        if (strcmp(E.identifier, 'Images:rgb2gray:invalidSizeForColormap'))
            % move on with origin image.
            Image = Image(:,:,1);
        else
            % otherwise rethrow the exception.
            rethrow(E);
        end
    end
end