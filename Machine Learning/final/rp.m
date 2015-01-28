function [ring_projection] = rp( Image, radius_step_size, degree_step_size)
%RP Calculates the ring projection of an image
    % gray scale the image
    if(size(Image, 3) > 1)
        Image = im2gray(Image);
    end
    % we have to inject some blank area into image so the whole image can
    % fit into our radius(R) zone
    % calc euclidean distance of image's size
    R = ceil(sqrt(sum(size(Image).^2, 2)) ./ 2);
    % expand the image to be confi. with ringing the image!
    padded_image = padarray(Image,[ceil(R - size(Image, 1) / 2) ceil(R - size(Image, 2) / 2)],'both');
    % calc the mid coordinate of expanded image
    mid = size(padded_image) ./ 2;
    % define the radius step size
    r_step = radius_step_size;
    % define the degree step size
    grad_step = degree_step_size;
    % expecting to #{r} rotate upto #{expected_x_size} times
    expected_x_size = ceil(R/r_step + 1);
    % expecting to #{theta} rotate upto #{expected_y_size} times [ no use ]
    %expected_y_size = ceil((R/r_step + 1) * floor((2*pi/grad_step + 1)));
    % create a ring projection 1-D array
    % in fact this would be the 1-D conversion of the 2-D image
    ring_projection = zeros(1, expected_x_size);
    % for every possible radius
    for r=r_step:r_step:R-1 % not that #{-1} is to prevent overflowing
        % calc. the index of current circle index in ring_projection matrix
        u = floor(r/r_step + 1);
        % reset cixel
        cixel = 0;
        % make circle rotation on current radius
        for theta=0:grad_step:2*pi
            % calc. the coords in polar coord. && 
            % sum the brightness of current pixcel with other brightnesses
            % of current circle 
            cixel = cixel + double(padded_image(ceil(mid(1) + r*cos(theta)), ceil(mid(2) + r*sin(theta))));
        end
        % set cixel value to current index
        ring_projection(u) = ((cixel));
    end
end