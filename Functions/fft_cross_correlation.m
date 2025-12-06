function [azimuth_shift, range_shift] = fft_cross_correlation(img1,img2)

    % evaluate the FFT of the images
    F1 = fft2(img1);
    F2 = fft2(img2);

    % multiple the reference spectrum by the complex
    % conjugate of the other spectrum (correlation)
    R = F1 .* conj(F2);

    % normalize
    R = R ./ abs(R);

    % find the index of the point that maximizes the
    % correlation
    c = abs(ifft2(R));
    c  = fftshift(c);  % center the peak
    
    [c_max, idx] = max(c(:));

    % find the shift in the x and y directions
    [idx_x, idx_y] = ind2sub(size(c),idx);

    % Compute shifts relative to center    
    azimuth_shift = idx_x - ceil(size(c,1)/2);    
    range_shift = idx_y - ceil(size(c,2)/2);

end

