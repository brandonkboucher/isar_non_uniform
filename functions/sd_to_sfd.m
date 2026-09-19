function sfd = sd_to_sfd(sd)

    %sfd = fftshift(fft2(ifftshift(sd)));
    sfd = (fft2(ifftshift(sd)));

end

