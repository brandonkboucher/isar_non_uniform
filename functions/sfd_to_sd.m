function sd = sfd_to_sd(sfd)

    sd = fftshift(ifft2(ifftshift(sfd)));
end

