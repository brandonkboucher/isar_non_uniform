function title_str = create_title(str,omega, omega_dot)
    %% omega
    if omega > 1 && omega ~= inf

        title_str = [str, ' ($\omega = \pi/',num2str(omega),'$ rad/s, '];
    
    elseif omega == inf

        title_str = [str, ' ($\omega = 0$ rad/s, '];
    
    else

        omega = 1 / omega;
        title_str = [str, ' ($\omega = ', num2str(omega),'\pi$ rad/s, '];
    
    end

    %% omega dot
    if omega_dot > 1 && omega_dot ~= inf

        title_str = [title_str, '$\dot{\omega} = \pi/', num2str(omega_dot),'$ rad/s/s)'];
    
    elseif omega_dot == inf

        title_str = [title_str, '$\dot{\omega} = 0$ rad/s/s)'];
    
    else

        omega_dot = 1 / omega_dot;
        title_str = [title_str, '$\dot{\omega} = ', num2str(omega_dot), '\pi$ rad/s/s)'];
    
    end
end

