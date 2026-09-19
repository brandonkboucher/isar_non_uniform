classdef Constants

    properties(Constant=true)
        c = 2.998e8; % [m/s]

        Hz2GHz = 1e-9;
        GHz2Hz = 1e9;
        MHz2Hz = 1e6;
        Hz2MHz = 1e-6;
        Hz2kHz = 1e-3;
        kHz2Hz = 1e3;

        MHz2GHz = 1e-3;
        GHz2MHz = 1e3;
        GHz2kHz = 1e6;
        kHz2GHz = 1e-6;
        kHz2MHz = 1e-3;
        MHz2kHz = 1e3;

        km2m = 1e3;
        m2km = 1e-3;

        deg2rad = pi/180;
        rad2deg = 180/pi;

        ns2s = 1e-9;
        s2ns = 1e9;
        us2s = 1e-6;
        s2us = 1e6;
        ms2s = 1e-3;
        s2ms = 1e3;

    end

end

