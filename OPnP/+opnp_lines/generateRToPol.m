function RtoPol = generateRToPol()
%1 a2 ab ac ad b2 bc bd c2 cd d2    
    RtoPol = zeros(11, 9);
    RtoPol(2, 1) = 1;
    RtoPol(6, 1) = 1;
    RtoPol(9, 1) = -1;
    RtoPol(11, 1) = -1;

    RtoPol(7, 2) = 2;
    RtoPol(5, 2) = -2;

    RtoPol(4, 3) = 2;
    RtoPol(8, 3) = 2;

    RtoPol(5, 4) = 2;
    RtoPol(7, 4) = 2;

    RtoPol(2, 5) = 1;
    RtoPol(6, 5) = -1;
    RtoPol(9, 5) = 1;
    RtoPol(11, 5) = -1;

    RtoPol(3, 6) = -2;
    RtoPol(10, 6) = 2;

    RtoPol(4, 7) = -2;
    RtoPol(8, 7) = 2;

    RtoPol(3, 8) = 2;
    RtoPol(10, 8) = 2;

    RtoPol(2, 9) = 1;
    RtoPol(6, 9) = -1;
    RtoPol(9, 9) = -1;
    RtoPol(11, 9) = 1;
end