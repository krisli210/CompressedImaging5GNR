function [azInd, rangeInd] = getDCM_Image()
    im = zeros(26, 128);
    

    % D 
    D_start = 8;
    im(9:17, D_start) = 1;
    im(9, D_start:D_start+12) = 1;
    im(17, D_start:D_start+12) = 1;
    
    curve_start = D_start+12;
    im(10, curve_start:curve_start+1) = 1;
    im(11, curve_start+2:curve_start+3) = 1;
    im(12, curve_start+4:curve_start+5) = 1;
    
    im(13, curve_start+4:curve_start+5) = 1;

    im(16, curve_start:curve_start+1) = 1;
    im(15, curve_start+2:curve_start+3) = 1;
    im(14, curve_start+4:curve_start+5) = 1;
    

    % C
    im(9:17, 58) = 1;
    im(9, 58:70) = 1;
    im(17, 58:70) = 1;
    
    % M 
    M_start = 103;
    im(9:17, M_start) = 1;
    %%% M diagonal
    im(9, M_start+1:M_start+2) = 1;
    im(10, M_start+3:M_start+4) = 1;
    im(11, M_start+5:M_start+6) = 1;
    im(12, M_start+7:M_start+8) = 1;
    im(11, M_start+9:M_start+10) = 1;
    im(10, M_start+11:M_start+12) = 1;
    im(9, M_start+13:M_start+14) = 1;
    im(9:17, M_start+15) = 1;

    [rangeInd, azInd] = find(im);
    
end