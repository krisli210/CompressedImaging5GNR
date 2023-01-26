function [dict] = generateDFT(arrDim, sceneDim)
    %ensure sceneDim is even to make life easier
    if (mod(sceneDim, 2) ~= 0 || mod(arrDim, 2) ~= 0) 
        error('Ensure sceneDim AND arrDim are even');
    end

    
    dft = (1/sqrt(sceneDim)) .* exp(-1j * 2*pi * (0:sceneDim-1).' * (0:sceneDim-1) ./ sceneDim);
    % -60 degrees at ind = (sceneDim/2) + 2

    % 0 degrees at ind = 1
    
    % 60 degrees at ind = sceneDim/2
    dict = [dft(:, sceneDim/2 + 2:end), dft(:, 1:sceneDim/2)];
    dict = fliplr(dict(1:arrDim, :));

    %(sceneDim/2) + 1 should be omitted
end