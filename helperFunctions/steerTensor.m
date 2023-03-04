function steered = steerTensor(tens, CB, I, mags)
    % Steer a tensor towards I indices of CB 
    % dim 1 of tens is assumed as space dimension

    sv = sum( mags.' .* CB(:, I), 2);
    sv = sv ./ norm(sv); 

    steered = tensorprod(tens, sv, 1);
end