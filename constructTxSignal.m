function x = constructTxSignal(prm, H_TX)
    s = zeros(prm.NumUsers, prm.Ns * prm.NumPackets);
    x = zeros(prm.NumBsElements, prm.Ns * prm.NumPackets);
    qpskmod = comm.QPSKModulator('BitInput',true);   
    
%     A_tx = (1/sqrt(prm.NumBsElements)) .* exp(-1j * 2*pi * (0:prm.NumBsElements-1).' * (0:prm.NumBsElements-1) ./ prm.NumBsElements);
    nDirs = size(H_TX, 2);
    for p = 1:prm.NumPackets
        cols = randperm(nDirs, prm.NumUsers); % DFT based precoding
        F = 1./sqrt(prm.NumUsers) * H_TX(:, cols); %power normalization here might be wrong
        % Above, assume an equal power allocation across users
    
        for u = 1:prm.NumUsers
            bits = randi(2, prm.M*prm.Ns, 1) - 1; 
            s(u, (p-1)*prm.Ns+1:p*prm.Ns) = qpskmod(bits);
        end
        x(:, (p-1)*prm.Ns+1:p*prm.Ns) = F * s(:, (p-1)*prm.Ns+1:p*prm.Ns);
        % x is NumBsElements x Ns
    end
end