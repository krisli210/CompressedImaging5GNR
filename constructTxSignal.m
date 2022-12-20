function x = constructTxSignal(prm)
    s = zeros(prm.NumUsers, prm.Ns * prm.NumPackets);
    x = zeros(prm.NumBsElements, prm.Ns * prm.NumPackets);
    qpskmod = comm.QPSKModulator('BitInput',true);   
    
    A_tx = (1/sqrt(prm.NumBsElements)) .* exp(-1j * 2*pi * (0:prm.NumBsElements-1).' * (0:prm.NumBsElements-1) ./ prm.NumBsElements);
    for p = 1:prm.NumPackets
        cols = randperm(prm.NumBsElements, prm.NumUsers); % DFT based precoding
        F = A_tx(:, cols);
    
        for u = 1:prm.NumUsers
            bits = randi(2, prm.M*prm.Ns, 1) - 1; 
            s(u, (p-1)*prm.Ns+1:p*prm.Ns) = qpskmod(bits);
        end
        x(:, (p-1)*prm.Ns+1:p*prm.Ns) = F * s(:, (p-1)*prm.Ns+1:p*prm.Ns);
        % x is NumBsElements x Ns
    end
end