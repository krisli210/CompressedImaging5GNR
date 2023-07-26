function W = getMF(txAnglesAtSubcarrierK, H_RX)
    W = H_RX(:, txAnglesAtSubcarrierK)';
end