function A = columnNorm(A)
    A = A./sqrt(sum(abs(A).^2));
end
