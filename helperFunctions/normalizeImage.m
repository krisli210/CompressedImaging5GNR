function normImage = normalizeImage(im)
    imAbs = abs(im).^2;
    normImage = (imAbs - min(imAbs, [], 'all')) / (max(imAbs, [], 'all') - min(imAbs, [], 'all'));
end