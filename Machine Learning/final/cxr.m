function cxrsig = cxr(image)
    a = 0.3;
    b = a;
    cxrsig = conv(cp(image, a, b), rp(image, a, b));
end

