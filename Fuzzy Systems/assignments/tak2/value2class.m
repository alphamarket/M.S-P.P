function z = value2class(y)
    z = ones(size(y));
    for i=1:size(z, 2)
        if(y(i) < 0)
            z(i) = 2;
        end
    end
end