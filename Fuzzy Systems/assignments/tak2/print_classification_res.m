function print_classification_res(xp, z)
    for i=1:size(z,1)
        fprintf('xp[%i] = [%.2f, %.2f] belongs to class G%i with %.2f confidence.\n', i, xp(i, 1), xp(i, 2), z(i, 1), z(i, 2));
    end
end