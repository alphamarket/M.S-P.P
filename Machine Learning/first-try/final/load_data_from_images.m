function [data] = load_data_from_images(address, label)
    imgs = getAllFiles(address);
    data = cell(size(imgs));
    resize_factor = [64 64];
    dispstat('','init');
    dispstat(sprintf('Processing %i found images...', size(imgs, 1)), 'keepthis');
    for img = 1:size(imgs, 1)
        img_adr = cell2mat(imgs(img));
        if isempty(regexp(img_adr, '.(jpg|jpeg|png|bmp)$', 'once'))
            continue
        end
        try
            dispstat(sprintf('[ %%%.2f ] processing "%s"...', 100 * img / size(imgs, 1), img_adr));
            data{img} = [cxr(imresize(imread(img_adr), resize_factor)) label];
        catch
            %delete(img_adr);
            % the file was not a valid image
            dispstat(sprintf('[ SKIPPED ] %s', img_adr), 'keepthis');
        end
    end
    data = cell2mat(data);
end