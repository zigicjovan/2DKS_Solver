function write_binary(data, filename)
    fid = fopen(filename, 'w');  % explicitly binary
    if fid == -1
        error('Could not open file for writing: %s', filename);
    end

    if ~isreal(data)
        tmp = zeros(numel(data)*2,1);
        tmp(1:2:end) = real(data(:));
        tmp(2:2:end) = imag(data(:));
        fwrite(fid, tmp, 'double');
    else
        fwrite(fid, data(:), 'double');
    end
    fclose(fid);
end
