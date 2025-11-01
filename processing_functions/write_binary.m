function write_binary(filename,data)
    fid = fopen(filename, 'w');  % explicitly binary
    if fid == -1
        error('Could not open file for writing: %s', filename);
    end

    % Write complex numbers interleaved (real, imag)
    if ~isreal(data)
        fwrite(fid, [real(data(:)) imag(data(:))], 'double');
    else
        fwrite(fid, data(:), 'double');
    end
    fclose(fid);
end