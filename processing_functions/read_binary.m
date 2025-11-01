function data = read_binary(filename, Nx, Ny, isComplex)
    if nargin < 4, isComplex = false; end
    fid = fopen(filename, 'r');
    if fid == -1, error('Could not open file: %s', filename); end
    raw = fread(fid, 'double');
    fclose(fid);

    nPerSnapshot = Nx*Ny;
    if isComplex
        data = complex(raw(1:2:end), raw(2:2:end));
    else
        data = raw;
    end

    nCols = numel(data)/nPerSnapshot;
    data = reshape(data, nPerSnapshot, nCols);
end
