function data = read_binary(filename, Nx, Ny, isComplex)
    if nargin < 4
        isComplex = false;  % default to real
    end

    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end

    raw = fread(fid, 'double');  % read all elements
    fclose(fid);

    nPerSnapshot = Nx*Ny;

    if isComplex
        % Complex data (interleaved)
        nCols = length(raw)/(2*nPerSnapshot);
        data = complex( raw(1:2:end), raw(2:2:end) );
    else
        % Real data
        nCols = length(raw)/nPerSnapshot;
        data = raw;
    end

    % Reshape to [Nx*Ny, nCols]
    data = reshape(data, nPerSnapshot, nCols);