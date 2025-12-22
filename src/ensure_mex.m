function ensure_mex(srcfile, outdir)

    if isdeployed
        return;
    end

    [~,name] = fileparts(srcfile);
    mexfile = fullfile(outdir, [name '.' mexext]);

    srcinfo = dir(srcfile);
    if isempty(srcinfo)
        error('Source file not found: %s', srcfile);
    end

    mexinfo = dir(mexfile);
    if isempty(mexinfo) || srcinfo.datenum > mexinfo.datenum
        fprintf('Compiling %s...\n', name);
        mex('-O', ...
            'CXXFLAGS=$CXXFLAGS -fopenmp', ...
            'LDFLAGS=$LDFLAGS -fopenmp', ...
            '-outdir', outdir, ...
            srcfile);
    end
end