function plot_file(fname, col, scale, labeltxt, varargin)
    data = load(fname);   % assumes whitespace-delimited numeric file
    % 
    seldata = data(:,col);
    if col == 5
        seldata = data(:,4) - data(:,2) - 2*(data(:,3) - data(:,2)) ; % lap
    elseif col == 6
        seldata = data(:,3) - data(:,2) ; % grad
    elseif col == 7
        seldata = data(:,4) - data(:,2) - 2*(data(:,3) - data(:,2)) ;
        seldata = (seldata - data(:,2)) ./ data(:,2) ; % lap/phi
        scale = 1;
    end
    
    plot(data(:,1), scale*seldata, ...
        'DisplayName', labeltxt, ...
        varargin{:});   % <-- pass styling through
end