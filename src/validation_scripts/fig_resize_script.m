%h = openfig(filename);
savename = 'fixedxy44.pdf';
set(h,'Visible','on')

hsurf = findobj(h,'Type','surface');
haxes = findobj(h,'Type','axes');

% take the main surface

X = get(hsurf(1),'XData');
Y = get(hsurf(1),'YData');
Z = get(hsurf(1),'ZData');

% shift by grid indices
fig2 = figure;
set(fig2,'Position',[100 100 900 750])
sx = 32;   % shift sx columns in x
sy = 24;   % shift sy rows in y
Zshift = circshift(Z,[sy sx]);
surfc(X,Y,Zshift);
shading('interp')
colormap(redblue)
view(2)
axis([0 sqrt(16)*1.02 0 sqrt(16)*1.02])
axis equal

% Copy styling / labels manually
title(haxes.Title.String,'Interpreter',haxes.Title.Interpreter);
subtitle(haxes.Subtitle.String,'Interpreter',haxes.Subtitle.Interpreter);
xlabel(haxes.XLabel.String,'Interpreter',haxes.XLabel.Interpreter);
ylabel(haxes.YLabel.String,'Interpreter',haxes.YLabel.Interpreter);
zlabel(haxes.ZLabel.String,'Interpreter',haxes.ZLabel.Interpreter);

lgd_old = findobj(h,'Type','legend');
if ~isempty(lgd_old)
    legend(lgd_old.String, 'Location', lgd_old.Location);
end

exportgraphics(fig2,savename)
