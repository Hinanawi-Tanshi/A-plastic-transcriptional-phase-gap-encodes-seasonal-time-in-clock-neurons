function exportKDEFigure(fig, outputBase, cfg)
% EXPORTKDEFIGURE Export editable SVG and optional PNG at the figure's size.
fig.PaperUnits = 'centimeters';
position = get(fig, 'Position'); units = get(fig, 'Units');
set(fig, 'Units', 'centimeters'); physical = get(fig, 'Position');
set(fig, 'Units', units, 'Position', position);
fig.PaperSize = physical(3:4);
fig.PaperPosition = [0 0 physical(3:4)];
fig.PaperPositionMode = 'manual';
fig.InvertHardcopy = 'off';
drawnow;
print(fig, char(outputBase + ".svg"), '-dsvg', '-vector');
if cfg.exportPNG
    print(fig, char(outputBase + ".png"), '-dpng', '-r200');
end
if ~cfg.showFigures, close(fig); end
end
