function plot_amplitude_bootstrap_panel(ax, wtResult, perResult, delta, colors)
%PLOT_AMPLITUDE_BOOTSTRAP_PANEL Compare bootstrap peak-to-trough amplitude.

arguments
    ax (1,1) matlab.graphics.axis.Axes
    wtResult struct
    perResult struct
    delta (1,1) double
    colors struct
end

hold(ax, "on");
wtAmplitude = wtResult.amplitudeBootstrap.amplitude;
perAmplitude = perResult.amplitudeBootstrap.amplitude;
perUpper = perResult.amplitudeBootstrap.amplitudeUpper95;

upperX = max([wtAmplitude; perAmplitude; delta; perUpper]);
if upperX <= 0
    upperX = 1;
end
edges = linspace(0, min(1, upperX * 1.08), 29);
if numel(unique(edges)) < 3
    edges = linspace(0, 1, 29);
end

wtHistogram = histogram(ax, wtAmplitude, edges, "Normalization", "probability", ...
    "FaceColor", colors.wt, "FaceAlpha", 0.42, "EdgeColor", "none");
perHistogram = histogram(ax, perAmplitude, edges, "Normalization", "probability", ...
    "FaceColor", colors.per0, "FaceAlpha", 0.52, "EdgeColor", "none");

deltaLine = xline(ax, delta, "--", "Color", [0.1, 0.1, 0.1], ...
    "LineWidth", 1.0);
upperLine = xline(ax, perUpper, ":", ...
    "Color", colors.per0, "LineWidth", 1.2);

xlim(ax, [0, edges(end)]);
xlabel(ax, "Bootstrap peak-to-trough amplitude");
ylabel(ax, "Bootstrap probability");
title(ax, sprintf("%s cells: amplitude equivalence", upper(wtResult.signal)));
legend(ax, [wtHistogram, perHistogram], ["Wild type", "per^0"], ...
    "Location", "best", "Box", "off");

equivalent = perUpper < delta;
if equivalent
    conclusion = sprintf("U_{95}=%.3g < \\Delta=%.3g\npractically arrhythmic", ...
        perUpper, delta);
else
    conclusion = sprintf("U_{95}=%.3g \\geq \\Delta=%.3g\nequivalence not shown", ...
        perUpper, delta);
end
text(ax, 0.98, 0.94, conclusion, "Units", "normalized", ...
    "HorizontalAlignment", "right", "VerticalAlignment", "top", ...
    "FontName", "Arial", "FontSize", 7);

deltaLine.HandleVisibility = "off";
upperLine.HandleVisibility = "off";
local_style_axes(ax);
end

function local_style_axes(ax)
ax.FontName = "Arial";
ax.FontSize = 7;
ax.LineWidth = 0.6;
ax.TickDir = "out";
ax.Box = "off";
ax.Color = "none";
end
