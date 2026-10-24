function plot_null_bootstrap_panel(ax, result, color)
%PLOT_NULL_BOOTSTRAP_PANEL Plot null-bootstrap deviance improvement.

arguments
    ax (1,1) matlab.graphics.axis.Axes
    result struct
    color (1,3) double
end

hold(ax, "on");
values = result.nullBootstrap.statistics;
observed = result.nullBootstrap.observedStatistic;
critical = prctile(values, 95);

histogram(ax, values, 28, "Normalization", "probability", ...
    "FaceColor", color, "FaceAlpha", 0.55, "EdgeColor", "none");
observedLine = xline(ax, observed, "-", ...
    "Color", [0.15, 0.15, 0.15], "LineWidth", 1.2);
criticalLine = xline(ax, critical, "--", ...
    "Color", [0.45, 0.45, 0.45], "LineWidth", 0.9);

upperX = max([values; observed; critical]);
if upperX <= 0
    upperX = 1;
end
xlim(ax, [0, upperX * 1.08]);
xlabel(ax, "\Delta deviance (flat - periodic)");
ylabel(ax, "Bootstrap probability");
title(ax, sprintf("%s %s: no-rhythm null", ...
    result.genotype, upper(result.signal)));

pValue = result.nullBootstrap.pValue;
annotationText = sprintf("T_{obs} = %.3g\nNull 95%% = %.3g\np = %.3g", ...
    observed, critical, pValue);
text(ax, 0.97, 0.93, annotationText, "Units", "normalized", ...
    "HorizontalAlignment", "right", "VerticalAlignment", "top", ...
    "FontName", "Arial", "FontSize", 7);

observedLine.HandleVisibility = "off";
criticalLine.HandleVisibility = "off";
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
