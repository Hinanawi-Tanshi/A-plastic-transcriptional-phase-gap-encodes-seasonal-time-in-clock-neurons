function plot_clock_timecourse_panel(ax, wtResult, perResult, cfg, colors, showLegend)
%PLOT_CLOCK_TIMECOURSE_PANEL Plot raw summaries and flexible periodic fits.

arguments
    ax (1,1) matlab.graphics.axis.Axes
    wtResult struct
    perResult struct
    cfg struct
    colors struct
    showLegend (1,1) logical = true
end

hold(ax, "on");

local_curve_band(ax, wtResult, colors.wt, 0.12);
local_curve_band(ax, perResult, colors.per0, 0.14);

wtSummary = summarize_clock_timecourse(wtResult.data);
perSummary = summarize_clock_timecourse(perResult.data);

if cfg.plotRawPoints
    local_raw_points(ax, wtResult.data, colors.wt, -0.10);
    local_raw_points(ax, perResult.data, colors.per0, 0.10);
end

errorbar(ax, wtSummary.ZT - 0.10, wtSummary.Mean, ...
    wtSummary.Mean - wtSummary.CI95Low, wtSummary.CI95High - wtSummary.Mean, ...
    "o", "Color", colors.wt, "MarkerFaceColor", "white", ...
    "MarkerSize", 3.5, "LineWidth", 0.8, "CapSize", 2);
errorbar(ax, perSummary.ZT + 0.10, perSummary.Mean, ...
    perSummary.Mean - perSummary.CI95Low, perSummary.CI95High - perSummary.Mean, ...
    "o", "Color", colors.per0, "MarkerFaceColor", "white", ...
    "MarkerSize", 3.5, "LineWidth", 0.8, "CapSize", 2);

wtLine = plot(ax, wtResult.fit.timeGrid, wtResult.fit.predicted, ...
    "-", "Color", colors.wt, "LineWidth", 1.4);
perLine = plot(ax, perResult.fit.timeGrid, perResult.fit.predicted, ...
    "-", "Color", colors.per0, "LineWidth", 1.4);

xlim(ax, [0, cfg.period]);
ylim(ax, [0, 1]);
xticks(ax, 0:4:cfg.period);
xlabel(ax, "ZT (h)");
ylabel(ax, "Signal fraction");
title(ax, sprintf("%s cells: observed data and %s fit", ...
    upper(wtResult.signal), strrep(cfg.modelName, "_", " ")));

if showLegend
    legend(ax, [wtLine, perLine], ["Wild type", "per^0"], ...
        "Location", "best", "Box", "off");
end

local_style_axes(ax);
end

function local_curve_band(ax, result, color, alphaValue)
grid = result.amplitudeBootstrap.timeGrid;
lower = result.amplitudeBootstrap.curveCI95Low;
upper = result.amplitudeBootstrap.curveCI95High;
fill(ax, [grid; flipud(grid)], [lower; flipud(upper)], color, ...
    "FaceAlpha", alphaValue, "EdgeColor", "none", "HandleVisibility", "off");
end

function local_raw_points(ax, data, color, xOffset)
timeLevels = unique(data.time, "sorted");
for timeIndex = 1:numel(timeLevels)
    values = data.signal(data.time == timeLevels(timeIndex));
    if numel(values) == 1
        jitter = 0;
    else
        jitter = linspace(-0.16, 0.16, numel(values))';
    end
    scatter(ax, timeLevels(timeIndex) + xOffset + jitter, values, 7, color, ...
        "filled", "MarkerFaceAlpha", 0.12, "MarkerEdgeAlpha", 0.12, ...
        "HandleVisibility", "off");
end
end

function local_style_axes(ax)
ax.FontName = "Arial";
ax.FontSize = 7;
ax.LineWidth = 0.6;
ax.TickDir = "out";
ax.Box = "off";
ax.Color = "none";
end
