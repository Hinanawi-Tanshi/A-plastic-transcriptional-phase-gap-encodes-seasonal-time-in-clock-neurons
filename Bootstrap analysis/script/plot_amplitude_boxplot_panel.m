function plot_amplitude_boxplot_panel(ax, results, adjustedPValues, colors, showScatter)
%PLOT_AMPLITUDE_BOXPLOT_PANEL Plot bootstrap amplitudes as points and boxplots.
%
% Comparisons shown:
%   WT M versus per0 M
%   WT E versus per0 E

arguments
    ax (1,1) matlab.graphics.axis.Axes
    results struct
    adjustedPValues (2,1) double
    colors struct
    showScatter (1,1) logical = true
end

hold(ax, "on");

amplitudes = {
    results.wt.m.amplitudeBootstrap.amplitude
    results.wt.e.amplitudeBootstrap.amplitude
    results.per0.m.amplitudeBootstrap.amplitude
    results.per0.e.amplitudeBootstrap.amplitude
};
positions = [1, 2, 4, 5];
boxColors = {
    colors.per0
    colors.wt
    colors.per0
    colors.wt
};

% Use a local stream so the visual jitter is reproducible without changing
% the random-number state used by the statistical analysis.
jitterStream = RandStream("mt19937ar", "Seed", 51827);
jitterHalfWidth = 0.17;

for groupIndex = 1:numel(amplitudes)
    values = amplitudes{groupIndex};
    values = values(isfinite(values));
    values = values(:);

    if showScatter
        jitter = (2 * rand(jitterStream, numel(values), 1) - 1) * jitterHalfWidth;
        scatter(ax, positions(groupIndex) + jitter, values, 6, ...
            "MarkerFaceColor", boxColors{groupIndex}, ...
            "MarkerEdgeColor", "none", ...
            "MarkerFaceAlpha", 0.14);
    end

    % Draw the boxplot after the points so its quartiles, median, and
    % whiskers remain legible over the bootstrap cloud.
    boxchart(ax, repmat(positions(groupIndex), numel(values), 1), values, ...
        "BoxFaceColor", boxColors{groupIndex}, ...
        "BoxFaceAlpha", 0.55, ...
        "BoxWidth", 0.58, ...
        "MarkerStyle", "none", ...
        "LineWidth", 0.8);
end

allValues = vertcat(amplitudes{:});
dataMaximum = max(allValues);
verticalStep = max(0.045, 0.07 * dataMaximum);
bracketHeight = max(0.012, 0.018 * dataMaximum);
firstY = dataMaximum + verticalStep;
secondY = firstY + 1.65 * verticalStep;

local_add_significance_bracket(ax, positions(1), positions(3), ...
    firstY, bracketHeight, local_significance_label(adjustedPValues(1)));
local_add_significance_bracket(ax, positions(2), positions(4), ...
    secondY, bracketHeight, local_significance_label(adjustedPValues(2)));

xlim(ax, [0.4, 5.6]);
ylim(ax, [0, secondY + 1.25 * verticalStep]);
xticks(ax, positions);
xticklabels(ax, {"WT M", "WT E", "per^0 M", "per^0 E"});
% ylabel(ax, "Bootstrap peak-to-trough amplitude");

ax.FontName = "Arial";
ax.FontSize = 8;
ax.LineWidth = 0.6;
ax.TickDir = "out";
ax.Box = "off";
ax.Color = "none";
end

function local_add_significance_bracket(ax, xLeft, xRight, y, height, label)
plot(ax, [xLeft, xLeft, xRight, xRight], ...
    [y, y + height, y + height, y], ...
    "-", "Color", [0.15, 0.15, 0.15], "LineWidth", 0.8);
text(ax, mean([xLeft, xRight]), y - 1.5 * height, label, ...
    "HorizontalAlignment", "center", ...
    "VerticalAlignment", "bottom", ...
    "FontName", "Arial", "FontSize", 8);
end

function label = local_significance_label(pValue)
if pValue < 1e-4
    label = "****";
elseif pValue < 1e-3
    label = "***";
elseif pValue < 1e-2
    label = "**";
elseif pValue < 0.05
    label = "*";
else
    label = "ns";
end
end
