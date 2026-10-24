function smoothed = compute_NW_activity(timeMinutes, activity, sigmaMinutes, targetMinutes)
% COMPUTE_NW_ACTIVITY Nadaraya-Watson regression for fractional activities.
% Source: current Figure 3 LP/Fig 3C/fig3c.m. No binarization or hard cutoff.
% Caller excludes globally missing source rows and fills remaining NaNs as 0.
% Local weight normalization differs from compute_KDE's fixed normalization.
smoothed = NaN(numel(targetMinutes),size(activity,2));
for i = 1:numel(targetMinutes)
    weights = exp(-((timeMinutes-targetMinutes(i)).^2)/(2*sigmaMinutes^2));
    total = sum(weights);
    if total > 1e-9
        smoothed(i,:) = sum(activity.*weights,1)/total;
    end
end
end
