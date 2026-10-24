function [MET_z, MET_p, z_near, z_far, boundaryZ_lower, ...
    boundaryZ_upper] = computeWaldMET(observedDiff, bootstrapSE, delta0, tailMethod)
% Wald-type minimum-effect test for the interval null:
%   H0: |Delta| <= delta0
%   H1: |Delta| >  delta0
% The least-favorable / nearest null boundary defines the Wald distance.

    if nargin < 4, tailMethod = "complement"; end
    absoluteDiff = abs(observedDiff);
    if ~isfinite(observedDiff) || ~isfinite(bootstrapSE) || bootstrapSE <= 0
        boundaryZ_lower = NaN;
        boundaryZ_upper = NaN;
        z_near = NaN;
        z_far = NaN;
        MET_z = NaN;
        MET_p = NaN;
        return;
    end

    boundaryZ_lower = (-delta0 - observedDiff) / bootstrapSE;
    boundaryZ_upper = ( delta0 - observedDiff) / bootstrapSE;
    z_near = (absoluteDiff - delta0) / bootstrapSE;
    z_far = (absoluteDiff + delta0) / bootstrapSE;

    if absoluteDiff <= delta0
        MET_z = 0;
        MET_p = 1;
    else
        MET_z = z_near;
        if tailMethod == "survival"
            MET_p = normcdf(-z_near) + normcdf(-z_far);
        else
            MET_p = 1 - normcdf(z_near) + normcdf(-z_far);
        end
        MET_p = min(max(MET_p, 0), 1);
    end
end
