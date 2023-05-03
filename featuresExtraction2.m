function [Higuchi, output_lnk, output_lnLk] = featuresExtraction2(segments, klin, kmax)
% Function for feature extraction
%
% INPUTS:
% segments: matrix of dimension NxM where N is the
% number of segments and M the length of the segment.
%
% OUTPUTS:
% features vector (one features for each segment)

if nargin<2
    kmax = 18;
    klin = 6;
end

        Higuchi = NaN(size(segments,1),1);

        for n = 1 : size(segments,1)
            [~,lnk,lnLk] = hfd(segments(n,:),kmax);

            % linear fit
            p = polyfit(lnk(1:klin),lnLk(1:klin),1);
            y = p(1)*lnk + p(2);

            % Higuchi Fractal Dimension
            Higuchi(n) = p(1);

            % Length
            output_lnLk(n,:) = lnLk;
            output_lnk(n,:) = -lnk;

        end

end