function [B] = findDampingBasis(xd)
    y1 = 1;
    y2 = -xd(1)/(xd(2)+eps);
    y = [y1;y2]; % y is perpendicular to xd
    B = [xd./norm(xd), y./norm(y)];
end

