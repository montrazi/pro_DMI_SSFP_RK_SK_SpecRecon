function matSpi=spiral_index(mat)

m = size(mat,1);
n = size(mat,2);

top = 1;
bottom = m;
left = 1;
right = n; %To keep track of the four directions
value = 1;%Let us strat from 1
k = 0;

matSpi = zeros(m*n,2);
while true
    %Top Row First
    if left>right
        break;
    end
    for i = left:right
        k = k +1;
        matSpi(k,1)= top;
        matSpi(k,2) = i;
        value = value + 1;
    end
    top = top + 1;
    %Then The RightMost Column
    if top>bottom
        break;
    end
    for i = top:bottom
        k = k +1;
        matSpi(k,1)= i;
        matSpi(k,2) = right;
        value = value + 1;
    end
    right = right - 1;
    %Then Bottom Row
    if left>right
        break;
    end
    for i = right:-1:left
        k = k +1;
        matSpi(k,1)= bottom;
        matSpi(k,2) = i;
        value = value + 1;
    end
    bottom = bottom - 1;
    %Then The Left Column
    if top>bottom
        break;
    end
    for i = bottom:-1:top
        k = k +1;
        matSpi(k,1)= i;
        matSpi(k,2) = left;
        value = value + 1;
    end
    left = left + 1;
end

 matSpi=flipud(matSpi);