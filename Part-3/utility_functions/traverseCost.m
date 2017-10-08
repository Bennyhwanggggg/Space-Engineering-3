function cost = traverseCost(cost1,cost2)
%TRAVERSECOST calculates the cost it takes to traverse from the centre of
%one cell to the centre of a neighbouring cell.

    cost = (cost1 + cost2)/2;

end