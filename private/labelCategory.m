function labelCategory(f,cat,pos,corner)
% LABELCATEGORY(F,CAT,POS,CORNER)
% Labels output from griddedSpots.m with spots denoting which spots belong
% to which categories. Input includes:
%   f: The figure axis handle.
%   cat: Structure containing category information, including .list and
%       .colour
%   pos: Locations of each of the spots in the grid
%   corner: {1}, 2, 3, or 4. Position to print label, starting in top left
%       moving clockwise
% C A Smith 2019

w = pos(1,3);
pos = pos(:,1:2);

col = cat.colour;

catl = cat.list(:)';

switch corner
    case 1
        x = +2;
        y = +2;
    case 2
        x = w-3;
        y = +2;
    case 3
        x = w-3;
        y = w-3;
    case 4
        x = +2;
        y = w-3;
end

for i = catl
    ipos = pos(i,:);
    plot(f,ipos(1)+x,ipos(2)+y,'o','MarkerFaceColor',col,'MarkerEdgeColor','w','MarkerSize',10);
end
    
    

