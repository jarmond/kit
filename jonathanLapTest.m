
% check that lap and lapMosek are doing roughly the same thing
% define a simple problem that both can be compared on
%%%%%%%%%%%%%

cc = [[-1 0.3 0.6 -1];
      [0.3,-1,0.9,-1];
      [0.6,0.9,-1,0.3];
      [-1,-1,0.3,-1]];
nlm = -1;
tic; [xMosek,yMosek] = lapMosek(cc, nlm, 1, 1); toc
tic; [xJV,yJV] = lap(cc,nlm,0,1,1); toc %note lap has extra argument

assert(all(xMosek==xJV),'Mosek and JV give the same solution');
