edges = [250. , 252.5, 255. , 257.5, 260. , 262.5, 265. , 267.5, 270. , ...
    272.5, 275. , 277.5, 280. , 282.5, 285. , 287.5, 290. , 292.5, ...
    295. , 297.5, 300. , 302.5, 305. , 307.5, 310. , 312.5, 315. , ...
    317.5, 320. , 322.5, 325. , 327.5, 330. , 332.5, 335. , 337.5, ...
    340. , 342.5, 345. , 347.5, 350. ]*3600/1000;

% 1hour
counts = [ 0.,  0.,  1.,  2.,  2.,  3.,  5.,  7., 16., 16., 22., 20., 29., ...
    26., 32., 41., 16., 27., 21., 15., 13., 11.,  4., 10.,  3.,  1., ...
    1.,  1.,  1.,  1.,  0.,  0.,  0.,  1.,  1.,  2.,  0.,  1.,  0., ...
    0.];

 h=histogram('BinEdges',edges,'BinCounts',counts);
myedges=(edges(2)-edges(1))/2 + edges; myedges=myedges(1:numel(counts));

%x = 0:0.1:10;
%y = gaussmf(x,[2 5]);
%plot(x,y)
% plot(myedges,normpdf(1035.89,sqrt(1035.89)));

