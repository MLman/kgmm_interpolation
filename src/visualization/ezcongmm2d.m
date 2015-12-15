function h = ezcongmm2d(gmm, range)
    h = ezcontour(@(x,y)log(pdf(gmm,[x y])),range(1,:),range(2,:));
end
