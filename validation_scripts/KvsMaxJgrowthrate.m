
maxJopt = load([pwd '/data/maxJopt.dat']);
newmaxJopt = [ maxJopt , NaN(size(maxJopt,1),3) ];

for i = 1:size(maxJopt,1)
    newmaxJopt(i, size(maxJopt,2) + 1) = log(newmaxJopt(i,4))/log(newmaxJopt(i,1));
    newmaxJopt(i, size(maxJopt,2) + 2) = newmaxJopt(i, size(maxJopt,2) + 1)/newmaxJopt(i,2);
    newmaxJopt(i, size(maxJopt,2) + 3) = log(newmaxJopt(i,4))/log(newmaxJopt(i,1));
end

