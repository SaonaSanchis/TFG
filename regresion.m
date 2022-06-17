%% Multiple Regression

 bl = [];
 bintl = [];
 rl = [];
 rintl= [];
 statsl = [];

for i=1:94
    y = subset(1:129,i);
    M = F_data(1:129,1:7);
    
    [b,bint,r,rint,stats] = regress(y,M);

    bl = [bl b];
    bintl = [bintl bint];
    rl = [rl r];
    rintl = [rintl rint];
    statsl = [statsl stats];

end
%% select top 10 variables for each dim
for i=1:7
    dim = bl(i,1:94);
    [sortvals, sortidx] = sort(abs(dim),'descend');
    dims_id{i} = sortidx(1:10);
    dims_vals{i} =  sortvals(1:10);
end

for i = 1:7
    j = dims_id{i};
    for n=1:10
        name = j(1,n);
        dims_var{i,n} = F_names(name);
    end
end