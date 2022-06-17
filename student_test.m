vh = [];
vp = [];
ALL = subset20;
for i = 1:601
    sum1=[];
    sum2=[];
    for j =1:129
        if ALL(j,603)==1
            sum1 = [sum1 ALL(j,i)];
        else
            sum2= [sum2 ALL(j,i)];

        end
    end
   [h,p] = ttest2(sum1,sum2,'Alpha',0.01);
   vh = [vh h];
   vp = [vp p];
end

%%
list_94 = [];
for i=1:601
    if vh(i)==1
        list_94=[list_94 i];
    end
end
%%

ALL = ALL(:,list_94);

ALL(isinf(ALL)|isnan(ALL)) = 0;

for i=1:94
   FEATURES{i,1}=ALL(:,i)';
end

subset_94 = data_to_paint(:,list_94);
F_names = subset_94.Properties.VariableNames;

for i = 1:94
    Features{i}.Feature = FEATURES{i};
    Features{i}.str = F_names(i);
end

Features = Features';