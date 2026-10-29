temp=whos;

I=find(arrayfun(@(x) strcmp(x.class,'gpuArray'),temp));

for i=1:length(I)
    evalin('base',[temp(I(i)).name ' = gather(' temp(I(i)).name ');'])
end