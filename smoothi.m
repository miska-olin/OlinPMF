function  out  = smoothi( in,span )

    temp=movmean(in,[0 span-1]);
    out=downsample(temp,span);

end

