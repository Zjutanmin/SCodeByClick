function Mresult = ShowResultDim(result,i)
Mresult= permute(result, [i, setdiff([1:ndims(result)], i)]);
Mresult= reshape(Mresult, size(Mresult,1), []);