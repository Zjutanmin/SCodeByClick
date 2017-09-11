function QueryFea = Mynomalize(QueryFea, dim, normweight)
if normweight == 1
    QueryFea = bsxfun(@times, QueryFea, 1./ sum(QueryFea,dim));
else if normweight == 2
        QueryFea = bsxfun(@times, QueryFea, 1./ sqrt(sum(QueryFea.^2,dim)));
    end
end