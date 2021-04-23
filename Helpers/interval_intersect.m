function result = interval_intersect(inter1 , inter2)
% calculates intersection between 2 closed intervals.
%   each argument is of the form: [a , b], and the interval represented is
%   all real values of x such that a<=x<=b.
% the empty set is represented by an empty array: [].


% example: result = interval_intersect([0 , 2] , [-0.5 , 1]) should return
% [1 , 2]

if inter1(2) < inter1(1) || inter2(2) < inter2(1)
    result = [];
    return
end

result = [max(inter1(1) , inter2(1)) , min(inter1(2) , inter2(2))];
if(result(1)>result(2))
    result = [];
end