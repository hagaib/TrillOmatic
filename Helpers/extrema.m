function [ ext , ind] = extrema( x , type , boundaryflag)
% Returns array of extremal values in "continuous" signal x
% using forward difference scheme for derivative
% Arguments:
%   x : signal
%   type : 1=maxima -1 =minima 0=both
%   boundaryflag : 1 = include boundary values. does not apply to type 0
%   cases fix: it does
%                  0 = exclude boundary values.
% Return Value:
%   ext: array of maxima values

    if(iscolumn(x)) , dim = 1; else dim=2; end
    
    dx = x(2:end) - x(1:end-1); %forward difference
    cross = dx(1:end-1) .* dx(2:end);
    ind = find(cross<=0)+1;
    
    if type==1 , ext_ind = find(dx(ind)<0);  
    elseif type==-1 , ext_ind = find(dx(ind)>0); 
    elseif type==0 , ext_ind = 1:length(ind);
    end;
    
%     ext = x(ind(ext_ind));

    ind = ind(ext_ind);
    ext = x(ind);
    

    if (boundaryflag) 
        if (isempty(ext))
            switch (type)
            case 1
                ext = max(x(1) , x(end)); 
            case -1
                ext = min(x(1) , x(end)); 
            case 0
                ext = x([1, end]);
            end
        else
            switch (type)
            case 1
                if(x(1)>x(2)) , ext = cat(dim , x(1) , ext); ind = cat(dim ,1, ind); end
                if(x(end)>x(end-1)) , ext = cat(dim , ext , x(end));  ind = cat(dim ,ind , length(x)); end
            case -1
                if(x(1)<x(2)) , ext = cat(dim , x(1) , ext); ind = cat(dim ,1, ind); end
                if(x(end)<x(end-1)) , ext = cat(dim , ext , x(end));  ind = cat(dim ,ind , length(x)); end
            case 0
                ext = cat(dim , x(1) , ext , x(end));
                ind = cat(dim ,1 , ind , length(x));
            end
        end
    end    

%     % monotonic signal
%     if isempty(ext)
%         switch (type)
%             case 1
%                 if x(1) > x(end) , ext = x(1); else ext = x(end); end
%             case -1
%                 if x(1) < x(end) , ext = x(1); else ext = x(end); end
%             case 0
%                 ext = x([1, end]);
%         end
% %         ext = x([1, end]);
%     end
end

