function w = polynomial_window(N , alpha, inorout)
% This function creates a fade in/out window defined in the article: 
%"Reconstruction Method for Missing or Damaged Long Portions in Audio Signal"
%by ISMO KAUPPINEN AND JYRKI KAUPPINEN

% N = window length
% alpha = control parameter. alpha = 1 is a linear window and alpha->infinity is a step function

    tt = (0:N-1)/(N-1);
    midpoint = floor(N/2);
    
    w = 1 - 0.5*(2*tt(1:midpoint)).^alpha;
    w(midpoint+1:N) = 0.5*(2 - 2*tt(midpoint+1:end)).^alpha;
    
    if (strcmp(inorout , 'in'))
        w = 1-w;
    end
    
    w = w';
end