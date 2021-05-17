function y  = homemadefft(c)
%   An implementation of the Fast Fourier Transform as described in
%   textbooks such as Introduction to Applied Mathematics by Strang and
%   Introduction to Linear Algebra (4th Ed.) by Strang.

    n = length(c);
    if n == 1
        y = c;
    else
        m = n / 2;
        w_n = exp(2*pi*1i / n);
        c_even = c(1:2:n-1);
        c_odd= c(2:2:n);
        y_even = homemadefft(c_even);
        y_odd = homemadefft(c_odd);

        w = (w_n.^(0:m-1)).';
        y_top = y_even + w.* y_odd;
        y_bottom = y_even - w.* y_odd;
        y = [y_top; y_bottom];
    end
end

