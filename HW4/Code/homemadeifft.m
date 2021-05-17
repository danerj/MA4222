function y  = homemadeifft(c)
%   An implementation of the inverse Fast Fourier Transform as described in
%   textbooks such as Introduction to Applied Mathematics by Strang and
%   Introduction to Linear Algebra (4th Ed.) by Strang.
%   There are two differences between this function and recursivefft. The
%   first is to replace w_n with it's complex conjugate omega_n. The second
%   is to divide y by 2 before returning y in the recursive cases. I do not
%   understand why dividing by 2 works but the inverse appears to only work
%   if this division is included. 

    n = length(c);
    if n == 1
        y = c;
    else
        m = n / 2;
        omega_n = exp(-2*pi*1i / n);
        c_even = c(1:2:n-1);
        c_odd= c(2:2:n);
        y_even = homemadeifft(c_even);
        y_odd = homemadeifft(c_odd);

        omega = (omega_n.^(0:m-1)).';
        y_top = y_even + omega.* y_odd;
        y_bottom = y_even - omega.* y_odd;
        y = [y_top; y_bottom] / 2;
    end
end

