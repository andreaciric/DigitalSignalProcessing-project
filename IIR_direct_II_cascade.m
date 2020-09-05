function [y] = IIR_direct_II_cascade(b, a, x)

N = length(a);
M = length(b);

Ns = floor((max(M, N)+1)/2);

y = zeros(1, length(x)); 

delay_line = zeros(max(M, N), Ns+1);

y_temp = zeros(1, Ns+2);
temp = zeros(1, Ns+1); 

for i = 1:length(x)
    y_temp(1) = x(i);
    for j = 1:Ns+1
        delay_line(:, j) = [temp(j) delay_line(1:end-1, j)'];
        temp(j) = y_temp(j)-a(j, 2:N)*delay_line(2:N, j);
        delay_line(:, j) = [temp(j) delay_line(2:end, j)'];
        y_temp(j+1) = b(j, :)*delay_line(1:M, j);
    end
    y(i) = y_temp(Ns+2);
end

end
