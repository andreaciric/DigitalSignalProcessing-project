function [y] = IIR_direct_II(b, a, x)

M = length(b);
N = length(a);

num = max(M, N);

delay_line = zeros(num, 1);

temp = zeros(1, length(x)); 
y = zeros(1, length(x));

b = b./a(1);
a = a./a(1);
    
for i = 1:length(x)
    delay_line = [temp(i); delay_line(1:end-1)];
    temp(i) = x(i)-a(2:N)*delay_line(2:N);
    delay_line = [temp(i); delay_line(2:end)];
    y(i) = b*delay_line(1:M);
end

end
