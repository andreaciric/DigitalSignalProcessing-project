function [y] = FI_IIR_direct_II_cascade(b, a, x)

N = length(a);
M = length(b);

Ns = floor((max(M, N)+1)/2);

fi_params = struct( 'OUT_SIGNAL_BITLENGTH', x.WordLength+b.WordLength, ...
                    'OUT_SIGNAL_FRAC', x.FractionLength+b.FractionLength);

y = fi(zeros(1, length(x)), true, fi_params.OUT_SIGNAL_BITLENGTH, fi_params.OUT_SIGNAL_FRAC, x.fimath);
delay_line = fi(zeros(max(M, N), Ns+1), x.Signed, x.WordLength, x.FractionLength, x.fimath);
y_temp = fi(zeros(1, Ns+2), true, fi_params.OUT_SIGNAL_BITLENGTH, fi_params.OUT_SIGNAL_FRAC, x.fimath);
temp = fi(zeros(1, Ns+1), true, fi_params.OUT_SIGNAL_BITLENGTH, fi_params.OUT_SIGNAL_FRAC, x.fimath);

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