function [y] = FI_IIR_direct_II(b, a, x)

N = length(a);

%changed parameters 3 times
fi_params = struct( 'OUT_SIGNAL_BITLENGTH', 27, ...
                    'OUT_SIGNAL_FRAC', 10);

FixedPointAttributes = fimath( 'RoundingMethod', 'Floor', 'OverflowAction', 'Wrap', ...
        'ProductMode', 'SpecifyPrecision', 'ProductWordLength', fi_params.OUT_SIGNAL_BITLENGTH , 'ProductFractionLength', fi_params.OUT_SIGNAL_FRAC, ...
        'SumMode', 'SpecifyPrecision', 'SumWordLength', fi_params.OUT_SIGNAL_BITLENGTH, 'SumFractionLength', fi_params.OUT_SIGNAL_FRAC ) ;
    
% delay_line = fi( zeros(N,1), x.Signed, x.WordLength, x.FractionLength, FixedPointAttributes); 
% y = fi(zeros(1,length(x)), true, fi_params.OUT_SIGNAL_BITLENGTH, fi_params.OUT_SIGNAL_FRAC, FixedPointAttributes);
% temp = fi(zeros(1, length(x)), true, fi_params.OUT_SIGNAL_BITLENGTH, fi_params.OUT_SIGNAL_FRAC,  FixedPointAttributes);

delay_line = fi( zeros(N, 1), x.Signed, x.WordLength, x.FractionLength, x.fimath);
y = fi(zeros(1, length(x)), true, fi_params.OUT_SIGNAL_BITLENGTH, fi_params.OUT_SIGNAL_FRAC, x.fimath);
temp = fi(zeros(1, length(x)), true, fi_params.OUT_SIGNAL_BITLENGTH, fi_params.OUT_SIGNAL_FRAC, x.fimath);
    
for i = 1:length(x)
    delay_line = [temp(i); delay_line(1:end-1)];
    temp(i) = x(i)-a(2:N)*delay_line(2:N);
    delay_line = [temp(i); delay_line(2:end)];
    y(i) = b*delay_line;
end

end
