function [y] = IIR_direct_II_cascade(b, a, x)

N = length(a);
M = length(b);

Ns = floor((max(M, N)+1)/2);

y = zeros(1, length(x)); 

delay_line = zeros(max(M, N), Ns+1);

y_temp = zeros(1, Ns+2);
w = zeros(1, Ns+1); % predstavlja medjurezultat koji se cuva unutar kaskade

% u spoljasnjoj petlji se prolazi kroz svaki odbirak ulaznog signala, dok se u
% unutrasnjoj petlji za svaki odbirak prolazi kroz sve kaskade
for k = 1:length(x)
    y_temp(1) = x(k);
    for i = 1:Ns+1
        delay_line(:, i) = [w(i) delay_line(1:end-1, i)'];
        w(i) = y_temp(i)-a(i, 2:N)*delay_line(2:N, i);
        delay_line(:, i) = [w(i) delay_line(2:end, i)'];
        y_temp(i+1) = b(i, :)*delay_line(1:M, i);
    end
    y(k) = y_temp(Ns+2);
end

end
