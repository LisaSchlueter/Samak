function g = gaussian(x, pos, sigma)
     g = exp( -0.5 .* ((x-pos)./sigma).^2 );
end