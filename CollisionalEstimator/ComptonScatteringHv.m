function hvPrimeReturn = ComptonScatteringHv(hv, theta)
hvPrimeReturn = hv / (1+(1-cos(theta))/0.511 * hv);
end

