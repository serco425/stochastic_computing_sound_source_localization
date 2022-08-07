function [prime_vec] = getPrimes(lower_than,N)
%getPrimes get the N primes lower than max

prime_vec = zeros(1,N);
found = 0;
for curNum = 1:lower_than
   if lower_than-curNum > 2 && found < N
       if isprime(lower_than-curNum)
           found = found+1;
           prime_vec(found) = lower_than-curNum;
       end
   else
       break;
   end
end

end

