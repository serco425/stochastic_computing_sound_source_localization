%vec1 = [0 1 1 0 0 0 1 1 0 1 0 1 0 1 1 1 1 0 0 0 ];
%vec2 = [1 0 1 1 1 1 1 1 0 1 0 1 0 1 1 1 1 1 1 1 ];

vec1 = [0 1 0 0 1 0 1 0];
vec2 = [0 0 1 0 0 0 1 0];

mat = [vec1;vec2];

N_bits = length(mat);

res1 = scTFFMuxMultiAdd(mat,N_bits, 1);
res2 = scTFFMuxMultiAdd(mat,N_bits, 0);

res1_binary = Unary2Binary(res1);
res2_binary = Unary2Binary(res2);