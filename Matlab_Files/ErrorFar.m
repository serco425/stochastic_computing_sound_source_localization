function [ Error ] = ErrorFar( phi, theta , D_vec, std )

%theta = elevation
%phi = azimuth

Doa = [cos(theta)*cos(phi) cos(theta)*sin(phi) sin(theta)];

tau_calc = -(Doa*std.Ind_Pair_Connecting_Vectors)/std.c;
tau_mesh = D_vec/std.fs;

Error = sqrt((tau_mesh-tau_calc)*transpose((tau_mesh-tau_calc)));

end

