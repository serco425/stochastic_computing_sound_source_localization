function [ stdvalues ] = generate_std( Array_Geo_input )
%% Generates Array Geometry related data, such as angles between microphone Pairs
%% Inputs: 
%%% 1) Array_Geo_input is (M x 5),  x,y,z,phi,polarangle(90-elevation)
%% Outputs:
%%% 1) stdvalues: 1x1 struct with data 
%%
Array_Angles = transpose(Array_Geo_input(:,4:5));
Array_Geo = Array_Geo_input(:,1:3);
Array_Geo = transpose(Array_Geo(:,1:3)); %3xM

stdvalues.m = length(Array_Geo(1,:));

XUnitVector = [1;0;0];
ind = 1;
for i1 = 1:stdvalues.m-1 
    for i2 = i1+1:stdvalues.m
       Pair_Vectors_All(:,ind) = Array_Geo(:,i1)-Array_Geo(:,i2);
       ind = ind +1;
    end
end

for i = 1:(stdvalues.m*(stdvalues.m -1))/2;
PairAngles(i) = acos(dot(XUnitVector,Pair_Vectors_All(:,i)) / norm(Pair_Vectors_All(:,i)));
DistPairs(i) = norm(Pair_Vectors_All(:,i)); %distances between microphone pairs
end

stdvalues.maxlagcircular7 = 3;
stdvalues.indpairs = (stdvalues.m*(stdvalues.m -1))/2;
stdvalues.Ind_Pair_Connecting_Vectors = Pair_Vectors_All;
stdvalues.IndPairDists = DistPairs;
stdvalues.IndPairAnglesToXAxis = PairAngles;
stdvalues.Array_Geo = Array_Geo;
stdvalues.Array_Angles = pi / 180 * (Array_Angles);
stdvalues.radiusmean = 0.033;

end
