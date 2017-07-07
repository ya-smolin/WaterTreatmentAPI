%Conventions
% !!!input
% C = C-NPi for i=1:m, substrate concentration [mg/l]
% T = ti for i=1:m, time [h]
% X = xi for i=1:m, Volatile Suspended Solids concentrations [mg/l]
% v = V, water treatment plant volume [l]
% ???output
% k1 = k-max, maximum specific growth rate [1/h]
% k2 = K-s, half-velocity concentration [mg/l]
% k3 = K-I, inhibition constant, the higher value of Ki means a less inhibitive substrate.[mg/l]
%================================
function data = loadData()
disp('load data');
%% 2ChlorPhenol
v = 2.25;
data{1}.C = [91.5 78 74.5 56 42.5 10 1 0.3]';
data{1}.T = [1e-10 26.4 48 69.6 81.6 112.8 124.8 144]';
data{1}.X = [205 203 186 188 189 184 210 187]';
data{1}.title='2ChlorPhenol-90';
data{1}.v = v;

data{2}.C = [80.1 69.5 65.2 55.5 22 0.2]';
data{2}.T = [0 14.88 24 57.6 96 115.2]';
data{2}.X = [196 195 199 213 248 248]';
data{2}.title='2ChlorPhenol-80';
data{2}.v = v;

data{3}.C = [54.5 54 50.8 39 30 9 0.01]';
data{3}.T = [0 6.96 24 54.96 69.6 96 110.4]';
data{3}.X = [186 183 197 198 218 227 227]';
data{3}.title='2ChlorPhenol-60';
data{3}.v = v;

data{4}.C = [36.3 35.5 31 25.1 13.8 4 0.01]';
data{4}.T = [0 0.167 1 2 3 3.5 4]'*24;
data{4}.X = [207 204 225 230 240 246 253]';
data{4}.title='2ChlorPhenol-40';
data{4}.v = v;

data{5}.C = [23.8 21.7 20 11 3 0.1]';
data{5}.T = [0 0.8 1 2.5 3.3 4]'*24;
data{5}.X = [234 237 249 264 263 262]';
data{5}.title='2ChlorPhenol-20';
data{5}.v = v;

data{6}.C = [24.9 20.5 18.8 14 7.7 0.5 0.01]';
data{6}.T = [0 18 24 43 55 72 91]';
data{6}.X = [197 174 183 189 190 202 220]';
data{6}.title='2ChlorPhenol-20';
data{6}.v = v;
%% series I 2NitroPhenol
vNitroPhenol = 0.15;
data{7}.title='2NitroPhenol-52-B1';
data{7}.C = [52 51 50 46 45 42 41 40 38.9 35.8 21 12.4 10.6 8.1 5.2 3.3 0.5]';
data{7}.T = [0 4 8 22 30 46 50 58 60 70 80 90 94 100 110 120 142]';
data{7}.X = repmat(675,length(data{7}.T),1);
data{7}.v = vNitroPhenol;

data{8}.title='2NitroPhenol-53.5-B2';
data{8}.C = [53.5 52 50.3 44.7 41.8 35 31 22 20 14.6 9 6.1 5 4 2 1.3 0.3]';
data{8}.T = data{7}.T;
data{8}.X = repmat(565,length(data{7}.T),1);
data{8}.v = vNitroPhenol;

data{8}.title='2NitroPhenol-54.5-BC';
data{8}.C = [54.5 53.6 53 48.7 45 38.7 35 28 25 15.8 11 9 8.4 7 3.5 1.3 2.1]';
data{8}.T = data{7}.T;
data{8}.X = repmat(625,length(data{7}.T),1);
data{8}.v = vNitroPhenol;
%% series II 2NitroPhenol
data{9}.title ='2NitroPhenol-44.8-B20';
data{9}.C = [44.8 44.1 27.5 6.2 3.7 1.1 3 3]';
data{9}.T = [0 2 24 48 74 94 120 144]'; 
data{9}.X = repmat(260,length(data{9}.T),1);
data{9}.v = vNitroPhenol;

data{10}.title = '2NitroPhenol-45-B21';
data{10}.C = [45 44 15 4.3 3.4 4.7 4.2 1.8]';
data{10}.T = data{9}.T; 
data{10}.X = repmat(545,length(data{9}.T),1);
data{10}.v = vNitroPhenol;

data{11}.title = '2NitroPhenol-140-B22';
data{11}.C = [140 137 67.6 6.2 4.3 3.1 2 2]';
data{11}.T = data{9}.T; 
data{11}.X = repmat(564,length(data{9}.T),1);
data{11}.v = vNitroPhenol;

data{12}.title  ='2NitroPhenol-142.7-B23';
data{12}.C = [142.7 141.8 98 9.6 5 1.5 1.6 2]';
data{12}.T = data{9}.T; 
data{12}.X = repmat(310,length(data{9}.T),1);
data{12}.v = vNitroPhenol;

data{13}.title = '2NitroPhenol-42.5-m24';
data{13}.C = [42.5 42 37.3 27.2 22.4 17.3 6.6 0.2]';
data{13}.T = data{9}.T; 
data{13}.X = repmat(17, length(data{9}.T),1);
data{13}.v = vNitroPhenol;

data{13}.title = '2NitroPhenol-41.1-m25';
data{13}.C = [41.1 40.7 37.7 8.1 3.8 2.6 0.6 0.8]';
data{13}.T = data{9}.T; 
data{13}.X = repmat(25,length(data{9}.T),1);
data{13}.v = vNitroPhenol;

%% III series 2NitroPhenol
Tobs = [0  22  45  70  93  117  138]';

dataP.title = '2NitroPhenol-51-b31';
dataP.C = [51  21.1  3.1  0.9]';
dataP.T = Tobs(1:length(dataP.C)); 
dataP.X = repmat(207, length(dataP.C), 1);
dataP.v = vNitroPhenol;
data{14} = dataP;

dataP.title = '2NitroPhenol-102-b32';
dataP.C = [102  75.7  22.3  1.9  1.8  2]';
dataP.T = Tobs(1:length(dataP.C)); 
dataP.X = repmat(207, length(dataP.C), 1);
dataP.v = vNitroPhenol;
data{15} = dataP;

dataP.title = '2NitroPhenol-200-b33';
dataP.C = [200  180  164  36.8  3.5  4  1.4]';
dataP.T = Tobs(1:length(dataP.C)); 
dataP.X = repmat(207, length(dataP.C), 1);
dataP.v = vNitroPhenol;
data{16} = dataP;

dataP.title = '2NitroPhenol-51.8-m31';
dataP.C = [51.8  31.4  5.7  0.9  0.4  1]';
dataP.T = Tobs(1:length(dataP.C)); 
dataP.X = repmat(47, length(dataP.C), 1);
dataP.v = vNitroPhenol;
data{17} = dataP;

dataP.title = '2NitroPhenol-103.5-m32';
dataP.C = [103.5  100.6  83.3  8.7  2.3  2]';
dataP.T = Tobs(1:length(dataP.C)); 
dataP.X = repmat(60, length(dataP.C), 1);
dataP.v = vNitroPhenol;
data{18} = dataP;

dataP.title = '2NitroPhenol-207.2-m33';
dataP.C = [207.2  206.6  204.6  187.2  159.8  16.1  1.9]';
dataP.T = Tobs(1:length(dataP.C)); 
dataP.X = repmat(87, length(dataP.C), 1);
dataP.v = vNitroPhenol;
data{18} = dataP;

end

