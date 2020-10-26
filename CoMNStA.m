
clear all;
close all;

ths = [0.2	0.25	0.15	0.4	0.175	0.25	0.275	0.325	0.3	0.2	0.275	0.2	0.35	0.4	0.4	0.25	0.275	0.2	0.3	0.25	0.2	0.275	0.325	0.275	0.125	0.25	0.2	0.2	0.175	0.15	0.3	0.275	0.175	0.175	0.25	0.2	0.225	0.25	0.175	0.275	0.3	0.25	0.275	0.3	0.325	0.175	0.175	0.225	0.35];

d = 0.1;

path = 'Images'; %path to images
fileList = dir(fullfile(path, '*.tif'));

ms = {};
Ng = {};
N2 = {};
N = {};
ccdf_woNg = {};
kk = 1;
for i = 1:length(fileList)
	clear c cb puntos L LL g
	display(fileList(i).name);

	c = imread(fullfile(path, fileList(i).name));
	[m n] = size(c);

	th = ths(i);
	th_d = [(th + d) th (th - d)];
	for k = 1:length(th_d)

		cb = im2bw(c, th_d(k));
		%cb = im2bw(c, th);
		%puntos = bwmorph(cb,'skel', Inf);

		%[L, num] = bwlabel(puntos, 8);
		[L, num] = bwlabel(cb, 8);
		LL = reshape(L, m*n, 1);

		for ii = 1:num
			g(ii) = length(find(LL == ii));
		end

		g = g(find(g > 100)); % denoising
		g = sort(g, 'descend');

		s = unique(g(2:length(g)));

		numer = [];
		denom = [];
		for ii = 1:length(s)
			Ns(ii) = length(find(g == s(ii)))/length(g(2:length(g)));
			numer(ii) = (s(ii)^2)*(Ns(ii));
			denom(ii) = s(ii)*Ns(ii);
		end

		ms{i}(k) = sum(numer)/sum(denom);
		Ng{i}(k) = g(1);
		N2{i}(k) = g(2);
		N{i}(k) = sum(g);
		ccdf_woNg{kk} = g(2:length(g));
		kk = kk + 1;
	end
end



%sort networks according to average size
avN = [];
sdN = [];
for k=1:length(N)
	avN = [avN mean(N{k})];
	sdN = [sdN std(N{k})];
end
[val id] = sort(avN); %id ranks mito nets
avN = avN(id);
sdN = sdN(id);



% CCDF
ccdf = {};
k = 1;
for i = 2:3:length(ccdf_woNg)
	ccdf{k} = ccdf_woNg{i};
	k = k + 1;
end

aux1 = [ccdf{id(1)} ccdf{id(2)} ccdf{id(3)} ccdf{id(4)} ccdf{id(5)} ccdf{id(6)} ccdf{id(7)} ccdf{id(8)} ccdf{id(9)} ccdf{id(10)}];
aux2 = [ccdf{id(11)} ccdf{id(12)} ccdf{id(13)} ccdf{id(14)} ccdf{id(15)} ccdf{id(16)} ccdf{id(17)} ccdf{id(18)} ccdf{id(19)} ccdf{id(20)}];
aux3 = [ccdf{id(21)} ccdf{id(22)} ccdf{id(23)} ccdf{id(24)} ccdf{id(25)} ccdf{id(26)} ccdf{id(27)} ccdf{id(28)} ccdf{id(29)} ccdf{id(30)}];
aux4 = [ccdf{id(31)} ccdf{id(32)} ccdf{id(33)} ccdf{id(34)} ccdf{id(35)} ccdf{id(36)} ccdf{id(37)} ccdf{id(38)} ccdf{id(39)} ccdf{id(40)}];
aux5 = [ccdf{id(41)} ccdf{id(42)} ccdf{id(43)} ccdf{id(44)} ccdf{id(45)} ccdf{id(46)} ccdf{id(47)} ccdf{id(48)} ccdf{id(49)}];

figure(1)
[bin, H] = cumhist2(aux1, max(aux1));
loglog(H, bin, 'o');
hold on
[bin, H] = cumhist2(aux2, max(aux2));
loglog(H, bin, 'o');
[bin, H] = cumhist2(aux3, max(aux3));
loglog(H, bin, 'o');
[bin, H] = cumhist2(aux4, max(aux4));
loglog(H, bin, 'o');
[bin, H] = cumhist2(aux5, max(aux5));
loglog(H, bin, 'o');



% N2
avN2 = [];
sdN2 = [];
for k=1:length(Ng)
	avN2 = [avN2 mean(N2{k})];
	sdN2 = [sdN2 std(N2{k})];
end

avN2 = avN2(id);
sdN2 = sdN2(id);

figure(2)
loglog(avN, avN2, 'o')



% susceptibility <s>
% split sets corresponding to different th
ccdf_1 = {}; % lower bound
cont = 1;
for i = 1:3:length(ccdf_woNg)
	ccdf_1{cont} = ccdf_woNg{i};
	cont = cont + 1;
end
% sort based on avN
ccdf1 = {};
for i = 1:length(ccdf_1)
	ccdf1{i} = ccdf_1{id(i)};
end

ccdf_2 = {}; % th*
cont = 1;
for i = 2:3:length(ccdf_woNg)
	ccdf_2{cont} = ccdf_woNg{i};
	cont = cont + 1;
end
% sort based on avN
ccdf2 = {};
for i = 1:length(ccdf_2)
	ccdf2{i} = ccdf_2{id(i)};
end

ccdf_3 = {}; % upper bound
cont = 1;
for i = 3:3:length(ccdf_woNg)
	ccdf_3{cont} = ccdf_woNg{i};
	cont = cont + 1;
end
% sort based on avN
ccdf3 = {};
for i = 1:length(ccdf_3)
	ccdf3{i} = ccdf_3{id(i)};
end

delta = 4;

sus = [];
meanN = [];
stdN = [];
for i = 1:(length(ccdf2) - delta)

	aux = [];
	for j = i:(i + delta)
		aux = [aux ccdf2{j}];
	end

	meanN(i) = mean(avN(i : (i + delta)));
	stdN(i) = std(avN(i : (i + delta)));

	s = unique(aux);
	numer = [];
	denom = [];
	Ns = [];
	for ii = 1:length(s)
		Ns(ii) = length(find(aux == s(ii)))/length(aux);
		numer(ii) = (s(ii)^2)*(Ns(ii));
		denom(ii) = s(ii)*Ns(ii);
	end
	sus(i) = sum(numer)/sum(denom);

end

% lower threshold (for error calculation)
sus_l = [];
for i = 1:(length(ccdf1) - delta)

	aux = [];
	for j = i:(i + delta)
		aux = [aux ccdf1{j}];
	end

	s = unique(aux);
	numer = [];
	denom = [];
	Ns = [];
	for ii = 1:length(s)
		Ns(ii) = length(find(aux == s(ii)))/length(aux);
		numer(ii) = (s(ii)^2)*(Ns(ii));
		denom(ii) = s(ii)*Ns(ii);
	end
	sus_l(i) = sum(numer)/sum(denom);

end

% upper threshold (for error calculation)
sus_u = [];
for i = 1:(length(ccdf3) - delta)

	aux = [];
	for j = i:(i + delta)
		aux = [aux ccdf3{j}];
	end

	s = unique(aux);
	numer = [];
	denom = [];
	Ns = [];
	for ii = 1:length(s)
		Ns(ii) = length(find(aux == s(ii)))/length(aux);
		numer(ii) = (s(ii)^2)*(Ns(ii));
		denom(ii) = s(ii)*Ns(ii);
	end
	sus_u(i) = sum(numer)/sum(denom);

end

av_sus = [];
sd_sus = [];
for i = 1:length(meanN)
	av_sus(i) = mean([sus_l(i) sus(i) sus_u(i)]);
	sd_sus(i) = std([sus_l(i) sus(i) sus_u(i)]);
end

figure(3)
loglog(meanN, av_sus, 'o')



% susceptibility X
%Ng
Ng1 = [];
Ng2 = [];
Ng3 = [];
for i = 1:length(Ng)
	Ng1(i) = Ng{i}(1);
	Ng2(i) = Ng{i}(2);
	Ng3(i) = Ng{i}(3);
end

%N
N1 = []; % upper bound
N2 = [];
N3 = []; % lower bound
for i = 1:length(N)
	N1(i) = N{i}(1);
	N2(i) = N{i}(2);
	N3(i) = N{i}(3);
end

Ng1sqN1 = (Ng1)./N1; % upper bound
Ng2sqN2 = (Ng2)./N2;
Ng3sqN3 = (Ng3)./N3; % lower bound

rN1 = N1(id);
rN2 = N2(id);
rN3 = N3(id);

rNg1sqN1 = Ng1sqN1(id);
rNg2sqN2 = Ng2sqN2(id);
rNg3sqN3 = Ng3sqN3(id);

delta = 10;

X = {};
avX = [];
sdX = [];

Nx = {};
avNx = [];
sdNx = [];
for i = 1:(length(rN1) - delta)
	X{i}(1) = mean( rN1(i: (i + delta)) ) * ( mean(rNg1sqN1(i: (i + delta)).^2) - mean(rNg1sqN1(i: (i + delta)))^2 );
	X{i}(2) = mean( rN2(i: (i + delta)) ) * ( mean(rNg2sqN2(i: (i + delta)).^2) - mean(rNg2sqN2(i: (i + delta)))^2 );
	X{i}(3) = mean( rN3(i: (i + delta)) ) * ( mean(rNg3sqN3(i: (i + delta)).^2) - mean(rNg3sqN3(i: (i + delta)))^2 );
	avX(i) = mean(X{i});
	sdX(i) = std(X{i});

	Nx{i}(1) = mean( rN1(i: (i + delta)) );
	Nx{i}(2) = mean( rN2(i: (i + delta)) );
	Nx{i}(3) = mean( rN3(i: (i + delta)) );
	avNx(i) = mean(Nx{i});
	sdNx(i) = std(Nx{i});
end

figure(4)
loglog(avNx, avX, 'o')

