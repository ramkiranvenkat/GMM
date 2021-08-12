clc
clear all
close all

for i=1:1000
	if (mod(i,3) == 0)
		data(i,:) = [2 3] + [randn(1) randn(1)]*[1 0.2;0.2 0.8];
	else
		data(i,:) = [15 16] + [randn(1) randn(1)]*[2.0 0.3;0.3 1.5];
	end
end

mean = [-1 2;16 19];
sig(:,:,1) = 0.1*eye(2);
sig(:,:,2) = 0.1*eye(2);
alpha = [0.5 0.5];
ds = length(data);
cc = 2;
D = 2;
for i = 1:100
	points(i,:) = 3*[cos(2*pi/100*i) sin(2*pi/100*i)];
end

figure,plot(data(:,1),data(:,2),'b*'),hold on,plot(mean(:,1),mean(:,2),'r+')

for itr=1:10
	%expectation
	llh = 0;
	for i = 1:100
		b(i,:,1) = mean(1,:) + points(i,:)*sig(:,:,1)^0.5;
		b(i,:,2) = mean(2,:) + points(i,:)*sig(:,:,2)^0.5;
	end
	plot(b(:,1,1),b(:,2,1),'r-'),plot(b(:,1,2),b(:,2,2),'r-'),grid minor

	for i=1:ds
		%gamma(i,:) = alpha.*[gaussianEval(data(i,:),mean(1,:),sig(:,:,1)) gaussianEval(data(i,:),mean(2,:),sig(:,:,2))];
		for j=1:cc
			gamma(i,j) = alpha(j)*gaussianEval(data(i,:),mean(j,:),sig(:,:,j));
		end
		gs = sum(gamma(i,:));
		gamma(i,:) = gamma(i,:)/gs;
		llh = llh + log(gs);
	end
	llh

	Nc = sum(gamma);
	%maximization

	%mean = [(sum([gamma(:,1).*data(:,1)   gamma(:,1).*data(:,2)]))/Nc(1);
	%	(sum([gamma(:,2).*data(:,1)   gamma(:,2).*data(:,2)]))/Nc(2)];

	for i=1:cc
		for j=1:D
			mean(i,j) = sum(gamma(:,i).*data(:,j))/Nc(i);
		end
	end

	sig(:,:,1) = [0 0;0 0];
	sig(:,:,2) = [0 0;0 0];
	for i=1:ds
		vec = data(i,:)-mean(1,:);
		sig(:,:,1) = sig(:,:,1) + gamma(i,1)*(vec'*vec);
		vec = data(i,:)-mean(2,:);
		sig(:,:,2) = sig(:,:,2) + gamma(i,2)*(vec'*vec);
	end
	sig(:,:,1) = diag(diag(sig(:,:,1)/Nc(1)));
	sig(:,:,2) = diag(diag(sig(:,:,2)/Nc(2)));

	alpha = Nc/sum(Nc);
end

