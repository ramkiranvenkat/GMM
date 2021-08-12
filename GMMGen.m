clc
clear all
close all

data = load('data.file');
dataInit = load('data.init');
dataLabels = load('data.labels');

mean = dataInit;
ds = length(data);
cc = 2;
D = size(data,2);
alpha = ones(1,cc)/cc;
for i = 1:cc
	sigL = zeros(D,D);
	count = 0;
	for j=1:ds
		if ((dataLabels(j)+1) == i)
			vec = data(j,:)-mean(i,:);
			sigL = sigL + vec'*vec;
			count = count + 1;
		end
	end
	sig(:,:,i) = sigL/count;
end

for itr=1:10
	%expectation
	llh = 0;
	for i=1:ds
		for j=1:cc
			gamma(i,j) = alpha(j)*gaussianEval(data(i,:),mean(j,:),sig(:,:,j));
		end
		gs = sum(gamma(i,:));
		gamma(i,:) = gamma(i,:)/gs;
		llh = llh + log(gs);
	end
	llhs(itr) = llh;
	Nc = sum(gamma);

	%maximization

	for i=1:cc
		for j=1:D
			mean(i,j) = sum(gamma(:,i).*data(:,j))/Nc(i);
		end
		sig(:,:,i) = zeros(D,D);
	end
	
	for i=1:ds
		for j=1:cc
			vec = data(i,:)-mean(j,:);
			sig(:,:,j) = sig(:,:,j) + gamma(i,j)*(vec'*vec)/Nc(j);
		end
	end

	alpha = Nc/sum(Nc);
end

save output.data alpha mean sig llhs

