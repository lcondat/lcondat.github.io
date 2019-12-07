% This is MATLAB code for computing the proximity operator 
% or projecting on the ball of several norms, which appear
% often in data science.
%
% Copyright Laurent Condat
% v1.0, May 11, 2018
% 
% If you find any bug, please contact me: 
% http://www.gipsa-lab.grenoble-inp.fr/~laurent.condat/



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% proximity operator of the l_1 norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% l1norm = @(x) lambda * sum(abs(x(:))),
% where lambda>=0,
% applied at the N-D array y.
% This is soft-thresholding elementwise.

function x = prox_l1 (y, lambda)
	x = y - max(min(y, lambda), -lambda);
end
% x = sign(y).*max(abs(y)-lambda,0); equivalent (slower)
% x = y - projball_linf(y,lambda); equivalent


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection on the l_infinity norm ball
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% { x : linfnorm(x)<=tau },
% where tau>=0 and
% linfnorm = @(x) max(abs(x(:))),
% of the N-D array y.
% This is just saturation elementwise.

function x = projball_linf (y, tau)
	x = max(min(y, tau), -tau);
end
% x = sign(y).*min(abs(y),tau); equivalent (slower)
% x = y./max(abs(y)/tau,1); works even with tau=0,
%	because max(NaN,1)=1 and y./Inf=0
% x = y.*(tau./max(abs(y),tau)); works if tau>0
% avoiding divisions whenever possible is better.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% proximity operator of the l_2 (Frobenius) norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% l2norm = @(x) lambda * sqrt(sum(x(:).^2)),
% where lambda>=0,
% applied at the N-D array y.

function x = prox_l2 (y, lambda)
	tmp = sqrt(sum(y(:).^2));
	if tmp<=lambda, x = zeros(size(y));
	else 
		x = y * ((tmp-lambda)/tmp);
	end
end
% x = y - projball_l2(y,lambda); equivalent
% an alternative, only if lambda>0:
%	tmp = max(sqrt(sum(y(:).^2)),lambda);
%	x = y * ((tmp-lambda)/tmp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection on the l_2 (Frobenius) norm ball
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% { x : l2norm(x)<=tau },
% where tau>=0 and
% l2norm = @(x) sqrt(sum(x(:).^2)),
% of the N-D array y.

function x = projball_l2 (y, tau)
	tmp = sqrt(sum(y(:).^2));
	if tmp<=tau, x = y;
	else 
		x = y * (tau/tmp);
	end
end
% x = y - prox_l2(y,tau); equivalent
% x = y * (tau/max(sqrt(sum(y(:).^2)),tau)); works if tau>0
% x = y / max(sqrt(sum(y(:).^2))/tau,1); works even with tau=0
%    because max(NaN,1)=1 and y/Inf=0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection on the l_1 norm ball
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% { x : l1norm(x)<=tau },
% where tau>=0 and
% l1norm = @(x) sum(abs(x(:))),
% of the N-D array y.

function x = projball_l1(y, tau)
	tmp = abs(y(:));
	if sum(tmp)<=tau, x = y; return; end
	lambda = max((cumsum(sort(tmp,1,'descend'))-tau)./(1:length(tmp))');
	x = y - max(min(y, lambda), -lambda);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% proximity operator of the l_infinity norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% linfnorm = @(x) lambda * max(abs(x(:))),
% where lambda>=0,
% applied at the N-D array y.

function x = prox_linf (y, lambda)
	tmp = abs(y(:));
	if sum(tmp)<=lambda, x = zeros(size(y)); return; end
	tau = max((cumsum(sort(tmp,1,'descend'))-lambda)./(1:length(tmp))');
	x = max(min(y, tau), -tau);
end
% x = y - projball_l1(y,lambda); equivalent


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% proximity operator of the l_1,2 norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% l12norm = @(x) lambda * sum(sqrt(sum(x.^2,2))),
% where lambda>0,
% applied at the matrix y.
% This is just prox_l2 on each row.

function x = prox_l12 (y, lambda)
	tmp = max(sqrt(sum(y.^2,2)), lambda);
	x = bsxfun(@times, y, (tmp-lambda)./tmp);
end
% x = bsxfun(@times, y, 1-lambda./max(sqrt(sum(y.^2,2)),lambda));
%	is equivalent but slightly less precise
% an alternative, which works even when lambda=0:
%	tmp = sqrt(sum(y.^2,2));
%	tmp2 = (tmp-lambda)./tmp;
%	tmp2(tmp<=lambda) = 0;
%	x = bsxfun(@times, y, tmp2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection on the l_infinity,2 norm ball
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% { x : linf2norm(x)<=tau },
% where tau>0 and
% linf2norm = @(x) max(sqrt(sum(x.^2,2))),
% of the matrix y.
% This is just projball_l2 on each row.

function x = projball_linf2 (y, tau)
	x = bsxfun(@times, y, tau./max(sqrt(sum(y.^2,2)), tau));
end
% x = y - prox_l12(y,lambda); equivalent
% x = bsxfun(@rdivide, y, max(sqrt(sum(y.^2,2))/tau,1)); 
%	works even with tau=0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection on the l_1,2 norm ball
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% { x : l12norm(x)<=tau },
% where tau>=0 and
% l12norm = @(x) sum(sqrt(sum(x.^2,2))),
% of the matrix y.

function x = projball_l12 (y, tau)
	tmp = sqrt(sum(y.^2,2));
	if sum(tmp)<=tau, x = y; return; end
	lambda = max((cumsum(sort(tmp,1,'descend'))-tau)./(1:length(tmp))');
	tmp = max(tmp, lambda);
	x = bsxfun(@times, y, (tmp-lambda)./tmp);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% proximity operator of the l_infinity,2 norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% linf2norm = @(x) lambda * max(sqrt(sum(x.^2,2))),
% where lambda>=0,
% applied at matrix y.

function x = prox_linf2 (y, lambda)
	tmp = sqrt(sum(y.^2,2));
	if sum(tmp)<=lambda, x = zeros(size(y)); return; end
	tau = max((cumsum(sort(tmp,1,'descend'))-lambda)./(1:length(tmp))');
	x = bsxfun(@times, y, tau./max(tmp, tau));
end
% x = y - projball_l12(y,lambda); equivalent


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection on the l_1,infinity norm ball
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% { x : l1infnorm(x)<=tau },
% where tau>0 and
% l1infnorm = @(x) sum(max(abs(x),[],2)),
% of the matrix y.
%
% This new algorithm is an extension of the sort-based 
% projection onto the l1 ball, see
% L. Condat, "Fast projection onto the simplex and the l1 ball",
% Mathematical Programming, vol. 158, no. 1, pp. 575-585, 2016. 
% The algorithm is noniterative and exact. Its average
% complexity, for y of size N x M, is O(NM.log(M)). 

function x = projball_l1inf(y, tau)
	x = abs(y);
	v = sum(max(x,[],2));
	if v<=tau, x = y; return; end
	x = sort(x,2,'descend');
	[N,M] = size(y);
	S = cumsum(x,2);
	idx = ones(N,1);
	theta = (v-tau)/N;
	mu = zeros(N,1);
	active = ones(N,1);
	thetaold = 0;
	while thetaold~=theta
		for n=1:N
			if active(n)
				j = idx(n);
				while (j<M)&&((S(n,j)-theta)/j<x(n,j+1)), j=j+1; end
				idx(n) = j;
				mu(n) = S(n,j)/j;
				if (j==M)&&(mu(n)-theta/j<=0), active(n)=0; mu(n)=0; end
			end
		end
		thetaold = theta;
		theta = (sum(mu)-tau)/sum(active./idx);
	end
	x = bsxfun(@min,abs(y),(mu-theta./idx).*active).*sign(y);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% proximity operator of the l_infinity,1 norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% linf1norm = @(x) lambda * max(sum(abs(x),2)),
% where lambda>0,
% applied at the matrix y.

function x = prox_linf1(y, lambda)
	x = abs(y);
	v = sum(max(x,[],2));
	if v<=lambda, x = zeros(size(y)); return; end
	x = sort(x,2,'descend');
	[N,M] = size(y);
	S = cumsum(x,2);
	idx = ones(N,1);
	theta = (v-lambda)/N;
	mu = zeros(N,1);
	active = ones(N,1);
	thetaold = 0;
	while thetaold~=theta
		for n=1:N
			if active(n)
				j = idx(n);
				while (j<M)&&((S(n,j)-theta)/j<x(n,j+1)), j=j+1; end
				idx(n) = j;
				mu(n) = S(n,j)/j;
				if (j==M)&&(mu(n)-theta/j<=0), active(n)=0; mu(n)=0; end
			end
		end
		thetaold = theta;
		theta = (sum(mu)-lambda)/sum(active./idx);
	end
	x = y - bsxfun(@min,abs(y),(mu-theta./idx).*active).*sign(y);
end
% x = y - projball_l1inf(y,lambda); equivalent


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% projection on the l_infinity,1 norm ball
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% { x : linf1norm(x)<=tau },
% where tau>=0 and
% linf1norm = @(x) max(sum(abs(x),2)),
% of the matrix y.
% This is just projball_l1 on each row.

function x = projball_linf1(y, tau)	
	x = y;
	mask = sum(abs(y),2)>tau;
	y = y(mask,:);
	lambda = max(bsxfun(@rdivide, cumsum(sort(abs(y),2,'descend'),2)-tau, ...
		(1:size(y,2))),[],2);
	x(mask,:) = y - bsxfun(@max, bsxfun(@min, y, lambda), -lambda);
end
% an alternative:
%	lambda = max(max(bsxfun(@rdivide, cumsum(sort(abs(y),2,'descend'),2)-tau, ...
%		(1:size(y,2))),[],2),0);
%	x = y - bsxfun(@max, bsxfun(@min, y, lambda), -lambda);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% proximity operator of the l_1,infinity norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% l1infnorm = @(x) lambda * sum(max(abs(x),[],2)),
% where lambda>=0,
% applied at the matrix y.
% This is just prox_linf on each row.

function x = prox_l1inf (y, lambda)
	x = zeros(size(y));
	mask = sum(abs(y),2)>lambda;
	y = y(mask,:);
	tau = max(bsxfun(@rdivide, cumsum(sort(abs(y),2,'descend'),2)-lambda, ...
		(1:size(y,2))),[],2);
	x(mask,:) = bsxfun(@max, bsxfun(@min, y, tau), -tau);
end
% x = y - projball_linf1(y,lambda); equivalent


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% proximity operator of the nuclear norm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nucnorm = @(x) lambda*sum(svd(x)),
% where lambda>0,
% applied at the matrix y of size N x 2.
% This is soft-thresholding of the singular values.

function x = prox_nucnorm(y, lambda)
	s = sum(y.^2,1);
	theta = atan2(2*dot(y(:,1),y(:,2)),s(1)-s(2))/2;
	c = cos(theta);
	s = sin(theta);
	v = [c -s; s c];
	x = y * v;
	tmp = max(sqrt(sum(x.^2,1)), lambda);
	x = bsxfun(@times, x, (tmp-lambda)./tmp) * v';
end
% see http://scipp.ucsc.edu/~haber/ph116A/diag2x2_11.pdf
