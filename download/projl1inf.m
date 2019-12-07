% Projection of the matrix Y onto the l1,infinity norm ball 
% {X : sum_n max_m |X(n,m)| <= c}
% for some given c>0
%
% Author: Laurent Condat
%
% Version: 1.0, Sept. 1, 2017
%
% This algorithm is new, to the author's knowledge. It is based
% on the same ideas as for projection onto the l1 ball, see
% L. Condat, "Fast projection onto the simplex and the l1 ball",
% Mathematical Programming, vol. 158, no. 1, pp. 575-585, 2016. 
%
% The algorithm is exact and terminates in finite time. Its
% average complexity, for Y of size N x M, is O(NM.log(M)). 
% Its worst case complexity, never found in practice, is
% O(NM.log(M) + N^2.M).


function X = projl1inf(Y, c)
	X = sort(abs(Y),2,'descend');
	v = sum(X(:,1));
	if v<=c, X=Y; return; end
	[N,M] = size(Y);
	S = cumsum(X,2);
	idx = ones(N,1);
	theta = (v-c)/N;
	mu = zeros(N,1);
	active = ones(N,1);
	thetaold = 0;
	while thetaold~=theta
		for n=1:N
			if active(n)
				j = idx(n);
				while (j<M)&&((S(n,j)-theta)/j<X(n,j+1)), j=j+1; end
				idx(n) = j;
				mu(n) = S(n,j)/j;
				if (j==M)&&(mu(n)-theta/j<=0), active(n)=0; mu(n)=0; end
			end
		end
		thetaold = theta;
		theta = (sum(mu)-c)/sum(active./idx);
	end
	X = bsxfun(@min,abs(Y),(mu-theta./idx).*active).*sign(Y);
end



