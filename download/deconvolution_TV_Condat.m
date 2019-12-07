% This Matlab file implements deconvolution regularized by total
% variation, as a simple example of application of the optimization 
% algorithm described in the article in French:
%
% L. Condat, "Un nouvel algorithme proximal pour l'optimisation convexe 
% non lisse", Proc. of GRETSI, Brest, France, Sept. 2013.
%
% The optimization algorithm is described in more details in 
% L. Condat, "A primal-dual splitting method for convex optimization 
% involving Lipschitzian, proximable and linear composite terms", 
% J. Optimization Theory and Applications, vol. 158, no. 2, pp. 460-479, 
% 2013
% and in
% L. Condat, "A Generic Proximal Algorithm for Convex Optimization - 
% Application to Total Variation Minimization," IEEE Signal Proc. 
% Letters, vol. 21, no. 8, pp. 1054-1057, Aug. 2014.
%
% This code has been written by Laurent Condat, CNRS research fellow in 
% the Dept. of Images and Signals of GIPSA-lab, a research center of the
% University of Grenoble-Alpes.
%
% For any comment or question, contact me at 
% laurent.condat@gipsa-lab.grenoble-inp.fr
% 
% Version 1.0, June 17, 2013.
%
% Tested on a Apple laptop with Mac OS 10.8 and Matlab R2011b.
%
% Replace 'monarch.tif' by any image of your choice below.


function main()
	blurlevel=5;
	noiselevel=3;
	lambda=0.02;
	nbiter=300;
	I=imread('monarch.tif');
	I=double(I);	
	filter = fspecial('gaussian',blurlevel*6+1,blurlevel);
	I=imfilter(I,filter,'symmetric')+noiselevel*randn(size(I));
	imwrite(I/255,'degraded.tif');
	tic
	J=deconvolution(I,filter,lambda,nbiter);
	toc
	imwrite(J/255,'restored.tif');
end


function Iout = deconvolution(Iin,filter,lambda,nbiter)
	sigma=lambda;
	tau=0.99/(0.5+8*sigma);
	[sizex,sizey]=size(Iin);
	Iout=Iin;
	Idual1=zeros(sizex,sizey);
	Idual2=Idual1;
	thewaitbar = waitbar(0,'Nb iterations'); 
	figure
	imshow(Iin);
	colormap gray
	axis image
	for iter=1:nbiter  
	    Iaux=Iout;
	    Iout=Iout-tau*(imfilter((imfilter(Iout,filter,'symmetric')-Iin),filter,'symmetric')+...
	    	[-Idual1(:,1),Idual1(:,1:end-1)-Idual1(:,2:end)]+...
	    	[-Idual2(1,:);Idual2(1:end-1,:)-Idual2(2:end,:)]);
	    Iout=min(255,max(Iout,0));
	    imshow(Iout/255);
	    Iaux=2*Iout-Iaux;
	    Idual1=Idual1+sigma*[Iaux(:,2:end)-Iaux(:,1:end-1), zeros(sizex,1)];
	    Idual2=Idual2+sigma*[Iaux(2:end,:)-Iaux(1:end-1,:); zeros(1,sizey)];
	    Iaux=max(1,sqrt(Idual1.^2+Idual2.^2)/lambda);
	    Idual1=Idual1./Iaux;
	    Idual2=Idual2./Iaux;
	    waitbar(iter/nbiter);
	end
	close(thewaitbar)
end

