% Author: Laurent Condat 
% homepage: http://www-m6.ma.tum.de/~condat/ 
% Last Modified  3/14/08
% Tested with Matlab 7.5.0
% Please, report to condat@ma.tum.de any bug or comment. 


iptsetpref('ImshowBorder','tight');


coltab1=[1 3 2;3 2 1;2 1 3]; 			% coltab1(i,j) gives a color different from i and j if i different from j
coltab2=[2 3;1 3;1 2];					% coltab2(i,:) gives the colors different from i

Nx=100; 								% number of horizontal pixels 
Ny=100;  								% number of vertical pixels 

Nbiter=100;								% number of iterations for the second step of the algorithm


[w1,w2]=meshgrid(-1/2:1/Nx:1/2-1/Nx,-1/2:1/Ny:1/2-1/Ny);  % domain for the plots of the spectra
w1=w1*2*pi;
w2=w2*2*pi;

Mosa=zeros(Ny,Nx);						% mosaic image

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Step 1        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


p=randperm(Ny*Nx)-1;					% permutation for fast dart throwing
for q=1:Ny*Nx							% main loop
	y=floor(p(q)/Nx)+1;					% coordinate of the pixel to be added
	x=mod(p(q),Nx)+1;					% coordinate of the pixel to be added
	tab=[1 1 1 1];						% color values possible for the pixel: tab(i+1)==1 if the color i is allowed
	if y==1, aux=Mosa(2,x); else  aux=Mosa(y-1,x);  end; 	  	% color of the first neighbor
	tab(aux+1)=0;
	if x==1, aux=Mosa(y,2); else  aux=Mosa(y,x-1);  end; 	  	% color of the second neighbor
	tab(aux+1)=0;
	if y==Ny, aux=Mosa(Ny-1,x); else  aux=Mosa(y+1,x);  end;   	% color of the third neighbor
	tab(aux+1)=0;
	if x==Nx, aux=Mosa(y,Nx-1); else  aux=Mosa(y,x+1);  end;   	% color of the fourth neighbor
	tab(aux+1)=0;
	tab(1)=0;
	if sum(tab)==0						% case where no color can be assigned under the minimum-distance condition
		Mosa(y,x)=floor(rand()*3.0)+1;	% the color is chosen randomly
	elseif sum(tab)==3					% case where no neighbor has been assigned a color yet
		Mosa(y,x)=floor(rand()*3.0)+1;	% the color is chosen randomly
	elseif sum(tab)==2					% we have the choice between two colors
		tab(1)=1;
		id=find(tab==0);
		Mosa(y,x)=coltab2(id-1,(rand()<0.5)+1);	% the color is chosen randomly
	else								% only one color is possible
		id=find(tab==1);	
		Mosa(y,x)=id-1;
	end;
end;
Mosa1=Mosa;								% we save the mosaic for eventual later use

figure									% we plot the mosaic image
imshow(Mosa,[1,3]);
colormap(eye(3));

figure									% we plot the spectrum of the blue-yellow chrominance channel
MosaF=Mosa;
MosaF(find(MosaF<3))=0;				
MosaF=(MosaF-1)/3.0;				
										% the amplitude spectrum is lowpass filtered
surf(w1,w2,filter2(ones(5,5),abs(fftshift(fft2(MosaF))))/sqrt(Ny*Nx));
shading flat
axis square
axis([-pi pi -pi pi 0 50 0 25]);
set(gca,'YTick',[-3 0 3])
set(gca,'XTick',[-3 0 3])
colormap (1-gray)
figure
caxis([0 sqrt(Ny*Nx)*25])
imshow(filter2(ones(5,5),abs(fftshift(fft2(MosaF)))),caxis)
colormap (1-gray)


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Step 2        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iter=1:Nbiter						
	total=0;							% number of modified pixels during one iteration
	p=randperm(Ny*Nx)-1;				% a new permutation is computed at each iteration
	for q=1:Ny*Nx						% loop over the pixels
		y=floor(p(q)/Nx)+1;
		x=mod(p(q),Nx)+1;
		tab=[0 0 0];					% tab(i) = number of neighbors having the color i
		if y==1, aux=Mosa(2,x); 	else  aux=Mosa(y-1,x);  end;
		tab(aux)=tab(aux)+1;
		if x==1, aux=Mosa(y,2); 	else  aux=Mosa(y,x-1);  end; 	  	
		tab(aux)=tab(aux)+1;
		if y==Ny, aux=Mosa(Ny-1,x); else  aux=Mosa(y+1,x);  end;  
		tab(aux)=tab(aux)+1;
		if x==Nx, aux=Mosa(y,Nx-1); else  aux=Mosa(y,x+1);  end;   
		tab(aux)=tab(aux)+1;
		if tab(Mosa(y,x))>0 			% the pixel has the same color as one neighbor
			if nnz(tab)==3				% no value can be asigned: we assign one of the less present
				Mosa(y,x)=coltab2(find(tab==2),(rand()<0.5)+1);
				% the following two alternatives give worse results
				% Mosa(y,x)=coltab2(Mosa(y,x),(rand()<0.5)+1);
				% Mosa(y,x)=floor(rand()*3.0)+1;
			elseif nnz(tab)==1			% we have the choice between two values
				Mosa(y,x)=coltab2(find(tab>0),(rand()<0.5)+1);
			else  
				Mosa(y,x)=find(tab==0);
			end;
			total=total+1;
		end;
	end;
	total
	iter
end;
Mosa2=Mosa;								% we save the mosaic for eventual later use

figure
imshow(Mosa,[1,3]);
colormap(eye(3));

figure							
MosaF=Mosa;
MosaF(find(MosaF<3))=0;				
MosaF=(MosaF-1)/3.0;				
								
surf(w1,w2,filter2(ones(5,5),abs(fftshift(fft2(MosaF))))/sqrt(Ny*Nx));
shading flat
axis square
axis([-pi pi -pi pi 0 50 0 25]);
set(gca,'YTick',[-3 0 3])
set(gca,'XTick',[-3 0 3])
colormap (1-gray)
figure
caxis([0 sqrt(Ny*Nx)*25])
imshow(filter2(ones(5,5),abs(fftshift(fft2(MosaF)))),caxis)
colormap (1-gray)


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         Step 3        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%


tab=zeros(3,3,3,3);		% used to check if the four neighbors are equal
tab(1,1,1,1)=1;
tab(2,2,2,2)=2;
tab(3,3,3,3)=3;
tab1=[10 8 6];			
tab2=[2  1 1];			% tab2(i) is the minimal color different from i
tab3=[3  3 2];			% tab3(i) is the maximal color different from i
for y=1:Ny
	for x=1:Nx
		if y==1, y1=2; else y1=y-1;  end;
		if x==1, x1=2; else x1=x-1; end;
		if y==Ny, y2=Ny-1; else y2=y+1; end;
		if x==Nx, x2=Nx-1; else x2=x+1; end;
		aux=tab(Mosa(y,x1),Mosa(y,x2),Mosa(y1,x),Mosa(y2,x));
		if aux>0		% the four neighbors are equal
			aux2=Mosa(y1,x1)+Mosa(y1,x2)+Mosa(y2,x1)+Mosa(y2,x2);
			if aux2>tab1(aux)		% we don't have two diagonal pixels for each of the two colors
				Mosa(y,x)=tab2(aux);
			elseif aux2<tab1(aux)	% we don't have two diagonal pixels for each of the two colors
				Mosa(y,x)=tab3(aux);
			end;
		end;
	end;
end;
Mosa3=Mosa;						% we save the mosaic for eventual later use

figure
imshow(Mosa,[1,3]);
colormap(eye(3));

figure							
MosaF=Mosa;
MosaF(find(MosaF<3))=0;				
MosaF=(MosaF-1)/3.0;				
								
surf(w1,w2,filter2(ones(5,5),abs(fftshift(fft2(MosaF))))/sqrt(Ny*Nx));
shading flat
axis square
axis([-pi pi -pi pi 0 50 0 25]);
set(gca,'YTick',[-3 0 3])
set(gca,'XTick',[-3 0 3])
colormap (1-gray)
figure
caxis([0 sqrt(Ny*Nx)*25])
imshow(filter2(ones(5,5),abs(fftshift(fft2(MosaF)))),caxis)
colormap (1-gray)



%imwrite(Mosa1,eye(3),'CFArandom1.tif','TIFF')
%imwrite(Mosa2,eye(3),'CFArandom2.tif','TIFF')
%imwrite(Mosa3,eye(3),'CFArandom3.tif','TIFF')


% * This program is free software; you can redistribute it and/or
% * modify it under the terms of the GNU General Public License
% * as published by the Free Software Foundation; either version 2
% * of the License, or (at your option) any later version. 
% *
% * This program is distributed in the hope that it will be useful,
% * but WITHOUT ANY WARRANTY; without even the implied warranty of
% * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% * GNU General Public License for more details.
% *
% * You can receive a copy of the GNU General Public License
% * by writing to the Free Software Foundation,
% * Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
