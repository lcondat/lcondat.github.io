% Author: Laurent Condat 
% homepage: http://www-m6.ma.tum.de/~condat/ 
% Last Modified  3/14/08
% Tested with Matlab 7.5.0
% Please, report to condat@ma.tum.de any bug or comment.


iptsetpref('ImshowBorder','tight');


coltab1=[1 3 2;3 2 1;2 1 3]; 		% coltab1(i,j) gives a color different from i and j if i different from j
coltab2=[2 3;1 3;1 2];				% coltab2(i,:) gives the colors different from i

Nx=100; 							% number of horizontal pixels 
Ny=100;  							% number of vertical pixels 

tiles=[1 2 3;1 3 2; 2 3 1; 2 1 3; 3 1 2;  3 2 1];	% the 6 tiles that tile the first row and column
allowed=[2 4;1 5;4 6;1 3;2 6;3 5];	% allowed(i) : tiles allowed as neighbors of tile i
tab=[1 3 5];						% tab(i) = a tile that begins with i


Mosa=zeros(Ny,Nx);
tile=1;
for x=1:Nx							% we generate the first row
	if mod(x,3)==1					% we choose a tile randomly
		tile=allowed(tile,floor(rand()*2.0)+1);	
	end;
	Mosa(1,x)=tiles(tile,mod(x-1,3)+1);
end
tile=tab(Mosa(1,1));				% the first vertical tile should be consistent with the pixel Mosa(1,1)
for y=2:3
	Mosa(y,1)=tiles(tile,mod(y-1,3)+1);
end
for y=4:Ny							% we generate the first column
	if mod(y,3)==1	
		tile=allowed(tile,floor(rand()*2.0)+1);
	end;
	Mosa(y,1)=tiles(tile,mod(y-1,3)+1);
end
for y=2:Ny							% we generate the rest of the mosaic in scanline order
	for x=2:Nx
		if Mosa(y-1,x)==Mosa(y,x-1)
			Mosa(y,x)=coltab1(Mosa(y,x-1),Mosa(y-1,x-1));
		else
			Mosa(y,x)=coltab1(Mosa(y,x-1),Mosa(y-1,x));
		end
	end
end

figure
imshow(Mosa,[1,3]);
colormap(eye(3));

figure							% we plot the spectrum of the blue-yellow chrominance channel
MosaF=Mosa;
MosaF(find(MosaF<3))=0;				
MosaF=(MosaF-1)/3;				
[w1,w2]=meshgrid(-1/2:1/Nx:1/2-1/Nx,-1/2:1/Ny:1/2-1/Ny);  % domain for the plots of the spectra
w1=w1*2*pi;
w2=w2*2*pi;
								% the amplitude spectrum is lowpass filtered
surf(w1,w2,filter2(ones(5,5),abs(fftshift(fft2(MosaF))))/sqrt(Ny*Nx));
shading flat
axis square
axis([-pi pi -pi pi 0 50 0 25]);
set(gca,'YTick',[-3 0 3])
set(gca,'XTick',[-3 0 3])
colormap (1-gray)
figure
MosaF2=1-filter2(ones(5,5),abs(fftshift(fft2(MosaF))))/sqrt(Ny*Nx)/25;
imshow(MosaF2,[0 1])


%imwrite(Mosa,eye(3),'CFAtile.tif','TIFF')


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
