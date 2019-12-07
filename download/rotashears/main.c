/*
Program for image rotation using shears, as described in the article:
L. Condat and Dimitri Van De Ville, “Fully reversible image rotation by 1-D filtering,” 
	Proc. of IEEE ICIP, Oct. 2008, San Diego, US.
(PDF available at the homepage of L. Condat).

Preliminary version. Written by Laurent Condat, 24 Feb. 2009.
Only the type II (centered) order 2 all-pass filter is implemented for the shifts.
Only works for grayscale images in TIFF format. If a color image is given, its red channel will be processed.
The rotation angle should be in the range [-Pi/4,Pi/4]. Not implemented else (just combine with trivial rotations of -90, 90, 180 degrees to obtain every angle). The program will deliver a result for an angle outside this range, but the shears will cut a big relevant part of the image.

Compile with gcc main.c -o main -ltiff -lm 
execute by ./main <your image in TIFF format>
Please report any bug or comment to laurent.condat@greyc.ensicaen.fr

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version. 

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You can receive a copy of the GNU General Public License
by writing to the Free Software Foundation,
Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>  /*for memcpy*/
#include <tiffio.h>


#define Pi  3.1415926535897932384626

const float ANGLE = Pi/5.0;		/*rotation angle in rad. */
const int EXTEND = 1; /*extend the image before the rotation for avoiding the cutting effects of the shears at the boundaries*/

int SWAP_BOOL = 0; 	/*0 for Windows and Mac Intel, 1 for Mac PowerPC. Controls the little/big endianness*/


/*Our image structure*/
typedef struct {
	int width, height;	
	float* dataF;		/*Buffer of floats*/
	uint32* dataI;		/*Buffer of ints (raw data of the TIFF file)*/
} Image;

#define 	EMPTY  {0,0,NULL,NULL} 


/*Check function for the mallocs*/
void* Check(void* pointer) {
	if (pointer == NULL) {
		printf("Memory allocation error\n");
		exit(1);
	}
	return pointer;
}


static void SwapBuffer (Image* Source, char reverse) {
	int    res;
	uint32* 	cursS;
	uint32* 	cursD;
	uint32		pixel;
	int		x, y;
	int		width = Source->width;
	int		height = Source->height;
	uint32*		Dest;
	
  	Dest = (uint32*)Check(_TIFFmalloc(width * height * sizeof(uint32)));
	cursS = Source->dataI;
	cursD = Dest + (height - 1) * width;
	for (y = 0; y < height; y++) {
		for (x = 0; x < width; x++) {
			pixel = *cursS;
			if (reverse)
				pixel =	((pixel & 0xFF000000) >> 24) |
					((pixel & 0x00FF0000) >>  8) |
					((pixel & 0x0000FF00) <<  8) |
					((pixel & 0x000000FF) << 24);
			*cursD = pixel;
			cursS++;
			cursD++;
		}
		cursD -= 2 * width;
	}
	_TIFFfree ((tdata_t)(Source->dataI));
	Source->dataI = Dest;
}


static void LoadPixels (TIFF* file, Image* Dest) {
  	int     res;
	
	res = TIFFGetField (file, TIFFTAG_IMAGEWIDTH, &(Dest->width));
  	if (res != 1) {
    		printf ("Tag IMAGEWIDTH was not found in tiff file\n");
    		exit(1);
  	}
	res = TIFFGetField (file, TIFFTAG_IMAGELENGTH, &(Dest->height));
  	if (res != 1) {
    		printf ("Tag IMAGEHEIGHT was not found in tiff file\n");
    		exit(1);
  	}
  	/* Allocate memory for the pixels and load them from file */
	Dest->dataI = (uint32*)Check(malloc(Dest->width * Dest->height * sizeof(uint32)));
  	res = TIFFReadRGBAImage (file, Dest->width, Dest->height, Dest->dataI, 0);
  	if (res != 1) {
    		printf ("Could not load image pixels into memory\n");
    		free ((tdata_t)(Dest->dataI));
    		exit(1);
  	}  
}


void LoadImage(Image* Dest, char* name) {
	TIFF*     	File;
  	File = TIFFOpen (name, "r");	
  	if (File == (TIFF*)NULL) {
    		printf ("Could not open file named %s for reading.\n", name);
    		exit (1);
  	}
	LoadPixels (File, Dest);	
	TIFFClose (File);			
	SwapBuffer (Dest, SWAP_BOOL); 
}


static void SavePixels (TIFF* file, Image* Source) {
  	int     res;
	int	width = Source->width;
	int	height = Source->height;
  	
  	res =   TIFFSetField (file, TIFFTAG_IMAGEWIDTH, width);
  	if (res == 1)
    		res = TIFFSetField (file, TIFFTAG_IMAGELENGTH, height);
  	if (res == 1)
    		res = TIFFSetField (file, TIFFTAG_BITSPERSAMPLE, (uint16)8);
  	if (res == 1)
    		res = TIFFSetField (file, TIFFTAG_PHOTOMETRIC, (uint16)2);
  	if (res == 1)
    		res = TIFFSetField (file, TIFFTAG_SAMPLESPERPIXEL, (uint16)4);
  	if (res == 1)
    		res = TIFFSetField (file, TIFFTAG_PLANARCONFIG, (uint16)PLANARCONFIG_CONTIG);
	if (res == 1) {
    		res = (int)TIFFWriteRawStrip (file, 0, Source->dataI, sizeof(uint32) * width * height);
		if (res == -1) res = 0;
	}	
  	if (res == 0) {
    		printf ("An error occured while trying to save TIFF file\n");
    		exit(1);
  	}
}


void SaveImage(Image* Source, char* name) {
	TIFF*     	File;
  	File = TIFFOpen (name, "w");	
   	if (File == (TIFF*)NULL) {
    		printf ("Could not open file named %s for writing.\n", name);
    		exit (1);
  	}
	SavePixels (File, Source);	
	TIFFClose (File);			
}


uint8 f2i(float f) {
	long i = floor(f+0.5);
	if (i>=255) return (uint8)255;
	if (i<=0) return (uint8)0;
	return (uint8)i;
}

uint8 truncint (int i) {
	if (i<=0) return 0;
	if (i>=255) return 255;
	return (uint8)i;
}


/*convert the buffer Image.dataF into the buffer Image.dataI for writing the image as a TIFF file.
The image will be written as a COLOR image, with three equal R,G,B channels.*/
void Float_to_int(Image *Source) {
	float* 	curs = Source->dataF;
	uint8*  pixd;
	uint8   aux;
	int 	i,j;
	int 	width = Source->width;
	int 	height = Source->height;
	
	Source->dataI=(uint32*)Check(_TIFFmalloc(width * height * sizeof(uint32)));
	pixd = (uint8*)(Source->dataI);
	for (i=0;i<width*height;i++) {
		aux = f2i(*curs++);
		*pixd++ = aux;
		*pixd++ = aux;
		*pixd++ = aux;
		*pixd++ = 255;			/*alpha channel*/
	}
}


/*convert the raw data of the TIFF file (uint8 pixel values) into a buffer of floats*/
void Int_to_float(Image* Source) {
	uint8* 	curs;
	float*  pixd;
	int 	width = Source->width;
	int 	height = Source->height;
	int 	i,j;
	
	Source->dataF=(float*)Check(malloc(width*height*sizeof(float)));
	pixd = Source->dataF;
	curs = (uint8*)(Source->dataI);
	for (i=0; i<height*width; i++) {
		*pixd++ = *curs;
		curs+=4;
	}
}


/* Grayscale Images assumed.
compute the MSE and PSNR between Source1 and Source2. If Diff!=NULL, the difference is put in Diff as a grayscale image
the value border allows not to take into account the first and last <border> rows and columns of the image for the computation of the error*/
void Difference(Image* Source1, Image* Source2, Image* Diff, int border) {
	uint8* 	pixs1 = (uint8*)(Source1->dataI);
	uint8* 	pixs2 = (uint8*)(Source2->dataI);
	uint8* 	pixd;
	int 	i,j,k;
	unsigned long 	erreur2 = 0;
	double  e1, e2, em, PSNR;
	int		diff;
	int    	width = Source1->width;
	int		height = Source1->height;
	
	if (width != Source2->width) {
		printf("Error: images should have the same width for computing their difference: %d!=%d\n",width,Source2->width);
		exit(1);
	}
	else if (height != Source2->height) {
		printf("Error: images should have the same height for computing their difference: %d!=%d\n",height,Source2->height);
		exit(1);
	}
	if (Diff != NULL) {
		Diff->dataI=(uint32*)Check(_TIFFmalloc (width*height*sizeof(uint32)));
		pixd = (uint8*)(Diff->dataI);
		Diff->width = width;
		Diff->height = height;
	}
	for (i=0; i<height; i++) {
		for (j=0; j<width; j++) {
			diff = *pixs2 - *pixs1;
			if ((i>=border) && (i<height-border) && (j>=border) && (j<width-border)) 
				erreur2 += diff*diff;
			if (Diff != NULL) {
				*(pixd++) = truncint(diff+128);
				*(pixd++) = truncint(diff+128);
				*(pixd++) = truncint(diff+128);
				*(pixd++) = 255;
			}
			pixs1 += 4; 
			pixs2 += 4;
		}
	}		
 	e2 = (double)erreur2 / ((width-2*border) * (height-2*border));
	PSNR = 10.0*log10(255.0*255.0/e2);
	printf("MSE :     %f\n",e2);
	printf("PSNR :     %f\n",PSNR);
}


/*crop the image by removing rows and columns at the borders*/
int Crop(Image* Source, Image* Dest, int newwidth, int newheight) {
	int 	width = Source->width;
	int 	height = Source->height;
	int		i,j,k,l;
	float* 	Buffer;
	float* 	pixd;
	float*	pixs;

	Check(Buffer = (float*)Check(malloc(newwidth*newheight*sizeof(float))));
	pixd = Buffer;
	pixs = Source->dataF;

	pixs += (height-newheight)/2*width;
	for (i=0; i<newheight; i++) {
		memcpy(pixd, pixs+(width-newwidth)/2, newwidth*sizeof(float));
		pixs += width;
		pixd += newwidth;
	}
	if (Dest == NULL) {
		Source->width = newwidth;
		Source->height = newheight;
		free(Source->dataF);
		Source->dataF = Buffer;
	} else {
		Dest->width = newwidth;
		Dest->height = newheight;
		Dest->dataF = Buffer;
	}
}


/*extend the size of the image by adding rows and columns at the borders. The new pixel values are computed by mirror symmetry*/
void Extend(Image* Source, Image* Dest, int borderx, int bordery) {
	int width = Source->width;
	int height = Source->height;
	int i,j;
	int newwidth = width+2*borderx;
	int newheight = height+2*bordery;
	float* Buffer=(float*)Check(malloc(newwidth*newheight*sizeof(float)));
	for (i=0; i<bordery; i++) 
		for (j=0; j<width; j++) 
			Buffer[i*newwidth+j+borderx]=Source->dataF[(bordery-i)*width+j];
	for (i=0; i<height; i++) 
		for (j=0; j<width; j++) 
			Buffer[(bordery+i)*newwidth+j+borderx]=Source->dataF[i*width+j];
	for (i=0; i<bordery; i++) 
		for (j=0; j<width; j++) 
			Buffer[(bordery+height+i)*newwidth+j+borderx]=Source->dataF[(height-2-i)*width+j];
	for (i=0; i<borderx; i++) 
		for (j=0; j<newheight; j++) {
			Buffer[j*newwidth+i]=Buffer[j*newwidth+2*borderx-i];
			Buffer[j*newwidth+borderx+width+i]=Buffer[j*newwidth+borderx+width-2-i];
		}		
	if (Dest==NULL) {
		Source->width = newwidth;
		Source->height = newheight;
		free(Source->dataF);
		Source->dataF = Buffer;
	} else {
		Dest->dataF = Buffer;
		Dest->width = newwidth;
		Dest->height = newheight;
	}
}

/*rotation using type 2 order 2 filter for the shifts (see article in IEEE Trans. Image Proc.) */
void Rotation_shears(Image* Source, Image* Dest, float angle) {
	int 	i, j, width, height;
	float* 	p; 
	float* 	Buffer;
	float 	centerX, centerY; 	/*rotation center*/ 
	float* 	output;
	float* 	p1;	/*the current pixel*/
	float* 	p2; /*the next pixel*/
	float* 	p3; /*the previous pixel*/
	float* 	p4; /*the last pixel*/
	int 	e;
	float   aux, val, tau;

	
	width=Source->width;
	height=Source->height;
	centerX=width/2.0;	/*the rotation center is put at the center of the image by default. This can be changed.*/
	centerY=height/2.0;
	if (Dest==NULL)
		p=Source->dataF;
	else {
		Buffer=(float*)Check(malloc(width*height*sizeof(float)));
		memcpy(Buffer, Source->dataF, width*height*sizeof(float));
		Dest->width=width;
		Dest->height=height;
		Dest->dataF = Buffer;
		p=Buffer;
	}
	output = p;
	for (i=0; i<height; i++) {		/*first horizontal shear = shifts on the rows*/
		tau=tan(angle/2)*(i-centerY);
		e=(int)floorf(tau+0.5);			/*integer part of the shift*/
		if ((tau>0)&&(tau-floorf(tau)==0.5)) e--;	/*for symmetry of the rounding operation if tau==0.5*/
		tau-=e;	/*remaining part of the shift*/		
		val = tau*tau;
		aux = 4.0 - val - sqrt(12.0-3.0*val);
		val = aux/(3.0*tau+2.0+val);						
		p3=output;
		p1=p3+1;
		p2=p1+1;
		p4=output+width-1;
		/*nothing is done on the first pixel value*/
		while (p1<p4) 
			*p1++ += (*p2++ - *p3++)*val;
		val = aux/(-3.0*tau+2.0+tau*tau);
		/*nothing is done on the last pixel value*/
		p3=p1;
		p1--;
		p2=p1-1;
		while (p1>output) 
			*p1-- += (*p2-- - *p3--)*val;
		*p1 += (*p1 - *p3)*val;		
		if (e<0) {
			for (j=0; j<width+e; j++)
				output[j]=output[j-e];
			for (; j<width; j++)
				output[j]=127.5;
		} else {
			for (j=width-1; j>=e; j--)
				output[j]=output[j-e];
			for (; j>=0; j--)
				output[j]=127.5;
		}
		output+=width;
	}
	output = p;
	for (i=0; i<width; i++) {		/*second vertical shear = shifts on the columns*/	
		tau=-sin(angle)*(i-centerX);
		e=(int)floorf(tau+0.5);			/*integer part of the shift*/
		if ((tau>0)&&(tau-floorf(tau)==0.5)) e--;	/*for symmetry of the rounding operation if tau==0.5*/
		tau-=e;	/*remaining part of the shift*/
		val = tau*tau;
		aux = 4.0 - val - sqrt(12.0-3.0*val);
		val = aux/(3.0*tau+2.0+val);						
		p3=output;
		p1=p3+width;
		p2=p1+width;
		p4=output+(height-1)*width;
		while (p1<p4) {
			*p1 += (*p2 - *p3)*val;
			p1+=width;
			p2+=width;
			p3+=width;
		}
		val = aux/(-3.0*tau+2.0+tau*tau);
		p3=p1;
		p1-=width;
		p2=p1-width;
		while (p1>output) {
			*p1 += (*p2 - *p3)*val;
			p1-=width;
			p2-=width;
			p3-=width;
		}
		*p1 += (*p1 - *p3)*val;
		if (e<0) {
			for (j=0; j<height+e; j++)
				output[j*width]=output[(j-e)*width];
			for (; j<height; j++)
				output[j*width]=127.5;
		} else {
			for (j=height-1; j>=e; j--)
				output[j*width]=output[(j-e)*width];
			for (; j>=0; j--)
				output[j*width]=127.5;
		}
		output++;
	}

	output = p;
	for (i=0; i<height; i++) {		/*third horizontal shear = shifts on the rows*/
		tau=tan(angle/2)*(i-centerY);
		e=(int)floorf(tau+0.5);			/*integer part of the shift*/
		if ((tau>0)&&(tau-floorf(tau)==0.5)) e--;	/*for symmetry of the rounding operation if tau==0.5*/
		tau-=e;	/*remaining part of the shift*/
		val = tau*tau;
		aux = 4.0 - val - sqrt(12.0-3.0*val);
		val = aux/(3.0*tau+2.0+val);						
		p3=output;
		p1=p3+1;
		p2=p1+1;
		p4=output+width-1;
		while (p1<p4) 
			*p1++ += (*p2++ - *p3++)*val;
		val = aux/(-3.0*tau+2.0+tau*tau);
		p3=p1;
		p1--;
		p2=p1-1;
		while (p1>output) 
			*p1-- += (*p2-- - *p3--)*val;
		*p1 += (*p1 - *p3)*val;
		if (e<0) {
			for (j=0; j<width+e; j++)
				output[j]=output[j-e];
			for (; j<width; j++)
				output[j]=127.5;
		} else {
			for (j=width-1; j>=e; j--)
				output[j]=output[j-e];
			for (; j>=0; j--)
				output[j]=127.5;
		}
		output+=width;
	}
}


/*this test program loads an image, rotate it, then rotate it back and computes the difference with the initial image to show that the latter is perfectly recovered, so that the rotation is truly lossless.*/
int main (int argc, char* argv[]) {
	Image	  	Source = EMPTY;
	Image 	  	Dest = EMPTY;
	Image 	  	Diff = EMPTY;
	int 		width, height;
	int			newsize;
	uint32* 	dataI;
	
	if (argc != 2) {	/*error in the command line*/
   		printf ("Bad number of parameters: should be \"%s <TIFF image filename> \"\n", argv[0]);
 		exit (1);
  	}
	
	LoadImage (&Source, argv[1]);
	SaveImage(&Source,"0_original.tif");
	width = Source.width;		
	height = Source.height;	
	printf ("Image size: width=%d, height=%d.\n",width,height);
	/*SaveImage(&original,"orig.tif");*/


	Int_to_float(&Source);		/*process only the first channel: a grayscale image is assumed. 
	Else, the red channel would be processed.*/
	
	if (EXTEND) {
		/*the image is first extended, so that the shears do not cut parts of the image. 
		The number of pixels added is not optimized: the minimal number of pixels it is necessary to add depends on the angle*/
		int		newsize = (int)ceilf(sqrt(width*width+height*height)/2.0)*2;
		dataI=Source.dataI;	/*we save the buffer. We use it at the end for comparison with the 2 times rotated image*/
		Extend(&Source,NULL, (int)ceilf(sqrt(width*width+height*height)/2.0)*2-width, 
			(int)ceilf(sqrt(width*width+height*height)/2.0)*2-height);
		Float_to_int(&Source);
		SaveImage(&Source,"1_extended.tif");
		_TIFFfree ((tdata_t)(Source.dataI));
		Source.dataI=dataI;
	}
	Rotation_shears(&Source, NULL, ANGLE);	/*first rotation. In-place*/
	dataI=Source.dataI;
	Float_to_int(&Source);
	SaveImage(&Source,"2_rotated.tif");
	
	/*the following code can be used to compute the PSNR between the initial and rotated image, 
	for a rotationnally invariant image (zoneplate)*/
	if (0) {
		if (!EXTEND) {
			Dest.width=width;
			Dest.height=height;
			Dest.dataI=dataI;
			Difference(&Dest, &Source, &Diff, 0);		
			SaveImage(&Diff, "2_diff.tif");
			_TIFFfree ((tdata_t)(Diff.dataI));
		}
	}
	
	_TIFFfree ((tdata_t)(Source.dataI));
	Source.dataI=dataI;
	
	if (EXTEND) {
		Crop(&Source, &Dest, width, height);
		Float_to_int(&Dest);			
		free(Dest.dataF);
		SaveImage(&Dest,"2_rotated_cropped.tif");
		_TIFFfree ((tdata_t)(Dest.dataI));
	}
	Rotation_shears(&Source, &Dest, -ANGLE);	/*second rotation*/
	free(Source.dataF);
	Float_to_int(&Dest);			
	SaveImage(&Dest,"3_rotated_back.tif");
	if (EXTEND) {
		_TIFFfree ((tdata_t)(Dest.dataI));
		Crop(&Dest, NULL, width, height);
		Float_to_int(&Dest);
		SaveImage(&Dest,"3_rotated_back_cropped.tif");
		Source.width=width;
		Source.height=height;
	}
	free(Dest.dataF);	
	Difference(&Dest, &Source, &Diff, 0);		/*compute the PSNR*/
	SaveImage(&Diff, "4_diff.tif");
	_TIFFfree ((tdata_t)(Diff.dataI));
	_TIFFfree ((tdata_t)(Dest.dataI));
	_TIFFfree ((tdata_t)(Source.dataI));
}



