#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tiff_utils.h"

extern int SWAP_BOOL;

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
  	File = TIFFOpen (name, "r");	/*ouverture du fichier TIFF*/
  	if (File == (TIFF*)NULL) {
    		printf ("Could not open file named %s for reading.\n", name);
    		exit (1);
  	}
	LoadPixels (File, Dest);	/*chargement de l'image TIFF*/
	TIFFClose (File);			/*fermeture du fichier TIFF*/
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
  	File = TIFFOpen (name, "w");	/*ouverture du fichier TIFF*/
  	if (File == (TIFF*)NULL) {
    		printf ("Could not open file named %s for writing.\n", name);
    		exit (1);
  	}
	SavePixels (File, Source);	/*chargement de l'image TIFF*/
	TIFFClose (File);			/*fermeture du fichier TIFF*/
}


