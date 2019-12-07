#ifndef _TIFF_UTILS_H
#define _TIFF_UTILS_H

#include <tiffio.h>

/*Our image structure*/
typedef struct {
	int width, height;	
	float* dataF;		/*Buffer of floats*/
	uint32* dataI;		/*Buffer of ints (raw data of the TIFF file)*/
} Image;

#define 	EMPTY  {0,0,NULL,NULL}   

void* Check(void* pointer);

void LoadImage(Image* Dest, char* name);

void SaveImage(Image* Source, char* name);

#endif
