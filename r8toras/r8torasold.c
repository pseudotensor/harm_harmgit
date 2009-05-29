
/* 
	convert "r8" file to sun raster file
	usage: r8toras nx ny palfile < infile > outfile

	cfg 8.22.98
*/

#define MAXFNAME  400		/* max length of filenames */

#include <stdio.h>
#include <math.h>
#include "rasterfile.h"

main(int argc, char *argv[])
{
	int nx, ny, nmap ;
	char dum ;
	char palfile[MAXFNAME] ;
	void decint(int *list, int num) ;
	struct rasterfile header ;
	FILE *fp ;

	if(argc < 4) {
		fprintf(stderr,"usage: r8toras nx ny palfile < infile > outfile\n") ;
		exit(0) ;
	}

	/* get args, open palette file */
	sscanf(argv[1],"%d",&nx) ;
	sscanf(argv[2],"%d",&ny) ;
	fp = fopen(argv[3],"r") ;
	if(fp == NULL) {
		fprintf(stderr,"trouble opening palette file.\n") ;
		exit(0) ;
	}

	/* set up header */
	header.ras_magic = RAS_MAGIC ;
	header.ras_width = nx ;
	header.ras_height = ny ;
	header.ras_depth = 8 ;
	header.ras_length = nx*ny ;
	header.ras_type = RT_STANDARD ;
	header.ras_maptype = RMT_EQUAL_RGB ;
	header.ras_maplength = 3*256 ;

	/* swap & write it */
	decint(&(header.ras_magic),8) ;
	fwrite(&(header.ras_magic), sizeof(int), 8, stdout) ;

	/* read & write the palette file */
	nmap = 0 ;
	while(fread(&dum, sizeof(char), 1, fp) == 1) {
		fwrite(&dum, sizeof(char), 1, stdout) ;
		nmap++ ;
	}
	decint(&(header.ras_maplength),1) ;
	if(nmap != header.ras_maplength) {
		fprintf(stderr,"err: header: %d, nmap: %d\n",
			header.ras_maplength, nmap) ;
	}

	/* now read & write the image */
	while(fread(&dum, sizeof(char), 1, stdin) == 1) {
		fwrite(&dum, sizeof(char), 1, stdout) ;
	}

	fprintf(stderr,"done\n") ;

	return(0) ;
}

void decint(int *lp, int n)
{
        unsigned int t;
        static unsigned long lmask = 0x00ff0000, rmask = 0x0000ff00;

        for(; n--; lp++) {
                t = *lp;
                *lp = (t >> 24) | (t << 24) | ((t << 8) & lmask) |
                        ((t >> 8) & rmask);
        }
}


