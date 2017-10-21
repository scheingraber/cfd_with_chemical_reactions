#ifndef PGM_T9W0L0U0
#define PGM_T9W0L0U0

/**
 * reads in a ASCII pgm-file and returns the colour information in a two-dimensional integer array.
 * At this, a boundary layer around the image is additionally stored and initialised with 0. 
 */
// implementation has to be inside the header, since the function is a template
template <typename T> T **read_pgm(const char *filename, int size[])
{
    FILE *input = NULL;
    char line[1024];
    int levels;
    int xsize, ysize;
    int i1, j1;
    T **pic = NULL;

    if ((input=fopen(filename,"rb"))==0)
    {
       char szBuff[80];
       sprintf( szBuff, "Can not read file %s !!!", filename );
       ERROR( szBuff );
    }

    /* check for the right "magic number" */
    if ( fread(line,1,3,input)!=3 )
    {
        fclose(input);
        ERROR("Error Wrong Magic field!");
    }

    /* skip the comments */
    do
    fgets(line,sizeof line,input);
    while(*line=='#');

    /* read the width and height */
    sscanf(line,"%d %d\n",&xsize,&ysize);

    // Pass size outside
    size[0] = xsize;
    size[1] = ysize;

    #ifdef DEBUG
    printf("Image size: %d x %d\n", xsize,ysize);
    #endif // DEBUG

    /* read # of gray levels */
    fgets(line,sizeof line,input);
    sscanf(line,"%d\n",&levels);

    /* allocate memory for image */
    pic = matrix<T>(0,xsize+2,0,ysize+2);

    #ifdef DEBUG
    printf("Image initialised...\n");
    #endif // DEBUG

    /* read pixel row by row */
    for(j1=1; j1 < ysize+1; j1++)
    {
        for (i1=1; i1 < xsize+1; i1++)
        {
            int byte;
            fscanf(input, "%d", &byte);

            if (byte==EOF)
            {
                fclose(input);
                ERROR("read failed");
            }
            else
            {
                pic[i1][ysize+1-j1] = byte;
                #ifdef DEBUG
                printf("%d,%d: %d\n", i1,ysize+1-j1,byte);
                #endif // DEBUG
            }
         }
    }
    for (i1 = 0; i1 < xsize+2; i1++)
    {
        pic[i1][0] = 0;
    }
    for (i1 = 0; i1 < xsize+2; i1++)
    {
        pic[i1][ysize+1] = 0;
    }
    for (j1 = 0; j1 < ysize+2; j1++)
    {
        pic[0][j1] = 0;
        pic[xsize+1][j1] = 0;
    }

    /* close file */
    fclose(input);

    return pic;
}


#endif /* end of include guard: PGM_T9W0L0U0 */
