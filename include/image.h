#ifndef IMAGE_H
#define IMAGE_H
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <string.h>
#include <cstdlib>
#include <math.h>
#include <time.h>

using namespace std;

typedef struct
{
    unsigned long   bfSize;
    unsigned short  bfReserved1;
    unsigned short  bfReserved2;
    unsigned long   bfOffBits;
}ClBitMapFileHeader;

typedef struct
{
    unsigned long   biSize;
    long            biWidth;
    long            biHeight;
    unsigned short  biPlanes;
    unsigned short  biBitCount;
    unsigned long   biCompression;
    unsigned long   biSizeImage;
    long            biXPelsPerMeter;
    long            biYPelsPerMeter;
    unsigned long   biClrUsed;
    unsigned long   biClrImportant;
}CliBitMapInfoHeader;

typedef struct
{
    unsigned char rgbBlue;
    unsigned char rgbGreen;
    unsigned char rgbRed;
    unsigned char rgbReserved;
}ClRgbQuad;

class image
{
    public:
        image();
        virtual ~image();
        virtual void readFile(string filename);
        virtual void saveFile();
        string getName();
        void getSize(int * img_size);
        void getPixels(unsigned char ** pixelbuff);
        void writePixels(unsigned char ** pixelbuff);
        void copy_img(image * pimg);
        void operator =(image & original);

        virtual void encrypt(string key,string append);
        virtual void decrypt(string key,string append);

        void pixel_loss(double lossrate);
        void saltPepper_noise(double noiserate);
        void differential_attack(int pos[2]);

    protected:
        string img_name;
        ClBitMapFileHeader fileHeader;
        CliBitMapInfoHeader infoHeader;
        ClRgbQuad rgbHeader;
        int width;
        int height;
        unsigned char ** pixels;

    private:

};

#endif // IMAGE_H
