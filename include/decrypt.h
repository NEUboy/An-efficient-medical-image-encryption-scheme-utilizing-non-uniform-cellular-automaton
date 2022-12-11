#ifndef DECRYPT_H_INCLUDED
#define DECRYPT_H_INCLUDED
#include <iostream>
#include <stdlib.h>
#include "image.h"

using namespace std;
static void key_to_bit(string key, unsigned char bit_buff[256]);
static double bit_to_val(unsigned char bit_buff[256],int start, int stop);
static unsigned int bit_to_int(unsigned char bit_buff[256],int start, int stop);

void derandNumCreate(double ** chaos_x,double ** chaos_y,int img_size[2],string key);
void depretreatment(image *img,double ** chaos_x);
void depermutate(image *img,double ** chaos_x);
void dediffusion(image *img,double ** chaos_y);


#endif // DECRYPT_H_INCLUDED
