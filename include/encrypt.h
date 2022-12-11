#ifndef ENCRYPT_H_INCLUDED
#define ENCRYPT_H_INCLUDED
#include <iostream>
#include<cmath>
#include <stdlib.h>
#include "image.h"

using namespace std;

static void key_to_bit(string key, unsigned char bit_buff[256]);
static double bit_to_val(unsigned char bit_buff[256],int start, int stop);
static unsigned int bit_to_int(unsigned char bit_buff[256],int start, int stop);


void randNumCreate(double ** chaos_x,double ** chaos_y,int img_size[2],string key);
void permutate(image *img, double ** chaos_x);
void diffusion(image *img, double ** chaos_y);

#endif // ENCRYPT_H_INCLUDED
