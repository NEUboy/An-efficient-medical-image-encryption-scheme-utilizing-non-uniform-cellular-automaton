#ifndef TEST_H
#define TEST_H
#include <iostream>
#include <iomanip>
#include <string>
#include <time.h>
#include "image.h"

using namespace std;
static void key_to_bit(string key, unsigned char bit_buff[256]);
static void bit_to_key(string &key, unsigned char bit_buff[256]);

class test
{
    public:
        test();
        virtual ~test();
        void pix_cov(image * img1);
        void NPCR_UACI(image *img1,image *img2);
        void MSE_PSNR(image *img1,image *img2);
        void diff_attack_test(image *img,string key);
        void encode_key_sensitivity(image *img,string key,int * id);
        void decode_key_sensitivity(image *img,string key,int * id);
        double block_entropy_test(image *img);
        void global_entropy_test(image *img);
        void entropy_test(image *img);
        void Chi_square_test(image &img);
        void robust_test(image *img,string key);

    protected:

    private:
};

#endif // TEST_H
