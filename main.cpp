#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "image.h"
#include "test.h"
#include "load.h"

using namespace std;

int main()
{
    vector<string> image_set;
    read_imageset(image_set,"fileset.txt");
    string key="EFC796D47FDF6A69AB7DF3DFAF3CE7AFDEAEFC6979757FC9DA69D93F4D76FC7F";
    test te;
    image * pimg1;
    for(int i=0;i<image_set.size();i++)
    {
        cout<<"****************************************";
        cout<<image_set[i];
        cout<<"****************************************"<<endl;

        load_img(&pimg1,image_set[i]);
        pimg1->readFile(image_set[i]);
        //te.entropy_test(pimg1);
        //te.pix_cov(pimg1);
        //pimg1->encrypt(key,"");
        //pimg1->saveFile();
        //pimg1->decrypt(key,"");
        //pimg1->saveFile();


        te.robust_test(pimg1,key);
        //te.diff_attack_test(pimg1,key);
        //int id[4];
        //te.encode_key_sensitivity(pimg1,key,id);
        //te.decode_key_sensitivity(pimg1,key,id);
        //te.pix_cov(img1);
        //img1.encrypt(key,"");
        //img1.saveFile();

        //img1.decrypt(key,"");
        //img1.saveFile();
        //img1.encrypt(key,"");
        //img1.saveFile();
        //te.pix_cov(pimg1);
        //te.Chi_square_test(img1);
        //te.entropy_test(pimg1);
    }
    return 0;
}

