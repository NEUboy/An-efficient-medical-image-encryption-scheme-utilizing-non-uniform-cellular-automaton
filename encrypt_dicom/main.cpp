#include <iostream>
#include "dicom_img.h"

using namespace std;

int main()
{
    string key="EFC796D47FDF6A69AB7DF3DFFF3CE7AFDEFEFC6979757FC9DA69D93F4D76FC7F";
    dicom_img* pimg1;
    pimg1=new dicom_img;
    pimg1->readFile("case6b_001.dcm");
    pimg1->encrypt(key,"");
    pimg1->saveFile();
    return 0;
}
