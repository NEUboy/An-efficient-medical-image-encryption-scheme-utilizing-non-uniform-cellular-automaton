#ifndef LOAD_H_INCLUDED
#define LOAD_H_INCLUDED
#include <iostream>
#include <string>
#include <vector>
#include "image.h"
#include "image_24.h"

void read_imageset(vector<string> &img_set,string fileset)
{
    string temp;
    ifstream in_file;
    in_file.open(fileset,ios::in);

    while(getline(in_file,temp))
    {
        img_set.push_back(temp);
    }
    in_file.close();
}

void load_img(image ** pimg,string filename)
{
    ifstream in_file;
    in_file.open(filename,ios::in|ios::binary);
    if(!in_file)
    {
        cout<<"can not open file "<<filename<<endl;
    }
    unsigned short  fileType;
    in_file.read((char*)&fileType,sizeof(unsigned short));
    if(fileType==0x4d42)
    {
        ClBitMapFileHeader filehead;
        CliBitMapInfoHeader infohead;
        in_file.read((char*)&filehead,sizeof(ClBitMapFileHeader));
        in_file.read((char*)&infohead,sizeof(CliBitMapInfoHeader));

        if(infohead.biBitCount==8)
        {
            * pimg = new image;
        }
        else if(infohead.biBitCount>=24)
        {
            * pimg =new image_24;
        }
    }
    in_file.close();
}

#endif // LOAD_H_INCLUDED
