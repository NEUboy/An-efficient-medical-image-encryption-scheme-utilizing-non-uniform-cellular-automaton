#include "image.h"
#include "encrypt.h"
#include "decrypt.h"

using namespace std;
image::image()
{
    //ctor
}

void image::readFile(string filename)
{
    img_name=filename;
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
        in_file.read((char*)&fileHeader,sizeof(ClBitMapFileHeader));
        in_file.read((char*)&infoHeader,sizeof(CliBitMapInfoHeader));
        for(int i=0;i<256;i++)
        {
            in_file.read((char*)&rgbHeader,sizeof(ClRgbQuad));
        }
    }
    height=infoHeader.biHeight;
    width=infoHeader.biHeight;

    pixels= new unsigned char* [height];
    for(int i=0;i<height;i++)
    {
        pixels[i]= new unsigned char [width];
        in_file.read((char*)pixels[i],width*sizeof(char));
    }
    in_file.close();
}

void image::saveFile()
{
    ofstream out_file;
    out_file.open(img_name,ios::out|ios::binary);
    if(!out_file)
    {
        cout<<"can not open file "<<img_name<<endl;
    }
    unsigned short filetype;
    filetype=0x4D42;
    out_file.write((char*)&filetype,sizeof(unsigned short));
    out_file.write((char*)&fileHeader,sizeof(ClBitMapFileHeader));
    out_file.write((char*)&infoHeader,sizeof(CliBitMapInfoHeader));

    for(int i=0;i<256;i++)
    {
        ClRgbQuad pBmpRgbQuad;
        pBmpRgbQuad.rgbBlue=i;
        pBmpRgbQuad.rgbGreen=i;
        pBmpRgbQuad.rgbRed=i;
        pBmpRgbQuad.rgbReserved=0;
        out_file.write((char*)&pBmpRgbQuad,sizeof(ClRgbQuad));
    }

    for(int i=0;i<height;i++)
    {
        out_file.write((char*)pixels[i],width*sizeof(char));
    }
    out_file.close();
}

string image::getName()
{
    return img_name;
}

void image::getSize(int * img_size)
{
    img_size[0]=height;
    img_size[1]=width;
}

void image::getPixels(unsigned char ** pixelbuff)
{
    for(int i=0;i<height;i++)
    {
        for(int j=0;j<width;j++)
        {
            pixelbuff[i][j]=pixels[i][j];
        }
    }
}

void image::writePixels(unsigned char ** pixelbuff)
{
    for(int i=0;i<height;i++)
    {
        memcpy(pixels[i],pixelbuff[i],width);
    }
}

void image::copy_img(image * pimg)
{
    img_name=pimg->img_name;
    memcpy(&fileHeader,&pimg->fileHeader,sizeof(ClBitMapFileHeader));
    memcpy(&infoHeader,&pimg->infoHeader,sizeof(CliBitMapInfoHeader));
    memcpy(&rgbHeader,&pimg->rgbHeader,sizeof(ClRgbQuad));
    width=pimg->width;
    height=pimg->height;
    pixels=new unsigned char * [height];
    for(int i=0;i<height;i++)
    {
        pixels[i]=new unsigned char [width];
        memcpy(pixels[i],pimg->pixels[i],width);
    }
}

void image::operator=(image & original)
{
    img_name=original.img_name;
    memcpy(&fileHeader,&original.fileHeader,sizeof(ClBitMapFileHeader));
    memcpy(&infoHeader,&original.infoHeader,sizeof(CliBitMapInfoHeader));
    memcpy(&rgbHeader,&original.rgbHeader,sizeof(ClRgbQuad));
    width=original.width;
    height=original.height;
    pixels=new unsigned char * [height];
    for(int i=0;i<height;i++)
    {
        pixels[i]=new unsigned char [width];
        memcpy(pixels[i],original.pixels[i],width);
    }
}

void image::encrypt(string key,string append)
{
    double ** chaos_x;
    double ** chaos_y;
    chaos_x= new double* [height];
    chaos_y= new double* [height];
    for(int i=0;i<height;i++)
    {
        chaos_x[i]= new double [width];
        chaos_y[i]= new double [width];
    }
    int img_size[2]={height,width};
    randNumCreate(chaos_x,chaos_y,img_size,key);
    permutate(this,chaos_x);
    diffusion(this,chaos_y);

    string bmp=".bmp";
    string en="_encrypted";
    int index=img_name.find(bmp);
    if(append.size()!=0)
        img_name.insert(index,en+=append);
    else
         img_name.insert(index,en);

    for(int i=0;i<height;i++)
    {
        delete []chaos_x[i];
        delete []chaos_y[i];
    }
    delete []chaos_x;
    delete []chaos_y;
}

void image::decrypt(string key,string append)
{
    double ** chaos_x;
    double ** chaos_y;
    chaos_x= new double* [height];
    chaos_y= new double* [height];
    for(int i=0;i<height;i++)
    {
        chaos_x[i]= new double [width];
        chaos_y[i]= new double [width];
    }
    int img_size[2]={height,width};

    derandNumCreate(chaos_x,chaos_y,img_size,key);
    dediffusion(this,chaos_y);
    depermutate(this,chaos_x);

    string en="_encrypted";
    string de="_decrypted";
    int index=img_name.find(en);
    img_name.erase(index,en.length());
    if(append.size()!=0)
        img_name.insert(index,de+=append);
    else
        img_name.insert(index,de);
    for(int i=0;i<height;i++)
    {
        delete []chaos_x[i];
        delete []chaos_y[i];
    }
    delete []chaos_x;
    delete []chaos_y;
}

void image::pixel_loss(double lossrate)
{
    int black_width=(int)(lossrate*width);
    int black_height=(int)(lossrate*height);

    for(int i=(width-black_width)/2; i<(width+black_width)/2;i++)
    {
        for(int j=(height-black_height)/2;j<(height+black_height)/2;j++)
        {
            pixels[i][j]=0;
        }
    }
    string bmp=".bmp";
    string loss="_loss";
    loss+=to_string(lossrate);
    int index=img_name.find(bmp);
    img_name.insert(index,loss);
}

void image::saltPepper_noise(double noiserate)
{
    int noiseNum=(int)(noiserate*width*height);
    srand(time(0));
    for(int i=0;i<noiseNum;i++)
    {
        int row=rand()%width;
        int col=rand()%height;
        int val=rand()%2;

        if(val==0)
            pixels[row][col]=0;
        else
            pixels[row][col]=255;
    }
    string bmp=".bmp";
    string noise="_noise";
    noise+=to_string(noiserate);
    int index=img_name.find(bmp);
    img_name.insert(index,noise);
}

void image::differential_attack(int pos[2])
{
    pixels[pos[0]][pos[1]]=pixels[pos[0]][pos[1]]^(0b00000001);
    string s_pos="(";
    s_pos.append(to_string(pos[0]));
    s_pos.append("_");
    s_pos.append(to_string(pos[1]));
    s_pos.append(")");
    string bmp=".bmp";
    int index=img_name.find(bmp);
    img_name.insert(index,s_pos);

}

image::~image()
{
    for(int i=0;i<height;i++)
    {
        delete pixels[i];
    }
    delete pixels;
    //dtor
}
