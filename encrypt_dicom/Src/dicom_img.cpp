#include "dicom_img.h"
#include "encrypt.h"

void head_init(file_head * p_head)
{
    p_head->bytenum=0;
    p_head->maxbyte=1*100;
    p_head->header=new unsigned char [100];
}

void write_file_head(file_head * p_head, unsigned char * buffer,int data_length)
{
    while(p_head->bytenum+data_length>p_head->maxbyte)
    {
        unsigned char * new_head;
        new_head = new unsigned char [(p_head->maxbyte)*2];
        memcpy(new_head, p_head->header,p_head->bytenum);
        delete []p_head->header;
        p_head->header=new_head;
        p_head->maxbyte*=2;
    }
    memcpy(&p_head->header[p_head->bytenum],buffer,data_length);
    p_head->bytenum+=data_length;
}

dicom_img::dicom_img()
{
    //ctor
}

dicom_img::~dicom_img()
{
    //dtor
}

void dicom_img ::readFile(string filename)
{
    img_name=filename;
    ifstream in_file;
    in_file.open(filename,ios::in|ios::binary);
    if(!in_file)
    {
        cout<<"can not open file "<<filename<<endl;
    }
    file_head temp_head;
    unsigned char introduce [128];

    head_init(&temp_head);
    in_file.read((char*)introduce,sizeof(char)*128);
    write_file_head(&temp_head,introduce,128);

    unsigned int tag=0;
    while(tag!=0x7fe0)
    {
        in_file.read((char*)&tag,sizeof(char)*2);
        write_file_head(&temp_head,(unsigned char *)&tag,2);
    }
    in_file.read((char*)&tag,sizeof(char)*2);
    write_file_head(&temp_head,(unsigned char *)&tag,2);
    unsigned int datasize=0;
    in_file.read((char*)&datasize,sizeof(char)*4);
    write_file_head(&temp_head,(unsigned char *)&datasize,4);

    pf_head= new unsigned char [temp_head.bytenum];
    memcpy(pf_head,temp_head.header,sizeof(char)*temp_head.bytenum);
    head_length=temp_head.bytenum;

    width=512;height=512;

    pixels= new unsigned int *[height];
    for(int i=0;i<height;i++)
    {
        pixels[i]= new unsigned int [width];
        for(int j=0;j<width;j++)
        {
            in_file.read((char*)&tag,sizeof(char)*2);
            pixels[i][j]=tag;
        }
    }
    in_file.close();
}

void dicom_img ::saveFile()
{
    ofstream out_file;
    out_file.open(img_name,ios::out|ios::binary);
    if(!out_file)
    {
        cout<<"can not open file "<<img_name<<endl;
    }
    out_file.write((char*)pf_head,head_length);

    for(int i=0;i<height;i++)
    {
        for(int j=0;j<width;j++)
        {
            unsigned int pixel=pixels[i][j];
            char temp=0;
            for(int k=0;k<2;k++)
            {
                temp=(pixel&0x000000ff);
                out_file.write((char*)&temp,sizeof(char));
                pixel=pixel>>8;
            }
        }
    }
    out_file.close();
}

string dicom_img::getName()
{
    return img_name;
}

void dicom_img::getSize(int * img_size)
{
    img_size[0]=height;
    img_size[1]=width;
}

void dicom_img::getPixels(unsigned int ** pixelbuff)
{
    for(int i=0;i<height;i++)
    {
        memcpy(pixelbuff[i],pixels[i],sizeof(unsigned int)*width);
    }
}

void dicom_img::writePixels(unsigned int ** pixelbuff)
{
    for(int i=0;i<height;i++)
    {
        memcpy(pixels[i],pixelbuff[i],sizeof(unsigned int)*width);
    }
}

void dicom_img::encrypt(string key,string append)
{
    double ** chaos_x;
    chaos_x= new double* [height];
    for(int i=0;i<height;i++)
    {
        chaos_x[i]= new double [width];
    }
    int img_size[2]={height,width};
    randNumCreate(chaos_x,img_size,key);

    permutate(this,chaos_x);
    diffusion(this,chaos_x);

    string dcm=".dcm";
    string en="_encrypted";
    int index=img_name.find(dcm);
    if(append.size()!=0)
        img_name.insert(index,en+=append);
    else
         img_name.insert(index,en);

    for(int i=0;i<height;i++)
    {
        delete []chaos_x[i];
    }
    delete []chaos_x;
}


