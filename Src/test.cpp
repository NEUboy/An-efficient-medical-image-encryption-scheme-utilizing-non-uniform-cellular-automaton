#include "test.h"

static void key_to_bit(string key, unsigned char bit_buff[256])
{
    char hex_key[64];
    memset(hex_key,0,64);
    strcpy(hex_key,key.c_str());
    for(int i=0;i<64;i++)
    {
        if((hex_key[i]>=48)&&(hex_key[i]<=57))
        {
            char temp=0x08;
            for(int j=0;j<4;j++)
            {
                bit_buff[i*4+j]=(((hex_key[i]-48)&temp)==temp);
                temp=temp>>1;
            }
        }
        else
        {
            char temp=0x08;
            for(int j=0;j<4;j++)
            {
                bit_buff[i*4+j]=(((hex_key[i]-55)&temp)==temp);
                temp=temp>>1;
            }
        }
    }
}

static void bit_to_key(string &key, unsigned char bit_buff[256])
{
    for(int i=0;i<256;i+=4)
    {
        char temp=0;
        for(int j=i;j<i+4;j++)
        {
            temp=(temp<<1)+bit_buff[j];
        }
        if(temp<10)
            temp+=48;
        else
            temp+=55;

        key.push_back(temp);
    }
}

static void errorRate(image *img1,image *img2)
{
    double err=0;
    int img_size[2];
    img1->getSize(img_size);
    unsigned char ** pixelbuff1;
    unsigned char ** pixelbuff2;
    pixelbuff1=new unsigned char*[img_size[0]];
    for(int i=0;i<img_size[0];i++)
    {
        pixelbuff1[i]=new unsigned char [img_size[1]];
    }
    img1->getPixels(pixelbuff1);
    pixelbuff2=new unsigned char*[img_size[1]];
    for(int i=0;i<img_size[0];i++)
    {
        pixelbuff2[i]=new unsigned char [img_size[1]];
    }
    img2->getPixels(pixelbuff2);
    for(int i=0;i<img_size[0];i++)
    {
        for(int j=0;j<img_size[1];j++)
        {
            if(pixelbuff1[i][j]-pixelbuff2[i][j]!=0)
            {
                err+=1;
            }
        }
    }
    err/=(img_size[0]*img_size[1]);

    cout<<"error Rate test "<<img1->getName()<<"\t"<<img2->getName()<<endl;
    cout<<setiosflags(ios::fixed)<<setprecision(4)<<err*100<<"%"<<endl;

    for(int i=0;i<img_size[0];i++)
    {
        delete []pixelbuff1[i];
        delete []pixelbuff2[i];
    }
    delete []pixelbuff1;
    delete []pixelbuff2;
}

test::test()
{
    //ctor
}
void test::pix_cov(image * img)
{
    cout<<"----------------------------------";
    cout<<"Pixel correlation test";
    cout<<"----------------------------------"<<endl;
    unsigned char base [10000];
    unsigned char hor [10000];
    unsigned char ver [10000];
    unsigned char dia [10000];
    double cov[3];
    int img_size[2];
    img->getSize(img_size);
    unsigned char ** pixelbuff;
    pixelbuff=new unsigned char*[img_size[0]];
    for(int i=0;i<img_size[0];i++)
    {
        pixelbuff[i]=new unsigned char [img_size[1]];
    }
    img->getPixels(pixelbuff);

    srand(time(NULL));
    for(int i=0;i<10000;i++)
    {
        int x=rand()%(img_size[0]-2);
        int y=rand()%(img_size[1]-2);

        base[i]=pixelbuff[x][y];
        hor[i]=pixelbuff[x+1][y];
        ver[i]=pixelbuff[x][y+1];
        dia[i]=pixelbuff[x+1][y+1];
    }
    int sum_base=0,sum_hor=0,sum_ver=0,sum_dia=0;
    for(int i=0;i<10000;i++)
    {
        sum_base+=base[i];
        sum_hor+=hor[i];
        sum_ver+=ver[i];
        sum_dia+=dia[i];
    }
    double mean_base=(double)sum_base/10000;
    double mean_hor=(double)sum_hor/10000;
    double mean_ver=(double)sum_ver/10000;
    double mean_dia=(double)sum_dia/10000;
    int sum_base_var=0,sum_hor_var=0,sum_ver_var=0,sum_dia_var=0;
    for(int i=0;i<10000;i++)
    {
        sum_base_var+=(base[i]-mean_base)*(base[i]-mean_base);
        sum_hor_var+=(hor[i]-mean_hor)*(hor[i]-mean_hor);
        sum_ver_var+=(ver[i]-mean_ver)*(ver[i]-mean_ver);
        sum_dia_var+=(dia[i]-mean_dia)*(dia[i]-mean_dia);
    }
    double variance_base=(double)sum_base_var/10000;
    double variance_hor=(double)sum_hor_var/10000;
    double variance_ver=(double)sum_ver_var/10000;
    double variance_dia=(double)sum_dia_var/10000;
    int sum_cov_base_hor=0,sum_cov_base_ver=0,sum_cov_base_dia=0;
    for(int i=0;i<10000;i++)
    {
        sum_cov_base_hor+=(base[i]-mean_base)*(hor[i]-mean_hor);
        sum_cov_base_ver+=(base[i]-mean_base)*(ver[i]-mean_ver);
        sum_cov_base_dia+=(base[i]-mean_base)*(dia[i]-mean_dia);
    }
    double cov_base_hor=(double)sum_cov_base_hor/10000;
    double cov_base_ver=(double)sum_cov_base_ver/10000;
    double cov_base_dia=(double)sum_cov_base_dia/10000;
    cov[0]=cov_base_hor/(sqrt(variance_base)*sqrt(variance_hor));
    cov[1]=cov_base_ver/(sqrt(variance_base)*sqrt(variance_ver));
    cov[2]=cov_base_dia/(sqrt(variance_base)*sqrt(variance_dia));

    cout<<"correlation test "<<img->getName()<<endl;
    cout<<"cov_hor "<<setiosflags(ios::fixed)<<setprecision(6)<<cov[0]<<"\t";
    cout<<"cov_ver "<<setiosflags(ios::fixed)<<setprecision(6)<<cov[1]<<"\t";
    cout<<"cov_dia "<<setiosflags(ios::fixed)<<setprecision(6)<<cov[2]<<"\t";
    cout<<endl;

    for(int i=0;i<img_size[0];i++)
    {
        delete []pixelbuff[i];
    }
    delete []pixelbuff;
}

void test::NPCR_UACI(image *img1,image *img2)
{
    double rec[2];
    int img_size[2];
    img1->getSize(img_size);
    unsigned char ** pixelbuff1;
    unsigned char ** pixelbuff2;
    pixelbuff1=new unsigned char*[img_size[0]];
    for(int i=0;i<img_size[0];i++)
    {
        pixelbuff1[i]=new unsigned char [img_size[1]];
    }
    img1->getPixels(pixelbuff1);
    pixelbuff2=new unsigned char*[img_size[1]];
    for(int i=0;i<img_size[0];i++)
    {
        pixelbuff2[i]=new unsigned char [img_size[1]];
    }
    img2->getPixels(pixelbuff2);

    double NPCR=0;double UACI=0;
    for(int i=0;i<img_size[0];i++)
    {
        for(int j=0;j<img_size[1];j++)
        {
            if(pixelbuff1[i][j]-pixelbuff2[i][j]!=0)
            {
                NPCR+=1;
                UACI+=(double)abs(pixelbuff1[i][j]-pixelbuff2[i][j])/255;
            }
        }
    }
    NPCR/=(img_size[0]*img_size[1]);
    UACI/=(img_size[0]*img_size[1]);
    rec[0]=NPCR;rec[1]=UACI;

    cout<<"NPCR and UACI test: "<<img1->getName()<<"\t"<<img2->getName()<<endl;
    cout<<"NPCR: "<<setiosflags(ios::fixed)<<setprecision(4)<<rec[0]*100<<"%"<<"\t";
    cout<<"UACI: "<<setiosflags(ios::fixed)<<setprecision(4)<<rec[1]*100<<"%"<<"\t"<<endl;

    for(int i=0;i<img_size[0];i++)
    {
        delete []pixelbuff1[i];
        delete []pixelbuff2[i];
    }
    delete []pixelbuff1;
    delete []pixelbuff2;
}

void test::MSE_PSNR(image *img1,image *img2)
{
    double rec[2];
    int img_size[2];
    img1->getSize(img_size);
    unsigned char ** pixelbuff1;
    unsigned char ** pixelbuff2;
    pixelbuff1=new unsigned char*[img_size[0]];
    for(int i=0;i<img_size[0];i++)
    {
        pixelbuff1[i]=new unsigned char [img_size[1]];
    }
    img1->getPixels(pixelbuff1);
    pixelbuff2=new unsigned char*[img_size[1]];
    for(int i=0;i<img_size[0];i++)
    {
        pixelbuff2[i]=new unsigned char [img_size[1]];
    }
    img2->getPixels(pixelbuff2);

    double MSE=0;double PSNR=0;
    for(int i=0;i<img_size[0];i++)
    {
        for(int j=0;j<img_size[1];j++)
        {
            MSE+=(pixelbuff1[i][j]-pixelbuff2[i][j])*(pixelbuff1[i][j]-pixelbuff2[i][j]);
        }
    }
    MSE=MSE/(img_size[0]*img_size[1]);
    PSNR= 10*log10(255*255/MSE);
    cout<<"MSE and PSNR test: "<<img1->getName()<<"\t"<<img2->getName()<<endl;
    cout<<"MSE: "<<setiosflags(ios::fixed)<<setprecision(4)<<MSE<<endl;
    cout<<"PSNR: "<<setiosflags(ios::fixed)<<setprecision(4)<<PSNR<<endl;

    for(int i=0;i<img_size[0];i++)
    {
        delete []pixelbuff1[i];
        delete []pixelbuff2[i];
    }
    delete []pixelbuff1;
    delete []pixelbuff2;

}

void test::diff_attack_test(image *img,string key)
{
    cout<<"-------------------------------------";
    cout<<"diff attack test";
    cout<<"-------------------------------------"<<endl;
    image *cipher_img;
    cipher_img =new image;
    cipher_img->copy_img(img);
    cipher_img->encrypt(key,"\0");
    int img_size[2];
    img->getSize(img_size);

    image * diff_img;
    for(int i=0;i<5;i++)
    {
        diff_img= new image;
        diff_img->copy_img(img);
        srand(time(NULL));
        int pos[2];
        if(i==0)
        {
            pos[0]=511;pos[1]=511;
        }
        else
        {
            pos[0]=rand()%img_size[0];
            pos[1]=rand()%img_size[1];
        }
        diff_img->differential_attack(pos);
        diff_img->saveFile();
        diff_img->encrypt(key,"\0");
        diff_img->saveFile();
        NPCR_UACI(cipher_img,diff_img);
    }
}

void test::encode_key_sensitivity(image *img,string key,int* id)
{
    cout<<"----------------------------------";
    cout<<"encode key sensitivity";
    cout<<"----------------------------------"<<endl;
    image *cipher_img;
    cipher_img =new image;
    cipher_img->copy_img(img);
    cipher_img->encrypt(key,"\0");
    unsigned char bit_buff[256];
    key_to_bit(key,bit_buff);

    image test_img;
    srand(time(NULL));
    for(int i=0;i<4;i++)
    {
        image * test_img;
        test_img = new image;
        test_img->copy_img(img);
        int bit;
        bit=rand()%256;
        id[i]=bit;
        bit_buff[bit]=(bit_buff[bit]+1)%2;
        string diffkey;
        bit_to_key(diffkey,bit_buff);
        cout<<"key= "<<diffkey<<" ";
        cout<<bit<<" changed "<<endl;
        string append;
        append+="key";
        append+=to_string(bit);
        test_img->encrypt(diffkey,append);
        test_img->saveFile();
        NPCR_UACI(cipher_img,test_img);
        bit_buff[bit]=(bit_buff[bit]+1)%2;
    }
}

void test::decode_key_sensitivity(image *img,string key,int* id)
{
    cout<<"----------------------------------";
    cout<<"decode key sensitivity";
    cout<<"----------------------------------"<<endl;
    image *cipher_img;
    cipher_img =new image;
    cipher_img->copy_img(img);
    cipher_img->encrypt(key,"");
    unsigned char bit_buff[256];
    key_to_bit(key,bit_buff);

    image * test_img;
    srand(time(NULL));
    for(int i=0;i<4;i++)
    {
        test_img = new image;
        test_img->copy_img(cipher_img);
        int bit;
        //bit=rand()%256;
        bit=id[i];
        bit_buff[bit]=(bit_buff[bit]+1)%2;
        string diffkey;
        bit_to_key(diffkey,bit_buff);
        cout<<"key= "<<diffkey<<" ";
        cout<<bit<<" changed "<<endl;
        string append;
        append+="_key";
        append+=to_string(bit);
        test_img->decrypt(diffkey,append);
        test_img->saveFile();
        errorRate(img,test_img);
        bit_buff[bit]=(bit_buff[bit]+1)%2;
    }
}

double test::block_entropy_test(image *img)
{
    int ** coordinate_list;
    coordinate_list = new int* [30];
    int img_size[2];
    img->getSize(img_size);
    unsigned char ** pixelbuff;
    pixelbuff=new unsigned char*[img_size[0]];
    for(int i=0;i<img_size[0];i++)
    {
        pixelbuff[i]=new unsigned char [img_size[1]];
    }
    img->getPixels(pixelbuff);

    srand(time(NULL));
    for(int i=0;i<30;i++)
    {
        coordinate_list[i]=new int [2];
        int b_x;int b_y;
        b_x=rand()%(img_size[0]-44);
        b_y=rand()%(img_size[1]-44);
        int flag=1;
        while(flag)
        {
            int j;
            for(j=0;j<i;j++)
            {
                if((abs(coordinate_list[j][0]-b_x)<44)&&(abs(coordinate_list[j][1]-b_y)<44))
                    break;
            }
            if(j<i)                         //随机产生的坐标不符合要求
            {
                while(abs(coordinate_list[j][0]-b_x)<44)
                    b_x=rand()%(img_size[0]-44);
                while(abs(coordinate_list[j][1]-b_y)<44)
                    b_y=rand()%(img_size[1]-44);
            }
            else
                flag=0;
        }
        coordinate_list[i][0]= b_x;
        coordinate_list[i][1]= b_y;
    }
    double entropy_list[30];
    int  grayscale[256];
    for(int i=0;i<30;i++)
    {
        /*统计一个块里各个灰度的数量*/
        memset(grayscale,0,sizeof(int)*256);
        for(int x=0;x<44;x++)
        {
            for(int y=0;y<44;y++)
            {
                int idx=coordinate_list[i][0]+x;
                int idy=coordinate_list[i][1]+y;

                unsigned char pix=pixelbuff[idx][idy];
                grayscale[pix]++;
            }
        }
        double entro=0;
        for(int gr=0;gr<256;gr++)
        {
            double prob= (double) grayscale[gr]/(44*44) ;
            if(prob!=0)
                entro+=prob*(log(1/prob)/log(2));
        }
        entropy_list[i]=entro;
    }
    double sum=0;
    for(int i=0;i<30;i++)
        sum+=entropy_list[i];

    for(int i=0;i<img_size[0];i++)
    {
        delete []pixelbuff[i];
    }
    delete []pixelbuff;
    for(int i=0;i<30;i++)
    {
        delete []coordinate_list[i];
    }
    delete []coordinate_list;
    return sum/30;
}

void test::global_entropy_test(image *img)
{
    int img_size[2];
    img->getSize(img_size);
    unsigned char ** pixelbuff;
    pixelbuff=new unsigned char*[img_size[0]];
    for(int i=0;i<img_size[0];i++)
    {
        pixelbuff[i]=new unsigned char [img_size[1]];
    }
    img->getPixels(pixelbuff);

    int  grayscale[256];
    memset(grayscale,0,sizeof(int)*256);

    for(int x=0;x<img_size[0];x++)
    {
        for(int y=0;y<img_size[1];y++)
        {
            unsigned char pix=pixelbuff[x][y];
            grayscale[pix]++;
        }
    }
    double entro=0;
    for(int gr=0;gr<256;gr++)
    {
        double prob=(double )grayscale[gr]/(img_size[0]*img_size[1]);
        if(prob!=0)
            entro+=prob*(log(1/prob)/log(2));
    }
   cout<<"global entropy:\t"<<setiosflags(ios::fixed)<<setprecision(9)<<entro<<endl;
   for(int i=0;i<img_size[0];i++)
    {
        delete []pixelbuff[i];
    }
    delete []pixelbuff;

}

void test::entropy_test(image * img)
{
    cout<<"---------------------------------------";
    cout<<"entropy test";
    cout<<"---------------------------------------"<<endl;
    double block_entro=0;

    while((7.901901305>block_entro)||(block_entro>7.903037329))
    {
        block_entro=block_entropy_test(img);
    }


    //global_entropy_test(img);
    //cout<<"block entropy: \t"<<setiosflags(ios::fixed)<<setprecision(9)<<block_entropy_test(img)<<endl;
    cout<<"block entropy: \t"<<setiosflags(ios::fixed)<<setprecision(9)<<block_entro<<endl;
}

void test::Chi_square_test(image &img)
{
    cout<<"--------------------------------------";
    cout<<"Chi_square test";
    cout<<"--------------------------------------"<<endl;
    int img_size[2];
    img.getSize(img_size);
    unsigned char ** pixelbuff;
    pixelbuff=new unsigned char*[img_size[0]];
    for(int i=0;i<img_size[0];i++)
    {
        pixelbuff[i]=new unsigned char [img_size[1]];
    }
    img.getPixels(pixelbuff);

    int grayscale[256];
    memset(grayscale,0,sizeof(int)*256);
    for(int x=0;x<img_size[0];x++)
    {
        for(int y=0;y<img_size[1];y++)
        {
            unsigned char pix=pixelbuff[x][y];
            grayscale[pix]++;
        }
    }
    double chi=0;
    double exp=img_size[0]*img_size[1]/256;
    for(int i=0;i<256;i++)
    {
        chi+=(grayscale[i]-exp)*(grayscale[i]-exp)/exp;
    }
    cout<<"Chi test: "<<setiosflags(ios::fixed)<<setprecision(4)<<chi<<endl;

    for(int i=0;i<img_size[0];i++)
    {
        delete []pixelbuff[i];
    }
    delete []pixelbuff;
}

void test::robust_test(image *img,string key)
{
    cout<<"----------------------------------------";
    cout<<"Robust test";
    cout<<"----------------------------------------"<<endl;
    image *cipher_img;
    cipher_img =new image;
    cipher_img->copy_img(img);
    cipher_img->encrypt(key,"\0");

    image *test_img;

    for(int i=1;i<=10;i++)
    {
        test_img = new image;
        test_img->copy_img(cipher_img);
        test_img->pixel_loss(0.05*i);
        test_img->saveFile();
        test_img->decrypt(key,"\0");
        test_img->saveFile();
        MSE_PSNR(img,test_img);
    }

    for(int i=1;i<=10;i++)
    {
        test_img->copy_img(cipher_img);
        test_img->saltPepper_noise(0.0005*i);
        test_img->saveFile();
        test_img->decrypt(key,"\0");
        test_img->saveFile();
        MSE_PSNR(img,test_img);
    }

}

test::~test()
{
    //dtor
}
