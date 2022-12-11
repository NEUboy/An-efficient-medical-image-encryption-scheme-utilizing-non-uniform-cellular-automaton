# include "encrypt.h"

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

static double bit_to_val(unsigned char bit_buff[256],int start, int stop)
{
    double val=0;
    for(int i=start;i<stop;i++)
    {
        val+=bit_buff[i]*pow(0.5,i-start+1);
    }
    return val;
}

static unsigned int bit_to_int(unsigned char bit_buff[256],int start, int stop)
{
    unsigned int rec=0;
    for(int i=start;i<stop;i++)
    {
        rec=(rec<<1)+bit_buff[i];
    }
    return rec;
}

static void trans_pixs_bit(unsigned char * bitbuff, unsigned int * pix)
{
    unsigned int temp=1;
    for(int i=0;i<16;i++)
    {
        bitbuff[i]=((*pix&temp)==temp);
        temp=temp<<1;
    }
}

static void trans_bit_pix(unsigned char * bitbuff, unsigned int *pix)
{
    unsigned int val=0;
    for(int i=15;i>=0;i--)
    {
        val=(val<<1)+bitbuff[i];
    }
    *pix=val;
}

static void swap_double(double *a, double *b)
{
    double temp;
    temp=*a;*a=*b; *b=temp;
}

static void swap_int(int *a,int*b)
{
    int temp;
    temp=*a;*a=*b; *b=temp;
}

static void swap_char(unsigned char *a, unsigned char *b)
{
    char temp;
    temp=*a;*a=*b; *b=temp;
}

void randNumCreate(double ** chaos_x,int img_size[2],string key)
{
    double pi=3.1415926536;
    unsigned char bit_buff[256];
    key_to_bit(key,bit_buff);
    double xi=bit_to_val(bit_buff,0,32)*bit_to_val(bit_buff,32,64);
    double yi=bit_to_val(bit_buff,64,96)*bit_to_val(bit_buff,96,128);
    double theta=bit_to_val(bit_buff,128,160)*bit_to_val(bit_buff,160,192);
    unsigned int alpha= bit_to_int(bit_buff,192,224);
    unsigned int beta = bit_to_int(bit_buff,224,256);
    unsigned int mu;

    mu=alpha+beta;
    theta=theta*mu-floor(theta*mu);

    for (int i=0;i<30;i++)
    {
        xi=sin(pi*((4*theta*xi*(1-xi))+(1-theta)*sin(pi*yi)));
        yi=sin(pi*((4*theta*yi*(1-yi))+(1-theta)*sin(pi*xi)));
    }
    for(int i=0;i<img_size[0];i++)
    {
        for(int j=0;j<img_size[1];j+=2)
        {
            xi=sin(pi*((4*theta*xi*(1-xi))+(1-theta)*sin(pi*yi)));
            yi=sin(pi*((4*theta*yi*(1-yi))+(1-theta)*sin(pi*xi)));

            chaos_x[i][j]=xi;
            chaos_x[i][j+1]=yi;
        }
    }
}

static void divide_into_subfigs(unsigned int ** pixelbuff, unsigned int ** sub_pixbuffs[4], double ** chaosbuff, double ** sub_chaosbuff[4], int img_size[2])
{
    for(int i=0;i<img_size[0]/2;i++)
    {
        for(int j=0;j<img_size[1]/2;j++)
        {
            sub_pixbuffs[0][i][j]=pixelbuff[i*2][j*2];
            sub_chaosbuff[0][i][j]=chaosbuff[i*2][j*2];

            sub_pixbuffs[1][i][j]=pixelbuff[i*2][j*2+1];
            sub_chaosbuff[1][i][j]=chaosbuff[i*2][j*2+1];

            sub_pixbuffs[2][i][j]=pixelbuff[i*2+1][j*2];
            sub_chaosbuff[2][i][j]=chaosbuff[i*2+1][j*2];

            sub_pixbuffs[3][i][j]=pixelbuff[i*2+1][j*2+1];
            sub_chaosbuff[3][i][j]=chaosbuff[i*2+1][j*2+1];
        }
    }
}

static void recover_pixelbuff(unsigned int ** pixelbuff, unsigned int ** sub_pixbuffs[4], int img_size[2])
{
    for(int i=0;i<img_size[0]/2;i++)
    {
        for(int j=0;j<img_size[1]/2;j++)
        {
            pixelbuff[i][j]=sub_pixbuffs[0][i][j];
            pixelbuff[i][j+img_size[1]/2]=sub_pixbuffs[1][i][j];
            pixelbuff[i+img_size[0]/2][j]=sub_pixbuffs[2][i][j];
            pixelbuff[i+img_size[0]/2][j+img_size[1]/2]=sub_pixbuffs[3][i][j];
        }
    }
}

static void q_sort(double ** sub_chaosbuff,int * order,int sub_height,int left, int right)
{
    if(left==right)
        return;
    double flag=sub_chaosbuff[left/sub_height][left%sub_height];
    int left_point=left+1;int right_point=right;
    while(right_point>left_point)
    {
        while(sub_chaosbuff[left_point/sub_height][left_point%sub_height]<=flag&&(left_point<right_point))
            left_point++;
        while(sub_chaosbuff[right_point/sub_height][right_point%sub_height]>=flag&&(left_point<right_point))
            right_point--;

        if(left_point<right_point)
        {
            swap_double(&sub_chaosbuff[left_point/sub_height][left_point%sub_height],&sub_chaosbuff[right_point/sub_height][right_point%sub_height]);
            swap_int(&order[left_point],&order[right_point]);
        }
    }
    int index=left_point-1;
    swap_double(&sub_chaosbuff[index/sub_height][index%sub_height],&sub_chaosbuff[left/sub_height][left%sub_height]);
    swap_int(&order[index],&order[left]);

    q_sort(sub_chaosbuff,order,sub_height,left,index);
    q_sort(sub_chaosbuff,order,sub_height,index+1,right);
}

static void shuffle_pixel(unsigned int ** sub_pixbuff, int * order,int sub_height,int sub_widght)
{
    unsigned int ** tempbuff;
    tempbuff = new unsigned int * [sub_height];
    for(int i=0;i<sub_height;i++)
    {
        tempbuff[i]= new unsigned int [sub_height];
    }

    for(int i=0;i<sub_height*sub_height;i++)
    {
        tempbuff[i/sub_height][i%sub_height]=sub_pixbuff[order[i]/sub_height][order[i]%sub_height];
    }
    for(int i=0;i<sub_height;i++)
    {
        for(int j=0;j<sub_widght;j++)
        {
            sub_pixbuff[i][j]=tempbuff[i][j];
        }
        delete [] tempbuff[i];
    }
    delete []tempbuff;
}

static void subfig_shuffle(unsigned int ** sub_pixbuffs[4],double ** sub_chaosbuff[4],int img_size[2])
{
    int sub_height=img_size[0]/2;
    int sub_widght=img_size[1]/2;

    for(int s=0;s<4;s++)
    {
        int * order;
        order = new int[sub_height*sub_widght];
        for (int i=0;i<sub_height*sub_widght;i++)
            order[i]=i;
        q_sort(sub_chaosbuff[s],order,sub_height,0,sub_height*sub_widght-1);
        shuffle_pixel(sub_pixbuffs[s],order,sub_height,sub_widght);
        delete []order ;
    }
}

static void difference_multiply(unsigned int ** sub_pixbuffs[4],int img_size[2])
{
    int sub_height=img_size[0]/2;
    int sub_widght=img_size[1]/2;

    for(int k=0;k<4;k++)
    {
        for(int i=0;i<sub_height;i++)
        {
            for(int j=0;j<sub_widght;j++)
            {
                sub_pixbuffs[(k+1)%4][i][j]=sub_pixbuffs[(k+1)%4][i][j]+sub_pixbuffs[k][i][j]%(1<<16);
            }
        }
    }
}

static unsigned int * load_pix(unsigned int ** sub_pixbuff,int mod,int x,int y,int height,int width)
{
    if(mod==1)
        return &sub_pixbuff[x][y];
    if(mod==2)
        return &sub_pixbuff[(height*3/2-x)%height][width-1-y];
    if(mod==3)
        return &sub_pixbuff[(height*3/2+x)%height][y];
    if(mod==4)
        return &sub_pixbuff[height-1-x][width-1-y];
}

static void CA_rule_combine(unsigned char* pix_bit, int length,int divide)
{
    unsigned char * new_bit;
    new_bit= new unsigned char [length];
    memset(new_bit,0,length);
    int left;int right;

    for(int i=0;i<length;i++)
    {
        left= i-1;right=i+1;
        if(i==0)
            left=length-1;
        if(i==length-1)
            right=0;

        char temp;
        temp=pix_bit[left]*4+pix_bit[i]*2+pix_bit[right];

        if(i<divide)
        {
            if (temp==0||temp==2||temp==5||temp==6)
                new_bit[i]=1;
            else
                new_bit[i]=0;
        }
        else
        {
            if(temp==1||temp==2||temp==3||temp==4)
                new_bit[i]=1;
            else
                new_bit[i]=0;
        }
    }
    memcpy(pix_bit,new_bit,16);
    delete []new_bit;
}

static void pix_balance(unsigned char ** bit_buff, double * chaos_buf)
{
    unsigned char order[4]={0,1,2,3};
    double temp_chaos[4];

    memcpy(temp_chaos,chaos_buf,4*sizeof(double));

    int round=8+(long long)((chaos_buf[0]+chaos_buf[1]+chaos_buf[2]+chaos_buf[3])*1.0e+14)%16;

    for(int i=0;i<4;i++)
    {
        for(int j=0;j<3-i;j++)
        {
            if(temp_chaos[j]>temp_chaos[j+1])
            {
                swap_double(&temp_chaos[j],&temp_chaos[j+1]);
                swap_char(&order[j],&order[j+1]);
            }
        }
    }
    unsigned char *temp_bitbuff[4];
    for(int t=0;t<round;t++)
    {
        // save a backups
        for (int i=0;i<4;i++)
        {
            temp_bitbuff[i]= new unsigned char [16];
            memcpy(temp_bitbuff[i],bit_buff[i],16);
        }
        // iterate for next step
        for (int i=1;i<4;i++)
        {
            int divide=((long long)(chaos_buf[order[i]]*1.0e+14)%16);
            CA_rule_combine(bit_buff[order[i]],16,divide);
        }
        // propose the last one
        for(int b=0;b<16;b++)
        {
            int sum=0;
            for(int i=0;i<4;i++)
            {
                sum+=bit_buff[order[i]][b];
            }
            bit_buff[order[3]][b]=sum%2;
        }
        // copy the others
        for(int i=0;i<3;i++)
        {
            memcpy(bit_buff[order[i]],temp_bitbuff[order[i+1]],16);
        }
        //free memory
        for(int i=0;i<4;i++)
        {
            delete [] temp_bitbuff[i];
        }
    }
}

static void subfig_balance(unsigned int ** sub_pixbuffs[4],double ** chaso_x,int img_size[2])
{
    int sub_height=img_size[0]/2;
    int sub_widght=img_size[1]/2;

    double * chaos_buf;
    chaos_buf=new double [4];
    for(int i=0;i<sub_height;i++)
    {
        for(int j=0;j<sub_widght;j++)
        {
            unsigned char ** bit_buff;
            bit_buff= new unsigned char * [4];

            for(int k=0;k<4;k++)
            {
                bit_buff[k]= new unsigned char [16];
                if(i+j>0)
                {
                    if(j>0)
                        *load_pix(sub_pixbuffs[k],k+1,i,j,sub_height,sub_widght)= *load_pix(sub_pixbuffs[k],k+1,i,j,sub_height,sub_widght)^
                                                                                  *load_pix(sub_pixbuffs[k],k+1,i,j-1,sub_height,sub_widght);
                    else
                        *load_pix(sub_pixbuffs[k],k+1,i,j,sub_height,sub_widght)= *load_pix(sub_pixbuffs[k],k+1,i,j,sub_height,sub_widght)^
                                                                                  *load_pix(sub_pixbuffs[k],k+1,i-1,sub_widght-1,sub_height,sub_widght);
                }
                trans_pixs_bit(bit_buff[k],load_pix(sub_pixbuffs[k],k+1,i,j,sub_height,sub_widght));
                chaos_buf[k]=chaso_x[(k/2)*sub_height+i][(k%2)*sub_widght+j];
            }
            pix_balance(bit_buff,chaos_buf);

            for(int k=0;k<4;k++)
            {
                trans_bit_pix(bit_buff[k],load_pix(sub_pixbuffs[k],k+1,i,j,sub_height,sub_widght));
                delete [] bit_buff[k];
            }
            delete [] bit_buff;
        }
    }
    delete []chaos_buf ;
}

void permutate(dicom_img * img, double ** chaos_x)
{
    unsigned int ** pixelbuff;
    unsigned int ** sub_figs[4];
    double ** sub_chaosbuff[4];
    int img_size[2];
    img->getSize(img_size);
    pixelbuff = new unsigned int * [img_size[0]];
    for(int i=0;i< img_size[0];i++)
    {
        pixelbuff[i]= new unsigned int [img_size[1]];
    }
    img->getPixels(pixelbuff);

    for(int s=0;s<4;s++)
    {
        sub_figs[s]=new unsigned int * [img_size[0]/2];
        sub_chaosbuff[s]=new double * [img_size[0]/2];
        for(int i=0;i<img_size[0]/2;i++)
        {
            sub_figs[s][i]=new unsigned int [img_size[1]/2];
            sub_chaosbuff[s][i]=new double [img_size[0]/2];
        }
    }
    divide_into_subfigs(pixelbuff,sub_figs,chaos_x,sub_chaosbuff,img_size);

    subfig_shuffle(sub_figs,sub_chaosbuff,img_size);
    for(int i=0;i<2;i++)
        difference_multiply(sub_figs,img_size);

    recover_pixelbuff(pixelbuff,sub_figs,img_size);
    subfig_balance(sub_figs,chaos_x,img_size);

    img->writePixels(pixelbuff);

    for(int i=0;i< img_size[0];i++)
    {
        delete []pixelbuff[i];
    }
    delete []pixelbuff;
    for(int s=0;s<4;s++)
    {
        for(int i=0;i<img_size[0]/2;i++)
        {
            delete []sub_figs[s][i];
            delete []sub_chaosbuff [s][i];
        }
        delete []sub_figs[s];
    }
}

void diffusion(dicom_img * img, double ** chaos_x)
{
    unsigned int ** pixelbuff;

    int img_size[2];
    img->getSize(img_size);

    pixelbuff = new unsigned int* [img_size[0]];
    for(int i=0;i< img_size[0];i++)
    {
        pixelbuff[i]= new unsigned int[img_size[1]];
    }
    img->getPixels(pixelbuff);

    unsigned char ** bit_buff;
    bit_buff= new unsigned char * [4];
    double * chaos_buf;
    chaos_buf=new double [4];

    for(int k=0;k<4;k++)
    {
        bit_buff [k]=new unsigned char [16];
        trans_pixs_bit(bit_buff [k],&pixelbuff[(img_size[0]-1)-(k%img_size[0])][k/img_size[0]]);
        chaos_buf[k]=chaos_x[(img_size[0]-1)-(k%img_size[0])][k/img_size[0]];
    }
    pix_balance(bit_buff,chaos_buf);
    for(int k=0;k<4;k++)
    {
        trans_bit_pix(bit_buff[k],&pixelbuff[k%img_size[0]][k/img_size[0]]);
        delete [] bit_buff[k];
    }
    for(int i=4;i<img_size[0]*img_size[1];i+=4)
    {
        for(int k=0;k<4;k++)
        {
            pixelbuff[(img_size[0]-1)-((i+k)%img_size[0])][(i+k)/img_size[0]]=pixelbuff[(img_size[0]-1)-((i+k)%img_size[0])][(i+k)/img_size[0]]
                                                                            ^pixelbuff[(img_size[0]-1)-((i+k-4)%img_size[0])][(i+k-4)/img_size[0]];
            bit_buff [k]=new unsigned char [16];
            trans_pixs_bit(bit_buff [k],&pixelbuff[(img_size[0]-1)-((i+k)%img_size[0])][(i+k)/img_size[0]]);
            chaos_buf[k]=chaos_x[(img_size[0]-1)-((i+k)%img_size[0])][(i+k)/img_size[0]];
        }
        pix_balance(bit_buff,chaos_buf);
        for(int k=0;k<4;k++)
        {
            trans_bit_pix(bit_buff[k],&pixelbuff[(img_size[0]-1)-((i+k)%img_size[0])][(i+k)/img_size[0]]);
            delete [] bit_buff[k];
        }
    }
        img->writePixels(pixelbuff);

    delete [] bit_buff;
    delete [] chaos_buf;
    for(int i=0;i< img_size[0];i++)
    {
        delete []pixelbuff[i];
    }
    delete []pixelbuff;
}
