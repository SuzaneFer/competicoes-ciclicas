#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <math.h>

#define Nx 512
#define Ny 512
#define pr 0.25
#define NG 2000
#define NF 2000/*numero de vezes do arquivo*/


const double pp[3][3]= {\
    {-1.0,1.0,0.0},\
    {0.0,-1.0,1.0},\
    {1.0,0.0,-1.0},\
};

void op(int, int *phi, float *d0,float *d1, float *d2, float *d3);
void ic( int *phi,double *d0,double *d1, double *d2, double *d3);

int main ()
{
    const gsl_rng_type * W;


    int i, j, n, ni, nj, m, teste, ativo, passivo, vizinho, counter, n_log=22;
    int l = NG/NF, aux=0, amostragem=1, n_amostra=2;
    int *phi;
    float *d0,*d1,*d2,*d3;
    double p;
    float k, pm;
    char nome[100];
    FILE *densidade1, *densidade2, *densidade3, *densidade0;


    sprintf(nome, "densidade0.dat");
    densidade0 = fopen(nome, "a");

    sprintf(nome, "densidade1.dat");
    densidade1 = fopen(nome, "a");

    sprintf(nome, "densidade2.dat");
    densidade2 = fopen(nome, "a");

    sprintf(nome, "densidade3.dat");
    densidade3 = fopen(nome, "a");



    d0 = (float *)malloc((NF)*sizeof(int));
    d1 = (double *)malloc((NF)*sizeof(int));
    d2 = (double *)malloc((NF)*sizeof(int));
    d3 = (double *)malloc((NF)*sizeof(int));
    pm=0.6;
    printf("PM = %f \n", pm);
    for(aux=0; aux<NF; aux++)
    {
        d0[aux] = 0;
        d1[aux] = 0;
        d2[aux] = 0;
        d3[aux] = 0;
    }
    while(amostragem<n_amostra)
    {

        gsl_rng *w ;
        gsl_rng_env_setup();
        W = gsl_rng_default;
        w = gsl_rng_alloc (W);
        phi = (int *)malloc((Nx * Ny)*sizeof(int));


        ic(phi,d0,d1,d2,d3);
        counter=1;


        for (n=0; n<NG; n++)
        {
            for(m=0; m<Nx*Ny; m++)
            {
                i= gsl_rng_uniform(w)*Nx;
                j= gsl_rng_uniform(w)*Ny;
                ativo= j*Nx+i;

                if(phi[ativo]!=0)
                {
                    vizinho=gsl_rng_uniform (w)*4;
                    switch(vizinho)
                    {
                    case 0:
                        passivo= j*Nx+(i+1+Nx)%Nx;
                        break;
                    case 1:
                        passivo= j*Nx+(i-1+Nx)%Nx;
                        break;
                    case 2:
                        passivo= ((j+1+Ny)%Ny)*Nx+i;
                        break;
                    case 3:
                        passivo= ((j-1+Ny)%Ny)*Nx+i;
                        break;
                    }

                    p=gsl_rng_uniform(w);

                    if(p<pm)
                    {
                        teste=phi[ativo];
                        phi[ativo]=phi[passivo];
                        phi[passivo]=teste;
                    }
                    else
                    {
                        if(p>=pm && p< (pm + pr))
                        {
                            if (phi[passivo]==0)
                            {
                                phi[passivo]=phi[ativo];
                            }
                        }
                        else
                        {
                            p=gsl_rng_uniform(w);

                            if(p<pp[phi[ativo]-1][phi[passivo]-1])
                            {
                                phi[passivo]=0;
                            }
                        }
                    }
                }
            }
            if(n>= (int)(exp(log(NG)*(n_log/100.0))))
            {
                ni = 0;
                for(nj= 0; nj<Nx*Ny; nj++)
                {
                    if(phi[nj]==0)
                    {
                        ni++;
                    }
                }

                n_log++;
            }
            if ((n%l)==0)
            {

                op(counter, phi,d0, d1, d2, d3);
                counter++;
            }
        }
        gsl_rng_free (w);
        free(phi);
        amostragem++;

    }

    for(aux=1; aux<NF; aux++)
    {
        fprintf(densidade0,"%d %.4f\n",aux,d0[aux]/(n_amostra-1));
        fprintf(densidade1,"%d %.4f\n",aux,d1[aux]/(n_amostra-1));
        fprintf(densidade2,"%d %.4f\n",aux,d2[aux]/(n_amostra-1));
        fprintf(densidade3,"%d %.4f\n",aux,d3[aux]/(n_amostra-1));
    }

return 0;
}

void op(int k, int *phi,float *d0, float *d1, float *d2, float *d3)
{
    int i, j, aux=0;

    double contd1,contd0,contd2, contd3;
    FILE *out,*cont1,*cont2,*cont3;
    char nome[100];
    contd1 = 0;
    contd0 = 0;
    contd2 = 0;
    contd3 = 0;
    sprintf(nome, "rps_saida_testematriz-%d.dat",k);
    out = fopen(nome, "w");

    sprintf(nome, "cont1-%d.dat",k);
    cont1 = fopen(nome, "w");

    sprintf(nome, "cont2-%d.dat",k);
    cont2 = fopen(nome, "w");

    sprintf(nome, "cont3-%d.dat",k);
    cont3 = fopen(nome, "w");


    for (j = 0; j<Nx; j++)
    {
        for (i=0; i<Ny; i++ )
        {
            fprintf(out, "%d ", phi[j*Nx + i] );
            if(phi[j*Nx + i]==1)
            {
                contd1 = contd1 +1;
                fprintf(cont1,"%d %d\n",i,j);
            }
            if(phi[j*Nx + i]==2)
            {
                contd2 = contd2 +1;
                fprintf(cont2,"%d %d\n",i,j);
            }
            if(phi[j*Nx + i]==3)
            {
                contd3 = contd3 +1;
                fprintf(cont3,"%d %d\n",i,j);
            }
            if(phi[j*Nx + i]==0)
            {
                contd0 = contd0 +1;
            }
        }

        fprintf(out,"\n" );
    }
    float den0=0, den1=0, den2=0, den3=0;

    den0 = contd0/(Nx*Nx);
    den1 = contd1/(Nx*Nx);
    den2 = contd2/(Nx*Nx);
    den3 = contd3/(Nx*Nx);

    d0[k] = d0[k] + den0;
    d1[k] = d1[k] + den1;
    d2[k] = d2[k] + den2;
    d3[k] = d3[k] + den3;

    fclose(out);
}

void ic(int *phi,double *d0,double *d1, double *d2, double *d3)
{
    const gsl_rng_type * W;
    gsl_rng *w ;
    gsl_rng_env_setup();
    gsl_rng_default_seed=time(0);
    W = gsl_rng_default;
    w = gsl_rng_alloc (W);

    int i, j, counter;

    for (i=0; i<Nx*Ny; i++ )
    {
        phi[i]=0;
    }
    for(i=1; i<4; i++)
    {
        counter=1;
        while(counter<Nx*Ny*0.25)
        {
            j= gsl_rng_uniform(w)*Nx*Ny;
            if(phi[j]==0)
            {
                phi[j]=i;
                counter++;
            }
        }
    }

    op(0,phi,d0,d1,d2,d3);
    gsl_rng_free (w);
}


