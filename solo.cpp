#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define	UNIT		2//�w�K�ɒ��ڂ������j���[������
#define	A_PLUS		0.100//STDP�w�K��A+
#define A_MINAS		0.105//STDP�w�K��A-

#define TAU_F		1.0//STDP�w�K�̌�������F
#define TAU_M		1.0//�j���[�����̌�������M
#define TAU_ALPHA	0.5//�X�p�C�N�̌������у�
#define TAU_A		100//�������x�̎��萔��A

#define	N_MAX		1000//�w�K���ɔ��Ύ�����ێ��ł����
#define V_TH		9.0//�j���[���������΂���臒l
#define V_MAX		100.0//�j���[�������Ύ��̓d��
#define I			41.0//���͂̓d��

#define	DT			0.0010//�����Q=�N�b�^�@�̍��ݕ�
#define TIME_MAX	1000.0//�w�K���s������
#define	PERIOD		2.50//�^����X�p�C�N��̎���
#define SEED		0//�����̎�

#define W_SPIKE		1.0//�j���[�����̏����������x
#define W_MAX		2.5//�������x�̍ő�l
#define W_MIN		-0.5//�������x�̍ŏ��l

#define DELAY		0.000//�����I�X�p�C�N��̒x��
#define NOISE		0.050//�m�C�Y�̕�

int	N_spike[UNIT];//�j���[�����̔��Ή�
int Nr_spike[UNIT];//�j���[�����̔��Ή�
int N_learn;//output�̔��Ή�
double t_spike[UNIT][N_MAX];//�j���[�����̔��Ύ���
double V_spike;//�j���[�����̓d��
double W_spike;//�j���[�����Ԃ̌������x

double f(double T,double V);//�j���[����i�̎�
void STDP(double T);// STDP�w�K

FILE *fp1;
FILE *fp2;
FILE *fp3;
FILE *fp4;
FILE *fp12;

void main()
{
	fp1 = fopen("stdp_0.txt","w+");//�t�@�C���I�[�v��
	if(fp1 == NULL){
		exit(1);
	}

	fp2 = fopen("spike1_0.txt","w+");//�t�@�C���I�[�v��
	if(fp2 == NULL){
		exit(1);
	}

	fp3 = fopen("spike2_0.txt","w+");//�t�@�C���I�[�v��
	if(fp3 == NULL){
		exit(1);
	}

	fp4 = fopen("output_0.txt","w+");//�t�@�C���I�[�v��
	if(fp4 == NULL){
		exit(1);
	}

	fp12 = fopen("count_0.txt","w+");//�t�@�C���I�[�v��
	if(fp12 == NULL){
		exit(1);
	}

	int		i;
	double noise;//�m�C�Y
    double	Time;//����
	double dt;//���ݕ�
	double k[4];//�����Q=�N�b�^�@�ł̌��z

	srand(SEED);

	Time = 0.0;
	noise = 0.0;
    dt = DT;
	W_spike	= W_SPIKE;
	V_spike	= 0.0;
	N_learn = 0;

	for(i=0; i<UNIT; i++){//������
		N_spike[i]		= 0;
		Nr_spike[i]		= 0;
		t_spike[i][0]	= 0.0;
	}

	printf("%9.3f %10f\r",Time, W_spike);
	fprintf(fp1,"%9.3f %10f\n",Time, W_spike);
	
    for(Time=dt; Time<TIME_MAX; Time+=dt){
		if(Time > PERIOD/2){//Pre���X�p�C�N����
			if((int)(Time/dt)%(int)(PERIOD/dt) == 0){
				t_spike[0][N_spike[0]] = Time;
				//printf("%9.3f\n",t_spike[0][N_spike[0]]);
				fprintf(fp2,"%9.3f\n",t_spike[0][N_spike[0]]);
				N_spike[0]++;
			}
		
			if((int)((Time-DELAY*1)/dt)%(int)(PERIOD/dt) == 0){//Post���X�p�C�N����
				noise = 2.0*NOISE*rand()/(RAND_MAX+0.0) - NOISE;//�m�C�Y
				t_spike[1][N_spike[1]] = Time + noise;//�X�p�C�N��������
				//printf("%9.3f %9.3f\n",noise,t_spike[1][N_spike[1]]);
				fprintf(fp3,"%9.3f %9.3f\n",noise,t_spike[1][N_spike[1]]);
				N_spike[1]++;
			}
		}
	}

	Time = 0.0;
	//printf("%9.3f %7.3f\n",Time,V_spike);
	fprintf(fp4,"%9.3f %7.3f\n",Time,V_spike);
	for(Time=dt; Time<TIME_MAX; Time+=dt){
		if(V_spike >= 0.0 && V_spike < V_TH){
			k[0] = f(Time, V_spike);
	        k[1] = f(Time+dt*0.5, V_spike+k[0]*dt*0.5);
    	    k[2] = f(Time+dt*0.5, V_spike+k[1]*dt*0.5);
       		k[3] = f(Time+dt, V_spike+k[2]*dt);
       		V_spike += dt*(k[0]+2.0*k[1]+2.0*k[2]+k[3])/6.0;
		
			if(V_spike >= V_TH && V_spike < V_MAX){
				V_spike = V_MAX;
				N_learn++;
			}
		}
		else if(V_spike == V_MAX){
			V_spike = 0.0;
		}

		if(Time == t_spike[0][Nr_spike[0]]){
			//printf("%9.3f\n",t_spike[0][Nr_spike[0]]);
			Nr_spike[0]++;
			STDP(Time);
		}
		
		else if(Time < t_spike[1][Nr_spike[1]] + dt && Time > t_spike[1][Nr_spike[1]] - dt && Time != t_spike[0][Nr_spike[0]-1]){
			//printf("%9.3f\n",t_spike[1][Nr_spike[1]]);
			Nr_spike[1]++;
			STDP(Time);
		}

		W_spike = W_spike + DT / TAU_A *(W_SPIKE - W_spike);//�������l��

		printf("%9.3f %10f\r",Time, W_spike);
		fprintf(fp1,"%9.3f %10f\n",Time, W_spike);
		//printf("%9.3f %7.3f\n",Time,V_spike);
		fprintf(fp4,"%9.3f %7.3f\n",Time,V_spike);
	}

	fprintf(fp12,"�j���[�����̔��Ή�%4d��\n",N_learn);

	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	fclose(fp12);
}

void STDP(double T)//STDP�w�K��
{
	int i;

	if(Nr_spike[0] > 0 && Nr_spike[1] > 0){//�ŏ��̓��͂œ����Ȃ��悤��
		if(T > t_spike[1][Nr_spike[1]-1] - DT && T < t_spike[1][Nr_spike[1]-1] + DT && T != t_spike[0][Nr_spike[0]-1]){
			for(i=0; i<Nr_spike[0]; i++){
				W_spike += A_PLUS*exp(-(T-t_spike[0][i])/TAU_F);
			}
		}
		if(T == t_spike[0][Nr_spike[0]-1]){
			for(i=0; i<Nr_spike[1]; i++){
				W_spike -= A_MINAS*exp(-(T-t_spike[1][i])/TAU_F);
			}
		}
		if(W_spike > W_MAX){
			W_spike = W_MAX;
		}
		else if(W_spike < W_MIN){
			W_spike = W_MIN;
		}
	}
}

double f(double T,double V)//output�̎�
{
	int i;
	double G;
	
	G = -V;
	for(i=0; i<Nr_spike[0]; i++){
		G += I*W_spike*((T-t_spike[0][i])/TAU_ALPHA)*exp(-(T-t_spike[0][i])/TAU_ALPHA);
	}

	for(i=0; i<Nr_spike[1]; i++){
		G += I*W_SPIKE*((T-t_spike[1][i])/TAU_ALPHA)*exp(-(T-t_spike[1][i])/TAU_ALPHA);
	}

	G = G/TAU_M;
	return G;
}