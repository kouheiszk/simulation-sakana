#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define	UNIT		3//学習に直接かかわるニューロン数
#define	A_PLUS		0.100//STDP学習のA+
#define A_MINAS		0.105//STDP学習のA-

#define TAU_F		1.0//STDP学習の時定数τF
#define TAU_M		1.0//ニューロンの時定数τM
#define TAU_ALPHA	0.1//スパイクの時定数τα
#define TAU_A		100//結合強度の時定数τA

#define	N_MAX		1000//学習時に発火時刻を保持できる個数
#define V_TH		9.0//ニューロンが発火する閾値
#define V_MAX		100.0//ニューロン発火時の電位
#define I			41.0//入力の電位

#define	DT			0.0010//ルンゲ=クッタ法の刻み幅
#define TIME_MAX	1000.0//学習を行う時間
#define	PERIOD		2.50//与えるスパイク列の周期
#define SEED		0//乱数の種

#define W_start		1.0//ニューロンの初期結合強度
#define W_MAX		2.5//結合強度の最大値
#define W_MIN		-0.5//結合強度の最小値
#define W_TH		1.8//発火結合強度

#define DELAY1		0.003//周期的スパイク列の遅れ
#define DELAY2		0.002//周期的スパイク列の遅れ
#define NOISE		0.050//ノイズの幅

int	N_pre_spike[UNIT];//ニューロンの発火回数
int Nr_pre_spike[UNIT];//ニューロンの発火回数
int	N_post_spike[UNIT];//ニューロンの発火回数
int Nr_post_spike[UNIT];//ニューロンの発火回数
int N_learn[UNIT];//outputの発火回数
int flag_learn[UNIT];//N_MAXを超えたときのフラグ
double t_pre_spike[UNIT][N_MAX];//ニューロンの発火時刻
double t_post_spike[UNIT][N_MAX];//ニューロンの発火時刻
double t_learn[UNIT][N_MAX];//outputの発火時刻
double V_learn[UNIT];//ニューロンの電位
double W_spike[UNIT];//ニューロン間の結合強度
double W_learn[UNIT][UNIT];//outputの結合強度

double f(double T,double V,int j);//ニューロンiの式
void STDP(double T,int j);// STDP学習

FILE *fp1;
FILE *fp2;
FILE *fp3;
FILE *fp4;
FILE *fp5;
FILE *fp6;
FILE *fp7;
FILE *fp8;
FILE *fp9;
FILE *fp10;
FILE *fp11;
FILE *fp12;
FILE *fp13;
FILE *fp14;
FILE *fp15;

void main()
{
	fp1 = fopen("stdp1_3_2_0.txt","w+");//ファイルオープン
	if(fp1 == NULL){
		exit(1);
	}

	fp2 = fopen("spike1_3_2_0.txt","w+");//ファイルオープン
	if(fp2 == NULL){
		exit(1);
	}

	fp3 = fopen("spike2_3_2_0.txt","w+");//ファイルオープン
	if(fp3 == NULL){
		exit(1);
	}

	fp4 = fopen("output1_3_2_0.txt","w+");//ファイルオープン
	if(fp4 == NULL){
		exit(1);
	}

	fp5 = fopen("raster1_3_2_0.txt","w+");//ファイルオープン
	if(fp5 == NULL){
		exit(1);
	}

	fp6 = fopen("stdp2_3_2_0.txt","w+");//ファイルオープン
	if(fp6 == NULL){
		exit(1);
	}

	fp7 = fopen("output2_3_2_0.txt","w+");//ファイルオープン
	if(fp7 == NULL){
		exit(1);
	}

	fp8 = fopen("raster2_3_2_0.txt","w+");//ファイルオープン
	if(fp8 == NULL){
		exit(1);
	}

	fp9 = fopen("stdp3_3_2_0.txt","w+");//ファイルオープン
	if(fp9 == NULL){
		exit(1);
	}

	fp10 = fopen("output3_3_2_0.txt","w+");//ファイルオープン
	if(fp10 == NULL){
		exit(1);
	}

	fp11 = fopen("raster3_3_2_0.txt","w+");//ファイルオープン
	if(fp11 == NULL){
		exit(1);
	}

	fp12 = fopen("count_3_2_0.txt","w+");//ファイルオープン
	if(fp12 == NULL){
		exit(1);
	}
	fp13 = fopen("raster'1_3_2_0.txt","w+");//ファイルオープン
	if(fp13 == NULL){
		exit(1);
	}

	fp14 = fopen("raster'2_3_2_0.txt","w+");//ファイルオープン
	if(fp14 == NULL){
		exit(1);
	}

	fp15 = fopen("raster'3_3_2_0.txt","w+");//ファイルオープン
	if(fp15 == NULL){
		exit(1);
	}

	int		i,j;
	double noise;//ノイズ
    double	Time;//時間
	double dt;//刻み幅
	double k[4];//ルンゲ=クッタ法での勾配

	Time = 0.0;
	noise = 0.0;
    dt = DT;

	for(i=0; i<UNIT; i++){//初期化
		for(j=0; j<UNIT; j++){
			W_learn[i][j] = W_start * 0.75;
		}
		N_pre_spike[i]		= 0;
		Nr_pre_spike[i]		= 0;
		N_learn[i]			= 0;
		flag_learn[i]		= 0;
		t_pre_spike[i][0]	= 0.0;
		N_post_spike[i]		= 0;
		Nr_post_spike[i]	= 0;
		t_post_spike[i][0]	= 0.0;
		t_learn[i][0]		= 0.0;
		W_spike[i]			= W_start;
		V_learn[i]			= 0.0;

	}
	
	srand(SEED);
	
	for(i=0; i<UNIT; i++){
		for(Time=dt; Time<TIME_MAX; Time+=dt){
			if(Time > PERIOD/2){//Pre側スパイク発生
				if((int)(Time/dt)%(int)(PERIOD/dt) == 0){
					t_pre_spike[i][N_pre_spike[i]] = Time;
					//printf("%9.3f\n",t_pre_spike[i][N_pre_spike[i]]);
					fprintf(fp2,"%9.3f\n",t_pre_spike[i][N_pre_spike[i]]);
					N_pre_spike[i]++;
				}
		
				if(i == 0){
					if((int)((Time-DELAY1*1)/dt)%(int)(PERIOD/dt) == 0){//Post側スパイク発生
						noise = 2.0*NOISE*rand()/(RAND_MAX+0.0) - NOISE;//ノイズ
						t_post_spike[i][N_post_spike[i]] = Time + noise;//スパイクが立つ時間
						//printf("%9.3f %9.3f\n",noise,t_post_spike[i][N_post_spike[i]]);
						fprintf(fp3,"%9.3f %9.3f\n",noise,t_post_spike[i][N_post_spike[i]]);
						N_post_spike[i]++;
					}
				}
				else{
					if((int)((Time-DELAY2*1)/dt)%(int)(PERIOD/dt) == 0){//Post側スパイク発生
						noise = 2.0*NOISE*rand()/(RAND_MAX+0.0) - NOISE;//ノイズ
						t_post_spike[i][N_post_spike[i]] = Time + noise;//スパイクが立つ時間
						//printf("%9.3f %9.3f\n",noise,t_post_spike[i][N_post_spike[i]]);
						fprintf(fp3,"%9.3f %9.3f\n",noise,t_post_spike[i][N_post_spike[i]]);
						N_post_spike[i]++;
					}
				}
			}
		}
		fprintf(fp2,"\n");
		fprintf(fp3,"\n");
	}

	Time = 0.0;

	for(i=0; i<UNIT; i++){
		if(i==0){
			printf("%9.3f %10f\r",Time, W_spike[i]);
			fprintf(fp1,"%9.3f %10f\n",Time, W_spike[i]);
			//printf("%9.3f %7.3f\n",Time,V_learn[i]);
			fprintf(fp4,"%9.3f %7.3f\n",Time,V_learn[i]);
		}
		else if(i==1){
			printf("%9.3f %10f\r",Time, W_spike[i]);
			fprintf(fp6,"%9.3f %10f\n",Time, W_spike[i]);
			//printf("%9.3f %7.3f\n",Time,V_learn[i]);
			fprintf(fp7,"%9.3f %7.3f\n",Time,V_learn[i]);
		}
		else{
			printf("%9.3f %10f\r",Time, W_spike[i]);
			fprintf(fp9,"%9.3f %10f\n",Time, W_spike[i]);
			//printf("%9.3f %7.3f\n",Time,V_learn[i]);
			fprintf(fp10,"%9.3f %7.3f\n",Time,V_learn[i]);
		}
	}
	for(Time=dt; Time<TIME_MAX; Time+=dt){
		for(i=0; i<UNIT; i++){
			if(V_learn[i] >= 0.0 && V_learn[i] < V_TH){
				k[0] = f(Time, V_learn[i],i);
				k[1] = f(Time+dt*0.5, V_learn[i]+k[0]*dt*0.5,i);
    			k[2] = f(Time+dt*0.5, V_learn[i]+k[1]*dt*0.5,i);
       			k[3] = f(Time+dt, V_learn[i]+k[2]*dt,i);
       			V_learn[i] += dt*(k[0]+2.0*k[1]+2.0*k[2]+k[3])/6.0;//outputのlif
		
				if(V_learn[i] >= V_TH && V_learn[i] < V_MAX){
					V_learn[i] = V_MAX;
					t_learn[i][N_learn[i]] = Time;
					N_learn[i]++;
				}
				if(N_learn[i]==N_MAX){
					N_learn[i] = 0;
					flag_learn[i]++;
				}
			}
			else if(V_learn[i] == V_MAX){
				V_learn[i] = 0.0;
			}

			if(Time == t_pre_spike[i][Nr_pre_spike[i]]){
				//printf("%9.3f\n",t_pre_spike[i][Nr_pre_spike[i]]);
				Nr_pre_spike[i]++;
				STDP(Time,i);
			}
		
			else if(Time < t_post_spike[i][Nr_post_spike[i]] + dt && Time > t_post_spike[i][Nr_post_spike[i]] - dt && Time != t_pre_spike[i][Nr_pre_spike[i]-1]){
				//printf("%9.3f\n",t_post_spike[i][Nr_post_spike[i]]);
				Nr_post_spike[i]++;
				STDP(Time,i);
			}

			W_spike[i] = W_spike[i] + DT / TAU_A *(W_start - W_spike[i]);//減衰を考慮

			if(i==0){
				printf("%9.3f %10f\r",Time, W_spike[i]);
				fprintf(fp1,"%9.3f %10f\n",Time, W_spike[i]);
				//printf("%9.3f %7.3f\n",Time,V_learn[i]);
				fprintf(fp4,"%9.3f %7.3f\n",Time,V_learn[i]);
			}
			else if(i==1){
				printf("%9.3f %10f\r",Time, W_spike[i]);
				fprintf(fp6,"%9.3f %10f\n",Time, W_spike[i]);
				//printf("%9.3f %7.3f\n",Time,V_learn[i]);
				fprintf(fp7,"%9.3f %7.3f\n",Time,V_learn[i]);
			}

			else{
				printf("%9.3f %10f\r",Time, W_spike[i]);
				fprintf(fp9,"%9.3f %10f\n",Time, W_spike[i]);
				//printf("%9.3f %7.3f\n",Time,V_learn[i]);
				fprintf(fp10,"%9.3f %7.3f\n",Time,V_learn[i]);
			}
			if(W_spike[i] >= W_TH){
				if(i==0){
					//printf("%9.3f %4.3f\n",Time,DELAY);
					fprintf(fp5,"%9.3f %4.3f\n",Time,DELAY1);
				}
				else if(i==1){
					//printf("%9.3f %4.3f\n",Time,DELAY);
					fprintf(fp8,"%9.3f %4.3f\n",Time,DELAY1);
				}
				else{
					//printf("%9.3f %4.3f\n",Time,DELAY);
					fprintf(fp11,"%9.3f %4.3f\n",Time,DELAY1);
				}
			}
		}
	}

	for(i=0; i<UNIT;i++){
		//printf("%d番目のニューロンの発火回数%4d回\n",i,N_learn[i]);
		fprintf(fp12,"%d番目のニューロンの発火回数%4d回\n",i + 1,N_learn[i]);
		for(j=0;j<N_learn[i];j++){
			if(i==0){
				//printf("%9.3f %d\n",t_learn[i][j],i);
				fprintf(fp13,"%9.3f %d\n",t_learn[i][j],i);
			}
			else if(i==1){
				//printf("%9.3f %d\n",t_learn[i][j],i);
				fprintf(fp14,"%9.3f %d\n",t_learn[i][j],i);
			}
			else{
				//printf("%9.3f %d\n",t_learn[i][j],i);
				fprintf(fp15,"%9.3f %d\n",t_learn[i][j],i);
			}
		}
	}

	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	fclose(fp5);
	fclose(fp6);
	fclose(fp7);
	fclose(fp8);
	fclose(fp9);
	fclose(fp10);
	fclose(fp11);
	fclose(fp12);
	fclose(fp13);
	fclose(fp14);
	fclose(fp15);

	printf("\a");//終了の合図
}

void STDP(double T, int j)//STDP学習則
{
	int i;

	if(Nr_pre_spike[j] > 0 && Nr_post_spike[j] > 0){//最初の入力で動かないように
		if(T > t_post_spike[j][Nr_post_spike[j]-1] - DT && T < t_post_spike[j][Nr_post_spike[j]-1] + DT && T != t_pre_spike[j][Nr_pre_spike[j]-1]){
			for(i=0; i<Nr_pre_spike[j]; i++){
				W_spike[j] += A_PLUS*exp(-(T-t_pre_spike[j][i])/TAU_F);
			}
		}
		if(T == t_pre_spike[j][Nr_pre_spike[j]-1]){
			for(i=0; i<Nr_post_spike[j]; i++){
				W_spike[j] -= A_MINAS*exp(-(T-t_post_spike[j][i])/TAU_F);
			}
		}
		if(W_spike[j] > W_MAX){//結合強度の上限
			W_spike[j] = W_MAX;
		}
		else if(W_spike[j] < W_MIN){//結合強度の下限
			W_spike[j] = W_MIN;
		}
	}
}

double f(double T,double V,int j)//outputの式
{
	int i,k;
	double G;
	
	G = -V;
	
	for(i=0; i<Nr_pre_spike[j]; i++){
		G += I*W_spike[j]*((T-t_pre_spike[j][i])/TAU_ALPHA)*exp(-(T-t_pre_spike[j][i])/TAU_ALPHA);
	}

	for(i=0; i<Nr_post_spike[j]; i++){
		G += I*W_start*((T-t_post_spike[j][i])/TAU_ALPHA)*exp(-(T-t_post_spike[j][i])/TAU_ALPHA);
	}

	for(k=0; k<UNIT;k++){
		if(k != j){
			if(flag_learn[k] == 0){
				for(i=0; i<N_learn[k]; i++){
					G += I*W_learn[j][k]*((T-t_learn[k][i])/TAU_ALPHA)*exp(-(T-t_learn[k][i])/TAU_ALPHA);
				}
			}
			else{
				for(i=0; i<N_MAX; i++){
					G += I*W_learn[j][k]*((T-t_learn[k][i])/TAU_ALPHA)*exp(-(T-t_learn[k][i])/TAU_ALPHA);
				}
			}
		}
	}

	G = G/TAU_M;
	return G;
}