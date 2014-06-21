#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define	UNIT		2//学習に直接かかわるニューロン数
#define	A_PLUS		0.100//STDP学習のA+
#define A_MINAS		0.105//STDP学習のA-

#define TAU_F		1.0//STDP学習の減衰率τF
#define TAU_M		1.0//ニューロンの減衰率τM
#define TAU_ALPHA	0.5//スパイクの減衰率τα
#define TAU_A		100//結合強度の時定数τA

#define	N_MAX		1000//学習時に発火時刻を保持できる個数
#define V_TH		9.0//ニューロンが発火する閾値
#define V_MAX		100.0//ニューロン発火時の電位
#define I			41.0//入力の電位

#define	DT			0.0010//ルンゲ=クッタ法の刻み幅
#define TIME_MAX	1000.0//学習を行う時間
#define	PERIOD		2.50//与えるスパイク列の周期
#define SEED		0//乱数の種

#define W_SPIKE		1.0//ニューロンの初期結合強度
#define W_MAX		2.5//結合強度の最大値
#define W_MIN		-0.5//結合強度の最小値

#define DELAY		0.000//周期的スパイク列の遅れ
#define NOISE		0.050//ノイズの幅

int	N_spike[UNIT];//ニューロンの発火回数
int Nr_spike[UNIT];//ニューロンの発火回数
int N_learn;//outputの発火回数
double t_spike[UNIT][N_MAX];//ニューロンの発火時刻
double V_spike;//ニューロンの電位
double W_spike;//ニューロン間の結合強度

double f(double T,double V);//ニューロンiの式
void STDP(double T);// STDP学習

FILE *fp1;
FILE *fp2;
FILE *fp3;
FILE *fp4;
FILE *fp12;

void main()
{
	fp1 = fopen("stdp_0.txt","w+");//ファイルオープン
	if(fp1 == NULL){
		exit(1);
	}

	fp2 = fopen("spike1_0.txt","w+");//ファイルオープン
	if(fp2 == NULL){
		exit(1);
	}

	fp3 = fopen("spike2_0.txt","w+");//ファイルオープン
	if(fp3 == NULL){
		exit(1);
	}

	fp4 = fopen("output_0.txt","w+");//ファイルオープン
	if(fp4 == NULL){
		exit(1);
	}

	fp12 = fopen("count_0.txt","w+");//ファイルオープン
	if(fp12 == NULL){
		exit(1);
	}

	int		i;
	double noise;//ノイズ
    double	Time;//時間
	double dt;//刻み幅
	double k[4];//ルンゲ=クッタ法での勾配

	srand(SEED);

	Time = 0.0;
	noise = 0.0;
    dt = DT;
	W_spike	= W_SPIKE;
	V_spike	= 0.0;
	N_learn = 0;

	for(i=0; i<UNIT; i++){//初期化
		N_spike[i]		= 0;
		Nr_spike[i]		= 0;
		t_spike[i][0]	= 0.0;
	}

	printf("%9.3f %10f\r",Time, W_spike);
	fprintf(fp1,"%9.3f %10f\n",Time, W_spike);
	
    for(Time=dt; Time<TIME_MAX; Time+=dt){
		if(Time > PERIOD/2){//Pre側スパイク発生
			if((int)(Time/dt)%(int)(PERIOD/dt) == 0){
				t_spike[0][N_spike[0]] = Time;
				//printf("%9.3f\n",t_spike[0][N_spike[0]]);
				fprintf(fp2,"%9.3f\n",t_spike[0][N_spike[0]]);
				N_spike[0]++;
			}
		
			if((int)((Time-DELAY*1)/dt)%(int)(PERIOD/dt) == 0){//Post側スパイク発生
				noise = 2.0*NOISE*rand()/(RAND_MAX+0.0) - NOISE;//ノイズ
				t_spike[1][N_spike[1]] = Time + noise;//スパイクが立つ時間
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

		W_spike = W_spike + DT / TAU_A *(W_SPIKE - W_spike);//減衰を考慮

		printf("%9.3f %10f\r",Time, W_spike);
		fprintf(fp1,"%9.3f %10f\n",Time, W_spike);
		//printf("%9.3f %7.3f\n",Time,V_spike);
		fprintf(fp4,"%9.3f %7.3f\n",Time,V_spike);
	}

	fprintf(fp12,"ニューロンの発火回数%4d回\n",N_learn);

	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	fclose(fp12);
}

void STDP(double T)//STDP学習則
{
	int i;

	if(Nr_spike[0] > 0 && Nr_spike[1] > 0){//最初の入力で動かないように
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

double f(double T,double V)//outputの式
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