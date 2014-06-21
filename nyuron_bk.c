#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// #define tane 7 //乱数の種
#define tm 1.0 //時定数
#define Vreset 0.0 //初期値
#define theta 9.0 // 閾値
#define h 0.001 //刻み幅
#define tmax 1000.0 // 時間
#define I0 41.0 //　外部電流
#define tauarufa 0.1
#define okure 0.005
//#define w1 2.5
#define w2 1.0
#define tauf 1.0
#define Apurasu 0.100
#define Amainasu 0.105
#define max 1000
#define taua 100

double func(double V,double I) {
	V = ((-1.0) * V + I) / tm; // V = dV/dtとする
	return V;
}

//ここでIを計算
double keisanI(double t, int M) {
	double I = 0.0;
	double ti = 0.0;
	double derutaW = 0.0;
	double derutat = - 0.01;
	int i = 1;

	if(M > 0) {
		if(t >= 2.5) {
			for(i = 1; i <= M; i++) {
				ti = 2.5*(double)i;
				I = I + I0 * (((t - ti) / tauarufa) * exp(-(t - ti) / tauarufa));
				//printf("t-ti:%lf\n",t-ti);
			}
			//printf("I:%lf\n",I);
		}
	}

	return I;
}

double keisanI2(double t, int N, double Itime[]) { //ここでI2の計算
	//ｔは現在時刻、Nはスパイク数の累計、M[]はスパイクが出る時間

	double I = 0.0;
	int i = 1;
	if(N >= 0) {
		for(i = 1; i <= N; i++) {
			if(t >= Itime[i]) {
				I = I + I0 * (((t - Itime[i]) / tauarufa) * exp(-(t - Itime[i]) / tauarufa));
			}
		}
		//printf("I:%lf\n",I);
	}
	if(I < 0.0) {
		//printf("I:%lf %lf %lf %d \n",I,t,Itime[i],N);
	}

	return I;
}

double runge(double V,double I) { 
	double k1, k2, k3, k4, x1, x2, x3;

	k1 = h * func(V, I);
	k2 = h * func(V + k1 / 2.0, I);
	k3 = h * func(V + k2 / 2.0, I);
	k4 = h * func(V + k3, I);
	V = V + 1.0 / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

	return V;
}

double Wkousin(double I1time ,double I2time, double W){
	double derutat=0.0;
	double derutaW=0.0;
	double Wgensui=0.0;

	derutat = I1time - I2time;

	Wgensui = - (derutat / taua) * (W - 1.0);

	if(derutat <= 0) {
		derutaW = Apurasu * exp(-(fabs(derutat) / tauf));
		W = W + derutaW + Wgensui;

	} else {
		derutaW = - Amainasu * exp(-(fabs(derutat) / tauf));
		W = W + derutaW + Wgensui;
	}

	//printf("%lf \n",derutaW);
	return W;
}

int seika(int tane) {
	int step = -1, step2 = -2, i = 0;
	int N = 1, N2 = 1, a = 0, j = 0, k = 0, s = 0, count = 0;
	
	double t,V,I1,I2,I =0.0;
	double w1=1.0; // シナプスの結合強度
	double Itime[max];
	double I1time=0.0;
	double I2time=0.0;
	double H_time[400]={0.0};
	double sum_f;

	FILE *my;

	/////////////////////////////////////////////////////////////////////////////////////
	// FILE *fp,*fp1,*fp2;
	
	V=0.0; 
	t=0.0;

	// my = fopen("result.csv", "w");

	/////////////////////////////////////////////////////////////////////////////////////
	// fp = fopen("7bth.txt","w"); //Vtの関係のファイル作成
	// fp1 = fopen("Ib.txt","w");  //Iの挙動のファイル作成
	// fp2 = fopen("+20.txt","w"); //w1の挙動のファイル作成

	srand(tane);//srand関数

	// I2のスパイクが出る時間が決まった
	for(i = 0; i < max; i++) {
		// ランダムの値を101で割った余り、つまり0から100のどれか
		a = rand() % 101;

		// I2のスパイクが出る時間、ノイズを50μsにした
		Itime[i] = (2.5 * i) + (okure + 0.001 * (a - 50));
		// printf("Itime %d %lf\n", i, Itime[i]);
	}

	// 時刻を進める
	for(t = 0.000; t <= tmax + h; t += h) {
		j = 0;
		k = 0;

		// printf("t:%lf Itime %lf\n",t,Itime[i]);
		// 0.001を2500回刻んだときNをプラス、Nはそれまでに出たスパイクの累計

		// I1について
		if(++step == (N * 2500)) {
		// if((int)(t * 1000) % 2500 == 0) {
			N += 1;
			I1time = t;
			j = 1;
		}

		// I2のスパイク数の更新と次のスパイク時間への更新
		if(Itime[N2] < t) {
			I2time = Itime[N2];
			//printf("t:%lf Itime %lf\n",t,Itime[i]);

			N2 += 1;
			k = 1;
		}

		// printf("M,step=%d,%d\n",N,step);
		I1 = keisanI(t, N - 1); // I1を計算
		I2 = keisanI2(t, N2 - 1, Itime); // I2を計算
		I = w1 * I1 + w2 * I2 + sum_f; // 入力電流の合計
		V = runge(V, I); // ルンゲクッタ法によるVの更新

		// I1もしくはI2のスパイクが出たさいにw1を更新する
		if(j == 1 || k == 1) {
			// この条件のもとでだけ、w1を変化するように条件を設定
			if(fabs(I1time - I2time) <= okure + 0.05) {
				w1 = Wkousin(I1time, I2time, w1);
			}
		}

		// 頭うち
		if(w1 >= 2.5) w1 = 2.5;

		//Vが閾値を超えたら50mVにいく
		if(V >= theta) {
			V = 50;
			H_time[count] = t;
			count++;
		}

		// fprintf(my, "%lf, %lf, %lf, %lf\n", t, V, I, w1);

		/////////////////////////////////////////////////////////////////////////////////////
		// fprintf(fp, "%lf %lf\n", t, V);
		// fprintf(fp1, "%lf %lf\n", t, I); // 時刻tと電流Iの電流を出力
		// fprintf(fp2, "%lf %lf\n", t, w1);

		//50mVに飛ばしたら次には0に戻す
		if(V >= theta) V = Vreset;
	}

	printf("%d\n", count);

	// fclose(my);

	/////////////////////////////////////////////////////////////////////////////////////
	// fclose(fp);
	// fclose(fp1);
	// fclose(fp2);

	return 0;
}

int main(void){
	seika(7);
}
