#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#define neuron 1
// #define tane 7 //乱数の種
#define tm 1.0 //時定数
#define Vreset 0.0 //初期値
#define theta 9.0 // 閾値
#define h 0.001 //刻み幅
#define tmax 1000.0 // 時間
#define I0 41.0 //　外部電流
#define tauarufa 0.1
//#define okure 0.005
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

//ここでI1を計算
double keisanI(double t, int M) {
	int i;
	double I = 0.0;
	double ti = 0.0;
	double derutaW = 0.0;
	double derutat = - 0.01;

	if(M > 0) {
		if(t >= 2.5) {
			for(i = 1; i <= M; i++) {
				ti = 2.5 * (double)i;
				I = I + I0 * (((t - ti) / tauarufa) * exp(-(t - ti) / tauarufa));
				//printf("t-ti:%lf\n",t-ti);
			}
			//printf("I:%lf\n",I);
		}
	}

	return I;
}

//ここでI2の計算
double keisanI2(double t, int N, double Itime[]) { 
	int i;
	double I = 0.0;

	if(N >= 0) {
		for(i = 1; i <= N; i++) {
			if(t >= Itime[i]) {
				I = I + I0 * (((t - Itime[i]) / tauarufa) * exp(-(t - Itime[i]) / tauarufa));
			}
		}
	}

	return I;
}

//ここでほかニューロンからの発火による入力信号を計算
double keisanHakka(double t, int N, double HAKKAtime[]) { 
	int i;
	double I = 0.0;

	if(N >= 0) {
		for(i = 1; i <= N; i++) {
			if(t >= HAKKAtime[i]) {
				I = I + I0 * (((t - HAKKAtime[i]) / tauarufa) * exp(-(t - HAKKAtime[i]) / tauarufa));
			}
		}
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

int seika(int tane, double okure) {
	int i, j, randam = 0;
	double t, I1, I2, HAKKA;

	int N1[neuron] = {1}, N2[neuron] = {1}, count[neuron] = {0};
	bool hakka1[neuron] = {0},  hakka2[neuron] = {0};
	
	double I[neuron] = {0.0}, V[neuron] = {0.0};
	double w1[neuron] = {1.0}; // シナプスの結合強度
	double Itime[neuron][max];
	double I1time[neuron] = {0.0};
	double I2time[neuron] = {0.0};
	double HAKKAtime[neuron][max] = {0.0};

	/////////////////////////////////////////////////////////////////////////////////////
	// FILE *fp,*fp1,*fp2;

	/////////////////////////////////////////////////////////////////////////////////////
	// fp = fopen("7bth.txt","w"); //Vtの関係のファイル作成
	// fp1 = fopen("Ib.txt","w");  //Iの挙動のファイル作成
	// fp2 = fopen("+20.txt","w"); //w1の挙動のファイル作成

	srand(tane);//srand関数

	// I2のスパイクが出る時間が決まった
	for(i = 0; i < neuron; i++) {
		for(j = 0; j < max; j++) {
			// ランダムの値を101で割った余り、つまり0から100のどれか
			randam = rand() % 101;

			// I2のスパイクが出る時間、ノイズを50μsにした
			Itime[i][j] = (2.5 * i) + (okure + 0.001 * (randam - 50));
		}
	}

	// 時刻を進める
	// t[ms]
	for(t = 0.000; t <= tmax + h; t += h) {
		// neuron個のニューロンを考える
		for(i = 0; i < neuron; i++) {
			hakka1[i] = false;
			hakka2[i] = false;

			// printf("t:%lf Itime %lf\n",t,Itime[i]);
			// 0.001を2500回刻んだときN1をプラス、N1はそれまでに出たスパイクの累計

			// I1について
			if((int)(t * 1000) % 2500 == 0) {
				N1[i]++;
				I1time[i] = t;
				hakka1[i] = true;
			}

			// I2のスパイク数の更新と次のスパイク時間への更新
			if(Itime[i][N2[i]] < t) {
				I2time[i] = Itime[i][N2[i]];
				N2[i]++;
				hakka2[i] = true;
			}

			// printf("M,step=%d,%d\n",N1,step);
			I1 = keisanI(t, N1[i] - 1); // I1を計算
			I2 = keisanI2(t, N2[i] - 1, Itime[i]); // I2を計算

			HAKKA = 0.0;
			for(j = 0; j < neuron - 1; j++)
			{
				HAKKA += keisanHakka(t, count[(i + j + 1) % neuron], HAKKAtime[(i + j + 1) % neuron]);
			}
			I[i] = w1[i] * I1 + w2 * I2 + 0.75 * HAKKA; // 入力電流の合計
			V[i] = runge(V[i], I[i]); // ルンゲクッタ法によるVの更新

			// I1もしくはI2のスパイクが出たさいにw1を更新する
			if(hakka1[i] || hakka2[i]) {
				// この条件のもとでだけ、w1を変化するように条件を設定
				if(fabs(I1time[i] - I2time[i]) <= okure + 0.05) {
					w1[i] = Wkousin(I1time[i], I2time[i], w1[i]);
				}
			}

			// 頭うち
			if(w1[i] >= 2.5) w1[i] = 2.5;

			//Vが閾値を超えたら50mVにいく
			if(V[i] >= theta) {
				V[i] = 50;
				HAKKAtime[i][count[i]] = t;
				count[i]++;
			}

			/////////////////////////////////////////////////////////////////////////////////////
			// fprintf(fp, "%lf %lf\n", t, V);
			// fprintf(fp1, "%lf %lf\n", t, I); // 時刻tと電流Iの電流を出力
			// fprintf(fp2, "%lf %lf\n", t, w1);

			//50mVに飛ばしたら次には0に戻す
			if(V[i] >= theta) V[i] = Vreset;
		}
	}

	printf("tane : %d\n", tane);
	for(i = 0; i < neuron; i++) {
		printf("%d : %d\n", i + 1, count[i]);
	}
	printf("----------\n");

	/////////////////////////////////////////////////////////////////////////////////////
	// fclose(fp);
	// fclose(fp1);
	// fclose(fp2);

	return 0;
}

int main(int argc, char *argv[]) {
	int tane;
	double okure;

	printf("tane : ");
	scanf("%d", &tane);
	printf("okure : ");
	scanf("%lf", &okure);
	
	printf("--start--\n");
	seika(tane, okure);

	return 0;
}

