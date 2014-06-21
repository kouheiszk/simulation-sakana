#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define UNIT_MAX            5 // ニューロンの最大数
#define DELAY               0.003 // PreとPostのスパイクの位相差[ms]
#define SEED_START          1 // srand関数の初期化のための数字
#define SEED_END            1 // srand関数の初期化のための数字

#define NOISE_MAX           50 // ノイズの最大時間[us]
#define T_MAX               1000.0 // 時間[ms]
#define DT                  0.001 // 刻み幅[ms]

#define V_RESET             0.0 // ニューロンの電位の初期値
#define V_THETA             9.0 // ニューロンが発火する閾値
#define V_MAX               20.0 // ニューロンが発火した時の電位

#define PRE_SPIKE_INTERVAL  2.5 // Pre側のスパイクの間隔 2.5[ms]
#define POST_SPIKE_INTERVAL 2.5 // Post側のスパイクの間隔 2.5[ms]

#define HAKKA_MAX           1000 // 発火数とスパイク数の理論的上限

#define TAU_F		        1.0 // STDP学習の減衰率τF[ms]
#define TAU_A		        100.0 // 減衰項における時定数τA[ms]
#define TAU_ALPHA           0.1 // スパイクの時定数 τα[ms]
#define TAU_M               1.0 // ニューロンの膜電位の時定数 1.0[ms]
#define W2                  1.0 // Post側スパイクの結合強度(抑制性) 1.0[一定]
#define I                   41.0 // EPSPの最大値 41.0[mV]
#define	A_PLUS		        0.100 // STDP学習のA+
#define A_MINUS		        0.105 // STDP学習のA-

#define W_START             1.0 // ニューロンの初期結合強度
#define W_MAX               1.778 // 結合強度の最大値
#define W_MIN               -0.50 // 結合強度の最小値

/**
 * α関数
 *
 * @param t 時刻
 * @param spike_time スパイクの入った時刻
 */
double alpha(double t, double spike_time) {
	return I * ((t - spike_time) / TAU_ALPHA) * exp(-(t - spike_time) / TAU_ALPHA);
}

/**
 * LIF
 *
 */

double func(double V, double W1, double t, double pre_spike_time[], double post_spike_time[]) {
	int i;

	V = -V;

	for(i = 0; pre_spike_time[i] <= t; i++) {
		V += W1 * alpha(t, pre_spike_time[i]);
	}

	for(i = 0; post_spike_time[i] <= t; i++) {
		V += W2 * alpha(t, post_spike_time[i]);
	}

	return V / TAU_M;
}

double LIF(double V, double W1, double t, double pre_spike_time[], double post_spike_time[]) {
	double k1, k2, k3, k4;

	k1 = DT * func(V,            W1, t,            pre_spike_time, post_spike_time);
	k2 = DT * func(V + k1 * 0.5, W1, t + DT * 0.5, pre_spike_time, post_spike_time);
	k3 = DT * func(V + k2 * 0.5, W1, t + DT * 0.5, pre_spike_time, post_spike_time);
	k4 = DT * func(V + k3,       W1, t + DT,       pre_spike_time, post_spike_time);

	// 先輩のプログラムのルンゲ・クッタはこんな感じ
	/*
	k1 = DT * func(V, W1, t, pre_spike_time, post_spike_time);
	k2 = DT * func(V + k1 * DT * 0.5, W1, t + DT * 0.5, pre_spike_time, post_spike_time);
	k3 = DT * func(V + k2 * DT * 0.5, W1, t + DT * 0.5, pre_spike_time, post_spike_time);
	k4 = DT * func(V + k3 * DT, W1, t + DT, pre_spike_time, post_spike_time);
	*/

	return V + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
}

/**
 * ネットワークLIF
 *
 */

double network_func(double V, double W1, double WIJ[UNIT_MAX][UNIT_MAX], double t, double pre_spike_time[], double post_spike_time[], int unit, double hakka_time[UNIT_MAX][HAKKA_MAX], int hakka_count[]) {
	int i, unit_interaction;

	V = -V;

	for(i = 0; pre_spike_time[i] <= t; i++) {
		V += W1 * alpha(t, pre_spike_time[i]);
	}

	for(i = 0; post_spike_time[i] <= t; i++) {
		V += W2 * alpha(t, post_spike_time[i]);
	}

	for(unit_interaction = 0; unit_interaction < UNIT_MAX; unit_interaction++) {
		// 発火していなければ足さない
		if(hakka_count[unit_interaction] == 0) continue;
		for(i = 0; hakka_time[unit_interaction][i] <= t && i < hakka_count[unit_interaction]; i++) {
			V += WIJ[unit][unit_interaction] * alpha(t, hakka_time[unit_interaction][i]);
		}
	}

	return V / TAU_M;
}

double NETWORK_LIF(double V, double W1, double WIJ[UNIT_MAX][UNIT_MAX], double t, double pre_spike_time[], double post_spike_time[], int unit, double hakka_time[UNIT_MAX][HAKKA_MAX], int hakka_count[]) {
	double k1, k2, k3, k4;

	k1 = DT * network_func(V,            W1, WIJ, t,            pre_spike_time, post_spike_time, unit, hakka_time, hakka_count);
	k2 = DT * network_func(V + k1 * 0.5, W1, WIJ, t + DT * 0.5, pre_spike_time, post_spike_time, unit, hakka_time, hakka_count);
	k3 = DT * network_func(V + k2 * 0.5, W1, WIJ, t + DT * 0.5, pre_spike_time, post_spike_time, unit, hakka_time, hakka_count);
	k4 = DT * network_func(V + k3,       W1, WIJ, t + DT,       pre_spike_time, post_spike_time, unit, hakka_time, hakka_count);

	return V + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
}

/**
 * STDP学習則
 *
 */
double STDP(double W, double pre_spike_time, double post_spike_time) {
	double dt = pre_spike_time - post_spike_time;
	double dw, W_new;

	if(dt <= 0) {
		dw = A_PLUS * exp(- fabs(dt) / TAU_F);
	} else {
		dw = - A_MINUS * exp(- fabs(dt) / TAU_F);
	}

	W_new = W + dw - dt * (W - W_START) / TAU_A;

	// 論文にはなかったが、先輩のプログラムにあった
	// おそらく結合強度の最大数を設定することで発火数を抑制するのが狙い
	// 先輩のプログラムは2.5だったけど、1.40から1.50の間の値にするとそれなりの結果になった
	if(W_new > W_MAX) {
		W_new = W_MAX;
	} else if(W_new < W_MIN) {
		W_new = W_MIN;
	}

	return W_new;
}

void calc(int seed) {
	int i, j, unit, unit_interaction;
	double t, noise;

	int pre_spike_count[UNIT_MAX] = {0}, post_spike_count[UNIT_MAX] = {0}; // スパイクがでた回数
	double pre_spike_time[UNIT_MAX][HAKKA_MAX] = {0.0}, post_spike_time[UNIT_MAX][HAKKA_MAX] = {0.0}; // スパイクが出る時刻
	double V[UNIT_MAX] = {0.0}; // ニューロンの電位
	double W1[UNIT_MAX] = {W_START}; // Pre側の結合強度
	double WIJ[UNIT_MAX][UNIT_MAX] = {W_START}; // ニューロン間の結合強度
	int hakka_count[UNIT_MAX] = {0}; // ニューロンの発火回数
	double hakka_time[UNIT_MAX][HAKKA_MAX] = {0}; // 発火した時刻

	// srand関数の初期化
	srand(seed);

	// スパイク時刻の初期化
	for(unit = 0; unit < UNIT_MAX; unit++) {
		// Pre側のスパイクが出るタイミングを設定する
		// Pre側のスパイクが出るのは、2.5[ms] に1回、つまり 2500[us] に1回
		for(i = 0; i <= (int) T_MAX / PRE_SPIKE_INTERVAL; i++) {
			t = i * PRE_SPIKE_INTERVAL;
			// Post側のスパイクが出る時刻がマイナスになってしまう可能性を考慮する
			if(t < POST_SPIKE_INTERVAL + DELAY - NOISE_MAX * 0.001) continue;

			pre_spike_time[unit][pre_spike_count[unit]] = t;
			pre_spike_count[unit]++;
		}

		// ちゃんと時刻が初期化できているか確認するためのもの
		// for(i = 0; i < pre_spike_count[unit]; i++) {
		// 	printf("%d %lf\n", i, pre_spike_time[unit][i]);
		// }

		// Post側のスパイクが出るタイミングを設定する
		// Post側のスパイクが出るのは、2.5[ms] に1回、つまり 2500[us] に1回
		// またPost側は位相差 DELAY[ms] を考慮し、さらにランダムにスパイクを出す時刻を変化させる
		for(i = 0; i <= (int) T_MAX / POST_SPIKE_INTERVAL; i++) {
			t = i * POST_SPIKE_INTERVAL;
			// Post側のスパイクが出る時刻がマイナスになってしまう可能性を考慮する
			if(t < POST_SPIKE_INTERVAL + DELAY - NOISE_MAX * 0.001) continue;

			// ノイズは -50[us] から 50[us] までの間
			noise = (NOISE_MAX - (rand() % (2 * NOISE_MAX + 1))) * 0.001;
			post_spike_time[unit][post_spike_count[unit]] = t + DELAY + noise;
			post_spike_count[unit]++;
		}

		// ちゃんと時刻が初期化できているか確認するためのもの
		// for(i = 0; i < post_spike_count[unit]; i++) {
		// 	printf("%d %lf\n", i, post_spike_time[unit][i]);
		// }

		pre_spike_count[unit] = 0;
		post_spike_count[unit] = 0;
	}

	// ニューロンネットワークによるニューロンの発火をしらべる
	for(t = 0; t < T_MAX; t += DT) {
		for(unit = 0; unit < UNIT_MAX; unit++) {
			// Pre側、Post側の結合強度を更新する
			// 更新のタイミングは、Pre側、Post側のどちらかまたは両方のスパイクが出たとき
			if(t >= pre_spike_time[unit][pre_spike_count[unit]] || t >= post_spike_time[unit][post_spike_count[unit]]) {
				// STDPの学習則に則って学習を行う
				W1[unit] = STDP(W1[unit], pre_spike_time[unit][pre_spike_count[unit]], post_spike_time[unit][post_spike_count[unit]]);

				// Pre側が発火したら、Pre側の発火数をカウントアップする
				if(t >= pre_spike_time[unit][pre_spike_count[unit]]) 
					pre_spike_count[unit]++;

				// Post側が発火したら、Post側の発火数をカウントアップする
				if(t >= post_spike_time[unit][post_spike_count[unit]]) 
					post_spike_count[unit]++;
			}

			// ニューロンの結合強度を更新する
			// ニューロンが複数個の時だけ STDP を計算する
			if(UNIT_MAX > 1) {
				for(unit_interaction = 0; unit_interaction < UNIT_MAX; unit_interaction++) {
					// 自分自身との結合強度は更新しない
					// 自分のニューロンもしくは、結合先のニューロンが発火していなければ更新しない
					if(unit_interaction == unit || hakka_count[unit] == 0 || hakka_count[unit_interaction] == 0) continue;

					// 更新のタイミングは、自分、もしくは結合しているニューロンが発火した場合
					// STDPの δt はfor文の実行間隔であるため、DT[ms] を用いる
					// 
					if(t >= hakka_time[unit][hakka_count[unit] - 1] || t >= hakka_time[unit_interaction][hakka_count[unit_interaction] - 1]) {
						WIJ[unit][unit_interaction] = STDP(WIJ[unit][unit_interaction], hakka_time[unit_interaction][hakka_count[unit_interaction] - 1], hakka_time[unit][hakka_count[unit] - 1]);
					}
				}
			}

			// ニューロンの電位を更新する
			// ニューロンの電位が閾値未満の場合、ニューロンの電位を更新する
			// ニューロンの電位が閾値以上のときは発火し終わったあとなのでもとに戻す
			if(V[unit] < V_THETA) {
				if(UNIT_MAX > 1) {
					V[unit] = NETWORK_LIF(V[unit], W1[unit], WIJ, t, pre_spike_time[unit], post_spike_time[unit], unit, hakka_time, hakka_count);
				} else {
					V[unit] = LIF(V[unit], W1[unit], t, pre_spike_time[unit], post_spike_time[unit]);
				}

				// ルンゲタックで更新後、ニューロンの電位が閾値を超えていたら、
				if(V[unit] >= V_THETA) {
					V[unit] = V_MAX;
					// 発火時刻を記録する
					hakka_time[unit][hakka_count[unit]] = t;
					hakka_count[unit]++;

					// ニューロンの発火時刻と発火回数を表示する
					// printf("Time: %9.3lf, Unit: %2d, %3d times\n", t, unit + 1, hakka_count[unit]);
				}
			} else if(V[unit] >= V_THETA) {
				V[unit] = V_RESET;
			}
		}
	}

	// 結果の出力（位相差と乱数の種）
	printf("Unit: %2d, Delay: %5.3lf, Seed: %2d\n", UNIT_MAX, DELAY, seed);
	// 結果の出力（発火回数）
	for(unit = 0; unit < UNIT_MAX; unit++) {
		printf("Unit%2d: %4d times\n", unit + 1, hakka_count[unit]);
	}
	printf("----------\n");
}

int main(void) {
	int i;
	for(i = SEED_START; i <= SEED_END; i++) calc(i);
	return EXIT_SUCCESS;
}
