#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "data.h"

const int HIGHEST = 3;
int taskperthr = 1;
int sizepernode;
int ITER = 1;

// global var
float preScore = FLT_MIN;
float maxScore[HIGHEST] = { FLT_MIN };

bool orders[NODE_N][NODE_N];
bool preOrders[NODE_N][NODE_N];
bool preGraph[NODE_N][NODE_N];
bool graph[NODE_N][NODE_N];
bool bestGraph[HIGHEST][NODE_N][NODE_N];
float *localscore, *D_localscore, *D_Score, *scores;
float *LG;
bool *D_parent;
int *D_resP, *parents;
int *D_data;

void initialize();			//initial orders and data
int genOrders();		//swap
int conCore(float score);			//discard new order or not
bool getparent(int *bit, int *pre, int posN, int *parent, int *parN, int time);	//get every possible set of parents for a node
void incr(int *bit, int n);	//binary code increases 1 each time
void incrS(int *bit, int n);	//STATE_N code increases 1 each time
bool getState(int parN, int *state, int time);	//get every possible combination of state for a parent set
float logGamma(int N);		// log and gamma
float findBestGraph();
void genScore();
int convert(int *parent, int parN);
void sortGraph();
void swap(int a, int b);
void Pre_logGamma();
int findindex(int *arr, int size);
int C(int n, int a);

void genScoreKernel(int sizepernode, float *D_localscore, int *D_data, float *D_LG);
void Dincr(int *bit, int n);
void DincrS(int *bit, int n);
bool D_getState(int parN, int *sta, int time);
void D_findComb(int *comb, int l, int n);
int D_findindex(int *arr, int size);
int D_C(int n, int a);


int main()
{
	int c = 0;
	clock_t total = 0, start;

	//cudaDeviceSynchronize();
	srand(time(NULL));

	start = clock();

	initialize();
	genScore();

	total += clock() - start;

	for (int i = 0; i < ITER; i++) {
		start = clock();

		for (int a = 0; a < NODE_N; a++) {
			for (int j = 0; j < NODE_N; j++) {
				orders[a][j] = preOrders[a][j];
			}
		}

		int tmp = rand() % 6;
		for (int j = 0; j < tmp; j++)
			genOrders();

		float score = findBestGraph();
		conCore(score);

		total += clock() - start;

		//store the top HIGHEST highest orders
		if (c < HIGHEST) {
			tmp = 1;
			for (int j = 0; j < c; j++) {
				if (maxScore[j] == preScore) {
					tmp = 0;
				}
			}
			if (tmp != 0) {
				maxScore[c] = preScore;
				for (int a = 0; a < NODE_N; a++) {
					for (int b = 0; b < NODE_N; b++) {
						bestGraph[c][a][b] = preGraph[a][b];
					}
				}
				c++;
			}

		} else if (c == HIGHEST) {
			sortGraph();
			c++;
		} else {

			tmp = 1;
			for (int j = 0; j < HIGHEST; j++) {
				if (maxScore[j] == preScore) {
					tmp = 0;
					break;
				}
			}
			if (tmp != 0 && preScore > maxScore[HIGHEST - 1]) {
				maxScore[HIGHEST - 1] = preScore;
				for (int a = 0; a < NODE_N; a++) {
					for (int b = 0; b < NODE_N; b++) {
						bestGraph[HIGHEST - 1][a][b] = preGraph[a][b];
					}
				}
				int b = HIGHEST - 1;
				for (int a = HIGHEST - 2; a >= 0; a--) {
					if (maxScore[b] > maxScore[a]) {
						swap(a, b);
						int tmpd = maxScore[a];
						maxScore[a] = maxScore[b];
						maxScore[b] = tmpd;
						b = a;
					}
				}
			}
		}
	}

	free(localscore);
	//cudaFree(D_localscore);
	//cudaFree(D_parent);

	free(scores);
	free(parents);
	//cudaFree(D_Score);
	//cudaFree(D_resP);

	printf("%d,%d,%d,%f\n", STATE_N, NODE_N, DATA_N, (float) total / CLOCKS_PER_SEC);
}

void sortGraph()
{
	float max = FLT_MIN;
	int maxi;

	for (int j = 0; j < HIGHEST - 1; j++) {
		max = maxScore[j];
		maxi = j;
		for (int i = j + 1; i < HIGHEST; i++) {
			if (maxScore[i] > max) {
				max = maxScore[i];
				maxi = i;
			}
		}

		swap(j, maxi);
		float tmp = maxScore[j];
		maxScore[j] = max;
		maxScore[maxi] = tmp;
	}
}

void swap(int a, int b)
{
	for (int i = 0; i < NODE_N; i++) {
		for (int j = 0; j < NODE_N; j++) {
			bool tmp = bestGraph[a][i][j];
			bestGraph[a][i][j] = bestGraph[b][i][j];
			bestGraph[b][i][j] = tmp;
		}
	}
}

void initialize()
{
	int i, j, tmp, a, b, r;
	bool tmpd;
	tmp = 1;
	for (i = 1; i <= 4; i++) {
		tmp += C(NODE_N - 1, i);
	}
	sizepernode = tmp;
	tmp *= NODE_N;

	localscore = (float *)malloc(tmp * sizeof(float));

	for (i = 0; i < tmp; i++)
		localscore[i] = 0;

	for (i = 0; i < NODE_N; i++) {
		for (j = 0; j < NODE_N; j++)
			orders[i][j] = 0;
	}
	for (i = 0; i < NODE_N; i++) {
		for (j = 0; j < i; j++)
			orders[i][j] = 1;
	}
	r = rand() % 10000;
	for (i = 0; i < r; i++) {
		a = rand() % NODE_N;
		b = rand() % NODE_N;
		for (j = 0; j < NODE_N; j++) {
			tmpd = orders[j][a];
			orders[j][a] = orders[j][b];
			orders[j][b] = tmpd;

		}

		for (j = 0; j < NODE_N; j++) {
			tmpd = orders[a][j];
			orders[a][j] = orders[b][j];
			orders[b][j] = tmpd;
		}
	}

	for (i = 0; i < NODE_N; i++) {
		for (j = 0; j < NODE_N; j++) {
			preOrders[i][j] = orders[i][j];
		}
	}

	D_data = (int *)malloc(DATA_N * NODE_N * sizeof(int));
	memcpy(D_data, data, DATA_N * NODE_N * sizeof(int));
}

 //generate ramdom order
int genOrders()
{

	int a, b, j;
	bool tmp;
	a = rand() % NODE_N;
	b = rand() % NODE_N;

	for (j = 0; j < NODE_N; j++) {
		tmp = orders[a][j];
		orders[a][j] = orders[b][j];
		orders[b][j] = tmp;
	}
	for (j = 0; j < NODE_N; j++) {
		tmp = orders[j][a];
		orders[j][a] = orders[j][b];
		orders[j][b] = tmp;
	}

	return 1;
}

//decide leave or discard an order
int conCore(float score)
{
	int i, j;
	float tmp;
	tmp = log((rand() % 100000) / 100000.0);
	if (tmp < (score - preScore)) {

		for (i = 0; i < NODE_N; i++) {
			for (j = 0; j < NODE_N; j++) {
				preOrders[i][j] = orders[i][j];
				preGraph[i][j] = graph[i][j];
			}

		}
		preScore = score;

		return 1;
	}

	return 0;

}

void genScore()
{
	Pre_logGamma();
	memset(localscore, 0.0, NODE_N * sizepernode * sizeof(float));

	genScoreKernel(sizepernode, localscore, D_data, LG);

	free(LG);

	scores = (float *)malloc((sizepernode / (256 * taskperthr) + 1) * sizeof(float));
	parents = (int *)malloc((sizepernode / (256 * taskperthr) + 1) * 4 * sizeof(int));
	D_Score = (float *)malloc((sizepernode / (256 * taskperthr) + 1) * sizeof(float));
	D_parent = (bool *)malloc(NODE_N * sizeof(bool));
	D_resP = (int *)malloc((sizepernode / (256 * taskperthr) + 1) * 4 * sizeof(int));

}

int convert(int *parent, int parN)
{
	int i, j, w = 1, tmp = 0;
	j = 0;
	for (i = 0; parN > 0 && i <= parent[parN - 1]; i++) {
		if (parent[j] == i) {
			j++;
			tmp += w;
		}
		w *= 2;
	}

	return tmp;
}

void Pre_logGamma()
{

	LG = (float *)malloc((DATA_N + 2) * sizeof(float));

	LG[1] = log(1.0);
	float i;
	for (i = 2; i <= DATA_N + 1; i++) {
		LG[(int)i] = LG[(int)i - 1] + log((float)i);
	}

}

void incr(int *bit, int n)
{

	bit[n]++;
	if (bit[n] >= 2) {
		bit[n] = 0;
		incr(bit, n + 1);
	}

	return;
}

void incrS(int *bit, int n)
{

	bit[n]++;
	if (bit[n] >= STATE_N) {
		bit[n] = 0;
		incr(bit, n + 1);
	}

	return;
}

bool getState(int parN, int *state, int time)
{
	int j = 1;

	j = pow(STATE_N, (float)parN) - 1;

	if (time > j)
		return false;

	if (time >= 1)
		incrS(state, 0);

	return true;

}

bool getparent(int *bit, int *pre, int posN, int *parent, int *parN, int time)
{
	int i, j = 1;

	*parN = 0;
	if (time == 0)
		return true;

	for (i = 0; i < posN; i++) {
		j = j * 2;
	}
	j--;

	if (time > j)
		return false;

	incr(bit, 0);

	for (i = 0; i < posN; i++) {
		if (bit[i] == 1) {
			parent[(*parN)++] = pre[i];
		}
	}

	return true;

}

float findBestGraph()
{
	float bestls = FLT_MIN;
	int bestparent[5];
	int bestpN, total;
	int node, index;
	int pre[NODE_N] = { 0 };
	int parent[NODE_N] = { 0 };
	int posN = 0, i, j, parN, tmp, k, l;
	float ls = FLT_MIN, score = 0;
	int blocknum;

	for (i = 0; i < NODE_N; i++)
		for (j = 0; j < NODE_N; j++)
			graph[i][j] = 0;

	for (node = 0; node < NODE_N; node++) {

		bestls = FLT_MIN;
		posN = 0;

		for (i = 0; i < NODE_N; i++) {
			if (orders[node][i] == 1) {
				pre[posN++] = i;
			}
		}

		if (posN >= 4) {
			for (i = 0; i < posN; i++) {
				for (j = i + 1; j < posN; j++) {
					for (k = j + 1; k < posN; k++) {
						for (l = k + 1; l < posN; l++) {
							parN = 4;
							if (pre[i] > node)
								parent[1] = pre[i];
							else
								parent[1] = pre[i] + 1;
							if (pre[j] > node)
								parent[2] = pre[j];
							else
								parent[2] = pre[j] + 1;
							if (pre[k] > node)
								parent[3] = pre[k];
							else
								parent[3] = pre[k] + 1;
							if (pre[l] > node)
								parent[4] = pre[l];
							else
								parent[4] = pre[l] + 1;

							index = findindex(parent, parN);
							index += sizepernode * node;
							ls = localscore[index];

							if (ls > bestls) {
								bestls = ls;
								bestpN = parN;
								for (tmp = 0; tmp < parN; tmp++)
									bestparent[tmp] = parent[tmp + 1];
							}

						}
					}
				}
			}
		}

		if (posN >= 3) {
			for (i = 0; i < posN; i++) {
				for (j = i + 1; j < posN; j++) {
					for (k = j + 1; k < posN; k++) {

						parN = 3;
						if (pre[i] > node)
							parent[1] = pre[i];
						else
							parent[1] = pre[i] + 1;
						if (pre[j] > node)
							parent[2] = pre[j];
						else
							parent[2] = pre[j] + 1;
						if (pre[k] > node)
							parent[3] = pre[k];
						else
							parent[3] = pre[k] + 1;

						index = findindex(parent, parN);
						index += sizepernode * node;
						ls = localscore[index];

						if (ls > bestls) {
							bestls = ls;
							bestpN = parN;
							for (tmp = 0; tmp < parN; tmp++)
								bestparent[tmp] = parent[tmp + 1];
						}

					}
				}
			}
		}

		if (posN >= 2) {
			for (i = 0; i < posN; i++) {
				for (j = i + 1; j < posN; j++) {

					parN = 2;
					if (pre[i] > node)
						parent[1] = pre[i];
					else
						parent[1] = pre[i] + 1;
					if (pre[j] > node)
						parent[2] = pre[j];
					else
						parent[2] = pre[j] + 1;

					index = findindex(parent, parN);
					index += sizepernode * node;
					ls = localscore[index];

					if (ls > bestls) {
						bestls = ls;
						bestpN = parN;
						for (tmp = 0; tmp < parN; tmp++)
							bestparent[tmp] = parent[tmp + 1];
					}

				}
			}
		}

		if (posN >= 1) {
			for (i = 0; i < posN; i++) {

				parN = 1;
				if (pre[i] > node)
					parent[1] = pre[i];
				else
					parent[1] = pre[i] + 1;

				index = findindex(parent, parN);
				index += sizepernode * node;
				ls = localscore[index];

				if (ls > bestls) {
					bestls = ls;
					bestpN = parN;
					for (tmp = 0; tmp < parN; tmp++)
						bestparent[tmp] = parent[tmp + 1];

				}
			}
		}

		parN = 0;
		index = sizepernode * node;

		ls = localscore[index];

		if (ls > bestls) {
			bestls = ls;
			bestpN = 0;
		}

		if (bestls > FLT_MIN) {

			for (i = 0; i < bestpN; i++) {
				if (bestparent[i] < node)
					graph[node][bestparent[i] - 1] = 1;
				else
					graph[node][bestparent[i]] = 1;
			}
			score += bestls;
		}

	}

	return score;
}

int findindex(int *arr, int size)
{				//reminder: arr[0] has to be 0 && size == array size-1 && index start from 0
	int i, j, index = 0;

	for (i = 1; i < size; i++) {
		index += C(NODE_N - 1, i);
	}

	for (i = 1; i <= size - 1; i++) {
		for (j = arr[i - 1] + 1; j <= arr[i] - 1; j++) {
			index += C(NODE_N - 1 - j, size - i);
		}
	}

	index += arr[size] - arr[size - 1];

	return index;

}

int C(int n, int a)
{
	int i, res = 1, atmp = a;

	for (i = 0; i < atmp; i++) {
		res *= n;
		n--;
	}

	for (i = 0; i < atmp; i++) {
		res /= a;
		a--;
	}

	return res;
}
void genScoreKernel(int sizepernode, float *D_localscore, int *D_data, float *D_LG)
{
	for (int id = 0; id < sizepernode; id++) {

		int node, index;
		bool flag;
		int parent[5] = { 0 };
		int pre[NODE_N] = { 0 };
		int state[5] = { 0 };
		int i, j, parN = 0, tmp, t;
		int t1 = 0, t2 = 0;
		float ls = 0;
		int Nij[STATE_N] = { 0 };

		D_findComb(parent, id, NODE_N - 1);

		for (i = 0; i < 4; i++) {
			if (parent[i] > 0)
				parN++;
		}

		for (node = 0; node < NODE_N; node++) {

			j = 1;
			for (i = 0; i < NODE_N; i++) {
				if (i != node)
					pre[j++] = i;

			}

			for (tmp = 0; tmp < parN; tmp++)
				state[tmp] = 0;

			index = sizepernode * node + id;

			t = 0;
			while (D_getState(parN, state, t++)) {	// for get state

				ls = 0;
				for (tmp = 0; tmp < STATE_N; tmp++)
					Nij[tmp] = 0;

				for (t1 = 0; t1 < DATA_N; t1++) {
					flag = true;
					for (t2 = 0; t2 < parN; t2++) {
						if (D_data[t1 * NODE_N + pre[parent[t2]]] != state[t2]) {
							flag = false;
							break;
						}
					}
					if (!flag)
						continue;

					Nij[D_data[t1 * NODE_N + node]]++;

				}

				tmp = STATE_N - 1;

				for (t1 = 0; t1 < STATE_N; t1++) {
					ls += D_LG[Nij[t1]];
					tmp += Nij[t1];
				}

				ls -= D_LG[tmp];
				ls += D_LG[STATE_N - 1];

				D_localscore[index] += ls;

			}

		}

	}
}

void Dincr(int *bit, int n)
{

	while (n <= NODE_N) {
		bit[n]++;
		if (bit[n] >= 2) {
			bit[n] = 0;
			n++;
		} else {
			break;
		}
	}

	return;
}

void DincrS(int *bit, int n)
{

	bit[n]++;
	if (bit[n] >= STATE_N) {
		bit[n] = 0;
		Dincr(bit, n + 1);
	}

	return;
}

bool D_getState(int parN, int *sta, int time)
{
	int i, j = 1;

	for (i = 0; i < parN; i++) {
		j *= STATE_N;
	}
	j--;
	if (time > j)
		return false;

	if (time >= 1)
		DincrS(sta, 0);

	return true;

}

void D_findComb(int *comb, int l, int n)
{
	const int len = 4;
	if (l == 0) {
		for (int i = 0; i < len; i++)
			comb[i] = -1;
		return;
	}
	int sum = 0;
	int k = 1;

	while (sum < l)
		sum += D_C(n, k++);
	l -= sum - D_C(n, --k);
	int low = 0;
	int pos = 0;
	while (k > 1) {
		sum = 0;
		int s = 1;
		while (sum < l)
			sum += D_C(n - s++, k - 1);
		l -= sum - D_C(n - (--s), --k);
		low += s;
		comb[pos++] = low;
		n -= s;
	}
	comb[pos] = low + l;
	for (int i = pos + 1; i < 4; i++)
		comb[i] = -1;
}

int D_findindex(int *arr, int size)
{				//reminder: arr[0] has to be 0 && size == array size-1 && index start from 0
	int i, j, index = 0;

	for (i = 1; i < size; i++) {
		index += D_C(NODE_N - 1, i);
	}

	for (i = 1; i <= size - 1; i++) {
		for (j = arr[i - 1] + 1; j <= arr[i] - 1; j++) {
			index += D_C(NODE_N - 1 - j, size - i);
		}
	}

	index += arr[size] - arr[size - 1];

	return index;

}

int D_C(int n, int a)
{
	int i, res = 1, atmp = a;

	for (i = 0; i < atmp; i++) {
		res *= n;
		n--;
	}

	for (i = 0; i < atmp; i++) {
		res /= a;
		a--;
	}

	return res;
}
