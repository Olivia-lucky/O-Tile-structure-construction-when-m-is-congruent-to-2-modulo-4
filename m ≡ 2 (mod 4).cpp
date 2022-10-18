#include "permuandcombi2.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>	
#include <new>
#include <set>
#include <iterator>

using namespace std;

// find the greatest common divisor of x1 and y1
int gcd(int x1, int y1)
{
	int z = y1;
	while (x1 % y1 != 0)
	{
		z = x1 % y1;
		x1 = y1;
		y1 = z;
	}
	return z;
}

void exgcd(int i_p, int i_m, int& i_v, int& x, int& y)
{
	if (!i_m) {
		i_v = i_p;
		x = 1;
		y = 0;
	}
	else {
		exgcd(i_m, i_p % i_m, i_v, y, x);
		y -= x * (i_p / i_m);
	}
}

// find modular multiplicative inverse of i_p under modulo i_m 
int inv(int i_p, int i_m)
{
	int i_v, x, y;
	exgcd(i_p, i_m, i_v, x, y);
	return i_v == 1 ? (x + i_m) % i_m : -1;
}

// romove duplicate elements from the vector v
vector<int> unique_element_in_vector(vector<int> v) {
	vector<int>::iterator vector_iterator;
	sort(v.begin(), v.end());
	vector_iterator = unique(v.begin(), v.end());
	if (vector_iterator != v.end()) {
		v.erase(vector_iterator, v.end());
	}
	return v;
}

// find the intersection of two vectors v1 and v2
vector<int> vectors_intersection(vector<int> v1, vector<int> v2) {
	vector<int> v;
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());
	set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v));
	return v;
}

// find the union of two vectors v1 and v2
vector<int> vectors_set_union(vector<int> v1, vector<int> v2) {
	vector<int> v;
	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());
	set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), back_inserter(v));
	return v;
}

// determine whether element v_e exists in vector v
bool is_element_in_vector(vector<int> v, int element) {
	vector<int>::iterator it;
	it = find(v.begin(), v.end(), element);
	if (it != v.end()) {
		return true;
	}
	else {
		return false;
	}
}

// modulo operation
inline int mod(int m_x, int m_m) {
	return (m_m + m_x % m_m) % m_m;
}

int m, s_max, grp, p, temp, k, number;

int** Tile;
int* selected;

int ept[10], none[10], c_s[4], combo[6][2], index[1000], MM[1000][1000], NN[1000][1000];

vector<vector<int> > M(1000), N(1000);
vector<int> Trans, Nums, SumM, SumN;

bool checkArr(int* a, int* l);
bool checkLen(int* l1, int b, int c, int d);
bool checkTileStructure(int* s_a1, int* s_l2, int s_i3, ofstream& s_fout);
int checkTilesOverlap();
bool checkSpecialRectangle();
void writeO_Tile(int* a1, int* l2, int i3, ofstream& fout);

int main() {
	int count1, count2, count3, count4, gap;

	cout << "This program will try to help you construct a set of O-tile strucutres, where each O-tile structure has m*m/4-m pieces of 2 * 2 tiles, m pieces of 1 * 2 tiles and m pieces of 2 * 1 tiles for m ¡Ô2 (mod 4) and 6 <= m <= 38." << endl;
	cout << "(1) m*m/4-m pieces of 2 * 2 tiles are divided into two cases: " << endl;
	cout << "case 1: m/2 pieces of 2 * 2 tiles {j, j+m/2} * {p*j, p*j+m/2}, 0 <= j <= m/2-1." << endl;
	cout << "case 2: (m*m-6*m)/4 pieces of 2 * 2 tiles M_{i,j} * N_{i,j}, where M_{i,j} = {j, j+l_{i}}, N_{i,j} = {p*j+a_{i}, p*j+a_{i}+l_{i}}, 0 <= i <= (m-10)/4, 0 <= j <= m-1." << endl;
	cout << "(2) m pieces of 1 * 2 tiles {j} * N_{j}, where N_{j} = {p*j+x, p*j+y}, 0 <= j <= m-1." << endl;
	cout << "(3) m pieces of 2 * 1 tiles M_{j} * {j}, where M_{j} = {p^{-1}*(j-z), p^{-1}*(j-w)}, 0 <= j <= m-1." << endl;
	cout << "And l_{i}, a_{i}, x, y, z and w correspond to len_{i}, arr_{i}, combo[i][0], combo[i][1], combo[5-i][0] and combo[5-i][1] in this program respectively." << endl;

	cout << "\nPlease enter m: ";
	cin >> m;
	s_max = m * m / 4 + m;
	grp = (m - 6) / 4;

	// create a file to store the constructed O-tile structrues
	stringstream sfile;
	sfile << "m" << m << "_s" << s_max << ".txt";
	string strfile = sfile.str();
	ofstream out(strfile.c_str());

	Tile = new int* [m];
	for (int i = 0; i < m; i++) {
		Tile[i] = new int[m];
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			Tile[i][j] = 0;
		}
	}

	selected = new int[m - 4];
	for (int i = 0; i < m - 4; i++) {
		selected[i] = 0;
	}

	selected[0] = 0;
	selected[1] = m / 2;

	if (m == 6) {  // the case of m = 6 is different from other values, since 6 * 6 grids only contain 2 * 2 tiles of the case 1
		p = 5;  // p = 1 is equivalent to p = 5
		
		for (int i = 0; i < m / 2; i++) {  // 2 * 2 tiles
			Tile[i][(p * i) % m] = Tile[i][(p * i + m / 2) % m] = Tile[i + m / 2][(p * i) % m] = Tile[i + m / 2][(p * i + m / 2) % m] = i + 1;
		}

		count3 = 0;

		for (int i = 0; i < m; i++) {
			count2 = 0;

			for (int j = 0; j < m - 4; j++) {
				if (i == selected[j]) {
					count2++;
					break;
				}
			}
			if (count2 == 0) {
				c_s[count3++] = i;
			}
		}

		count4 = 0;

		for (int i = 0; i < 3; i++) {
			for (int j = i + 1; j < 4; j++) {
				combo[count4][0] = c_s[i];
				combo[count4][1] = c_s[j];
				count4++;
			}
		}

		for (int i = 0; i < 6; i++) {
			gap = (inv(p, m) * (combo[5 - i][1] - combo[5 - i][0])) % m;
			if ((combo[i][1] - combo[i][0]) != m / 2 && gap != m / 2) {
				for (int j = 0; j < m; j++) {  // 1 * 2 tiles 
					Tile[j][(p * j + combo[i][0]) % m] = Tile[j][(p * j + combo[i][1]) % m] = (j + 1) + m * m / 4 - m;
				}

				for (int j = 0; j < m; j++) {  // 2 * 1 tiles 
					Tile[mod((inv(p, m) * (j - combo[5 - i][0])), m)][j] = (j + 1) + m * m / 4;
					Tile[mod((inv(p, m) * (j - combo[5 - i][1])), m)][j] = (j + 1) + m * m / 4;
				}

				// check whether the constructed tile structure is an O-tile structure
				if (!checkTileStructure(ept, none, i, out)) return 0;

				for (int ii = 0; ii < m; ii++) {
					for (int jj = 0; jj < 4; jj++) {
						Tile[ii][(p * ii + c_s[jj]) % m] = 0;
					}
				}
			}
		}
	}
	else {
		for (p = 3; p < m - 2; p++) {  // p != 1, since when p = 1, 2 * 2 tiles will overlap, the same as p = m - 1
			if (gcd(p, m) == 1) {
				for (int i = 0; i < m / 2; i++) {  // 2 * 2 tiles(case 1)
					Tile[i][(p * i) % m] = Tile[i][(p * i + m / 2) % m] = Tile[i + m / 2][(p * i) % m] = Tile[i + m / 2][(p * i + m / 2) % m] = (i + 1) + grp * m;
				}
				permuandcombi2 combi(1, m / 2, grp);
				do {  // iteratively pick numbers len_0 < len_1 < ... < len_{grp-1} from the set {1, 2, ... ,m / 2 - 1}
					permuandcombi2 permu(1, m, grp);
					do {  // iteratively pick distinct numbers {arr_0, arr_1, ... , arr_{grp-1}} from the set {1, ... ,m - 1} \ {m / 2}
						if (checkArr(permu.arr, combi.len)) { 
						count1 = 2;
						for (int i = 0; i < grp; i++) {  // 2 * 2 tiles(case 2)
							for (int j = 0; j < m; j++) {
								Tile[j][(p * j + permu.arr[i]) % m] = (j + 1) + i * m;
								Tile[j][(p * j + permu.arr[i] + combi.len[i]) % m] = (j + 1) + i * m;
								Tile[(j + combi.len[i]) % m][(p * j + permu.arr[i]) % m] = (j + 1) + i * m;
								Tile[(j + combi.len[i]) % m][(p * j + permu.arr[i] + combi.len[i]) % m] = (j + 1) + i * m;
							}
						}

						for (int i = 0; i < grp; i++) {
							selected[count1++] = permu.arr[i];
							selected[count1++] = (permu.arr[i] + combi.len[i]) % m;
							selected[count1++] = mod((permu.arr[i] - p * combi.len[i]), m);
							selected[count1++] = mod((permu.arr[i] - (p - 1) * combi.len[i]), m);
						}

						count3 = 0;

						for (int i = 0; i < m; i++) {
							count2 = 0;

							for (int j = 0; j < m - 4; j++) {
								if (i == selected[j]) {
									count2++;
									break;
								}
							}
							if (count2 == 0) {
								c_s[count3++] = i;
							}
						}

						count4 = 0;

						for (int i = 0; i < 3; i++) {
							for (int j = i + 1; j < 4; j++) {
								combo[count4][0] = c_s[i];
								combo[count4][1] = c_s[j];
								count4++;
							}
						}

						for (int i = 0; i < 6; i++) {
							gap = (inv(p, m) * (combo[5 - i][1] - combo[5 - i][0])) % m;
							if (checkLen(combi.len, combo[i][0], combo[i][1], gap)) {
								for (int j = 0; j < m; j++) {  // 1 * 2 tiles 
									Tile[j][(p * j + combo[i][0]) % m] = Tile[j][(p * j + combo[i][1]) % m] = (j + 1) + m * m / 4 - m;
								}

								for (int j = 0; j < m; j++) {  // 2 * 1 tiles 
									Tile[mod((inv(p, m) * (j - combo[5 - i][0])), m)][j] = (j + 1) + m * m / 4;
									Tile[mod((inv(p, m) * (j - combo[5 - i][1])), m)][j] = (j + 1) + m * m / 4;
								}

								// check whether the constructed tile structure is an O-tile structure
								if (!checkTileStructure(permu.arr, combi.len, i, out)) return 0;

								for (int ii = 0; ii < m; ii++) {
									for (int jj = 0; jj < 4; jj++) {
										Tile[ii][(p * ii + c_s[jj]) % m] = 0;
									}
								}
							}
						}

						for (int i = 0; i < m / 2; i++) {
							for (int j = 0; j < m; j++) {
								if (j != (p * i) % m && j != (p * i + m / 2) % m) {
									Tile[i][j] = 0;
									Tile[i + m / 2][j] = 0;
								}
							}
						}
					}
				} while (permu.next());
			} while (combi.next_len());

			for (int i = 0; i < m; i++) {
				for (int j = 0; j < m; j++) {
					Tile[i][j] = 0;
				}
			}
		}
	}
}

delete[]selected;

for (int i = 0; i < m; i++) {
	delete[]Tile[i];
}
delete[]Tile;

cout << "\nDone constructing! A total of " << number << " O-tile structures have been constructed." << endl;
cout << "\nOutput written to file \"m" << m << "_s" << s_max << " .txt\"." << endl;

return 0;
}

/*
   This function makes sure the m - 4 values in S(of the paper) are distinct modulo m.
*/
bool checkArr(int* a, int* l) {
	if (grp == 1) {
		if (a[0] == mod((a[0] - p * l[0]), m) || a[0] == mod((a[0] - (p - 1) * l[0]), m) || (a[0] + l[0]) % m == mod((a[0] - p * l[0]), m) || (a[0] + l[0]) % m == mod((a[0] - (p - 1) * l[0]), m)) {
			return false;
		}
		if (a[0] == m / 2 || (a[0] + l[0]) % m == 0 || (a[0] + l[0]) % m == m / 2 || mod((a[0] - p * l[0]), m) == 0 || mod((a[0] - p * l[0]), m) == m / 2 || mod((a[0] - (p - 1) * l[0]), m) == 0 || mod((a[0] - (p - 1) * l[0]), m) == m / 2) {
			return false;
		}
	}
	else {
		int n1, n2, n3, n4;
		for (int i1 = 0; i1 < grp - 1; i1++) {
			n1 = mod((a[i1] - p * l[i1]), m);
			n2 = mod((a[i1] - (p - 1) * l[i1]), m);
			if (a[i1] == n1 || a[i1] == n2 || (a[i1] + l[i1]) % m == n1 || (a[i1] + l[i1]) % m == n2) {
				return false;
			}
			if (a[i1] == m / 2 || (a[i1] + l[i1]) % m == 0 || (a[i1] + l[i1]) % m == m / 2 || n1 == 0 || n1 == m / 2 || n2 == 0 || n2 == m / 2) {
				return false;
			}
			for (int j1 = i1 + 1; j1 < grp; j1++) {
				n3 = mod((a[j1] - p * l[j1]), m);
				n4 = mod((a[j1] - (p - 1) * l[j1]), m);
				if (a[i1] == a[j1] || a[i1] == (a[j1] + l[j1]) % m || a[i1] == n3 || a[i1] == n4) {
					return false;
				}
				if ((a[i1] + l[i1]) % m == a[j1] || (a[i1] + l[i1]) % m == (a[j1] + l[j1]) % m || (a[i1] + l[i1]) % m == n3 || (a[i1] + l[i1]) % m == n4) {
					return false;
				}
				if (n1 == a[j1] || n1 == (a[j1] + l[j1]) % m || n1 == n3 || n1 == n4) {
					return false;
				}
				if (n2 == a[j1] || n2 == (a[j1] + l[j1]) % m || n2 == n3 || n2 == n4) {
					return false;
				}
				if (i1 == grp - 2 && j1 == grp - 1) {
					if (a[j1] == n3 || a[j1] == n4 || (a[j1] + l[j1]) % m == n3 || (a[j1] + l[j1]) % m == n4) {
						return false;
					}
					if (a[j1] == m / 2 || (a[j1] + l[j1]) % m == 0 || (a[j1] + l[j1]) % m == m / 2 || n3 == 0 || n3 == m / 2 || n4 == 0 || n4 == m / 2) {
						return false;
					}
				}
			}
		}
	}
	return true;
}

/*
   This function avoids 2 * 2 tiles and 1 * 2 tiles(or 2 * 1 tiles) forming special rectangles.
*/
bool checkLen(int* l1, int b, int c, int d) {
	for (int i2 = 0; i2 < grp; i2++) {
		if (l1[i2] == c - b || (l1[i2] + (c - b)) == m || l1[i2] == d || (l1[i2] + d) == m || m / 2 == c - b || m / 2 == d) {
			return false;
		}
	}
	return true;
}

/*
   This function checks whether the constructed tile structure is an O-tile structure.
*/
bool checkTileStructure(int* s_a1, int* s_l2, int s_i3, ofstream& s_fout) {
	for (int v = 0; v < m; v++) {
		for (int w = 0; w < m; w++) {
			MM[Tile[v][w] - 1][index[Tile[v][w] - 1]] = v;
			NN[Tile[v][w] - 1][index[Tile[v][w] - 1]] = w;
			index[Tile[v][w] - 1]++;
		}
	}

	M = vector<vector<int>>(s_max, vector<int>(1));
	N = vector<vector<int>>(s_max, vector<int>(1));

	for (int v = 0; v < s_max; v++) {
		M[v].assign(index[v], 0);
		N[v].assign(index[v], 0);
	}

	for (int v = 0; v < s_max; v++) {
		for (int w = 0; w < index[v]; w++) {
			M[v][w] = MM[v][w];
			N[v][w] = NN[v][w];
		}
	}

	for (int v = 0; v < s_max; v++) {
		M[v] = unique_element_in_vector(M[v]);
		N[v] = unique_element_in_vector(N[v]);
	}

	if (checkTilesOverlap()) {
		cout << "\nSorry, the construction is wrong!" << endl;
		cout << "Because when m = " << m << ", tile " << checkTilesOverlap() << " does not meet the construction requirements£¡" << endl;
		cout << "\nThe wrong constructed tile structure is" << endl;
		for (int vv = 0; vv < m; vv++) {
			for (int ww = 0; ww < m; ww++) {
				cout << Tile[vv][ww] << " ";
			}
			cout << endl;
		}
		return false;
	}

	if (checkSpecialRectangle()) {
		writeO_Tile(s_a1, s_l2, s_i3, s_fout);

		number++;
		cout << "\nO-tile structure " << number << " has been constructed." << endl;
	}

	for (int v = 0; v < s_max; v++) {
		for (int w = 0; w < index[v]; w++) {
			MM[v][w] = NN[v][w] = 0;
		}
	}

	for (int v = 0; v < s_max; v++) {
		index[v] = 0;
	}

	M.clear();
	M.shrink_to_fit();
	N.clear();
	N.shrink_to_fit();

	Trans.clear();
	Trans.shrink_to_fit();
	Nums.clear();
	Nums.shrink_to_fit();
	SumM.clear();
	SumM.shrink_to_fit();
	SumN.clear();
	SumN.shrink_to_fit();

	return true;
}

/*
   This function checks whether the constructed tile structure is indeed a tile structure.
*/
int checkTilesOverlap() {
	for (int v1 = 0; v1 < s_max; v1++) {
		if (v1 < m * m / 4 - m) {
			if (M[v1].size() != 2 || N[v1].size() != 2 || index[v1] != 4) {
				return v1 + 1;
			}
		}
		else if (v1 >= m * m / 4) {
			if (M[v1].size() != 2 || N[v1].size() != 1  || index[v1] != 2) {
				return v1 + 1;
			}
		}
		else {
			if (M[v1].size() != 1 || N[v1].size() != 2  || index[v1] != 2) {
				return v1 + 1;
			}
		}
	}
	return 0;
}

/*
   This function checks whether the constructed tile structure contains other special rectangles in addition to the whole grid itself.
*/
bool checkSpecialRectangle() {
	for (int v2 = 0; v2 < s_max - 1; v2++) {
		for (int w2 = v2 + 1; w2 < s_max; w2++) {
			Trans.emplace_back(v2 + 1);
			Trans.emplace_back(w2 + 1);

			while (!Trans.empty()) {
				temp = Trans.back();
				Trans.pop_back();
				Nums.emplace_back(temp);
				SumM = vectors_set_union(SumM, M[temp - 1]);
				SumN = vectors_set_union(SumN, N[temp - 1]);

				if (SumM.size() == m && SumN.size() == m) {
					goto nextpair;
				}

				for (vector<int>::iterator it = SumM.begin(); it < SumM.end(); it++) {
					for (vector<int>::iterator vit = SumN.begin(); vit < SumN.end(); vit++) {
						k = Tile[(*it)][(*vit)] - 1;
						if ((k < v2) || (k != v2 && k < w2)) {
							goto nextpair;
						}
						if ((is_element_in_vector(Nums, k + 1) == false) && (is_element_in_vector(Trans, k + 1) == false)) {
							Trans.emplace_back(k + 1);
						}
					}
				}
			}

			if (SumM.size() < m || SumN.size() < m) {
				return false;
			}

		nextpair:
			Trans.clear();
			Trans.shrink_to_fit();
			Nums.clear();
			Nums.shrink_to_fit();
			SumM.clear();
			SumM.shrink_to_fit();
			SumN.clear();
			SumN.shrink_to_fit();
		}
	}
	return true;
}

/*
   This function writes the O-tile structures that have been constructed to the output file.
*/
void writeO_Tile(int* a1, int* l2, int i3, ofstream& fout) {
	fout << "The value of m is: " << m << endl;
	fout << "The value of p is: " << p << endl;
	if (m == 6) {
		fout << "The side length of 2 * 2 tiles is: {" << m / 2 << "}" << endl;
	}
	else if (m == 10) {
		fout << "The side lengths of 2 * 2 tiles are: {" << l2[0] << ", " << m / 2 << "}" << endl;
	}
	else if (m > 10) {
		for (int v3 = 0; v3 < grp; v3++) {
			if (v3 == 0) {
				fout << "The side lengths of 2 * 2 tiles are: {" << l2[v3] << ", ";
			}
			else if (v3 == grp - 1) {
				fout << l2[v3] << ", " << m / 2 << "}";
			}
			else {
				fout << l2[v3] << ", ";
			}
		}
		fout << endl;
	}
	if (m == 10) {
		fout << "The value of arr is: {" << a1[0] << "}" << endl;
	}
	else if (m > 10) {
		for (int v3 = 0; v3 < grp; v3++) {
			if (v3 == 0) {
				fout << "The value of arr is: {" << a1[v3] << ", ";
			}
			else if (v3 == grp - 1) {
				fout << a1[v3] << "}";
			}
			else {
				fout << a1[v3] << ", ";
			}
		}
		fout << endl;
	}
	fout << "The column indices of 1 * 2 tiles are: {" << p << "i+" << combo[i3][0] << ", " << p << "i+" << combo[i3][1] << "} for 0<=i<=m-1." << endl;
	for (int v3 = 0; v3 < m; v3++) {
		for (int w3 = 0; w3 < m; w3++) {
			fout << Tile[v3][w3] << " ";
		}
		fout << endl;
	}
	fout << endl;
}



