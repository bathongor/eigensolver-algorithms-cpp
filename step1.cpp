#include <iostream>
#include <cstdlib>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace std;

struct Pair {
	int rowIndex;
	int colIndex;
	double value;

	void setPair(int row, int col, double val){
		this->rowIndex = row;
		this->colIndex = col;
		this->value = val;
	}
};

double norm(vector<double> &a){
	double sumSqr = 0;
	for(int i = 0; i < a.size(); i++){
		sumSqr += a[i]*a[i];
	}
	return sqrt(sumSqr);
}

double dotProduct(vector<double> a, vector<double> b){
    double tmp = 0;
    for(int i = 0; i < a.size(); i++){
        tmp += a[i]*b[i];
    }
    return tmp;
}

void printVector(string name, vector<double> x){
	cout << name << ": ";
	for(int i = 0; i < x.size(); i++){	
		cout << " " << x[i];
	}
	cout << '\n';
}

void printMatrix(string name, vector<vector<double>> &A){
	cout << name << ": " << '\n';
	for(int i = 0; i < A.size(); i++){
		for(int j = 0; j < A[i].size(); j++){	
			cout << " " << A[i][j];
		}
		cout << '\n';
	}
	cout << '\n';
}

bool sortHelper(Pair a, Pair b){
	return a.rowIndex < b.rowIndex;
}

vector<double> matrixVectorProduct(vector<int> &I, vector<int> &J, vector<double> &V, vector<double> x, bool isSymmetric){
	vector<double> y(x.size(), 0);
	for(int i = 0; i < I.size()-1; i++){
		int i1 = I[i]-1;
		int i2 = I[i+1]-1;
		for(int k = i1; k < i2; k++){
			y[i] += V[k]*x[J[k]-1];
			if(isSymmetric && J[k] != i+1){
				y[J[k]-1] += V[k]*x[i];
			}
		}
	}

	return y;
}

void getKrylov(vector<int> &I, vector<int> &J, vector<double> &V, vector<vector<double>> &v, int j, vector<vector<double>> &h, int m, bool isSymmetric){
	vector<double> w = matrixVectorProduct(I, J, V, v[j], isSymmetric);
	vector<double> h_colj;
	for(int i = 0; i <= j; i++){
		h_colj.push_back(inner_product(w.begin(), w.end(), v[i].begin(), 0.0));
		for(int k = 0; k < w.size(); k++){
			w[k] = w[k] - h_colj[i]* v[i][k];
		}
	}
	h_colj.push_back(norm(w));
	vector<double> v_next;
	for(int l = 0; l < w.size(); l++){
		v_next.push_back(w[l]/h_colj[j+1]);
	}
	h_colj.resize(m+1, 0);
	h[j] = h_colj;
	v.push_back(v_next);
}

vector<double> normalize(vector<double> v){
    double v_norm = norm(v);
    for(int i = 0; i < v.size(); i++){
        v[i] = v[i] / v_norm;
    }
    return v;
}

double powerIteration(vector<int> &I, vector<int> &J, vector<double> &V, int n, bool isSymmetric){
    vector<double> q, lambda, z(n, (double) 1/sqrt(n)), q_prev(n, (double) 1/sqrt(n));
    lambda.push_back(0);
    int k = 1;
    while(true){
        // printVector("z[k]", matrixVectorProduct(I, J, V, q[k-1], isSymmetric));
        z = matrixVectorProduct(I, J, V, q_prev, isSymmetric);
        lambda.push_back(dotProduct(q_prev, z));
        if(fabs(lambda[k] - lambda[k-1]) < 0.00000001)break;
        q_prev = normalize(z);
        k++;
    };
    cout << "lambda: " << lambda[k] << " iterations: " << k;
    return lambda[k];
}

int main(){
	
	// Open the file:
	ifstream fin("nos6.mtx");
	int M, N, L;
	while(fin.peek() == '%') fin.ignore(2048, '\n');
	fin >> M >> N >> L;
	Pair csrMatrix[L+1], originalMatrix[L+1];

	for(int i = 0; i < L; i++){
		int m, n;
		double data;
		fin >> m >> n >> data;
		csrMatrix[i].setPair(m, n, data);
		originalMatrix[i].setPair(m, n, data);
	}

	fin.close();
	sort(csrMatrix, csrMatrix+L, sortHelper);
	bool isSymmetric = true;

	// Initialize your matrix;
	vector<int> I, J; // initialize row and col index array for CSR
	vector<double> V; // initialize values array that contain CSR matrix values

	// Assemble CSR row, col and V based on the symmetry
	int prevRow = -1;
	for(int i = 0; i < L; i++){
		if((isSymmetric == true && csrMatrix[i].colIndex <= csrMatrix[i].rowIndex) || isSymmetric == false){
			J.push_back(csrMatrix[i].colIndex);
			V.push_back(csrMatrix[i].value);
			if(prevRow != csrMatrix[i].rowIndex){
				I.push_back(J.size());
			}
			prevRow = csrMatrix[i].rowIndex;	
		}
	} 
	I.push_back(J.size()+1);
    powerIteration(I, J, V, N, isSymmetric);
	return 0;
}