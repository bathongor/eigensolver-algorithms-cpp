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

    void printPair(){
        cout << this->rowIndex << " " << this->colIndex << " " << this->value;
    }
};

struct CSR {
    vector<int> I, J;
    vector<double> V;

    void setCSR(vector<int> I, vector<int> J, vector<double> V){
        this->I = I;
        this->J = J;
        this->V = V;
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

vector<double> scaleVector(double s, vector<double> v){
    for(int i = 0; i < v.size(); i++){
        v[i] = s*v[i];
    }
    return v;
}

vector<double> subtractVectors(vector<double> a, vector<double> b){
    vector<double> c;
    for(int i = 0; i < a.size(); i++){
        c.push_back(a[i] - b[i]);
    }
    return c;
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

void printCSR(string name, CSR A){
    cout << name << ": " << '\n';
    for(int i = 0; i < A.V.size(); i++){
        cout << A.I[i] << " " << A.J[i] << " " << A.V[i] << '\n';
    }
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

CSR assembleCSR(vector<Pair> csrMatrix, bool isSymmetric){
    // Initialize your matrix;
	vector<int> I, J; // initialize row and col index array for CSR
	vector<double> V; // initialize values array that contain CSR matrix values

	// Assemble CSR row, col and V based on the symmetry
	int prevRow = -1;
    // cout << "csrMatrix: " << csrMatrix.size();
	for(int i = 0; i < csrMatrix.size(); i++){
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
    CSR temp;
    temp.setCSR(I, J, V);
    return temp;
}

vector<double> normalize(vector<double> v){
    double v_norm = norm(v);
    for(int i = 0; i < v.size(); i++){
        v[i] = v[i] / v_norm;
    }
    return v;
}

double powerIteration(vector<int> &I, vector<int> &J, vector<double> &V, int n, bool isSymmetric, double tol){
    vector<double> q, lambda, z(n, (double) 1/sqrt(n)), q_prev(n, (double) 1/sqrt(n));
    lambda.push_back(0);
    int k = 1;
    while(true){
        // printVector("z[k]", matrixVectorProduct(I, J, V, q[k-1], isSymmetric));
        z = matrixVectorProduct(I, J, V, q_prev, isSymmetric);
        lambda.push_back(dotProduct(q_prev, z));
        if(fabs(lambda[k] - lambda[k-1]) < tol)break;
        q_prev = normalize(z);
        k++;
    };
    cout << "lambda: " << lambda[k] << " iterations: " << k << '\n';
    return lambda[k];
}

CSR lanczos(vector<int> &I, vector<int> &J, vector<double> &V, int n, int m, bool isSymmetric){
    vector<vector<double>> v(m+2, vector<double> (n, 0.0));
    vector<double> alpha(m+1, 0.0), beta(m+1, 0.0);
    vector<double> v1(n, (double) 1/sqrt(n)), w;
    v[1] = v1;
    for(int j = 1; j <= m; j++){
        w = matrixVectorProduct(I, J, V, v[j], isSymmetric);
        w = subtractVectors(w, scaleVector(beta[j-1], v[j-1]));
        alpha[j] = dotProduct(v[j], w);
        w = subtractVectors(w, scaleVector(alpha[j], v[j]));
        beta[j] = norm(w);
        v[j+1] = scaleVector((double) 1/beta[j], w);
    }

    int count = 0;
    vector<Pair> TriDiag(m*m);
    for(int i = 1; i <= m; i++){
        for(int j = 1; j <= m; j++){
            if(i == j) {
                TriDiag[count].setPair(i, j, alpha[i]);
                count++;
            }
            else if(i == j+1){
                TriDiag[count].setPair(i, j, beta[i]);
                count++;
            }             
        }
    }
    TriDiag.resize(count);
    CSR Tm = assembleCSR(TriDiag, isSymmetric);
    return Tm;
}

int main(){
	// Open the file:
	ifstream fin("s3rmt3m3.mtx");
	int M, N, L;
	while(fin.peek() == '%') fin.ignore(2048, '\n');
	fin >> M >> N >> L;
	vector<Pair> csrMatrix(L);

	for(int i = 0; i < L; i++){
		int m, n;
		double data;
		fin >> m >> n >> data;
		csrMatrix[i].setPair(m, n, data);
	}

	fin.close();
	sort(csrMatrix.begin(), csrMatrix.begin()+L, sortHelper);
	bool isSymmetric = true;
    CSR lgMatrix = assembleCSR(csrMatrix, isSymmetric);
    CSR smMatrix = lanczos(lgMatrix.I, lgMatrix.J, lgMatrix.V, N, 100, isSymmetric);
    powerIteration(smMatrix.I, smMatrix.J, smMatrix.V, N, isSymmetric, 0.0000000001);
	return 0;
}