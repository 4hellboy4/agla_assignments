#include <iostream>
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iomanip>
#include <set>

using namespace std;

class Matrix {
public:
    int n;
    int m;
    vector<vector<double>> matrix;

    Matrix(int n, int m) {
        this->n = n;
        this->m = m;
        for (int i = 0; i < n; ++i) {
            vector<double> temp(m);
            for (int j = 0; j < m; ++j) {
                temp[j] = 0;
            }
            matrix.emplace_back(temp);
        }
    }

    Matrix(const Matrix& matrix) {
        this->n = matrix.n;
        this->m = matrix.m;
        this->matrix = matrix.matrix;
    }

    void identifyMatrix(int tempN, int tempM) {
        for (int i = 0; i < tempN; ++i) {
            for (int j = 0; j < tempM; ++j) {
                cin >> matrix[i][j];
            }
        }
    }

    int getN() {
        return this->n;
    }

    int getM() {
        return this->m;
    }

    vector<vector<double>> getMatrix() {
        return this->matrix;
    }

    void setNumber(int tempN, int tempM, double number) {
        this->matrix[tempN][tempM] = number;
    }

    virtual void outputMatrix() {
        for (vector<double> arr : getMatrix()) {
            for (double element : arr) {
                cout << fixed << setprecision(4) << element << " ";
            }
            cout << endl;
        }
    }

    Matrix operator+(const Matrix& other) const  {
        if (n != other.n || m != other.m) {
            throw invalid_argument("Error: the dimensional problem occurred");
        }

        Matrix result(n, m);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                result.matrix[i][j] = matrix[i][j] + other.matrix[i][j];
            }
        }

        return result;
    }

    Matrix operator-(const Matrix& other) const  {
        if (n != other.n || m != other.m) {
            throw invalid_argument("Error: the dimensional problem occurred");
        }

        Matrix result(n, m);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                result.matrix[i][j] = matrix[i][j] - other.matrix[i][j];
            }
        }

        return result;
    }

    Matrix operator*(const Matrix& other) const  {
        if (m != other.n) {
            throw invalid_argument("Error: the dimensional problem occurred");
        }

        Matrix result(n, other.m);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < other.m; ++j) {
                double sum = 0;
                for (int k = 0; k < m; ++k) {
                    sum += matrix[i][k] * other.matrix[k][j];
                }
                result.setNumber(i, j, sum);
            }
        }

        return result;
    }

    Matrix transpose() {
        Matrix result(m, n);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                result.matrix[j][i] = matrix[i][j];
            }
        }

        return result;
    }
};

class SquareMatrix : public Matrix {
public:
    SquareMatrix(int n) : Matrix(n, n) {};
    SquareMatrix(const Matrix& matrix) : Matrix(matrix) {};
};

class IdentityMatrix : public Matrix {
public:
    IdentityMatrix(int n) : Matrix(n, n) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                matrix[i][j] = 0;
            }
        }
        for (int i = 0; i < n; ++i) {
            matrix[i][i] = 1;
        }
    }
    IdentityMatrix(const Matrix& matrix) : Matrix(matrix) {};
};


class EliminationMatrix : public Matrix {
public:
    EliminationMatrix(const Matrix& A, int row, int col) : Matrix(A.matrix.size(), A.matrix[0].size()) {
        for (int i = 0; i < A.matrix.size(); ++i) {
            for (int j = 0; j < A.matrix[0].size(); ++j) {
                if (i == row && j == col) {
                    double factor = -A.matrix[i][col] / A.matrix[col][col];
                    matrix[i][j] = factor;
                } else {
                    matrix[i][j] = (i == j) ? 1 : 0;
                }
            }
        }
    }


};

class PermutationMatrix : public Matrix {
public:
    PermutationMatrix(int n, int i1, int i2) : Matrix(n, n) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == i1) {
                    matrix[i][j] = (j == i2) ? 1 : 0;
                } else if (i == i2) {
                    matrix[i][j] = (j == i1) ? 1 : 0;
                } else {
                    matrix[i][j] = (i == j) ? 1 : 0;
                }
            }
        }
    }
};

double findDeterminant(Matrix A) {
    double det = 1.0;
    int cnt = 1;
    int n = A.n;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i > j) {
                int pivot_row = j;

                for (int k = i; k < n; ++k) {
                    if (abs(A.matrix[k][j]) > abs(A.matrix[pivot_row][j])) {
                        pivot_row = k;
                    }
                }

                if (pivot_row != j) {
                    PermutationMatrix tempPermutation(n, j, pivot_row);
                    A = tempPermutation * A;
                    det *= -1;
                }

                if (A.matrix[i][j] != 0) {
                    EliminationMatrix tempElimination(A, i, j);
                    A = SquareMatrix(tempElimination * A);
                }
            }
        }
    }


    for (int i = 0; i < n; ++i) {
        det *= A.matrix[i][i];
    }

    return det;

}

class ColumnSpace : public Matrix {
public:

    ColumnSpace(int n) : Matrix(n, n) {};

    ColumnSpace(const Matrix& matrix) : Matrix(matrix) {};

    void identifyMatrix(int Bn) {
        for (int i = 0; i < n; ++i) {
            cin >> matrix[i][0];
        }
    }

    void outputMatrix() {
        for (int i = 0; i < n; ++i) {
            cout << fixed << setprecision(2) << matrix[i][0] << " ";
        }
        cout << endl;
    }
};

class InverseMatrix : public Matrix {
public:
    InverseMatrix(int n) : Matrix(n, n) {}

    InverseMatrix(const Matrix& matrix) : Matrix(matrix) {}
};

Matrix getInverseMatrix(const Matrix& matrix) {

    SquareMatrix A(matrix);

    IdentityMatrix I(matrix.n);

    int n = matrix.n;
    int cnt = 1;

    set<int> pivots;

    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            if (i > j) {
                if (pivots.count(j) != 1) {
                    int pivot_row = j;

                    for (int k = j+1; k < n; ++k) {
                        if (abs(A.matrix[k][j]) > abs(A.matrix[pivot_row][j])) {
                            pivot_row = k;
                        }
                    }
                    if (pivot_row != j) {
                        PermutationMatrix tempPermutation(n, j, pivot_row);
                        A = tempPermutation * A;
                        I = tempPermutation * I;

                        cnt++;
                    }
                    pivots.insert(j);
                }

                if (A.matrix[i][j] != 0) {
                    EliminationMatrix tempElimination(A, i, j);
                    A = SquareMatrix(tempElimination * A);
                    I = tempElimination * I;
                }
            }
        }
    }

    for (int j = n-1; j >= 0; --j) {
        for (int i = n-1; i >= 0; --i) {
            if (i < j) {
                if (A.matrix[i][j] != 0) {
                    EliminationMatrix tempElimination(A, i, j);
                    A = SquareMatrix(tempElimination * A);
                    I = tempElimination * I;
                }
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        if (A.matrix[i][i] != 1) {
            for (int j = 0; j < n; ++j) {
                I.matrix[i][j] /= A.matrix[i][i];
            }
            A.matrix[i][i] = 1;
        }
    }

    return I;
}

bool isJacobiApplicable(const Matrix& A) {
    int n = A.n;

    for (int i = 0; i < n; ++i) {
        double diagonal = abs(A.matrix[i][i]);
        double sum = 0;

        for (int j = 0; j < n; ++j) {
            if (j != i) {
                sum += abs(A.matrix[i][j]);
            }
        }

        if (diagonal <= sum) {
            return false;
        }
    }

    return true;
}


int main() {
    int aSize, bSize;
    double eps;
    cin >> aSize;

    Matrix A(aSize, aSize);

    A.identifyMatrix(aSize,aSize);

    cin >> bSize;

    Matrix b(bSize, 1);

    b.identifyMatrix(bSize, 1);

    cin >> eps;

    if (!isJacobiApplicable(A)) {
        cout << "The method is not applicable"<<endl;
        return 0;
    }

    Matrix D(aSize, aSize);
    Matrix L(aSize, aSize);
    Matrix U(aSize, aSize);
    for (int i = 0; i < aSize; ++i) {
        D.matrix[i][i] = A.matrix[i][i];
    }
    for (int i = 0; i < aSize; ++i) {
        for (int j = 0; j < aSize; ++j) {
            if (i < j) {
                U.matrix[i][j] = A.matrix[i][j]*(-1);
            }
        }
    }

    for (int i = 0; i < aSize; ++i) {
        for (int j = 0; j < aSize; ++j) {
            if (i > j) {
                L.matrix[i][j] = A.matrix[i][j]*(-1);
            }
        }
    }
    Matrix B(aSize, aSize);
    Matrix C(aSize, aSize);
    IdentityMatrix I(aSize);

    Matrix alpha(aSize, aSize);
    Matrix beta(bSize, 1);

    Matrix dInverse(aSize, aSize);

    dInverse = getInverseMatrix(D);

    alpha = dInverse*(L+U);
    cout << "alpha:" << endl;

    alpha.outputMatrix();

    beta = dInverse*b;

    cout << "beta:" << endl;

    beta.outputMatrix();

    for (int i = 0; i < aSize; ++i) {
        for (int j = 0; j < aSize; ++j) {
            if (i >= j && alpha.matrix[i][j] != 0) {
                B.matrix[i][j] = alpha.matrix[i][j];
            }
        }
    }

    cout << "B:" << endl;
    B.outputMatrix();

    for (int i = 0; i < aSize; ++i) {
        for (int j = 0; j < aSize; ++j) {
            if (i < j) {
                C.matrix[i][j] = alpha.matrix[i][j];
            }
        }
    }

    cout << "C:" << endl;
    C.outputMatrix();

    Matrix I_B(aSize, aSize);
    I_B = I - B;

    cout << "I-B:" << endl;
    I_B.outputMatrix();

    Matrix I_B_1(aSize,aSize);
    I_B_1 = getInverseMatrix(I_B);

    cout << "(I-B)_-1:" << endl;
    I_B_1.outputMatrix();

    Matrix currentX(bSize, 1);
    currentX = beta;

    int cnt = 1;
    while (true) {
        Matrix xNew(bSize, 1);

        xNew = I_B_1*C*currentX + I_B_1*beta;
        cout << "x(" << cnt << "):" << endl;
        xNew.outputMatrix();

        Matrix tempErrorMatrix(bSize, 1);
        tempErrorMatrix = xNew - currentX;
        double res = 0;
        for (int i = 0; i < bSize; ++i) {
            res += pow(tempErrorMatrix.matrix[i][0],2);
        }
        res = sqrt(res);
        cout << "e: ";
        cout << fixed << setprecision(4) << res << endl;

        currentX = xNew;

        cnt++;
        if (res < eps) {
            cout << "x~:"<< endl;
            xNew.outputMatrix();
            break;
        }
    }
}
