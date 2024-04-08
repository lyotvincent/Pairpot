#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <unordered_map>
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
namespace py = pybind11;
using namespace py;
class matCoo {
public:
    std::vector<double> v;
    std::vector<int> row, col;
    int n, m;
    int totalElements;
    void createMat(int n0, int m0);
    void matTimes(double c);
    void append(int n0, int m0, double val);
    matCoo& operator=(const matCoo& other);
    matCoo()
    {
        row.clear(), col.clear(), v.clear();
        n = 0, m = 0;
        totalElements = 0;
    }
    matCoo(int n0, int m0)
    {
        row.clear(), col.clear(), v.clear();
        n = n0, m = n0;
        totalElements = 0;
    }

};
class mat {
public:
    int n, m;
    std::vector<std::vector<double>> v;
    void createMat(int n0, int m0);
    void matTimes(double c);
    double findDiff(mat);
    mat()
    {
        m = 1, n = 1;
        createMat(1, 1);
    }
    mat(int n0, int m0)
    {
        m = m0, n = n0;
        createMat(n0, m0);
    }
    mat& operator=(const mat& other);
    mat(const mat& other) {
        n = other.n;
        m = other.m;
        createMat(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                v[i][j] = other.v[i][j];
            }
        }
    }
    double getval(int x,int y){return v[x][y];}
    void editval(int x,int y,double val){v[x][y]=val;}
    void setneg(){v.clear();v.resize(n, std::vector<double>(m, -1.0));}
    void editval2(int x,int y){for(int i=0;i<m;i++)v[x][i]=0.0;v[x][y]=1.0;}
};
bool isZero(double a)
{
    return abs(a) < 0.0000000001 ? true : false;
}

void matCoo::createMat(int n0, int m0)
{
    n = n0;
    m = m0;
    totalElements=0;
    row.clear();
    col.clear();
    v.clear();
    return;
}
void matCoo::matTimes(double c)
{
    long long num = col.size();
    if (c == 0)
    {
        row.clear(), col.clear(), v.clear();
        n = 0, m = 0;
        totalElements = 0;
        return;
    }
    for (int i = 0; i < num; i++)
    {
        v[i] *= c;
    }
    return;
}
void matCoo::append(int n0, int m0, double val)
{
    if (isZero(val))  return;
    row.push_back(n0);
    col.push_back(m0);
    v.push_back(val);
    totalElements++;
    return;
}
matCoo& matCoo::operator=(const matCoo& other)
{
    this->v = other.v;
    this->col = other.col;
    this->row = other.row;
    this->m = other.m;
    this->n = other.n;
    this->totalElements = other.totalElements;
    return *this;
}

void mat::createMat(int n0, int m0)
{
    n = n0, m = m0;
    v.clear();
    v.resize(n, std::vector<double>(m, 0));
    return;
}
void mat::matTimes(double c)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            v[i][j] *= c;
    return;
}
double mat::findDiff(mat other)
{
    double curr = 0.0;
    if (n != other.n || m != other.m)  return -1;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
        {
            curr += abs(v[i][j] - other.v[i][j]);
        }
    return curr;
}
mat& mat::operator=(const mat& other)
{
    if (this != &other)
    {
        createMat(other.n, other.m);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                v[i][j] = other.v[i][j];
    }

    return *this;
}
void matMultiply(matCoo &x1, mat &x2,mat &res) {
    int n = x1.n, m = x2.m;
    res.createMat(n, m);
    std::sort(x1.row.begin(), x1.row.end());
    for (int i = 0; i < n; i++) {
        long long p = std::lower_bound(x1.row.begin(), x1.row.end(), i) - x1.row.begin();
        for (int j = 0; j < m; j++) {
            double sum = 0.0;
            int q = p;
            while (q < x1.row.size() && x1.row[q] == i) {
                sum += x1.v[q] * x2.v[x1.col[q]][j];
                q++;
            }
            res.v[i][j] = sum;
        }
    }
    return;
}

void matMultiply(mat &x1, mat &x2,mat &res)
{
    int n = x1.n, m = x2.m, k = x1.m;
    res.createMat(n, m);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            for (int t = 0; t < k; t++)
                res.v[i][j] += x1.v[i][t] * x2.v[t][j];
}
void labelPropagation(matCoo& X, mat& y_label, mat& y_pred,double alpha=0.5, int max_iter=1000) 
{
    int n_samples = X.n;
    int n_classes = y_label.m;
    double diff2=0,diff1=0,diff=0;
    // Initialize
    mat Y = y_label;

    // Compute similarity matrix
    matCoo W(n_samples, n_samples);
    for (int i = 0; i < X.totalElements; ++i) {
        int row = X.row[i];
        int col = X.col[i];
        double val = X.v[i];
        double dist = val * val;
        double similarity = exp(-alpha * dist);
        W.row.push_back(row);
        W.col.push_back(col);
        W.v.push_back(similarity);
    }
    W.totalElements = int(W.row.size());
    mat Y_old;
    mat Y_new;
    // LPA
    for (int iter = 0; iter < max_iter; ++iter) {
        //cout << iter << endl;
        Y_old = Y_new;
        matMultiply(W, Y,Y_new);
        for (int i = 0; i < n_samples; ++i) {
            double row_sum = 0.0;
            
            for (int j = 0; j < n_classes; ++j) {
                row_sum += Y_new.v[i][j];
            }
            for (int j = 0; j < n_classes; ++j) {
                Y_new.v[i][j] /= row_sum;
            }
            bool has_prior = false;
            for (int j = 0; j < n_classes; ++j) {
                if (y_label.v[i][j] != -1) {
                    Y_new.v[i][j] = y_label.v[i][j];
                    has_prior = true;
                    break;
                }
            }
            if (!has_prior) {
                for (int j = 0; j < n_classes; ++j) {
                    Y.v[i][j] = Y_new.v[i][j];
                }
            }
        }
        diff2 = diff1;
        diff1 = diff;
        diff = Y_old.findDiff(Y_new);
        
        diff = diff / n_samples;
        if (iter == 0)  continue;
        if (diff<1e-3) {
            break;
        }
        if (abs(diff - diff2) < 1e-3 && abs(diff2 - diff1) < 1e-3)
        {
            //break;
        }
    }

    // result
    y_pred.createMat(n_samples, 1);
    for (int i = 0; i < n_samples; ++i) {
        int max_index = 0;
        double max_value = Y.v[i][0];
        for (int j = 1; j < n_classes; ++j) {
            if (Y.v[i][j] > max_value) {
                max_value = Y.v[i][j];
                max_index = j;
            }
        }
        y_pred.v[i][0] = max_index;
    }
    // ofstream out1("res.txt");
    // for (int i = 0; i < y_pred.n; ++i) 
    // {
    //     out1 <<y_pred.v[i][0] << endl;
    // }
}

// int main()
// {
    // srand(time(0));
    // clock_t start = clock();
    // ifstream in1("out.txt");
    // ifstream in2("out2.txt");
    // // ofstream out1("res.txt");
    // int n0,m0,m1;
    // in1>>n0;
    // in2>>m0>>m1;
    // matCoo X(m0,m0);
    // mat y_label(m0,m1);
    // int x,y;
    // double z;
    // ios::sync_with_stdio(false);
    // for(int i=0;i<n0;i++)
    // {
    //     in1>>x>>y>>z;
    //     X.append(x,y,z);
    // }

    // for(int i=0;i<m0;i++)
    // {
    //     in2>>x;
    //     y=rand()%10;
    //     if(!y)
    //       for(int r=0;r<m1;r++)  y_label.v[i][r]=(r==x)?1:0;
    //     else
    //       for(int r=0;r<m1;r++)  y_label.v[i][r]=-1;
    // }
    // mat res;
    // mat y_pred;
    // labelPropagation(X, y_label,y_pred);
    // clock_t end = clock();
    // double elapsed = double(end - start) / CLOCKS_PER_SEC;
    // std::cout << elapsed << std::endl;
    // getchar();
// }
PYBIND11_MODULE(label_propagation, m) {
    py::class_<matCoo>(m, "matCoo")
        .def(py::init<int, int>())
        .def("createMat", &matCoo::createMat)
        .def("matTimes", &matCoo::matTimes)
        .def("append", &matCoo::append);

    py::class_<mat>(m, "mat")
        .def(py::init<int, int>())
        .def("createMat", &mat::createMat)
        .def("matTimes", &mat::matTimes)
        .def("findDiff", &mat::findDiff)
        .def("editval",&mat::editval)
        .def("getval",&mat::getval)
        .def("setneg",&mat::setneg)
        .def("editval2",&mat::editval2);
    m.def("labelPropagation", &labelPropagation);
}
