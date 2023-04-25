// TRS2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <math.h>
#include <vector>
#include <tuple>
#include <iomanip>
namespace var9 
{
    int k0 = 1;
    double f(double t, double x)
    {
        return t * x * x - t * t;
    }
    int alpha0 = 0;
    int betta0 = 1;
    int alpha1 = 1;
    int betta1 = 1;
    double phi(double x)
    {
        return x;
    }
    double psi0(double t)
    {
        return 1;
    }
    double psi1(double t)
    {
        return 3*t*t/2+2;
    }
    double l = 1.0;
};
using namespace std;
using namespace var9;
void PrintMatrix(vector<vector<double>> Matrix)
{
    cout << fixed << std::setprecision(4);
    cout << "-------------------------------------------------------------" << endl;
    for (int i = 0; i < Matrix.size(); i++)
    {
        for (int j = 0; j < Matrix[i].size(); j++)
        {
            cout << setw(5) << Matrix[i][j] << "  ";
        }
        cout << endl;
    }
    cout << "-------------------------------------------------------------" << endl;
}
vector<vector<double>> ExplicitSchemeMethod(double tau, double h)
{
    int n_big = int(l / h)+1;
    vector <double> X;
    vector <double> T;
    vector <vector<double>> U; // for U t is the first arg, x is second
    for (int i = 0; i < n_big; i++)
    {
        X.push_back(i * h);
        T.push_back(i * tau);
    }
    for (int m = 0; m < n_big; m++)
    {
        vector<double> temp;
        U.push_back(temp);
        for (int n = 0; n < n_big; n++)
        {
            U[m].push_back(0.0);
        }
    }
    for (int i = 0; i < n_big; i++)
    {
        U[0][i] = phi(X[i]);
    }
    for (int j = 1; j < n_big; j++)
    {
        for (int i = 1; i < n_big - 1; i++)
        {
            U[j][i] = tau / (h * h) * (U[j - 1][i - 1] - 2 * U[j - 1][i] + U[j - 1][i + 1]) + tau * (T[j] * X[i] * X[i]
                + T[j] * T[j]) + U[j - 1][i];
        }
        
    }
    for (int j = 1; j < n_big; j++)
    {
        U[j][0] = -h + U[j][1];
        U[j][n_big - 1] = (h * (3 / 2 * T[j] * T[j] + 2) + U[j][n_big - 2]) / ( h + 1);
    }
    
    PrintMatrix(U);
    return U;
}
vector<vector<double>> ImplicitSchemeMethod(double tau, double h)
{
    int n_big = int(l / h) + 1;
    vector <double> X;
    vector <double> T;
    vector <vector<double>> U; // for U t is the first arg, x is second
    
    for (int i = 0; i < n_big; i++)
    {
        X.push_back(i * h);
        T.push_back(i * tau);
    }
    for (int m = 0; m < n_big; m++)
    {
        vector<double> temp(n_big, 0.);
        U.push_back(temp);
    }
    for (int i = 0; i < n_big; i++)
    {
        U[0][i] = phi(X[i]);
    }
    for (int n = 1; n < n_big; n++)
    {
        vector <double> altha(n_big,0), betta(n_big,0);
        altha[0] = 1; betta[0] = -h * psi0(0);
        for (int i = 1; i <= n_big - 2; i++)
        {
            double a = -tau / (h * h);
            double b = 1 + 2 * tau / (h * h);
            double c = -tau / (h * h);
            double z = U[n][i] + tau * (T[n] * X[i] * X[i] + T[n] * T[n]);
            altha[i] = -a / (b + c * altha[i - 1]);
            betta[i] = (z - c * betta[i - 1]) / (b + c * altha[i - 1]);
        }
        U[n][n_big - 1] = ( (3 * T[n] * T[n] / 2 + 2) * h + betta[n_big - 2] ) / (1 + h - altha[n_big - 2]);
        for (int i = n_big - 2; i >= 0; i--)
        {
            U[n][i] = altha[i] * U[n][i + 1] + betta[i];
        }
    }
    PrintMatrix(U);
    return U;
}
int main()
{
    // x from 0 to 1, t is more than 0
    double h = 0.1; // so then N is 10, and T=tau*N=50 ms
    double tau = 0.005;
    ExplicitSchemeMethod(tau, h);
    ImplicitSchemeMethod(tau, h);
}

// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
