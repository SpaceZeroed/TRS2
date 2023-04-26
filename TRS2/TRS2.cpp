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
    double k(double u)
    {
        return sin(u);
    }
    double F(double u)
    {
        return u;
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
    int n_big = int(l / h) + 1;
    vector <double> X(n_big);
    vector <double> T(n_big);
    vector <vector<double>> U; // for U t is the first arg, x is second
    for (int i = 0; i < n_big; i++)
    {
        X[i] = i * h;
        T[i] = i * tau;
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
    for (int j = 1; j < n_big; j++) 
    {
        for (int i = 1; i <= n_big - 2; i++) 
        {
            U[j][i] = tau / (h * h) * (U[j - 1][i - 1] - 2 * U[j - 1][i] + U[j - 1][i + 1]) + tau * f(T[j], X[i]) + U[j - 1][i];
            
        }
        U[j][0] = -h + U[j][1];
        U[j][n_big - 1] = ( h * (3 / 2 * T[j] * T[j] + 2) + U[j][n_big - 2] ) / ( h + 1);
    }

    PrintMatrix(U);
    return U;
}
vector<vector<double>> ImplicitSchemeMethod(double tau, double h)
{
    int n_big = int(l / h) + 1;
    vector <double> X(n_big);
    vector <double> T(n_big);
    vector <vector<double>> U; // for U t is the first arg, x is second
    for (int i = 0; i < n_big; i++)
    {
        X[i] = i * h;
        T[i] = i * tau;
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
        altha[0] = 1; betta[0] = -h * psi0(T[n]);
        for (int i = 1; i <= n_big - 2; i++)
        {
            double a = -tau / (h * h);
            double b = 1 + 2 * tau / (h * h);
            double c = -tau / (h * h);
            double z = U[n - 1][i] + tau * f(T[n], X[i]);
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
vector<vector<double>> KrankNicholsonScheme(double tau, double h)
{
    int n_big = int(l / h) + 1;
    vector <double> X(n_big);
    vector <double> T(n_big);
    vector <vector<double>> U; // for U t is the first arg, x is second
    for (int i = 0; i < n_big; i++)
    {
        X[i] = i * h;
        T[i] = i * tau;
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
        vector <double> altha(n_big, 0), betta(n_big, 0);
        altha[0] = 1; betta[0] = -h * psi0(T[n]);
        for (int i = 1; i <= n_big - 2; i++)
        {
            double a = -tau / (2 * h * h);
            double b = 1 +  tau / ( h * h);
            double c = -tau / (2 * h * h);
            double z = U[n - 1][i] + tau * (U[n-1][i-1] - 2* U[n-1][i] + U[n-1][i+1]) / (2 * h * h ) +  tau * f(T[n], X[i]);
            altha[i] = -a / (b + c * altha[i - 1]);
            betta[i] = (z - c * betta[i - 1]) / (b + c * altha[i - 1]);
        }
        U[n][n_big - 1] = ((3 * T[n] * T[n] / 2 + 2) * h + betta[n_big - 2]) / (1 + h - altha[n_big - 2]);
        for (int i = n_big - 2; i >= 0; i--)
        {
            U[n][i] = altha[i] * U[n][i + 1] + betta[i];
        }
    }
    PrintMatrix(U);
    return U;
}
vector<vector<double>> ConservativeScheme(double tau, double h)
{
    int n_big = int(l / h) + 1;
    vector <double> X(n_big);
    vector <double> T(n_big);
    vector <vector<double>> U; // for U t is the first arg, x is second
    for (int i = 0; i < n_big; i++)
    {
        X[i] = i * h;
        T[i] = i * tau;
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
        vector <double> altha(n_big, 0), betta(n_big, 0);
        altha[0] = 1; betta[0] = -h * psi0(T[n]);
        for (int i = 1; i <= n_big - 2; i++)
        {
            double a = -tau * k(U[n-1][i] / 2. + U[n-1][i-1] /  2.) / (h * h);
            double b = 1 + tau * k(U[n - 1][i] / 2. + U[n - 1][i - 1] / 2.) / (h * h) + tau *  k(U[n - 1][i + 1] / 2. + U[n - 1][i] / 2.) / h /h ;
            double c = -tau * k(U[n - 1][i+1] / 2. + U[n - 1][i] / 2.) / ( h * h);
            double z = U[n - 1][i] +  tau * F(U[n-1][i]) * f(T[n], X[i]) / h;
            altha[i] = -a / (b + c * altha[i - 1]);
            betta[i] = (z - c * betta[i - 1]) / (b + c * altha[i - 1]);
		}
        U[n][n_big - 1] = ((3 * T[n] * T[n] / 2 + 2) * h + betta[n_big - 2]) / (1 + h - altha[n_big - 2]);
        for (int i = n_big - 2; i >= 0; i--)
        {
            U[n][i] = altha[i] * U[n][i + 1] + betta[i];
        }
    }
    PrintMatrix(U);
    return U;
}
vector<vector<double>> SpesialConservativeScheme(double tau, double h)
{
    int n_big = int(l / h) + 1;
    vector <double> X(n_big);
    vector <double> T(n_big);
    vector <vector<double>> U; // for U t is the first arg, x is second
    for (int i = 0; i < n_big; i++)
    {
        X[i] = i * h;
        T[i] = i * tau;
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
        double diff = 0;
        vector <double> altha(n_big, 0), betta(n_big, 0);
        altha[0] = 1; betta[0] = -h * psi0(T[n]);
        for (int i = 1; i <= n_big - 2; i++)
        {
            double a = -tau * k(U[n - 1][i] / 2. + U[n - 1][i - 1] / 2.) / (h * h);
            double b = 1 + tau * k(U[n - 1][i] / 2. + U[n - 1][i - 1] / 2.) / (h * h) + tau * k(U[n - 1][i + 1] / 2. + U[n - 1][i] / 2.) / h / h;
            double c = -tau * k(U[n - 1][i + 1] / 2. + U[n - 1][i] / 2.) / (h * h);
            double z = U[n - 1][i] + tau * F(U[n - 1][i]) * f(T[n], X[i]) / h;
            altha[i] = -a / (b + c * altha[i - 1]);
            betta[i] = (z - c * betta[i - 1]) / (b + c * altha[i - 1]);
        }
        U[n][n_big - 1] = ((3 * T[n] * T[n] / 2 + 2) * h + betta[n_big - 2]) / (1 + h - altha[n_big - 2]);
        for (int i = n_big - 2; i >= 0; i--)
        {
            U[n][i] = altha[i] * U[n][i + 1] + betta[i];
        }

        int Q = 0, M;
        vector<double> prev;
        do
        {
            M = 0;
            Q++;
            prev.assign(U[n].begin(), U[n].end());
            altha[0] = 1; betta[0] = -h * psi0(T[n]);
            for (int i = 1; i <= n_big - 2; i++)
            {
                double A = -tau * k(U[n][i] / 2. + U[n][i - 1] / 2.) / (h * h);
                double B = 1 + tau * k(U[n][i] / 2. + U[n][i - 1] / 2.) / (h * h) + tau * k(U[n][i + 1] / 2. + U[n][i] / 2.) / h / h;
                double C = -tau * k(U[n][i + 1] / 2. + U[n][i] / 2.) / (h * h);
                double Z = U[n - 1][i] + tau * F(U[n - 1][i]) * f(T[n], X[i]) / h;
                altha[i] = -A / (B + C * altha[i - 1]);
                betta[i] = (Z - C * betta[i - 1]) / (B + C * altha[i - 1]);

            }
            U[n][n_big - 1] = ((3 * T[n] * T[n] / 2 + 2) * h + betta[n_big - 2]) / (1 + h - altha[n_big - 2]);
            for (int i = n_big - 2; i >= 0; i--)
            {
                U[n][i] = altha[i] * U[n][i + 1] + betta[i];
            }

            for (int s = 1; s < n_big - 1; s++)
            {
                diff = abs(k(U[n][s - 1] / 2. + U[n][s - 1] / 2.) - k(prev[s - 1] / 2. + prev[s - 1] / 2.));
                cout << diff << endl;
                if (diff > M) {
                    M = diff;
                }
            }
            prev.clear();
            cout << Q << " " << M << endl;
        } while (M > 1e-15);
    }  

    PrintMatrix(U);
    return U;
}
int main()
{
    // x from 0 to 1, t is more than 0
    double h = 0.1; // so then N is 10, and T=tau*N=50 ms
    double tau = 0.005;
    //ex1
    ExplicitSchemeMethod(tau, h);
    //ex2
    //ImplicitSchemeMethod(tau, h);
    //KrankNicholsonScheme(tau, h);
    //ex3
    ConservativeScheme(tau, h);
    SpesialConservativeScheme(tau, h);
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
