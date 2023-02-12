#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

void print(int n, double** tab)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << tab[i][j];
            cout << "\t\t";
        }
        cout << endl;
    }
}

double** copy(int n, double** tab, double** tab_copy)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            tab_copy[i][j] = tab[i][j];
        }
    }
    return tab_copy;
}

double** check_LU(int n, double** tab)
{
    double** tab2;
    double** tab4;
    tab2 = new double* [n];
    tab4 = new double* [n];
    for (int i = 0; i < n; i++)
    {
        tab2[i] = new double[n + 1];
        tab4[i] = new double[n + 1];
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                tab2[i][j] = 1;
            }
            else if (i > j)
            {
                tab2[i][j] = tab[i][j];
                tab[i][j] = 0;
            }
            else
            {
                tab2[i][j] = 0;
            }
        }
    }

    cout << "\nMacierz L" << endl;
    print(n, tab2);
    cout << "\nMacierz U" << endl;
    print(n, tab);

    for (int i = 0; i < n; i++) //mno¿enie macierzy
    {
        for (int j = 0; j < n; j++)
        {
            double s = 0;
            for (int k = 0; k < n; k++)
            {
                s += tab2[i][k] * tab[k][j];
                tab4[i][j] = s;
            }
        }
    }
    tab = copy(n, tab4, tab);
    for (int i = 0; i < n; i++)
    {
        delete[] tab2[i];
        delete[] tab4[i];
    }
    return tab;
}

void LU_error(int n, double** tab_ob, double** tab_pier)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            tab_ob[i][j] = abs(tab_ob[i][j] - tab_pier[i][j]);
        }
    }

    double max = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (tab_ob[i][j] > max)
            {
                max = tab_ob[i][j];
            }
        }
    }
    cout << setprecision(6) << scientific;
    cout << "\nBlad: " << max << endl;
}

bool gauss(const int n, double** tab, double* r)
{
    double** tab3;
    double** tab4;
    tab3 = new double* [n];
    tab4 = new double* [n];
    for (int i = 0; i < n; i++)
    {
        tab3[i] = new double[n + 1];
        tab4[i] = new double[n + 1];
    }
    tab3 = copy(n, tab, tab3);

    for (int k = 0; k < n - 1; k++) //dekompozycja LU
    {
        for (int i = k + 1; i < n; i++)
        {
            double aux = tab[i][k] / tab[k][k];
            for (int j = k + 1; j <= n; j++)
            {
                tab[i][j] = tab[i][j] - tab[k][j] * aux;
            }
            tab[i][k] = aux;
        }
    }

    cout << "\nMacierz po dekompozycji LU" << endl;
    print(n, tab);
    tab4 = copy(n, tab, tab4);
    tab = check_LU(n, tab);

    cout << "\nSprawdzenie - L*U = A" << endl;
    print(n, tab);
    LU_error(n, tab, tab3);

    double wek3[4];
    for (int i = 0; i < n; i++)
    {
        wek3[i] = r[i];
    }

    for (int k = 0; k < n - 1; k++)//eliminacja w przod
    {
        for (int i = k + 1; i < n; i++)
        {
            r[i] = r[i] - r[k] * tab4[i][k];
        }
    }
    r[n - 1] = r[n - 1] / tab4[n - 1][n - 1];//podstawienie wstecz
    for (int i = n - 2; i >= 0; i--)
    {
        double s = 0;
        for (int j = i + 1; j < n; j++)
        {
            s = s + tab4[i][j] * r[j];
        }
        r[i] = (r[i] - s) / tab4[i][i];
    }

    cout << "\nEliminacja w przod + podstawienie wstecz" << endl;
    cout << "\nWynik" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << setprecision(6) << scientific;
        cout << r[i] << "\t";
    }

    double wek2[4];
    for (int i = 0; i < n; i++) //A*x
    {
        double s = 0;
        for (int j = 0; j < n; j++)
        {
            s += tab3[i][j] * r[j];
        }
        wek2[i] = s;
    }

    for (int i = 0; i < n; i++)//|a*x-r|
    {
        r[i] = abs(wek2[i] - wek3[i]);
    }
    double max = 0;
    for (int i = 0; i < n; i++)
    {
        if (r[i] > max)
        {
            max = r[i];
        }
    }
    cout << setprecision(6) << scientific;
    cout << "\n\nBlad: " << max << endl;

    for (int i = 0; i < n; i++)
    {
        delete[] tab3[i];
        delete[] tab4[i];
    }
    return true;
}

int main()
{
    const int n = 3;
    double** tab;
    double r[n] = { 0.7351, 0.9036, 0.7427 };
    tab = new double* [n];
    ifstream file;
    file.open("dane.txt", ios::in);
    for (int i = 0; i < n; i++)
    {
        tab[i] = new double[n + 1];
    }
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {

            file >> tab[i][j];
        }
    }
    file.close();

    cout << "\nMacierz pierwotna A" << endl;
    cout << setprecision(2) << fixed;
    print(n, tab);

    cout << "\nWektor r" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << r[i];
        cout << "\t";
    }
    cout << "\n\nDekompozycja LU" << endl;
    gauss(n, tab, r);
    for (int i = 0; i < n; i++)
    {
        delete[] tab[i];
    }
}