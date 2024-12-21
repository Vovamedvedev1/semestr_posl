#include <iostream>
#include <iomanip>
#include <windows.h>
#include <sstream>
#include <fstream>
#include <string>
#include <time.h>
#include <mpi.h>
#include <vector>
#include <map>
using namespace std;

void permutation(vector<int>& vc)
{
    int i, j;
    for (i = vc.size() - 2; i >= 1; i--)
    {
        if (vc[i] < vc[i + 1])
            break;
    }
    for (j = vc.size() - 2; j > i; j--)
    {
        if (vc[j] > vc[i])
            break;
    }
    if (i == -1 || j == -1)
        return;
    int temp;
    temp = vc[i];
    vc[i] = vc[j];
    vc[j] = temp;
    int k = vc.size() - 2;
    for (int m = i + 1; m <= (k + i + 1) / 2; m++)
    {
        temp = vc[m];
        vc[m] = vc[k - (m - i - 1)];
        vc[k - (m - i - 1)] = temp;
    }
}

int notUsed(vector<int>& used, long long blockNum) 
{
    int j;
    long long pos = 0;
    for (j = 1; j < (int)used.size(); j++) {
        if (!used[j]) pos++;
        if (blockNum == pos)
            break;
    }
    return j;
}


void permutation_by_num(vector<int>& re, int n, long long num, long long countMarshruts)
{
    vector<int>used(n + 1, 0);
    int temp = n;
    countMarshruts = countMarshruts / temp;
    for (int i = 0; i < n; i++)
    {
        long long blockNum = (num - 1) / countMarshruts + 1;
        int j = notUsed(used, blockNum);
        re[i + 1] = j;
        used[j] = 1;
        num = (num - 1) % countMarshruts + 1;
        temp--;
        if (temp != 0)
            countMarshruts = countMarshruts / temp;
    }
}

void solver(int& res, vector<vector<int>>& re_res, int**& a, int n, long long startPermutation, long long endPermutation, long long countMarshruts)
{
    res = 0;
    vector<int>re(n + 1, 0);
    re[0] = 0;
    re[n] = 0;
    permutation_by_num(re, n - 1, startPermutation + 1, countMarshruts);
    for (int i = 1; i < (int)re.size(); i++)
        res += a[re[i - 1]][re[i]];
    re_res.push_back(vector<int>(re));
    int res_temp = 0;
    while (startPermutation != endPermutation)
    {
        permutation(re);
        res_temp = 0;
        for (int i = 1; i < (int)re.size(); i++)
            res_temp += a[re[i - 1]][re[i]];
        if (res_temp < res)
        {
            res = res_temp;
            re_res.clear();
            re_res.push_back(re);
        }
        else if (res_temp == res)
        {
            re_res.push_back(re);
        }
        startPermutation++;
    }
}

int main(int argc, char** argv)
{
    setlocale(LC_ALL, "Russian");
    int rank, size;
    double start = 0, stop = 0;
    MPI_Status stat;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n;
    string file = "ex1.txt";
    int** a{nullptr};
    start = MPI_Wtime();
    ifstream ofs(file);
    if (ofs.is_open())
    {
        ofs >> n;
        a = new int* [n];
        a[0] = new int[n * n];
        for (int i = 1; i < n; i++)
            a[i] = a[i - 1] + n;
        ofs.clear();
        ofs.seekg(0, ios::beg);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++)
            {
                ofs >> a[i][j];
            }
        }

    }
    else
    {
        cout << "----Ошибка открытия файла----" << endl;
        ofs.close();
        MPI_Abort(MPI_COMM_WORLD, 1);
        MPI_Finalize();
        return 0;
    }
    ofs.close();
    vector<vector<int>>re_res(0);
    long long countMarshruts = 1;
    for (int i = 1; i <= n - 1; i++)
        countMarshruts *= i;
    long long startPermutation = 0;
    long long endPermutation = startPermutation + countMarshruts;
    int res = 0;
    solver(res, re_res, a, n, startPermutation, endPermutation, countMarshruts);
    ofstream ifs("output.txt");
    if (ifs.is_open())
    {
        for (int i = 0; i < (int)re_res.size(); i++)
        {
            for (int j = 0; j < (int)re_res[i].size(); j++)
                ifs << re_res[i][j] << " ";
            ifs << endl;
        }
    }
    else
    {
        cout << "ERROR!!!" << endl;
        MPI_Finalize();
        return 0;
    }
    ifs.close();
    stop = MPI_Wtime();
    cout << setprecision(10) << fixed << "TIME:  " << stop - start << endl;
    cout << "Длина оптимального пути:  " << res << "   " << endl;
    for (int i = 0; i < (int)re_res.size(); i++)
    {
        for (int j = 0; j < (int)re_res[i].size(); j++)
            cout << re_res[i][j] << " ";
        cout << endl;
    }
    delete[]a[0];
    delete[]a;
    MPI_Finalize();
    return 0;
}
