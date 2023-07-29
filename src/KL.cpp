#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>
#include <cmath>
#include <queue>
#include <sys/timeb.h>

using namespace std;

// 添加零元素起始编号（即最大处理 10000 个任务图）
#define ZEROBEGINNUM 10000
// 迭代 5 次，求出最小总代价的划分作为最终总代价的划分
#define EPOCH 5
// 当有小于 20 个任务图时，以临界矩阵形式显示，不小于时，以 【头节点->尾节点 权值】形式显示
#define MATRIXORARC 20

struct KL1_2
{
    set<int> A;
    set<int> B;
    double *D;
    int EN;
    int n;
};

struct Dendrogram
{
    double threshold;
    int blockNum;
    set<set<int>> block;
};

int countZeroNum(KL1_2 kl)
{
    int zeroNum = 0;
    set<int>::iterator it = kl.A.begin();
    for (it; it != kl.A.end(); it++)
    {
        if (*it >= ZEROBEGINNUM)
        {
            zeroNum++;
        }
    }
    for (it = kl.B.begin(); it != kl.B.end(); it++)
    {
        if (*it >= ZEROBEGINNUM)
        {
            zeroNum++;
        }
    }
    return zeroNum;
}

int calD(KL1_2 kl, const int i, double **am)
{
    // int zeroNum = countZeroNum(kl);
    if (i < 0 || i > kl.n)
    {
        return 0;
    }
    double E = 0, I = 0, D;
    // i in set A
    if (kl.A.count(i))
    {
        set<int>::iterator it = kl.A.find(i);
        // exterior cost
        set<int>::iterator it1 = kl.B.begin();
        for (it1; it1 != kl.B.end(); it1++)
        {
            // 防止越界
            if (*it1 >= kl.n || *it1 < 0)
            {
                E += 0;
            }
            else
            {
                E += am[i][*it1];
            }
        }
        // interior cost
        set<int>::iterator it2 = kl.A.begin();
        for (it2; it2 != kl.A.end(); it2++)
        {
            if (it2 != it && *it2 < kl.n && *it2 >= 0)
            {
                I += am[i][*it2];
            }
        }
        // i in set B
    }
    else
    {
        set<int>::iterator it = kl.B.find(i);
        // exterior cost
        set<int>::iterator it1 = kl.A.begin();
        for (it1; it1 != kl.A.end(); it1++)
        {
            if (*it1 >= kl.n || *it1 < 0)
            {
                E += 0;
            }
            else
            {
                E += am[i][*it1];
            }
        }
        // interior cost
        set<int>::iterator it2 = kl.B.begin();
        for (it2; it2 != kl.B.end(); it2++)
        {
            if (it2 != it && *it2 >= 0 && *it2 < kl.n)
            {
                I += am[i][*it2];
            }
        }
    }
    return E - I;
}

void createNumDAG(const string graphName, KL1_2 &kl, double **&am)
{
    /*
    M: The number of modules
    EN: Max of the communication costs
    n: The number of tasks
    m: The number of arches
    u: The precursor of arch
    v: The successor of arch
    */
    int n = 0, m, k;
    int u, v;
    double c;
    ifstream infile(graphName);

    while (!infile.eof())
    {
        infile >> u >> v >> c;
        kl.A.insert(u);
        kl.A.insert(v);
    }
    n = kl.A.size();
    kl.n = n;
    kl.EN = n;

    infile.clear(ios::goodbit);
    infile.seekg(0);

    am = new double *[n];
    for (int i = 0; i < n; i++)
    {
        am[i] = new double[n]();
    }
    while (!infile.eof())
    {

        infile >> u >> v >> c;
        am[u][v] = c;
        am[v][u] = c;
    }
    infile.close();

    // print adjacency matrix

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (n < MATRIXORARC)
            {
                printf("%.1f   ", am[i][j]);
            }
            else
            {
                if (am[i][j])
                {
                    cout << i << "->" << j << "  " << am[i][j] << endl;
                }
            }
        }
        if (n < MATRIXORARC)
        {
            cout << endl;
        }
    }
    cout << endl;

    kl.D = new double[n]();
    for (int i = 0; i < n; i++)
    {
        kl.D[i] = calD(kl, i, am);
    }
}

// 计算并返回集合 A、B 中最大增益，及交换对 a、b
double maxProceeds(KL1_2 kl, double **am, int &a, int &b)
{
    double G = -__FLT_MAX__;
    int a1 = -1, b1 = -1;
    set<int>::iterator it1 = kl.A.begin();
    for (it1; it1 != kl.A.end(); it1++)
    {
        set<int>::iterator it2 = kl.B.begin();
        for (it2; it2 != kl.B.end(); it2++)
        {
            // 对增加的额外零元素不进行增益计算
            if (*it1 >= kl.n || *it2 >= kl.n || *it1 < 0 || *it2 < 0)
            {
                break;
            }

            double Da = calD(kl, *it1, am);
            double Db = calD(kl, *it2, am);
            double g = Da + Db - 2 * am[*it1][*it2];
            if (g > G)
            {
                G = g;
                a1 = *it1;
                b1 = *it2;
            }
        }
    }
    a = a1;
    b = b1;
    return G;
}

void copy(KL1_2 kl, KL1_2 &klCopy)
{
    int zeroNum = 0;
    set<int>::iterator ita = kl.A.begin();
    for (ita; ita != kl.A.end(); ita++)
    {
        klCopy.A.insert(*ita);
        if (*ita >= ZEROBEGINNUM)
        {
            zeroNum++;
        }
    }
    set<int>::iterator itb = kl.B.begin();
    for (itb; itb != kl.B.end(); itb++)
    {
        klCopy.B.insert(*itb);
        if (*itb >= ZEROBEGINNUM)
        {
            zeroNum++;
        }
    }
    klCopy.D = new double[kl.n - zeroNum]();
    for (int i = 0; i < kl.n - zeroNum; i++)
    {
        klCopy.D[i] = kl.D[i];
    }
    klCopy.n = kl.n;
    klCopy.EN = kl.EN;
}

void klPartition(KL1_2 &kl, double **am)
{
    KL1_2 klCopy;
    copy(kl, klCopy);
    int a[kl.n / 2], b[kl.n / 2];
    double g[kl.n / 2];
    int zeroNum = countZeroNum(kl);

    for (int i = 0; i < kl.n / 2; i++)
    {
        g[i] = maxProceeds(klCopy, am, a[i], b[i]);
        // a,b 值为-1 表示出现零元素
        if (a[i] != -1 || b[i] != -1)
        {
            klCopy.A.erase(a[i]);
            klCopy.B.erase(b[i]);
        }

        for (int j = 0; j < kl.n - zeroNum; j++)
        {
            klCopy.D[j] = calD(kl, j, am);
        }
    }

    double maxG = 0;
    int k;
    for (int i = 0; i < kl.n / 2; i++)
    {
        double G = 0;
        int j = i;
        for (j; j >= 0; j--)
        {
            G += g[j];
        }
        if (G > maxG)
        {
            maxG = G;
            k = i;
        }
    }
    // cout << "g[i]:\t";
    // for (int i = 0; i < kl.EN; i++)
    // {
    //     cout << g[i] << "\t";
    // }
    // cout << endl;

    if (maxG > 0)
    {
        for (int i = 0; i <= k; i++)
        {
            // a、b 不为零元素时，交换a、b
            if (a[i] != -1 && b[i] != -1)
            {
                kl.A.erase(a[i]);
                kl.A.insert(b[i]);
                kl.B.erase(b[i]);
                kl.B.insert(a[i]);
            }
        }
        klPartition(kl, am);
    }
    else
    {

        // cout << "final partition:" << endl;
        // print(kl, am);
    }
}

void addZero(set<int> &s, double **&am, int setNum, int zeroNum)
{
    int i = 0;
    double **newAm = new double *[setNum + zeroNum];
    double **temp = am;
    for (int i = 0; i < setNum + zeroNum; i++)
    {
        newAm[i] = new double[setNum + zeroNum]();
    }
    for (int i = 0; i < setNum; i++)
    {
        for (int j = 0; j < setNum; j++)
        {
            newAm[i][j] = am[i][j];
        }
    }

    while (zeroNum > 0)
    {
        s.insert(ZEROBEGINNUM + i);
        i++;
        zeroNum--;
    }

    am = newAm;

    for (int i = 0; i < setNum; i++)
    {
        delete[] temp[i];
    }
    delete[] temp;
}

void deleteZero(KL1_2 &kl)
{
    int zeroNum = 0;
    set<int>::iterator it, temp;
    for (it = kl.A.begin(); it != kl.A.end(); it)
    {
        if (*it < 0 || *it >= ZEROBEGINNUM)
        {
            temp = ++it;
            kl.A.erase(--it);
            it = temp;
            zeroNum++;
        }
        else
        {
            it++;
        }
    }
    for (it = kl.B.begin(); it != kl.B.end(); it)
    {
        if (*it < 0 || *it >= ZEROBEGINNUM)
        {
            temp = ++it;
            kl.B.erase(--it);
            it = temp;
            zeroNum++;
        }
        else
        {
            it++;
        }
    }
    kl.n -= zeroNum;
}

// 此时数据全部在集合A中，通过平衡，将元素平衡的划分到A和B集合中，最终A集合元素大于等于B集合元素
void setBalancez(KL1_2 &kl, int kModule)
{
    int i, n = kl.n;
    int en = kl.EN;
    vector<int> shuffled;

    set<int>::iterator ita = kl.A.begin();
    for (ita; ita != kl.A.end(); ita++)
    {
        if (*ita < ZEROBEGINNUM)
        {
            shuffled.push_back(*ita);
        }
    }
    set<int>::iterator itb = kl.B.begin();
    for (itb; itb != kl.B.end(); itb++)
    {
        if (*itb < ZEROBEGINNUM)
        {
            shuffled.push_back(*itb);
        }
    }
    // 产生毫秒级随机种子
    struct timeb timeSeed;
    ftime(&timeSeed);
    srand(timeSeed.time * 1000 + timeSeed.millitm);
    // 将打乱A集合中元素分布，尽量扩大“探索”范围，以避免避免局部最优情况
    shuffle(shuffled.begin(), shuffled.end(), default_random_engine((unsigned int)rand()));
    kl.A.clear();

    if (kModule <= 1)
    {
        return;
    }
    else if (kModule == 2)
    {
        // 保证A集合元素不少于B集合元素
        if ((n % 2 == 0 && en < n / 2) || (n % 2 == 1 && en < (n + 1) / 2))
        {
            en = n - en;
        }
        for (i = 0; i < en; i++)
        {
            kl.A.insert(shuffled[i]);
        }
        for (i = en; i < n; i++)
        {
            kl.B.insert(shuffled[i]);
        }
    }
    else
    {
        int temp = n / 2;
        if (n % 2 == 1)
        {
            temp++;
        }
        if (temp <= en)
        {
            temp = en;
        }
        for (i = 0; i < temp; i++)
        {
            kl.A.insert(shuffled[i]);
        }
        for (i = temp; i < n; i++)
        {
            kl.B.insert(shuffled[i]);
        }
    }
}

void setSplit(KL1_2 &kl, double **am, KL1_2 &k1, KL1_2 &k2)
{
    // set<int>::iterator ita = kl.A.begin();
    // for (ita; ita != kl.A.end(); ita++)
    // {
    //     k1.A.insert(*ita);
    // }
    int zeroNum = 0;
    set<int>::iterator it;
    if (kl.B.empty())
    {
        return;
    }
    else
    {
        it = kl.B.begin();
    }

    for (it; it != kl.B.end(); it++)
    {
        k2.A.insert(*it);
        if (*it >= ZEROBEGINNUM)
        {
            zeroNum++;
        }
    }
    k1.B.clear();
    k2.n = k2.A.size() - zeroNum;

    zeroNum = 0;
    it = k1.A.begin();
    for (auto e : kl.A)
    {
        if (e >= ZEROBEGINNUM)
        {
            zeroNum++;
        }
    }
    k1.n = k1.A.size() - zeroNum;

    
    k1.D = new double[kl.n / 2]();
    for (int i = 0; i < kl.n / 2; i++)
    {
        k1.D[i] = calD(k1, i, am);
    }
    
    k2.D = new double[k2.n / 2]();
    for (int i = 0; i < k2.n / 2; i++)
    {
        k2.D[i] = calD(k2, i, am);
    }

}

void mmmsaa(Dendrogram &dendrogram, double **&am, const int M, const int EN)
{
    while (dendrogram.threshold >= 0 && dendrogram.blockNum != M)
    {
        set<set<int>>::iterator it1 = dendrogram.block.begin();
        while (it1 != dendrogram.block.end())
        {
            set<int> k1 = *it1;
            it1++;
            set<set<int>>::iterator it2 = it1;
            while (it2 != dendrogram.block.end())
            {
                set<int> k2 = *it2;
                it2++;
                double largestCost = 0;
                for (set<int>::iterator it3 = k1.begin(); it3 != k1.end(); it3++)
                {
                    for (set<int>::iterator it4 = k2.begin(); it4 != k2.end(); it4++)
                    {
                        if (am[*it3][*it4] > largestCost)
                        {
                            largestCost = am[*it3][*it4];
                        }
                    }
                }
                if (largestCost >= dendrogram.threshold && k1.size() + k2.size() <= EN)
                {
                    set<int> k12(k1);
                    k12.insert(k2.begin(), k2.end());

                    if (*it1 == k1 || *it1 == k2)
                    {
                        it1++;
                    }
                    dendrogram.block.erase(dendrogram.block.find(k1));
                    k1 = k12;
                    if (*it2 == k2 || *it2 == k1)
                    {
                        it2++;
                    }
                    dendrogram.block.erase(dendrogram.block.find(k2));

                    dendrogram.block.insert(k12);
                    dendrogram.blockNum--;
                    // cout << dendrogram.threshold << "\t" << dendrogram.blockNum << "\t";
                    // printSet(dendrogram, am);
                    // cout << endl;
                }
            }
        }
        dendrogram.threshold--;
    }
}

void creatDendrogram(Dendrogram &dendrogram, KL1_2 *kl, double **am, int blockNum)
{
    set<int> s;
    int zeroNum = 0;
    double threshold = am[0][0];
    for (int i = 0; i < blockNum / 2; i++)
    {
        for (auto e : kl[i].A)
        {
            s.insert(e);
        }
        if(!s.empty()){
            dendrogram.block.insert(s);
        }
        s.clear();

        for (auto e : kl[i].B)
        {
            s.insert(e);
        }
        if(!s.empty()){
            dendrogram.block.insert(s);
        }
        s.clear();

        zeroNum = countZeroNum(kl[i]);
        for (int i = 0; i < kl[i].n - zeroNum; i++)
        {
            for (int j = 0; j < kl[i].n - zeroNum; j++)
            {
                if (threshold < am[i][j])
                {
                    threshold = am[i][j];
                }
            }
        }
    }
    dendrogram.blockNum = dendrogram.block.size();

    dendrogram.threshold = ++threshold;
}

void printSet(Dendrogram dendrogram, double **am)
{
    set<set<int>>::iterator it = dendrogram.block.begin();
    int i = 0;
    while (it != dendrogram.block.end())
    {
        set<int> subset = *it;
        set<int>::iterator ele = subset.begin();
        cout << "set " << i << " : ";
        for (ele; ele != subset.end(); ele++)
        {
            cout << *ele << "  ";
        }
        cout << endl;
        ;
        it++;
        i++;
    }
}

Dendrogram klPart(KL1_2 &kl, double **&am, const int moduleNum)
{
    // 两个集合元素随机初始化，A集合包含集合最多包含元素数，B集合包含剩余元素
    Dendrogram dendrogram;
    int n = kl.n;
    int zeroNum = 0;
    int en = kl.EN;

    // 默认集合 A 优先集合 B “填满”
    if (moduleNum == 2 && n % 2 == 0)
    {
        setBalancez(kl, moduleNum);
        // 1优化2划分
        klPartition(kl, am);
        creatDendrogram(dendrogram, &kl, am, moduleNum);
    }
    else if (moduleNum == 2 && n % 2 == 1)
    {
        setBalancez(kl, moduleNum);
        // 奇数个元素集合的划分
        zeroNum = 1;
        addZero(kl.B, am, kl.n, zeroNum);
        kl.n += zeroNum;
        klPartition(kl, am);
        deleteZero(kl);
        creatDendrogram(dendrogram, &kl, am, moduleNum);
    }
    else if (moduleNum == 2)
    {
        setBalancez(kl, moduleNum);
        // 模块尺寸不相等划分
        addZero(kl.B, am, kl.n, 2 * en - n);
        klPartition(kl, am);
        deleteZero(kl);
        creatDendrogram(dendrogram, &kl, am, moduleNum);
    }
    else
    {
        // kl多式划分
        int fL = floor(log2(moduleNum));
        int cL = ceil(log2(moduleNum));

        if (fL == cL)
        {
            setBalancez(kl, moduleNum);
            // case 1
            addZero(kl.B, am, n, moduleNum * en - n);
            // 每个 kl 包含 A、B 两个集合
            KL1_2 *klRst = new KL1_2[moduleNum / 2];
            for (int i = 0; i < moduleNum / 2; i++)
            {
                klRst[i].A.clear();
                klRst[i].B.clear();
                klRst[i].D = NULL;
                klRst[i].EN = kl.EN;
                klRst[i].n = 0;
            }
            copy(kl, klRst[0]);
            int h = 0;
            int maxLast = 0;
            queue<int> pending;
            pending.push(0);
            // 经过 h 次划分，得到 moduleNum/2 个kl，即 moduleNum 个集合
            while (!pending.empty() && h < log2(moduleNum / 2))
            {
                int setOne = pending.front();
                int setTwo = setOne + moduleNum / (2 * (pow(2, h + 1)));
                pending.pop();

                setSplit(klRst[setOne], am, klRst[setOne], klRst[setTwo]);
                setBalancez(klRst[setOne], moduleNum);
                setBalancez(klRst[setTwo], moduleNum);
                klPartition(klRst[setOne], am);
                klPartition(klRst[setTwo], am);
                if (setOne != setTwo)
                {
                    pending.push(setOne);
                    pending.push(setTwo);
                }

                if (maxLast == setOne)
                {
                    h++;
                    maxLast = setTwo;
                }
            }
            for(int i = 0; i < moduleNum/2 ; i++){
                deleteZero(klRst[i]);
            }
            creatDendrogram(dendrogram, klRst, am, moduleNum);
        }
        else
        {
            // case 2
            dendrogram = klPart(kl, am, pow(2, cL));
            // printSet(dendrogram, am);
            mmmsaa(dendrogram, am, moduleNum, en);
        }
    }
    return dendrogram;
}

double totalCost(Dendrogram dendrogram, double **am)
{
    double cost = 0, costGroup = 0, mean = 0, var = 0;
    set<set<int>>::iterator it1 = dendrogram.block.begin();
    while (it1 != dendrogram.block.end())
    {
        set<int> k1 = *it1;
        it1++;
        set<set<int>>::iterator it2 = it1;
        while (it2 != dendrogram.block.end())
        {
            set<int> k2 = *it2;
            it2++;

            for (set<int>::iterator it3 = k1.begin(); it3 != k1.end(); it3++)
            {
                for (set<int>::iterator it4 = k2.begin(); it4 != k2.end(); it4++)
                {
                    cost += am[*it3][*it4];
                }
            }
        }
    }
    return cost;
}

// 计算总代价、均值、方差
void calCMV(Dendrogram dendrogram, double **am)
{
    // 计算总通信代价
    double cost = 0, costGroup = 0, mean = 0, var = 0;
    set<set<int>>::iterator it1 = dendrogram.block.begin();
    while (it1 != dendrogram.block.end())
    {
        set<int> k1 = *it1;
        it1++;
        set<set<int>>::iterator it2 = it1;
        while (it2 != dendrogram.block.end())
        {
            set<int> k2 = *it2;
            it2++;

            for (set<int>::iterator it3 = k1.begin(); it3 != k1.end(); it3++)
            {
                for (set<int>::iterator it4 = k2.begin(); it4 != k2.end(); it4++)
                {
                    cost += am[*it3][*it4];
                }
            }
        }
    }

    mean = cost / dendrogram.block.size();

    // 计算方差
    it1 = dendrogram.block.begin();
    while (it1 != dendrogram.block.end())
    {
        cout << "set(" << *(*it1).begin() << ")~~";
        set<int> k1 = *it1;
        it1++;
        set<set<int>>::iterator it2 = it1;
        while (it2 != dendrogram.block.end())
        {
            cout << "~set(" << *(*it2).begin() << ")  ";
            set<int> k2 = *it2;
            it2++;
            costGroup = 0;

            for (set<int>::iterator it3 = k1.begin(); it3 != k1.end(); it3++)
            {
                for (set<int>::iterator it4 = k2.begin(); it4 != k2.end(); it4++)
                {
                    costGroup += am[*it3][*it4];
                }
            }
            cout << "[group cost: " << costGroup << "]\t";
        }
        cout << endl;
        var += pow(fabs(mean - costGroup), 2);
    }
    var /= dendrogram.block.size();
    cout << endl;
    cout << "cost: " << cost << endl;
    cout << "mean: " << mean << endl;
    cout <<  "var: " << var << endl;
}

void minCostSelect(KL1_2 kl, double **am, int k, int en)
{
    KL1_2 klCopy;
    kl.EN = en;
    copy(kl, klCopy);

    Dendrogram d[EPOCH];
    double minCost = INT_MAX, cost = 0;
    int minD = 0;
    // cout << "epoch: " << EPOCH << " cost: ";
    for (int i = 0; i < EPOCH; i++)
    {
        d[i] = klPart(kl, am, k);
        cost = totalCost(d[i], am);
        if (cost < minCost)
        {
            minCost = cost;
            minD = i;
        }
        // cout << cost << "  ";
        copy(klCopy, kl);
    }
    // cout << endl;
    printSet(d[minD], am);
    cout << endl;
    calCMV(d[minD], am);
}



int main()
{
    KL1_2 kl;
    double **am;
    Dendrogram dendrogram;


    // 7-3
    // createNumDAG("data/7-3.csv", kl, am);
    // minCostSelect(kl, am, 3, 4);

    // 21
    // createNumDAG("data/t21.csv", kl, am);
    // 21-5
    // minCostSelect(kl, am, 5, 5);
    // 21-6
    // minCostSelect(kl, am, 6, 4);
    // 21-7
    // minCostSelect(kl, am, 7, 4);

    // 33
    // createNumDAG("data/t33.csv", kl, am);
    // 33-5
    // minCostSelect(kl, am, 5, 7);
    // 33-6
    // minCostSelect(kl, am, 6, 6);
    // 33-7
    // minCostSelect(kl, am, 7, 6);

    // 51
    createNumDAG("data/t51.csv", kl, am);
    // 51-5
    // minCostSelect(kl, am, 5, 11);
    // 51-6
    // minCostSelect(kl, am, 6, 9);
    // 51-7
    minCostSelect(kl, am, 7, 9);

    // printSet(dendrogram, am);
    // calCMV(dendrogram, am);
    // totalCost(dendrogram, am);

    return 0;
}