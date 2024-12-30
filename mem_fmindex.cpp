//For running:     
//     clang++ -std=c++14 -o fm mem_fmindex.cpp
//    ./fm

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include <algorithm>
#include <random>
#include <cstdint>
#include <chrono>
using namespace std;

// defining block size and integer size
static int blockSize = 92;
static int intSize = 32;

// implementing on a bitvector is rank1(i) = count the number of 1s in B[1...i]
int rank1(unsigned int bitvector, int i) {
    unsigned int rank = 0;
    while (i > 0) {
        rank += bitvector & 1;
        bitvector = bitvector >> 1;
        i--;
    }
    return rank;
}

//implementing on a bitvector is rank0(i) = count the number of 0s in B[1...i]
int rank0(unsigned int bitvector, int i) {
    unsigned int rank = i;
    while (i > 0) {
        rank -= bitvector & 1;
        bitvector = bitvector >> 1;
        i--;
    }
    return rank;
}

// Radix sort implementaion


inline bool leq(int a1, int a2, int b1, int b2)
{
    return(a1 < b1 || a1 == b1 && a2 <= b2);
}

inline bool leq(int a1, int a2, int a3, int b1, int b2, int b3)
{
    return(a1 < b1 || a1 == b1 && leq(a2, a3, b2, b3));
}

// stably sort a[0..n-1] to b[0..n-1] with keys in 0..K from r
static void radixPass(int* a, int* b, int* r, int n, int K)
{// count occurrences
    int* c = new int[K + 1]; // counter array
    for (int i = 0; i <= K; i++)
        c[i] = 0; // reset counters
    for (int i = 0; i < n; i++)
        c[r[a[i]]]++; // count occurrences
    for (int i = 0, sum = 0; i <= K; i++) // exclusive prefix sums
    {
        int t = c[i];
        c[i] = sum;
        sum += t;
    }
    for (int i = 0; i < n; i++)
        b[c[r[a[i]]]++] = a[i]; // sort
    delete[] c;
}

// Suffix array
void suffixArray(int* s, int* SA, int n, int K)
{
    int n0 = (n + 2) / 3, n1 = (n + 1) / 3, n2 = n / 3, n02 = n0 + n2;
    int* s12 = new int[n02 + 3]; s12[n02] = s12[n02 + 1] = s12[n02 + 2] = 0;
    int* SA12 = new int[n02 + 3]; SA12[n02] = SA12[n02 + 1] = SA12[n02 + 2] = 0;
    int* s0 = new int[n0];
    int* SA0 = new int[n0];

    for (int i = 0, j = 0; i < n + (n0 - n1); i++)
        if (i % 3 != 0)
            s12[j++] = i;
    radixPass(s12, SA12, s + 2, n02, K);
    radixPass(SA12, s12, s + 1, n02, K);
    radixPass(s12, SA12, s, n02, K);
    int name = 0, c0 = -1, c1 = -1, c2 = -1;
    for (int i = 0; i < n02; i++)
    {
        if (s[SA12[i]] != c0 || s[SA12[i] + 1] != c1 || s[SA12[i] + 2] != c2)
        {
            name++;
            c0 = s[SA12[i]];
            c1 = s[SA12[i] + 1];
            c2 = s[SA12[i] + 2];
        }
        if (SA12[i] % 3 == 1)
        {
            s12[SA12[i] / 3] = name;
        }
        else
        {
            s12[SA12[i] / 3 + n0] = name;
        }
    }
    if (name < n02)
    {
        suffixArray(s12, SA12, n02, name);
        for (int i = 0; i < n02; i++)
            s12[SA12[i]] = i + 1;
    }
    else
    {
        for (int i = 0; i < n02; i++)
        SA12[s12[i] - 1] = i;
    }
    for (int i = 0, j = 0; i < n02; i++)
        if (SA12[i] < n0)
            s0[j++] = 3 * SA12[i];
    radixPass(s0, SA0, s, n0, K);
    for (int p = 0, t = n0 - n1, k = 0; k < n; k++)
    {
        #define GetI() (SA12[t] < n0 ? SA12[t]*3+1: (SA12[t] - n0) * 3 + 2)
        int i = GetI();
        int j = SA0[p];
        if (SA12[t] < n0 ? leq(s[i], s12[SA12[t] + n0], s[j], s12[j / 3]) : leq(s[i], s[i + 1], s12[SA12[t] - n0 + 1], s[j], s[j + 1], s12[j / 3 + n0]))
        {
            SA[k] = i; t++;
            if (t == n02)
                for (k++; p < n0; p++, k++) SA[k] = SA0[p];
        }
        else
        {
            SA[k] = j; p++;
            if (p == n0)
                for (k++; t < n02; t++, k++)
                    SA[k] = GetI();
        }
    }
    delete[] s12; delete[] SA12; delete[] SA0; delete[] s0;
}

//Building Burrows-Wheeler transform (BWT) from Suffix Array
int BurrowsWheelerTransform(int* s, int* BurrowsWheelerTransform, int* SA, int n)
{
    int j;
    int bwtresult = std::numeric_limits<int>::min();

    //int bwtresult = INT_MIN ;
    for (int i = 0; i < n; i++)
    {
        j = SA[i] - 1;
        if (j < 0)
            j = j + n;
        BurrowsWheelerTransform[i] = s[j];
        bwtresult = max(BurrowsWheelerTransform[i], bwtresult);
    }
    return bwtresult;
}

//Reading text file
std::string readTXTFile(const char* filename)
{
    std::FILE* fp = std::fopen(filename, "rb");
    if (fp)
    {
        std::string filecontent;
        std::fseek(fp, 0, SEEK_END);
        filecontent.resize(std::ftell(fp));
        std::rewind(fp);
        std::fread(&filecontent[0], 1, filecontent.size(), fp);
        std::fclose(fp);
        return(filecontent);
    }
    throw(errno);
}

// Creating wavelet tree
class waveletTree {
    public:

    int minval, maxval;

    waveletTree* left, * right;

    vector<unsigned int> packedArray;
    vector<int> indexArray;

    waveletTree(int* beg, int* endp, int low, int high)
    {
        unsigned int currentInteger = 0;
        int bitPosition = 0;
        int oneCount = 0;

        minval = low, maxval = high;
        if (beg >= endp)
            return;

        packedArray.reserve(ceil((endp - beg + 1) / blockSize));
        if (maxval == minval) {
            for (auto j = beg; j < endp; j++)
            {
                bitPosition = ((j - beg) % intSize);
                currentInteger += pow(2, bitPosition);
                oneCount += 1;
                if (bitPosition == (intSize - 1)) {
                    packedArray.push_back(currentInteger);
                    currentInteger = 0;
                }
                if (bitPosition == (blockSize - 1)) {
                    indexArray.push_back(oneCount + (indexArray.empty() ? 0 : indexArray.back()));
                    oneCount = 0;
                }
            }
            if (bitPosition != (intSize - 1)) {
                packedArray.push_back(currentInteger);
            }
            if (bitPosition != (blockSize - 1)) {
                indexArray.push_back(oneCount + (indexArray.empty() ? 0 : indexArray.back()));
            }
            return;
        }

        int mid = (minval + maxval) / 2;

        auto lessThanMid = [mid](int x) {
            return x <= mid;
        };

        for (auto j = beg; j < endp; j++)
        {
            int countLessThanMid = !lessThanMid(*j);
            bitPosition = ((j - beg) % intSize);
            currentInteger += countLessThanMid * pow(2, bitPosition);
            oneCount += countLessThanMid;
            if (bitPosition == (intSize - 1)) {
                packedArray.push_back(currentInteger);
                currentInteger = 0;
            }
            if ((bitPosition + 1) % blockSize == 0) {
                indexArray.push_back(oneCount + (indexArray.empty() ? 0 : indexArray.back()));
                oneCount = 0;
            }
        }
        if (bitPosition != (intSize - 1)) {
            packedArray.push_back(currentInteger);
        }
        if (bitPosition != (blockSize - 1)) {
            indexArray.push_back(oneCount + (indexArray.empty() ? 0 : indexArray.back()));
        }

        auto pivotVal = stable_partition(beg, endp, lessThanMid);

        left = new waveletTree(beg, pivotVal, minval, mid);
        right = new waveletTree(pivotVal, endp, mid + 1, maxval);
    }

    // Printing packed array
    void print()
    {
        for (int i = 0; i < packedArray.size(); i++)
            cout << packedArray.at(i) << ", ";
        cout << "\n";
    }

    int rank(int charVal, int index)
    {
        auto lessThanMid = [this](int x) {
            return x <= (minval + maxval) / 2;
        };

        if (minval == maxval)
        {
            return index;
        }
        else
        {
            int intindex = (int) index / intSize;
            int blockindex = (int)index / blockSize;
            int bitnumber = index % intSize;
            unsigned int bitvector = packedArray[intindex] >> ((int)bitnumber/blockSize)*blockSize;
            int i = bitnumber % blockSize;
            if (bool(lessThanMid(charVal)))
            {
                int result = rank0(bitvector, i);
                if (blockindex > 0) {
                    result += blockSize*(blockindex) - indexArray[blockindex - 1];
                }

                return (*left).rank(charVal, result);
            }
            else
            {
                int result = rank1(bitvector, i);
                if (blockindex > 0) {
                    result += indexArray[blockindex - 1];
                }
                return (*right).rank(charVal, result);
            }
        }
    }
};

 // implementing backward pattern matching using FM-index
int patternMatching(int patternlength, int pattern[], int* charactertable, int startpoint, int endpoint, waveletTree wtree)
{
    int count = 1;
    for (int i = patternlength - 1; i >= 0; i--)
    {
        int charpos = pattern[i];
        startpoint = charactertable[charpos] + wtree.rank(charpos, startpoint - 1) + 1;
        endpoint = charactertable[charpos] + wtree.rank(charpos, endpoint);
    }
    return endpoint - startpoint + 1;
}

void computeCharacterRank(int* charactertable, int* SA, int* s, int n)
{
    int count = 1;
    for (int i = 1; i < n; i++)
    {
        if (s[SA[i]] == s[SA[i - 1]])
            count += 1;
        else
        {
            charactertable[s[SA[i]]] = charactertable[s[SA[i - 1]]] + count;
            count = 1;
        }
    }
}

void printTotalSpace(int textLength, int alphabetSize) {
    size_t inputTextSpace = sizeof(char) * textLength;
    size_t textArraySpace = sizeof(int) * (textLength + 3);
    size_t suffixArraySpace = sizeof(int) * (textLength + 1);
    size_t bwtSpace = sizeof(int) * textLength;
    size_t waveletTreeSpace = textLength * log2(alphabetSize) / 8; // Approximation

    size_t totalSpace = inputTextSpace + textArraySpace + suffixArraySpace + bwtSpace + waveletTreeSpace;

    cout << "Space Breakdown (in MB):\n";
    cout << "Input Text: " << inputTextSpace / (1024.0 * 1024.0) << " MB\n";
    cout << "Text Array: " << textArraySpace / (1024.0 * 1024.0) << " MB\n";
    cout << "Suffix Array: " << suffixArraySpace / (1024.0 * 1024.0) << " MB\n";
    cout << "BWT: " << bwtSpace / (1024.0 * 1024.0) << " MB\n";
    cout << "Wavelet Tree: " << waveletTreeSpace / (1024.0 * 1024.0) << " MB\n";
    cout << "Total Space: " << totalSpace / (1024.0 * 1024.0) << " MB\n";
}


int main()
{
    cout << "Reading text file \n";
    string input = readTXTFile("dna_50mb.txt");
    int textLength = input.size();

    int n = input.size();
    int* s = new int[n + 3];

    for (int i = 0; i < n; i++)
        s[i] = (int)input[i];

    cout << "Creating suffix array\n";
    namespace sc = std::chrono;
    auto timenow = sc::system_clock::now();
    auto timeSinceEpoch = timenow.time_since_epoch();
    auto timemillisec = sc::duration_cast<sc::milliseconds>(timeSinceEpoch);
    long starttime = timemillisec.count();
    int* SA = new int[n + 1];
    int K = 128;
    int* C = new int[K];
    std::fill(C, C + K, 0);
    suffixArray(s, SA, n, K);
    timenow = sc::system_clock::now();
    timeSinceEpoch = timenow.time_since_epoch();
    timemillisec = sc::duration_cast<sc::milliseconds>(timeSinceEpoch);
    long endtime = timemillisec.count();
    cout << "Suffix Array Creation Time: " << endtime-starttime << " milliseconds.\n";

    int startpoint = 1, endpoint = n;;
    cout << "Building Burrows-Wheeler Transform (BWT) from suffix array : \n";
    int* bwt = new int[n];

    int maximum = BurrowsWheelerTransform(s, bwt, SA, n);

    cout << "Constructing Wavelet Tree\n";
    timenow = sc::system_clock::now();
    timeSinceEpoch = timenow.time_since_epoch();
    timemillisec = sc::duration_cast<sc::milliseconds>(timeSinceEpoch);
    starttime = timemillisec.count();
    waveletTree wt = waveletTree(bwt, bwt + n, 64, maximum);
    timenow = sc::system_clock::now();
    timeSinceEpoch = timenow.time_since_epoch();
    timemillisec = sc::duration_cast<sc::milliseconds>(timeSinceEpoch);
    endtime = timemillisec.count();
    cout << "Wavelet Tree Creation Time: " << endtime-starttime << " milliseconds.\n";

    computeCharacterRank(C, SA, s, n);
    cout << "Backward Pattern Matching\n";
    string pattern;
    cout << "Enter Pattern: ";
    cin >> pattern;

    int P[pattern.length()];
    for (int i = 0; i < pattern.length(); i++)
    {
        P[i] = (int)pattern[i];
    }

    timenow = sc::system_clock::now();
    timeSinceEpoch = timenow.time_since_epoch();
    timemillisec = sc::duration_cast<sc::milliseconds>(timeSinceEpoch);
    starttime = timemillisec.count();

    int flag = 0;
    flag = patternMatching(pattern.length(), P, C, startpoint, endpoint, wt);

    if(flag == 0)
    {
        cout << "No Pattern Matches.\n";
    }
    else
    {
        cout << flag << " Patterns found.\n";
    }

    timenow = sc::system_clock::now();
    timeSinceEpoch = timenow.time_since_epoch();
    timemillisec = sc::duration_cast<sc::milliseconds>(timeSinceEpoch);
    endtime = timemillisec.count();

    cout << "Pattern Matching Time: " << endtime-starttime << " milliseconds.\n";

    delete s; delete SA; delete bwt;

    // Calculate and print total space usage
    int alphabetSize = 128; // ASCII size
    printTotalSpace(textLength, alphabetSize);

    return 0;

}
