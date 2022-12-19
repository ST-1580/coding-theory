#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

int n, k;
vector<vector<bool>> g;

vector<bool> encode(const bool v[k]) {
    vector<bool> res;

    for (int j = 0; j < n; j++) {
        res.push_back(false);
        for (int i = 0; i < k; i++) {
            res[j] = res[j] ^ (v[i] & g[i][j]);
        }
    }

    return res;
}

vector<bool> decode(double v[n]) {

}

void simulate(double noise_lvl, int num_of_operations, int max_errors) {

}

int main() {
    ifstream in("/Users/st1580/CLionProjects/professional/input.txt");
    ofstream out("/Users/st1580/CLionProjects/professional/output.txt");

    ios_base::sync_with_stdio(false);
    in.tie(0);
    out.tie(0);

    in >> n >> k;

    vector<bool> empty;
    for (int i = 0; i < k; i++) {
        g.push_back(empty);
        for (int j = 0; j < n; j++) {
            int x;
            in >> x;
            g[i].push_back(x == 1);
        }
    }


    string s;
    while (in >> s) {
        if (s[0] == 'E') {
            bool v[k];
            for (int i = 0; i < k; i++) {
                int x;
                in >> x;
                v[i] = x == 1;
            }

            for (bool val : encode(v)) {
                out << (val ? 1 : 0) << " ";
            }

        } else if (s[0] == 'D') {
            double v[n];
            for (int i = 0; i < n; i++) {
                double x;
                in >> x;
                v[i] = x;
            }

            for (bool val : decode(v)) {
                out << (val ? 1 : 0) << " ";
            }

        } else {
            double noise_lvl;
            int num_of_operations, max_errors;

            in >> noise_lvl >> num_of_operations >> max_errors;

            simulate(noise_lvl, num_of_operations, max_errors);
        }
    }

    in.close();
    out.close();

    return 0;
}
