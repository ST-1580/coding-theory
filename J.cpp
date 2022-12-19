#include <bits/stdc++.h>
//#include <iostream>
//#include <fstream>
//#include <vector>
//#include <algorithm>
//#include <random>
//#include <cmath>

using namespace std;

int n, p, d, k, rr;
vector<vector<int>> c;
vector<int> a;
vector<int> a_rev;
vector<pair<int, vector<bool>>> m; // index, mask
vector<bool> g;

void construct_cyclic_classes() {
    bool used[n];

    for (int i = 0; i < n; i++) {
        used[i] = false;
    }
    used[0] = true;

    int i = 1;
    while (i < n) {
        if (i >= d) {
            return;
        }

        vector<int> curr;

        int j = i;
        while (!used[j]) {
            curr.push_back(j);
            used[j] = true;
            j = (2 * j) % n;
        }

        if (!curr.empty()) {
            c.push_back(curr);
        }
        i++;
    }
}

void construct_alphas() {
    a.push_back(1);
    a_rev.push_back(-1);

    for (int i = 1; i <= n; i++) {
        int next = (a[i - 1] << 1);
        if (next > n) {
            next = (next ^ p) & n;
        }

        a.push_back(next);
        a_rev.push_back(-1);
    }

    for (int i = 0; i <= n; i++) {
        a_rev[a[i]] = i;
    }
}

void construct_min_polynomes() {
    construct_alphas();

    for (int i = 0; i < c.size(); i++) {
        int max_power = c[i].size();
        vector<int> coefs[max_power + 1]; // степени a при x^i
        vector<int> empty;
        for (int j = 0; j <= max_power; j++) {
            coefs[j] = empty;
        }
        coefs[0].push_back(c[i][0]);
        coefs[1].push_back(0);

        for (int j = 1; j < max_power; j++) {
            // умножаем на (x + a^c[i][j]);

            int curr_max_pow = -1;
            vector<vector<int>> mul_x;
            for (int i = 0; i <= max_power; i++) {
                mul_x.push_back(empty);
            }

            // умножили на x
            for (int q = max_power; q >= 0; q--) {
                if (curr_max_pow == -1 && !coefs[q].empty()) {
                    curr_max_pow = q;
                }
                if (q != 0) {
                    mul_x[q] = coefs[q - 1];
                }
            }

            // умножили на a^c[i][j]
            for (int q = 0; q <= curr_max_pow; q++) {
                for (int ii = 0; ii < coefs[q].size(); ii++) {
                    coefs[q][ii] += c[i][j];
                    coefs[q][ii] %= n;
                }
            }

            // сложили два результата
            for (int q = 0; q <= max_power; q++) {
                for (int ii = 0; ii < mul_x[q].size(); ii++) {
                    coefs[q].push_back(mul_x[q][ii]);
                }
            }
        }

        vector<bool> mask;
        int index = c[i][0];

        for (int j = 0; j <= max_power; j++) {
            // получаем коефф при x^j

            int now = 0;
            for (int q = 0; q < coefs[j].size(); q++) {
                now = now ^ a[coefs[j][q]];
            }

            assert(now == 1 || now == 0);

            mask.push_back(now == 1);
        }

        m.emplace_back(index, mask);
    }

}

void construct_g_polynome() {
    construct_cyclic_classes();
    construct_min_polynomes();

    // умножение многочленов в столбик
    g.push_back(true);

    for (int i = 0; i < m.size(); i++) {
        if (m[i].first >= d) {
            break;
        }
        vector<bool> mask = m[i].second;

        vector<bool> curr;
        for (int ii = 0; ii < mask.size() + g.size() - 1; ii++) {
            curr.push_back(false);
        }

        for (int im = 0; im < mask.size(); im++) {
            if (mask[im]) {
                for (int ig = 0; ig < g.size(); ig++) {
                    curr[im + ig] = curr[im + ig] ^ g[ig];
                }
            }
        }

        g = curr;
    }
}

vector<bool> encode(const vector<bool>& v) {
    vector<bool> res;
    for (int i = 0; i < rr; i++) {
        res.push_back(false);
    }
    for (int i = 0; i < k; i++) {
        res.push_back(v[i]);
    }

    // деление в столбик
    vector<bool> reminder = res;
    for (int i = res.size() - 1; i >= rr; i--) {
        if (reminder[i]) {
            for (int j = 0; j <= rr; j++) {
                reminder[i - j] = reminder[i - j] ^ g[rr - j];
            }
        }
    }

    for (int i = 0; i < rr; i++) {
        res[i] = res[i] ^ reminder[i];
    }

    return res;
}

vector<int> calc_syndrome(const vector<bool>& v) {
    vector<int> s;
    for (int j = 1; j <= d - 1; j++) {

        int curr = 0;
        for (int i = 0; i < n; i++) {
            if (v[i]) {
                curr = curr ^ a[(i * j) % n];
            }
        }

        s.push_back(curr);
    }

    return s;
}

int calc_in_module(int i, int j) {
    if (i * j == 0) {
        return 0;
    }
    return a[(a_rev[i] + a_rev[j]) % n];
}

vector<bool> decode(const vector<bool>& v) {
    vector<int> s = calc_syndrome(v);

    int r = 1;
    int mm = 0;
    int L = 0;
    vector<int> locators;
    locators.push_back(1);
    vector<int> b;
    b.push_back(1);

    while (r <= d - 1) {
        int delta = 0;
        for (int j = 0; j <= L; j++) {
            delta = delta ^ calc_in_module(locators[j], s[r - j - 1]);
        }

        if (delta != 0) {
            vector<int> t = locators;
            for (int i = locators.size(); i < r - mm + b.size(); i++) {
                t.push_back(0);
            }
            for (int i = 0; i < b.size(); i++) {
                t[r - mm + i] = t[r - mm + i] ^ calc_in_module(delta, b[i]);
            }

            if (2 * L <= r - 1) {
                b.clear();
                int delta_rev = a[(n - a_rev[delta]) % n];
                for (int i = 0; i < locators.size(); i++) {
                    b.push_back(calc_in_module(delta_rev, locators[i]));
                }

                L = r - L;
                mm = r;
            }

            locators = t;
        }

        r++;
    }

    if (L != locators.size() - 1) {
        return v;
    }

    // ищем все a^(-i) == 0
    vector<int> err_in;
    for (int i = 0; i < n; i++) {
        int curr = locators[0];
        for (int j = 1; j < locators.size(); j++) {
            curr = curr ^ calc_in_module(locators[j], a[(i * j) % n]);
        }

        if (curr == 0) {
            err_in.push_back(i);
        }
    }

    //a^(-i) = a^j => (n - i) % n = j
    vector<bool> res = v;
    for (int i = 0; i < err_in.size(); i++) {
        res[(n - err_in[i]) % n] = !res[(n - err_in[i]) % n];
    }

    return res;
}

double get_rand() {
    return (double) rand() / (RAND_MAX);
}

vector<bool> gen_word() {
    vector<bool> res;

    for (int i = 0; i < rr; i++) {
        int rand_gen = round(get_rand());
        res.push_back(rand_gen == 1);
    }

    return res;
}

double simulate(double noise_lvl, int num_of_operations, int max_errors) {
    int curr_operations = 0;
    int curr_errors = 0;

    random_device rd{};
    mt19937 gen{rd()};
    normal_distribution<> nd{0, sqrt(0.5 * pow(10,  -noise_lvl / 10.0) * ((double) n / k))};

    while (curr_operations < num_of_operations && curr_errors < max_errors) {
        curr_operations++;

        vector<bool> word = gen_word();
        vector<bool> encoded_word = encode(word);

        vector<bool> converted_word;
        for (int i = 0; i < n; i++) {
            converted_word.push_back(get_rand() <= noise_lvl ? !encoded_word[i] : encoded_word[i]);
        }

        vector<bool> decoded_word = decode(converted_word);

        for (int i = 0; i < encoded_word.size(); i++) {
            if (encoded_word[i] != decoded_word[i]) {
                curr_errors++;
                break;
            }
        }

    }

    return (double) curr_errors / curr_operations;
}

int main() {
//    ifstream in("/Users/st1580/CLionProjects/professional/input.txt");
//    ofstream out("/Users/st1580/CLionProjects/professional/output.txt");
    ifstream in("input.txt");
    ofstream out("output.txt");

    ios_base::sync_with_stdio(false);
    in.tie(0);
    out.tie(0);

    in >> n >> p >> d;

    construct_g_polynome();

    rr = g.size() - 1;
    k = n - rr;

    out << k << "\n";
    for (int i = 0; i <= rr; i++) {
        out << (g[i] ? 1 : 0) << " ";
    }

    string s;
    while (in >> s) {
        out << "\n";

        if (s[0] == 'E') {
            vector<bool> v;
            for (int i = 0; i < k; i++) {
                int x;
                in >> x;
                v.push_back(x == 1);
            }

            for (bool val : encode(v)) {
                out << (val ? 1 : 0) << " ";
            }

        } else if (s[0] == 'D') {
            vector<bool> v;
            for (int i = 0; i < n; i++) {
                int x;
                in >> x;
                v.push_back(x == 1);
            }

            for (bool val : decode(v)) {
                out << (val ? 1 : 0) << " ";
            }

        } else {
            double noise_lvl;
            int num_of_operations, max_errors;

            in >> noise_lvl >> num_of_operations >> max_errors;
            out << setprecision(15) << simulate(noise_lvl, num_of_operations, max_errors);
        }
    }

    in.close();
    out.close();

    return 0;
}
