//#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>

using namespace std;

int n, p, d, k, rr;
vector<int> a;
vector<int> a_rev;
vector<int> g;

void construct_alphas() {
    a.push_back(1);
    a_rev.push_back(0);

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


void construct_g_polynome() {
    construct_alphas();

    vector<int> coefs[d]; // степени a при x^i
    vector<int> empty;
    for (int j = 0; j < d; j++) {
        coefs[j] = empty;
    }
    coefs[0].push_back(1);
    coefs[1].push_back(0);

    for (int j = 2; j < d; j++) {
        // умножаем на (x + a^j);

        int curr_max_pow = -1;
        vector<vector<int>> mul_x;
        for (int i = 0; i < d; i++) {
            mul_x.push_back(empty);
        }

        // умножили на x
        for (int q = d - 1; q >= 0; q--) {
            if (curr_max_pow == -1 && !coefs[q].empty()) {
                curr_max_pow = q;
            }
            if (q != 0) {
                mul_x[q] = coefs[q - 1];
            }
        }

        // умножили на a^j
        for (int q = 0; q <= curr_max_pow; q++) {
            for (int ii = 0; ii < coefs[q].size(); ii++) {
                coefs[q][ii] += j;
                coefs[q][ii] %= n;
            }
        }

        // сложили два результата
        for (int q = 0; q < d; q++) {
            for (int ii = 0; ii < mul_x[q].size(); ii++) {
                coefs[q].push_back(mul_x[q][ii]);
            }
        }

        for (int q = 0; q < d; q++) {

            int now = 0;
            for (int ii = 0; ii < coefs[q].size(); ii++) {
                now = now ^ a[coefs[q][ii]];
            }

            coefs[q].clear();
            if (now != 0) {
                coefs[q].push_back(a_rev[now]);
            }
        }
    }

    for (int j = 0; j < d; j++) {
        // получаем коефф при x^j

        int now = 0;
        for (int q = 0; q < coefs[j].size(); q++) {
            now = now ^ a[coefs[j][q]];
        }

        g.push_back(now);
    }
}

int mul_in_module(int i, int j) {
    if (i * j == 0) {
        return 0;
    }
    return a[(a_rev[i] + a_rev[j]) % n];
}

int div_in_module(int up, int down) {
    if (up * down == 0) {
        return 0;
    }
    return a[(a_rev[up] - a_rev[down] + n) % n];
}

vector<int> clean_polynome(vector<int>& v) {
    int cnt = 0;
    for (int i = v.size() - 1; i >= 0; i--) {
        if (v[i] != 0) {
            break;
        }
        if (i != 0) {
            cnt++;
        }
    }
    if (cnt != 0) {
        v.erase(v.end() - cnt, v.end());
    }

    return v;
}

pair<vector<int>, vector<int>> div(const vector<int>& up, const vector<int>& down) {
    vector<int> reminder = up;
    vector<int> res;

    for (int i = (int) up.size() - 1; i >= (int) down.size() - 1; i--) {
        if (reminder[i] > 0) {
            int curr = div_in_module(reminder[i], down[(int) down.size() - 1]);
            for (int j = 0; j < down.size(); j++) {
                reminder[i - j] = reminder[i - j] ^ mul_in_module(curr, down[(int) down.size() - 1 - j]);
            }
            res.push_back(curr);
        } else {
            res.push_back(0);
        }

        if (i == 0) {
            break;
        }
    }

    reverse(res.begin(), res.end());

    return {clean_polynome(res), clean_polynome(reminder)};
}

vector<int> mul(const vector<int>& aa, const vector<int>& b) {
    vector<int> res;
    for (int i = 0; i < aa.size() + b.size() - 1; i++) {
        res.push_back(0);
    }

    for (int i = 0; i < aa.size(); i++) {
        if (aa[i] > 0) {
           for (int j = 0; j < b.size(); j++) {
               res[i + j] = res[i + j] ^ mul_in_module(aa[i], b[j]);
           }
        }
    }

    return clean_polynome(res);
}

vector<int> encode(const vector<int>& v) {
    vector<int> res;
    for (int i = 0; i < rr; i++) {
        res.push_back(0);
    }
    for (int i = 0; i < k; i++) {
        res.push_back(v[i]);
    }

    vector<int> reminder = div(res, g).second;

    for (int i = 0; i < rr; i++) {
        int reminder_val = i >= reminder.size() ? 0 : reminder[i];
        res[i] = res[i] ^ reminder_val;
    }

    return res;
}

vector<int> calc_syndrome(const vector<int>& v) {
    vector<int> s;
    for (int j = 1; j < d; j++) {

        int curr = 0;
        for (int i = 0; i < n; i++) {
            curr = curr ^ mul_in_module(v[i], a[(i * j) % n]);
        }

        s.push_back(curr);
    }

    return s;
}

vector<int> decode(const vector<int>& v) {
    vector<int> s = calc_syndrome(v);

    vector<int> r_pre_last;
    vector<int> r_last = s;
    vector<int> a_pre_last;
    vector<int> a_last;

    for (int i = 0; i < rr; i++) {
        r_pre_last.push_back(0);
    }
    r_pre_last.push_back(1);

    a_last.push_back(1);
    a_pre_last.push_back(0);

    while (r_last.size() - 1 >= rr / 2) {
        auto res = div(r_pre_last, r_last);
        vector<int> q = res.first;
        vector<int> r_now = res.second;
        vector<int> a_now;

        vector<int> second = mul(q, a_last);
        for (int i = 0; i < max(second.size(), a_pre_last.size()); i++) {
            int a_val = i >= a_pre_last.size() ? 0 : a_pre_last[i];
            int second_val = i >= second.size() ? 0 : second[i];
            a_now.push_back(a_val ^ second_val);
        }

        r_pre_last = r_last;
        r_last = r_now;
        a_pre_last = a_last;
        a_last = a_now;
    }

    vector<int> val;
    for (int i = 0; i < a_last.size(); i++) {
        if (a_last[i] != 0) {
            val.push_back(a_last[i]);
            break;
        }
    }

    vector<int> locators = div(a_last, val).first;

    vector<int> x_2t;
    for (int i = 0; i < d - 1; i++) {
        x_2t.push_back(0);
    }
    x_2t.push_back(1);

    vector<int> omega_without_module = mul(s, locators);
    vector<int> omega = div(omega_without_module, x_2t).second;

    // ищем все a^(-i) == 0
    vector<int> err_in;
    for (int i = 0; i < n; i++) {
        int curr = locators[0];
        for (int j = 1; j < locators.size(); j++) {
            curr = curr ^ mul_in_module(locators[j], a[(i * j) % n]);
        }

        if (curr == 0) {
            err_in.push_back((n - i) % n);
        }
    }

    //aлгоритм Форни
    vector<int> res = v;
    for (int i = 0; i < err_in.size(); i++) {
        int x_rev = div_in_module(1, a[err_in[i]]);

        int up = 0;
        int now_pow_val = 1;
        for (int j = 0; j < omega.size(); j++) {
            up = up ^ mul_in_module(omega[j], now_pow_val);
            now_pow_val = mul_in_module(now_pow_val, x_rev);
        }
        up = mul_in_module(x_rev, up);

        int down = 1;
        for (int j = 0; j < err_in.size(); j++) {
            if (err_in[i] != err_in[j]) {
                down = mul_in_module(down, (1 ^ mul_in_module(a[err_in[j]], x_rev)));
            }
        }

        res[err_in[i]] = res[err_in[i]] ^ div_in_module(up, down);
    }

    return res;
}

double get_rand() {
    return (double) rand() / (RAND_MAX);
}

int get_rand_int_to_n() {
    return (int) round(get_rand() * n) % n;
}

vector<int> gen_word() {
    vector<int> res;

    for (int i = 0; i < k; i++) {
        int rand_gen = get_rand_int_to_n();
        res.push_back(rand_gen);
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

        vector<int> word = gen_word();
        vector<int> encoded_word = encode(word);

        vector<int> converted_word;
        for (int i = 0; i < n; i++) {
            converted_word.push_back(get_rand() <= noise_lvl ? encoded_word[i] ^ a[get_rand_int_to_n()] : encoded_word[i]);
        }

        vector<int> decoded_word = decode(converted_word);

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
    ifstream in("/Users/st1580/CLionProjects/professional/input.txt");
    ofstream out("/Users/st1580/CLionProjects/professional/output.txt");
//    ifstream in("input.txt");
//    ofstream out("output.txt");

    ios_base::sync_with_stdio(false);
    in.tie(0);
    out.tie(0);

    in >> n >> p >> d;

    construct_g_polynome();

    rr = g.size() - 1;
    k = n - rr;

    out << k << "\n";
    for (int i = 0; i <= rr; i++) {
        out << g[i] << " ";
    }

    string s;
    while (in >> s) {
        out << "\n";

        if (s[0] == 'E') {
            vector<int> v;
            for (int i = 0; i < k; i++) {
                int x;
                in >> x;
                v.push_back(x);
            }

            for (int val : encode(v)) {
                out << val << " ";
            }

        } else if (s[0] == 'D') {
            vector<int> v;
            for (int i = 0; i < n; i++) {
                int x;
                in >> x;
                v.push_back(x);
            }

            for (int val : decode(v)) {
                out << val << " ";
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
