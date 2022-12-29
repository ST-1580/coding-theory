#include <bits/stdc++.h>

using namespace std;

int n, p, d, k, rr;
vector<int> a;      // предподсчитанные alphas
vector<int> a_rev;  // вектор, a_rev[i] = x => a^x = i
vector<int> g;      // порождающий многочлен

// строим альфы
void construct_alphas() {
    a.push_back(1);
    a_rev.push_back(0);

    // a[i] = (a[i - 1] * x) % p
    for (int i = 1; i <= n; i++) {
        int next = (a[i - 1] << 1);
        if (next > n) {
            // взятие по модулю
            // поскольку n = 2^m - 1, то его битовая маска - все единицы
            // поэтому можно "обрубить" ненужные биты при помощи логического и
            next = (next ^ p) & n;
        }

        // заполняем a_rev
        a.push_back(next);
        a_rev.push_back(-1);
    }

    for (int i = 0; i <= n; i++) {
        a_rev[a[i]] = i;
    }
}


void construct_g_polynome() {
    construct_alphas();

    // степени aльф при x^i
    vector<int> coefs[d];
    // вспомогательный пустой массив
    vector<int> empty;
    for (int j = 0; j < d; j++) {
        coefs[j] = empty;
    }

    // проиницилизируем коэффициенты первым множителем (x + a^1)
    coefs[0].push_back(1);
    coefs[1].push_back(0);

    for (int j = 2; j < d; j++) {
        // умножаем на (x + a^j)
        // будем делать это в три этапа
        // 1 - умножим на х
        // 2 - умножим на a^j
        // 3 - сложим полученные значения

        // текущая максимальная степень
        int curr_max_pow = -1;

        // вспомогательный вектор для умножения на х
        vector<vector<int>> mul_x;
        for (int i = 0; i < d; i++) {
            mul_x.push_back(empty);
        }

        // умножим на x, для этого просто сдвинем все коэффициенты
        for (int q = d - 1; q >= 0; q--) {
            if (curr_max_pow == -1 && !coefs[q].empty()) {
                curr_max_pow = q;
            }
            if (q != 0) {
                mul_x[q] = coefs[q - 1];
            }
        }

        // умножим на a^j
        for (int q = 0; q <= curr_max_pow; q++) {
            for (int ii = 0; ii < coefs[q].size(); ii++) {
                coefs[q][ii] += j;
                coefs[q][ii] %= n;
            }
        }

        // сложим два результата
        for (int q = 0; q < d; q++) {
            for (int ii = 0; ii < mul_x[q].size(); ii++) {
                coefs[q].push_back(mul_x[q][ii]);
            }
        }

        // для ускорения работы сократим степени альф
        for (int q = 0; q < d; q++) {
            // текущий результат
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

    // собираем порождающий многочлен
    for (int j = 0; j < d; j++) {

        // получаем коефф при x^j
        int now = 0;
        for (int q = 0; q < coefs[j].size(); q++) {
            now = now ^ a[coefs[j][q]];
        }

        g.push_back(now);
    }
}

// вспомогательная функция для умножения в поле
int mul_in_module(int i, int j) {
    if (i * j == 0) {
        return 0;
    }
    return a[(a_rev[i] + a_rev[j]) % n];
}

// вспомогательная функция для деления в поле
// up - числитель
// down - знаменатель
int div_in_module(int up, int down) {
    if (up * down == 0) {
        return 0;
    }
    return a[(a_rev[up] - a_rev[down] + n) % n];
}

// вспомогательная функция для стирания ведущих нулей в коэффициентах многочлена
vector<int> clean_polynome(vector<int>& v) {
    // число нулей
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

// деление многочленов в столбик
// вернет результат и остаток
// up - числитель
// down - знаменатель
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

// умножение многочленов в столбик
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
    // res = v * x^rr + reminder
    vector<int> res;
    for (int i = 0; i < rr; i++) {
        res.push_back(0);
    }
    for (int i = 0; i < k; i++) {
        res.push_back(v[i]);
    }

    // поиск остатка
    vector<int> reminder = div(res, g).second;

    // прибавление остатка
    for (int i = 0; i < rr; i++) {
        int reminder_val = i >= reminder.size() ? 0 : reminder[i];
        res[i] = res[i] ^ reminder_val;
    }

    return res;
}

// подсчет синдрома
vector<int> calc_syndrome(const vector<int>& v) {
    vector<int> s;
    for (int j = 1; j < d; j++) {
        // s_j = v(a^j)
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

    // инициализация
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

    // алгоритм Сугиямы (Евклида)
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

    // поиск коэффициента не равного 0
    vector<int> val;
    for (int i = 0; i < a_last.size(); i++) {
        if (a_last[i] != 0) {
            val.push_back(a_last[i]);
            break;
        }
    }

    // нахождение многочлена локаторов
    vector<int> locators = div(a_last, val).first;

    // вектор x^2t
    vector<int> x_2t;
    for (int i = 0; i < d - 1; i++) {
        x_2t.push_back(0);
    }
    x_2t.push_back(1);

    // нахождение омеги
    vector<int> omega_without_module = mul(s, locators);
    vector<int> omega = div(omega_without_module, x_2t).second;

    // поиск ошибок: ищем все a^(-i) == 0
    vector<int> err_in;
    for (int i = 0; i < n; i++) {
        int curr = locators[0];
        for (int j = 1; j < locators.size(); j++) {
            curr = curr ^ mul_in_module(locators[j], a[(i * j) % n]);
        }

        if (curr == 0) {
            // необходимо перевести -i в j
            // a^(-i) = a^j => (n - i) % n = j
            err_in.push_back((n - i) % n);
        }
    }

    // aлгоритм Форни
    vector<int> res = v;
    for (int i = 0; i < err_in.size(); i++) {
        // инвертируем
        int x_rev = div_in_module(1, a[err_in[i]]);

        // числитель
        int up = 0;
        int now_pow_val = 1;
        for (int j = 0; j < omega.size(); j++) {
            up = up ^ mul_in_module(omega[j], now_pow_val);
            now_pow_val = mul_in_module(now_pow_val, x_rev);
        }
        up = mul_in_module(x_rev, up);

        // знаменатель
        int down = 1;
        for (int j = 0; j < err_in.size(); j++) {
            if (err_in[i] != err_in[j]) {
                down = mul_in_module(down, (1 ^ mul_in_module(a[err_in[j]], x_rev)));
            }
        }

        // исправляем ошибки
        res[err_in[i]] = res[err_in[i]] ^ div_in_module(up, down);
    }

    return res;
}

// генерация рандома от 0 до 1
double get_rand() {
    return (double) rand() / (RAND_MAX);
}

// генерация рандома от 0 до n
int get_rand_int_to_n() {
    return (int) round(get_rand() * n) % n;
}

// генерация вектора для encode
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
            // случайное изменение бита
            converted_word.push_back(get_rand() <= noise_lvl ? encoded_word[i] ^ a[get_rand_int_to_n()] : encoded_word[i]);
        }

        vector<int> decoded_word = decode(converted_word);

        // проверка
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