#include <iostream>
#include <vector>

class phi {};

double eps = 1e-16;

// как мне как-то проще определить как считать произвоидные для j, k. Можно сделать вектор из
// phi

// i - это номер уравнения по сути (phi в произведении, которое не под дифференциалом), j и k -
// это phi которые под дифференциалом.
double skalar_of_differential(int i, int j, int k, int M) {
    if (i == 0) {
        if (k == 0 && j == 0)
            return -2. / 3;
        else if ((k == 0 && j == 1) || (k == 1 && j == 0))
            return 1. / 6;
        else if (k == 1 && j == 1)
            return 1. / 3;
    }

    if (i == M) {
        if (k == M && j == M)
            return 2. / 3;
        else if ((k == M && j == M - 1) || (k == M - 1 && j == M))
            return -1. / 6;
        else if (k == M - 1 && j == M - 1)
            return -1. / 3;
    }

    if (j == i - 1) {
        if (k == i - 1)
            return -1. / 3;
        else if (k == i)
            return -1. / 6;
        else if (k == i + 1)
            return 0;
        else
            return 0;
    } else if (j == i) {
        if (k == i - 1)
            return -1. / 6;
        else if (k == i)
            return 0;
        else if (k == i + 1)
            return 1. / 6;
        else
            return 0;
    } else if (j == i + 1) {
        if (k == i - 1)
            return 0;
        else if (k == i)
            return 1. / 6;
        else if (k == i + 1)
            return 1. / 3;
        else
            return 0;
    } else
        return 0;
}

double skalar_without_differential(int i, int j, double h, int M) {
    if (j == i + 1)
        if (i == M)
            return 0;
        else
            return (1. / 6) * h;

    else if (j == i)
        if (i == 0 || i == M)
            return (1. / 3) * h;
        else
            return (2. / 3) * h;

    else if (j == i - 1)
        if (i == 0)
            return 0;
        else
            return (1. / 6) * h;

    else
        return 0;
}

class func {
   public:
    std::vector<double> c_i;
    double m_tau, m_h;
    int m_M;  // это количество точек на 1 больше.

    func(double tau, double h, int M) {
        double default_value = 1;
        m_tau = tau;
        m_h = h;
        m_M = M;
        for (int i = 0; i <= m_M; i++) {
            c_i.push_back(default_value);
        }
    }

    func(std::vector<double> vec, double tau, double h, int M) {
        if (static_cast<int>(vec.size()) != M + 1) {
            std::cout << "Это какой-то неправильный вектор. Он имеет размер " << vec.size()
                      << ", а должен " << M + 1 << " \n";
        }
        m_tau = tau;
        m_h = h;
        m_M = M;
        for (int i = 0; i <= m_M; i++) {
            c_i.push_back(vec[i]);
        }
    }

    ~func() = default;

    double approx(double x) {
        int index = static_cast<int>(x / m_h);
        if (index == m_M) return c_i[m_M];
        double buffer = x - index * m_h;
        return c_i[index] * (1 - buffer) + c_i[index + 1] * buffer;
    }
};

struct equation {
    int num_of_eq = -1;
    double below_diag;
    double onthe_diag;
    double above_diag;
    double right_part;

    int N, M;
    double tau, h;

    void normalize_eq() {
        if (!(std::abs(below_diag) < eps)) {
            std::cout << "Не удалось нормировать уравнение " << num_of_eq
                      << ". Коэффициент под диагональю равен" << below_diag << "\n";
            return;
        } else {
            above_diag = above_diag / onthe_diag;
            right_part = right_part / onthe_diag;
            onthe_diag = 1.;
        }
    }

    void init_equation(const func& u /* значения прошлой функции в точках интерполяции */,
                       int num /* номер уравнения от 0 до M*/, int M /* количество отрезков */) {
        auto function_of_coefficient = [&u, M](int i, int j) {
            if (i == 0)
                return u.c_i[i] * skalar_of_differential(i, j, i, M) +
                       u.c_i[i + 1] * skalar_of_differential(i, j, i + 1, M);

            else if (i == M)
                return u.c_i[i - 1] * skalar_of_differential(i, j, i - 1, M) +
                       u.c_i[i] * skalar_of_differential(i, j, i, M);

            else
                return u.c_i[i - 1] * skalar_of_differential(i, j, i - 1, M) +
                       u.c_i[i] * skalar_of_differential(i, j, i, M) +
                       u.c_i[i + 1] * skalar_of_differential(i, j, i + 1, M);
        };

        if (num == 0) {
            num_of_eq = num;
            below_diag = 0;
            onthe_diag = (1. / tau) * skalar_without_differential(num, num, h, M) +
                         function_of_coefficient(num, num);
            above_diag = (1. / tau) * skalar_without_differential(num, num + 1, h, M) +
                         function_of_coefficient(num, num + 1);
            right_part =
                (1. / tau) * (skalar_without_differential(num, num, h, M) * u.c_i[num] +
                              skalar_without_differential(num, num + 1, h, M) * u.c_i[num + 1]);
            ;
        } else if (num == M) {
            num_of_eq = num;
            below_diag = (1. / tau) * skalar_without_differential(num, num - 1, h, M) +
                         function_of_coefficient(num, num - 1);
            onthe_diag = (1. / tau) * skalar_without_differential(num, num, h, M) +
                         function_of_coefficient(num, num);
            above_diag = 0;
            right_part =
                (1. / tau) * (skalar_without_differential(num, num - 1, h, M) * u.c_i[num - 1] +
                              skalar_without_differential(num, num, h, M) * u.c_i[num]);
        } else {
            num_of_eq = num;
            below_diag = (1. / tau) * skalar_without_differential(num, num - 1, h, M) +
                         function_of_coefficient(num, num - 1);
            onthe_diag = (1. / tau) * skalar_without_differential(num, num, h, M) +
                         function_of_coefficient(num, num);
            above_diag = (1. / tau) * skalar_without_differential(num, num + 1, h, M) +
                         function_of_coefficient(num, num + 1);
            right_part =
                (1. / tau) * (skalar_without_differential(num, num - 1, h, M) * u.c_i[num - 1] +
                              skalar_without_differential(num, num, h, M) * u.c_i[num] +
                              skalar_without_differential(num, num + 1, h, M) * u.c_i[num + 1]);
        }
    }
};

void eq_below_is_eq_below_minus_eq_above(equation& eq2 /* это уравнение внизу */,
                                         equation& eq1 /* это уравнение сверху */) {
    if (eq1.num_of_eq == -1 || eq1.num_of_eq == -1) {
        std::cout << "Одно из уравнений неинициализировано:  \n  eq1 = " << eq1.num_of_eq
                  << ", eq2 = " << eq2.num_of_eq << "\n";
        return;
    }
    if (eq1.num_of_eq - eq2.num_of_eq != 1 && eq1.num_of_eq - eq2.num_of_eq != -1) {
        if (eq1.num_of_eq == eq2.num_of_eq)
            std::cout << "Ты вычитаешь уравнение из самого себя\n";
        else
            std::cout << "Вероятно есть ошибка, потому что вычитается не предыдущее и не следующее "
                         "уравнение\n";
        return;
    } else {
        if (std::abs(eq1.below_diag) < eps) {
            eq2.below_diag -= eq1.onthe_diag;
            eq2.onthe_diag -= eq1.above_diag;
            eq2.right_part -= eq1.right_part;
            return;
        } else {
            std::cout << "Коэффициент под диагональью уравнения " << eq1.num_of_eq
                      << " ненулевой. Он равен " << eq1.below_diag << "\n";
            return;
        }
    }
}

func solver_of_matrix(std::vector<equation>& matrix, double tau, double h,
                      int M /* это число отрезков, точек же на одну больше */) {
    if (static_cast<int>(matrix.size()) != M + 1)
        std::cout << "Короче, неправильно заданы уравнения (их должно быть M + 1)\n";
    /* номер последнего уравнения равен M */

    for (int i = 0; i < M; i++) {
        matrix[i].normalize_eq();
        eq_below_is_eq_below_minus_eq_above(matrix[i + 1], matrix[i]);
    }

    matrix[M].normalize_eq();

    for (int i = M; i > 0; i++) {
        matrix[i].normalize_eq();
        eq_below_is_eq_below_minus_eq_above(matrix[i - 1], matrix[i]);
    }

    matrix[0].normalize_eq();

    std::vector<double> next_approx;
    for (int i = 0; i <= M; i++) {
        next_approx.push_back(matrix[i].right_part);
    }

    return func(next_approx, tau, h, M);
}

int main() {
    try {
        std::cout << "Начало курсовой программы" << std::endl;

        // задание матрицы. Нужно задать векторы для c_j, c_j_prev (это с прошлого шага, его
        // сохраняем), но их мы ищем просто.

        double X = 1., T = 1.;   // длина сегментов по пространству и по времени соответственно
        int N = 1000, M = 1000;  // M - это шаги по пространству, N - это шаги по времени (точек на
                                 // одну больше, так как считаем с 0)
        double tau = T / N, h = X / M;

        func u(tau, h, N);
        std::vector<equation> matrix;
        // задаем матрицу, она должна быть трехдиагональной. Надо создать структуру из векторов
        // коэффициентов уравнений Имеется n+1 уравнений.

    } catch (const std::exception& e) {
        std::cerr << "Стандартное исключение: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Неизвестный тип исключения!" << std::endl;
    }
    return 0;
}