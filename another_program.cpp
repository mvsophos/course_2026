#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

class phi {};

// надо нормальную функцию написать. Может быть вообще сделать сначала чисто метод конечных
// элементов на один элемент, то есть на кусочно-постоянных функциях

struct common_params {
    double tau, h;
    int N, M;
    // M - это число отрезков, а точек на одну больше

    common_params(int _N, double _tau, int _M, double _h) : N(_N), tau(_tau), M(_M), h(_h) {}

    ~common_params() = default;
};

double eps = 1e-12;

int function_type = 0;

double accuracy_right_part(double t, double x) {
    switch (function_type) {
        case 0:
            return std::exp(t) * (1.5 + std::cos(3 * M_PI * x)) *
                   (1 - 6 * M_PI * std::sin(3 * M_PI * x));
    }
    return 0;
}

double accuracy_u(double t, double x) {
    return std::exp(t) * (1.5 + std::cos(3 * M_PI * x));
}

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

// Интеграл по Симпсону. Это нужно для рассчета правой части. Она берет интеграл
double simpson_integral(double (*accuracy_func)(double /* time */, double /* spaceX */), int i,
                        double time, const common_params& params) {
    double value = 0;
    if (i != params.M)  // берем интеграл с \phi^{+}
        value += accuracy_func(time, i * params.h) + 2 * accuracy_func(time, (i + 0.5) * params.h);
    if (i != 0)  // берем интеграл с \phi^{-}
        value += 2 * accuracy_func(time, (i - 0.5) * params.h) + accuracy_func(time, i * params.h);
    return (params.h / 6) * value;
}

class func {
   public:
    std::vector<double> c_i;
    const common_params& m_params;

    func(const common_params& parameters) : m_params(parameters) {
        double default_value = 1;
        for (int i = 0; i <= m_params.M; i++) {
            c_i.push_back(default_value);
        }
    }

    func(std::vector<double> vec, const common_params& parameters) : m_params(parameters) {
        if (static_cast<int>(vec.size()) != m_params.M + 1) {
            std::cout << "Это какой-то неправильный вектор. Он имеет размер " << vec.size()
                      << ", а должен " << m_params.M + 1 << " \n";
        }
        for (int i = 0; i <= m_params.M; i++) {
            c_i.push_back(vec[i]);
        }
    }

    ~func() = default;

    // Конструктор перемещения
    func(func&& other) noexcept : c_i(std::move(other.c_i)), m_params(other.m_params) {}

    // Оператор перемещающего присваивания
    func& operator=(func&& other) noexcept {
        if (this != &other) {
            // Ссылка на параметры должна быть на тот же объект
            assert(&m_params == &other.m_params);
            c_i = std::move(other.c_i);
        }
        return *this;
    }

    double approx(double x) {
        int index = static_cast<int>(x / m_params.h);
        if (index == m_params.M) return c_i[m_params.M];
        double buffer = x - index * m_params.h;
        return c_i[index] * (1 - buffer) + c_i[index + 1] * buffer;
    }
};

struct equation {
    int num_of_eq = -1;
    double below_diag;
    double onthe_diag;
    double above_diag;
    double right_part;

    const common_params& m_params;

    equation(int number, const common_params& parameters)
        : num_of_eq(number), m_params(parameters) {}

    ~equation() = default;

    void normalize_eq(bool is_forward) {
        if (is_forward) {
            if (!(std::abs(below_diag) < eps * std::abs(onthe_diag))) {
                std::cout << "Не удалось нормировать уравнение " << num_of_eq
                          << ". Коэффициент под диагональю равен " << below_diag << "\n";
                throw -10;
                return;
            } else {
                above_diag = above_diag / onthe_diag;
                right_part = right_part / onthe_diag;
                onthe_diag = 1.;
            }
        } else {
            if (!(std::abs(above_diag) < eps * std::abs(onthe_diag)) ||
                !(std::abs(below_diag) < eps * std::abs(onthe_diag))) {
                std::cout << "Не удалось нормировать уравнение " << num_of_eq
                          << ". Коэффициент под диагональю равен " << above_diag << "\n";
                throw -11;
                return;
            } else {
                // above_diag = above_diag / onthe_diag;
                right_part = right_part / onthe_diag;
                onthe_diag = 1.;
            }
        }
    }

    void init_equation(const func& u /* значения прошлой функции в точках интерполяции */,
                       int num /* номер уравнения от 0 до M*/, int time_step, bool is_default) {
        auto function_of_coefficient = [&u, this](int i, int j) {
            if (i == 0)
                return u.c_i[i] * skalar_of_differential(i, j, i, m_params.M) +
                       u.c_i[i + 1] * skalar_of_differential(i, j, i + 1, m_params.M);

            else if (i == m_params.M)
                return u.c_i[i - 1] * skalar_of_differential(i, j, i - 1, m_params.M) +
                       u.c_i[i] * skalar_of_differential(i, j, i, m_params.M);

            else
                return u.c_i[i - 1] * skalar_of_differential(i, j, i - 1, m_params.M) +
                       u.c_i[i] * skalar_of_differential(i, j, i, m_params.M) +
                       u.c_i[i + 1] * skalar_of_differential(i, j, i + 1, m_params.M);
        };

        if (num == 0) {
            num_of_eq = num;
            below_diag = 0;
            onthe_diag = (1. / m_params.tau) *
                             skalar_without_differential(num, num, m_params.h, m_params.M) +
                         function_of_coefficient(num, num);
            above_diag = (1. / m_params.tau) *
                             skalar_without_differential(num, num + 1, m_params.h, m_params.M) +
                         function_of_coefficient(num, num + 1);
            right_part =
                (1. / m_params.tau) *
                (skalar_without_differential(num, num, m_params.h, m_params.M) * u.c_i[num] +
                 skalar_without_differential(num, num + 1, m_params.h, m_params.M) *
                     u.c_i[num + 1]);
            ;
        } else if (num == m_params.M) {
            num_of_eq = num;
            below_diag = (1. / m_params.tau) *
                             skalar_without_differential(num, num - 1, m_params.h, m_params.M) +
                         function_of_coefficient(num, num - 1);
            onthe_diag = (1. / m_params.tau) *
                             skalar_without_differential(num, num, m_params.h, m_params.M) +
                         function_of_coefficient(num, num);
            above_diag = 0;
            right_part =
                (1. / m_params.tau) *
                (skalar_without_differential(num, num - 1, m_params.h, m_params.M) *
                     u.c_i[num - 1] +
                 skalar_without_differential(num, num, m_params.h, m_params.M) * u.c_i[num]);
        } else {
            num_of_eq = num;
            below_diag = (1. / m_params.tau) *
                             skalar_without_differential(num, num - 1, m_params.h, m_params.M) +
                         function_of_coefficient(num, num - 1);
            onthe_diag = (1. / m_params.tau) *
                             skalar_without_differential(num, num, m_params.h, m_params.M) +
                         function_of_coefficient(num, num);
            above_diag = (1. / m_params.tau) *
                             skalar_without_differential(num, num + 1, m_params.h, m_params.M) +
                         function_of_coefficient(num, num + 1);
            right_part =
                (1. / m_params.tau) *
                (skalar_without_differential(num, num - 1, m_params.h, m_params.M) *
                     u.c_i[num - 1] +
                 skalar_without_differential(num, num, m_params.h, m_params.M) * u.c_i[num] +
                 skalar_without_differential(num, num + 1, m_params.h, m_params.M) *
                     u.c_i[num + 1]);
        }
        if (is_default) {
            right_part +=
                simpson_integral(accuracy_right_part, num, time_step * m_params.tau, m_params);
            // right_part += accuracy_right_part(time_step * m_params.tau, num * m_params.h);
        }
    }
};

std::vector<double> init_vector_of_function_u_value(int time_step,
                                                    const common_params& parameters) {
    std::vector<double> c_i;
    for (int i = 0; i <= parameters.M; i++) {
        c_i.push_back(accuracy_u(time_step * parameters.tau, i * parameters.h));
    }
    return c_i;
}

std::vector<double> init_vector_of_function_right_part_value(int time_step,
                                                             const common_params& parameters) {
    std::vector<double> c_i;
    for (int i = 0; i <= parameters.M; i++) {
        c_i.push_back(
            simpson_integral(accuracy_right_part, i, time_step * parameters.tau, parameters));
        // c_i.push_back(accuracy_right_part(time_step * parameters.tau, i * parameters.h));
    }
    return c_i;
}

// Проблема в том, что эта функция не работает, если надо вычитать в обратную сторону. Надо ввести
// ей параметр стороны. Сейчас это работает вниз (как бы прямой ход Гаусса) уравнение первого
// аргумента.
// (eq2) -= уравнение второго аргумента (eq1)
void eq_below_is_eq_below_minus_eq_above(equation& eq2 /* это уравнение внизу */,
                                         equation& eq1 /* это уравнение сверху */,
                                         bool is_forward) {
    if (eq1.num_of_eq == -1 || eq2.num_of_eq == -1) {
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
        if (is_forward) {
            // надо вычитать не само уравнение, а домноженное на какое-то число. По смыслу, так
            // чтобы обнулить коэффициент поддиагональю.
            double K = eq2.below_diag / eq1.onthe_diag;
            if (std::abs(eq1.below_diag) < eps * std::abs(eq1.onthe_diag)) {
                eq2.below_diag -= K * eq1.onthe_diag;
                eq2.onthe_diag -= K * eq1.above_diag;
                eq2.right_part -= K * eq1.right_part;
                return;
            } else {
                std::cout << "Коэффициент под диагональью уравнения " << eq1.num_of_eq
                          << " ненулевой. Он равен " << eq1.below_diag << "\n";
                throw -2;
                return;
            }
        } else {
            double K = eq2.above_diag / eq1.onthe_diag;
            if (std::abs(eq1.above_diag) < eps * std::abs(eq1.onthe_diag)) {
                eq2.above_diag -= K * eq1.onthe_diag;
                // eq2.onthe_diag -= K * eq1.below_diag;
                eq2.right_part -= K * eq1.right_part;
                return;
            } else {
                std::cout << "Коэффициент под диагональью уравнения " << eq1.num_of_eq
                          << " ненулевой. Он равен " << eq1.below_diag << "\n";
                throw -3;
                return;
            }
        }
    }
}

void print_vector(std::vector<double> func) {
    for (auto c_i : func) std::cout << c_i << "\n";
}

void print_diff_vector(std::vector<double> func, std::vector<double> accuracy_func) {
    for (int i = 0; i < func.size(); i++) {
        std::cout << func[i] - accuracy_func[i] << "  ";
    }
    std::cout << "\n";
}

func solver_of_matrix(std::vector<equation>& matrix, const common_params& parameters) {
    const int& M = parameters.M;
    const double &h = parameters.h, &tau = parameters.tau;

    if (static_cast<int>(matrix.size()) != M + 1)
        std::cout << "Короче, неправильно заданы уравнения (их должно быть M + 1)\n";
    /* номер последнего уравнения равен M */

    for (int i = 0; i < M; i++) {
        matrix[i].normalize_eq(true);
        eq_below_is_eq_below_minus_eq_above(matrix[i + 1], matrix[i], true);
    }

    // matrix[M].normalize_eq();

    for (int i = M; i > 0; i--) {
        matrix[i].normalize_eq(false);
        eq_below_is_eq_below_minus_eq_above(matrix[i - 1], matrix[i], false);
    }

    matrix[0].normalize_eq(false);

    std::vector<double> next_approx;
    for (int i = 0; i <= M; i++) {
        next_approx.push_back(matrix[i].right_part);
    }

    return func(next_approx, parameters);
}

int main() {
    try {
        std::cout << "Начало курсовой программы" << std::endl;

        // задание матрицы. Нужно задать векторы для c_j, c_j_prev (это с прошлого шага, его
        // сохраняем), но их мы ищем просто.

        double X = 1., T = 1.;   // длина сегментов по пространству и по времени соответственно
        int N = 1000, M = 1000;  // M - это шаги по пространству, N - это шаги по времени (точек на
        // одну больше, так как считаем с 0)

        std::cout << "Введи аргументы в таком порядке:  N шагов по времени,  T отрезок времени,  M "
                     "шагов по пространству,  X отрезок пространства.\n";
        std::cin >> N >> T >> M >> X;

        double tau = T / N, h = X / M;

        common_params params(N, tau, M, h);

        printf("h=%lf, tau=%lf; M=%d, N=%d\n", params.h, params.tau, params.M, params.N);

        bool _default = true;  // если истинно, то значит считаем для конкретного точного u
        func u(init_vector_of_function_u_value(0, params), params);

        // задаем матрицу, она должна быть трехдиагональной. Надо создать структуру из векторов
        // коэффициентов уравнений Имеется n+1 уравнений.
        for (int time_step = 0; time_step <= N; ++time_step) {
            std::vector<equation> matrix;

            // инициализируем матрицу
            for (int i = 0; i <= params.M; i++) {
                equation eq(i, params);
                eq.init_equation(u, i, time_step, _default);
                matrix.push_back(eq);
            }

            // тут решаем матрицу и заменяем приближение на правую часть матрицы
            u = solver_of_matrix(matrix, params);
            // тут запишем решение задачи в какой-то файл. Наверное, в виде трехмерной поверхности
            // (время / пространство / значение функции)
            if (time_step == 0) {
                std::cout << "Первое уравнение \n";
                print_diff_vector(u.c_i,
                                  func(init_vector_of_function_u_value(N, params), params).c_i);
                std::cout << "\n\n";
            }
        }
        print_diff_vector(u.c_i, func(init_vector_of_function_u_value(N, params), params).c_i);

    } catch (int e) {
        switch (e) {
            case -10:
                std::cerr << "Проблема в том, что не удалось нормировать уравнение с номером выше "
                             "(при прямом ходу) \n";
            case -11:
                std::cerr << "Проблема в том, что не удалось нормировать уравнение с номером выше "
                             "(при обратном ходу) \n";
            case -2:
                std::cerr << "Проблема в том, что коэффициент под диагональю не равен нулю, хотя "
                             "должен. Это возникло в прямом ходе \n";
            case -3:
                std::cerr << "Проблема в том, что коэффициент под диагональю не равен нулю, хотя "
                             "должен. Это возникло в обратном ходе \n";
        };
    } catch (const std::exception& e) {
        std::cerr << "Стандартное исключение: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Неизвестный тип исключения!" << std::endl;
    }
    return 0;
}