#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

struct common_params {
    double tau, h;
    int N, M;
    // M - это число отрезков, а точек на одну больше

    common_params(int _N, double _tau, int _M, double _h) : N(_N), tau(_tau), M(_M), h(_h) {}

    ~common_params() = default;

    std::vector<double> get_vector_OX() {
        std::vector<double> x;
        for (int i = 0; i <= M; i++) x.push_back(i * h);
        return x;
    }
};

// РЕШАЕМ УРАВНЕНИЕ БЮРГЕРСА ВИДА:       du/dt  +  0.5 * d(u^2)/dx  =  0

void plot_two(const std::vector<double>& x, const std::vector<double>& y1,
              const std::vector<double>& y2, const common_params& params, double time = -1) {
    std::ofstream dataFile("data.txt");
    for (size_t i = 0; i < x.size(); ++i) {
        dataFile << x[i] << " " << y1[i] << " " << y2[i] << std::endl;
    }
    dataFile.close();

    if (time < 0) time = params.tau * params.N;

    char filename[256];
    snprintf(filename, sizeof(filename), "graph_another_N=%d_M=%d_T=%lf_X=%lf.png", params.N,
             params.M, time, params.h * params.M);

    char command[512];
    snprintf(command, sizeof(command),
             "gnuplot -e \"set terminal png size 1200,800; set output '%s'; plot 'data.txt' "
             "using 1:2 with lines linewidth 3 title 'Numerical', 'data.txt' using 1:3 with lines "
             "linewidth 3 title 'Exact'\"",
             filename);
    system(command);

    char open_command[512];
    snprintf(open_command, sizeof(open_command), "xdg-open '%s' 2>/dev/null &", filename);
    system(open_command);
}

void plot_error_log(const std::vector<double>& x, const std::vector<double>& y1,
                    const std::vector<double>& y2, const common_params& params, double time = -1) {
    std::ofstream dataFile("data.txt");
    for (size_t i = 0; i < x.size(); ++i) {
        double error = std::abs(y1[i] - y2[i]);
        // Защита от log(0)
        if (error < 0.0000001) error = 0.0000001;
        dataFile << x[i] << " " << error << std::endl;
    }
    dataFile.close();

    if (time < 0) time = params.tau * params.N;

    char filename[256];
    snprintf(filename, sizeof(filename), "error_log_another_N=%d_M=%d_T=%lf_X=%lf.png", params.N,
             params.M, time, params.h * params.M);

    char command[512];
    snprintf(command, sizeof(command),
             "gnuplot -e \"set terminal png size 1200,800; "
             "set output '%s'; "
             "set xlabel 'x'; "
             "set ylabel '|Numerical - Exact|'; "
             "set title 'Absolute Error (log scale)'; "
             "set grid; "
             "set logscale y; "
             "plot 'data.txt' using 1:2 with lines linewidth 3 title 'Error'\"",
             filename);

    system(command);

    char open_command[512];
    snprintf(open_command, sizeof(open_command), "xdg-open '%s' 2>/dev/null &", filename);
    system(open_command);
}

void plot_error(const std::vector<double>& x, const std::vector<double>& y1,
                const std::vector<double>& y2, const common_params& params, double time = -1) {
    std::ofstream dataFile("data.txt");
    for (size_t i = 0; i < x.size(); ++i) {
        double error = std::abs(y1[i] - y2[i]);
        dataFile << x[i] << " " << error << std::endl;
    }
    dataFile.close();

    if (time < 0) time = params.tau * params.N;

    char filename[128];
    snprintf(filename, sizeof(filename), "error_another_N=%d_M=%d_T=%lf_X=%lf.png", params.N,
             params.M, time, params.h * params.M);

    char command[512];
    snprintf(command, sizeof(command),
             "gnuplot -e \"set terminal png size 1200,800; "
             "set output '%s'; "
             "set xlabel 'x'; "
             "set ylabel '|Numerical - Exact|'; "
             "set title 'Absolute Error'; "
             "set grid; "
             "plot 'data.txt' using 1:2 with lines linewidth 3 title 'Error'\"",
             filename);

    system(command);

    char open_command[512];
    snprintf(open_command, sizeof(open_command), "xdg-open '%s' 2>/dev/null &", filename);
    system(open_command);
}

double eps = 1e-12;

int function_type = 0;

double accuracy_right_part(double t, double x) {
    switch (function_type) {
        case 0:
            return std::exp(t) * (1.5 + std::cos(3 * M_PI * x)) *
                   (1 - 3 * M_PI * std::exp(t) * std::sin(3 * M_PI * x));
        case 1:
            return (x > 0.2 && x < 0.4) ? 1 : 0;
        case 2:
            return 1 + x + t;
        case 3:
            return x * x * (1 + 2 * t * t * x);
    }
    return 0;
}

double accuracy_u(double t, double x) {
    switch (function_type) {
        case 0:
            return std::exp(t) * (1.5 + std::cos(3 * M_PI * x));
        case 1:
            return 0;
        case 2:
            return x + t;
        case 3:
            return t * x * x;
    }
    return 0;
}

// как мне как-то проще определить как считать произвоидные для j, k. Можно сделать вектор из
// phi

// i - это номер уравнения по сути (phi в произведении, которое не под дифференциалом), j и k -
// это phi которые под дифференциалом.

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



// --- Ограничитель наклонов (TVBM minmod) ---
double minmod(double a, double b, double c) {
    if (a > 0 && b > 0 && c > 0) return std::min(a, std::min(b, c));
    if (a < 0 && b < 0 && c < 0) return std::max(a, std::max(b, c)); // все отрицательные -> ближайшее к нулю
    return 0.0;
}

void apply_slope_limiter(func& u, const common_params& p) {
    const double h = p.h;
    const double M_tvb = 100.0;   // оценка второй производной (для тестовой задачи)
    std::vector<double> limited = u.c_i; // копия

    for (int i = 1; i < p.M; ++i) {
        double uL = u.c_i[i-1], uC = u.c_i[i], uR = u.c_i[i+1];

        // центральная разность и односторонние
        double d_center = (uR - uL) / 2.0;
        double d_left   = uC - uL;
        double d_right  = uR - uC;

        // TVBM minmod
        double d_limited;
        if (std::abs(d_center) <= M_tvb * h * h) {
            d_limited = d_center; // сохраняем точность в экстремумах
        } else {
            d_limited = minmod(d_center, d_left, d_right);
        }

        // Восстанавливаем значение в узле i так, чтобы наклон на левом и правом элементах
        // соответствовал ограниченной центральной разности.
        // Мы просто заменяем u_i на (u_{i+1} + u_{i-1})/2 - (d_limited * h)/2? 
        // Несколько сложнее - сохраним среднее значение на элементе, но для непрерывного Галёркина
        // проще всего скорректировать u_i, чтобы подавить осцилляции.
        // Используем подход: u_i^{new} = u_i, если ограниченный наклон совпадает с центральным,
        // иначе сглаживаем.
        if (std::abs(d_limited - d_center) > 1e-12) {
            // ограничитель сработал — заменяем u_i на (u_{i-1} + u_{i+1})/2 
            // (обнуляем вторую производную)
            limited[i] = (uL + uR) / 2.0;
        }
    }

    u.c_i.swap(limited);
}



// Ошибка в норме L∞ (maximum norm)
double compute_L_inf_error(const func& numerical, double time) {
    double max_error = 0.0;
    const common_params& p = numerical.m_params;

    for (int i = 0; i <= p.M; i++) {
        double exact = accuracy_u(time, i * p.h);
        double error = std::abs(numerical.c_i[i] - exact);
        if (error > max_error) max_error = error;
    }
    return max_error;
}

// Ошибка в норме L1 (интеграл модуля ошибки)
// Используем составную формулу трапеций
double compute_L1_error(const func& numerical, double time) {
    const common_params& p = numerical.m_params;
    double sum = 0.0;

    for (int i = 0; i <= p.M; i++) {
        double exact = accuracy_u(time, i * p.h);
        double error = std::abs(numerical.c_i[i] - exact);

        if (i == 0 || i == p.M)
            sum += error / 2.0;
        else
            sum += error;
    }
    return sum * p.h;
}

// Ошибка в норме L2 (среднеквадратичная)
// Интеграл: sqrt(\int(u_num - u_exact)^2 dx)
// Используем точную квадратуру для кусочно-линейных функций
double compute_L2_error(const func& numerical, double time) {
    const common_params& p = numerical.m_params;
    double integral_sq = 0.0;

    // Интегрируем по каждому элементу точно
    for (int i = 0; i < p.M; i++) {
        // Значения ошибки в узлах i и i+1
        double uL_num = numerical.c_i[i];
        double uR_num = numerical.c_i[i + 1];
        double uL_ex = accuracy_u(time, i * p.h);
        double uR_ex = accuracy_u(time, (i + 1) * p.h);

        double eL = uL_num - uL_ex;
        double eR = uR_num - uR_ex;

        // Точный интеграл от (eL*(1-ξ) + eR*ξ)² по элементу длины h
        // ∫₀¹ (eL(1-ξ) + eRξ)² h dξ = h*(eL² + eL*eR + eR²)/3
        integral_sq += (p.h / 3.0) * (eL * eL + eL * eR + eR * eR);
    }

    return std::sqrt(integral_sq);
}

// Удобная функция для вывода всех норм
void print_all_errors(const func& numerical, double time) {
    double L_inf = compute_L_inf_error(numerical, time);
    double L1 = compute_L1_error(numerical, time);
    double L2 = compute_L2_error(numerical, time);

    std::cout << " (t = " << time << "): "
              << "L_inf = " << L_inf << ", "
              << "L1 = " << L1 << ", "
              << "L2 = " << L2 << "\n";
}

struct equation {
    int num_of_eq = -1;
    double below_diag = 0;
    double onthe_diag = 0;
    double above_diag = 0;
    double right_part = 0;

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
                          << ". Коэффициент над диагональю равен " << above_diag << "\n";
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
                       int num /* номер уравнения от 0 до M*/, int time_step,
                       bool is_default /* правая часть не равна нулю */) {
        double h = m_params.h;
        double tau = m_params.tau;
        int M = m_params.M;

        num_of_eq = num;
        below_diag = 0;
        onthe_diag = 0;
        above_diag = 0;
        right_part = 0;

        // мб лучше задать точное решение на границе, вместо коэффициентов на диагонали.
        // Это будет лучше
        if (num == 0 || num == M) {
            right_part = accuracy_u((time_step + 1) * tau, num * h);
            onthe_diag = 1;
            return;
        }

        /* if (num == 0) onthe_diag -= 0.5 * u.c_i[num] * (tau);
        if (num == M) onthe_diag += 0.5 * u.c_i[num] * (tau); */

        if (num != 0) below_diag = h / 6;
        if (num != M) above_diag = h / 6;

        if (num != 0) {
            below_diag += 0 + ((-1. / 6) * u.c_i[num - 1] + (-1. / 12) * u.c_i[num]) * (tau);
            onthe_diag += h / 3 + (-1. / 12) * (u.c_i[num - 1] + 2 * u.c_i[num]) * (tau);
            above_diag += 0;
        }
        if (num != M) {
            below_diag += 0;
            onthe_diag += h / 3 + (1. / 12) * (2 * u.c_i[num] + u.c_i[num + 1]) * (tau);
            above_diag += 0 + ((1. / 12) * u.c_i[num] + (1. / 6) * u.c_i[num + 1]) * (tau);
        }

        if (is_default)
            right_part += tau * simpson_integral(accuracy_right_part, num,
                                                 (time_step + 1) * m_params.tau, m_params);

        if (num != 0) right_part += (h / 6) * (2 * u.c_i[num] + u.c_i[num - 1]);
        if (num != M) right_part += (h / 6) * (2 * u.c_i[num] + u.c_i[num + 1]);
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
                std::cout << "Коэффициент под диагональю уравнения " << eq1.num_of_eq
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
                std::cout << "Коэффициент над диагональю уравнения " << eq1.num_of_eq
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

    matrix[M].right_part /= matrix[M].onthe_diag;
    matrix[M].onthe_diag = 1.;

    for (int i = M; i > 0; i--) {
        matrix[i].right_part /= matrix[i].onthe_diag;
        matrix[i].onthe_diag = 1.;
        eq_below_is_eq_below_minus_eq_above(matrix[i - 1], matrix[i], false);
    }

    matrix[0].right_part /= matrix[0].onthe_diag;

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
                     "шагов по пространству,  X отрезок пространства,  TYPE - тип функции.\n";
        std::cin >> N >> T >> M >> X >> function_type;

        double tau = T / N, h = X / M;

        common_params params(N, tau, M, h);

        printf("h=%lf, tau=%lf; M=%d, N=%d\n", params.h, params.tau, params.M, params.N);

        bool _default = true;  // если истинно, то значит считаем для конкретного точного u и правая
                               // часть не равна 0
        func u(init_vector_of_function_u_value(0, params), params);

        // задаем матрицу, она должна быть трехдиагональной. Надо создать структуру из векторов
        // коэффициентов уравнений Имеется n+1 уравнений.
        for (int time_step = 0; time_step < N; ++time_step) {
            std::vector<equation> matrix;

            // инициализируем матрицу
            for (int i = 0; i <= params.M; i++) {
                equation eq(i, params);
                eq.init_equation(u, i, time_step, _default);
                matrix.push_back(eq);
            }

            // тут решаем матрицу и заменяем приближение на правую часть матрицы
            u = solver_of_matrix(matrix, params);
            //apply_slope_limiter(u, params);

            // тут запишем решение задачи в какой-то файл. Наверное, в виде трехмерной поверхности
            // (время / пространство / значение функции)
            if (time_step == 0 && false) {
                std::cout << "\n Первое уравнение \n";
                print_diff_vector(u.c_i,
                                  func(init_vector_of_function_u_value(T, params), params).c_i);
                std::cout << "\n\n";
            }
        }

        std::cout << "\nОшибки\n";
        print_all_errors(u, N * params.tau);
        print_all_errors(u, T);

        // Сохраняем данные для построения графиков
        std::ofstream dataFile("numerical_solution.txt");
        dataFile << "# x_center   numerical   exact" << std::endl;
        for (int i = 0; i < M; i++) {
            double x_center = (i + 0.5) * h;
            dataFile << x_center << " " << u.c_i[i] << " "
                     << func(init_vector_of_function_u_value(T, params), params).c_i[i]
                     << std::endl;
        }
        dataFile.close();

        plot_two(params.get_vector_OX(), u.c_i,
                 func(init_vector_of_function_u_value(N, params), params).c_i, params);

        plot_error_log(params.get_vector_OX(), u.c_i,
                       func(init_vector_of_function_u_value(N, params), params).c_i, params);

        plot_error(params.get_vector_OX(), u.c_i,
                   func(init_vector_of_function_u_value(N, params), params).c_i, params);

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