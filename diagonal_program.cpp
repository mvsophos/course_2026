#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

// ПРОБЛЕМА В ТОМ, ЧТО РЕШЕНИЕ НЕ МЕНЯЕТСЯ, ПОТОМУ-ЧТО ПРОИЗВОДНАЯ ОТ КВАДРАТА РАВНА 0 ПО СУТИ. ОТ
// ЭТОГО ОСТАЕТСЯ ТОЛЬКО ПРАВАЯ ЧАСТЬ И ПРОИЗВОДНАЯ ПО t
// это уже неклассический метод конечных элементов

void plot_two(const std::vector<double>& x, const std::vector<double>& y1,
              const std::vector<double>& y2) {
    std::ofstream dataFile("data.txt");
    for (size_t i = 0; i < x.size(); ++i) {
        dataFile << x[i] << " " << y1[i] << " " << y2[i] << std::endl;
    }
    dataFile.close();

    system(
        "gnuplot -e \"set terminal png size 1200,800; set output 'graph.png'; plot 'data.txt' "
        "using 1:2 with lines linewidth 3 title 'Numerical', 'data.txt' using 1:3 with lines "
        "linewidth 3 title 'Exact'\"");
    system("xdg-open graph.png 2>/dev/null &");
}

void plot_error_log(const std::vector<double>& x, const std::vector<double>& y1,
                    const std::vector<double>& y2) {
    std::ofstream dataFile("data.txt");
    for (size_t i = 0; i < x.size(); ++i) {
        double error = std::abs(y1[i] - y2[i]);
        // Защита от log(0)
        if (error < 1e-300) error = 1e-300;
        dataFile << x[i] << " " << error << std::endl;
    }
    dataFile.close();

    system(
        "gnuplot -e \"set terminal png size 1200,800; "
        "set output 'graph_error_log_graph.png'; "
        "set xlabel 'x'; "
        "set ylabel '|Numerical - Exact|'; "
        "set title 'Absolute Error (log scale)'; "
        "set grid; "
        "set logscale y; "
        "plot 'data.txt' using 1:2 with lines linewidth 3 title 'Error'\"");
    system("xdg-open graph_error_log_graph.png 2>/dev/null &");
}

void plot_error(const std::vector<double>& x, const std::vector<double>& y1,
                const std::vector<double>& y2) {
    std::ofstream dataFile("data.txt");
    for (size_t i = 0; i < x.size(); ++i) {
        double error = std::abs(y1[i] - y2[i]);
        // Защита от log(0)
        if (error < 1e-300) error = 1e-300;
        dataFile << x[i] << " " << error << std::endl;
    }
    dataFile.close();

    system(
        "gnuplot -e \"set terminal png size 1200,800; "
        "set output 'graph_error_graph.png'; "
        "set xlabel 'x'; "
        "set ylabel '|Numerical - Exact|'; "
        "set title 'Absolute Error (log scale)'; "
        "set grid; "
        "plot 'data.txt' using 1:2 with lines linewidth 3 title 'Error'\"");
    system("xdg-open graph_error_graph.png 2>/dev/null &");
}

struct common_params {
    double tau, h;
    int N, M;
    // M - это число отрезков, а точек на одну больше

    common_params(int _N, double _tau, int _M, double _h) : N(_N), tau(_tau), M(_M), h(_h) {}

    ~common_params() = default;

    std::vector<double> get_vector_OX() {
        std::vector<double> x;
        for (int i = 0; i < M; i++) x.push_back((i + 0.5) * h);
        return x;
    }
};

double eps = 1e-12;

int function_type = 0;

// неправильно решение задал. Это не решение уравнения
double accuracy_right_part(double t, double x) {
    switch (function_type) {
        case 0:
            return std::exp(t) * (1.5 + std::cos(3 * M_PI * x)) *
                   (1 - 3 * M_PI * std::exp(t) * std::sin(3 * M_PI * x));
        case 1:
            return 0;
    }
    return 0;
}

double accuracy_u(double t, double x) {
    switch (function_type) {
        case 0:
            return std::exp(t) * (1.5 + std::cos(3 * M_PI * x));
        case 1:
            return (x > 0.2 && x < 0.4) ? 1 : 0;
    }
    return 0;
}

class func {
   public:
    std::vector<double> c_i;
    const common_params& m_params;

    func(const common_params& parameters) : m_params(parameters) {
        double default_value = 1;
        for (int i = 0; i < m_params.M; i++) {
            c_i.push_back(default_value);
        }
    }

    func(std::vector<double> vec, const common_params& parameters) : m_params(parameters) {
        if (static_cast<int>(vec.size()) != m_params.M) {
            std::cout << "Это какой-то неправильный вектор. Он имеет размер " << vec.size()
                      << ", а должен " << m_params.M << " \n";
        }
        for (int i = 0; i < m_params.M; i++) {
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

//

std::vector<double> init_vector_of_function_u_value(int time_step,
                                                    const common_params& parameters) {
    std::vector<double> c_i;
    for (int i = 0; i < parameters.M; i++) {
        c_i.push_back(accuracy_u(time_step * parameters.tau, (i + 0.5) * parameters.h));
    }
    return c_i;
}

void print_vector(std::vector<double> func) {
    for (auto c_i : func) std::cout << c_i << "\n";
}

void print_diff_vector(std::vector<double> func, std::vector<double> accuracy_func) {
    if (func.size() > 200) return;
    for (int i = 0; i < func.size(); i++) {
        std::cout << func[i] - accuracy_func[i] << "  ";
    }
    std::cout << "\n";
}

//

double flow_from(double a, double b) {
    //return 0;
    return 0.25 * (a * a + b * b) - 0.5 * std::max(std::abs(a), std::abs(b)) * (b - a);
}

int main() {
    try {
        std::cout << "Начало курсовой программы" << std::endl;

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

        func u_current(params);

        // задаем матрицу, она должна быть трехдиагональной. Надо создать структуру из векторов
        // коэффициентов уравнений Имеется n+1 уравнений.
        double flux = 0;
        for (int time_step = 0; time_step <= N; ++time_step) {
            // инициализируем матрицу
            for (int i = 0; i < params.M; i++) {
                if (i == 0)
                    // flux = flow_from(u.c_i[0], u.c_i[1]) - u.c_i[0] * u.c_i[0] / 2;
                    flux = flow_from(u.c_i[0], u.c_i[1]) -
                           flow_from(accuracy_u(time_step * params.tau, 0.0), u.c_i[0]);
                else if (i == params.M - 1)
                    flux = flow_from(u.c_i[M - 1], accuracy_u(time_step * params.tau, X)) -
                           flow_from(u.c_i[M - 2], u.c_i[M - 1]);
                else
                    flux = flow_from(u.c_i[i], u.c_i[i + 1]) - flow_from(u.c_i[i - 1], u.c_i[i]);

                // сразу считаем уравнение. Это чисто здесь только решение.
                u_current.c_i[i] =
                    params.tau *
                        (accuracy_right_part(time_step * params.tau, (i + 1) * params.h) +
                         4 * accuracy_right_part(time_step * params.tau, (i + 0.5) * params.h) +
                         accuracy_right_part(time_step * params.tau, i * params.h)) /
                        6 +
                    u.c_i[i] - flux * params.tau / params.h;
            }
            // тут мы как бы уже решили
            for (int o = 0; o < params.M; o++) {
                u.c_i[o] = u_current.c_i[o];
                if (std::isnan(u_current.c_i[o]) || std::isinf(u_current.c_i[o])) {
                    std::cerr << "Ошибка в решении на шаге " << time_step << "\n";
                    // вероятнее всего из-за того что нарушается условие CFL. Должно быть меньше 1.
                    throw -100;
                }
            }

            // тут запишем решение задачи в какой-то файл. Наверное, в виде трехмерной поверхности
            // (время / пространство / значение функции)
            if (time_step < 2 && u.c_i.size() <= 250) {
                std::cout << "Приближение " << time_step + 1 << "\n";
                print_diff_vector(u.c_i,
                                  func(init_vector_of_function_u_value(0, params), params).c_i);
                std::cout << "\n\n";
            }
        }
        std::cout << "Конечное приближение (разница между точным и найденным) \n";
        print_diff_vector(u.c_i, func(init_vector_of_function_u_value(N, params), params).c_i);

        //
        // Вычисляем точное решение в центрах ячеек
        std::vector<double> exact_solution(M);
        double t_final = T;  // конечное время

        for (int i = 0; i < M; i++) {
            double x_center = (i + 0.5) * h;
            exact_solution[i] = accuracy_u(t_final, x_center);
        }

        // Считаем ошибки
        double L1 = 0.0, L2 = 0.0, Linf = 0.0;
        for (int i = 0; i < M; i++) {
            double err = std::abs(u.c_i[i] - exact_solution[i]);
            L1 += err * h;
            L2 += err * err * h;
            Linf = std::max(Linf, err);
        }
        L2 = std::sqrt(L2);

        /* std::cout << "L1-норма ошибки:    " << L1 << std::endl;
        std::cout << "L2-норма ошибки:    " << L2 << std::endl;
        std::cout << "L_inf-норма ошибки: " << Linf << std::endl; */

        std::cout << " & " << L1 << " & " << L2 << " & " << Linf << "\n";

        // Сохраняем данные для построения графиков
        std::ofstream dataFile("numerical_solution.txt");
        dataFile << "# x_center   numerical   exact" << std::endl;
        for (int i = 0; i < M; i++) {
            double x_center = (i + 0.5) * h;
            dataFile << x_center << " " << u.c_i[i] << " " << exact_solution[i] << std::endl;
        }
        dataFile.close();

        //

        std::cout << "Максимальный элемент: " << *std::max_element(u.c_i.begin(), u.c_i.end())
                  << "\n";
        std::cout << "Минимальный  элемент: " << *std::min_element(u.c_i.begin(), u.c_i.end())
                  << "\n";

        plot_two(params.get_vector_OX(), u.c_i,
                 func(init_vector_of_function_u_value(N, params), params).c_i);

        plot_error_log(params.get_vector_OX(), u.c_i,
                       func(init_vector_of_function_u_value(N, params), params).c_i);

        plot_error(params.get_vector_OX(), u.c_i,
                   func(init_vector_of_function_u_value(N, params), params).c_i);

        // std::system("xdg-open graph.png");

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
            case -100:
                std::cerr << "Проблема в том, что решение плохое, там NaN или Inf\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Стандартное исключение: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Неизвестный тип исключения!" << std::endl;
    }
    return 0;
}