#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// Класс для решения задачи методом конечных элементов
class FEMSolver {
   private:
    int N;                         // количество элементов
    double h;                      // шаг сетки
    std::vector<double> x_nodes;   // узлы сетки
    std::vector<double> solution;  // решение в узлах

    // Параметры задачи (можно менять под свою постановку)
    double k_func(double x) const {
        return 1.0 + x;  // k(x) = 1 + x
    }

    double p_func(double x) const {
        return 1.0;  // p(x) = 1
    }

    double f_func(double x) const {
        return 1.0;  // f(x) = 1
    }

   public:
    // Конструктор: инициализация сетки
    FEMSolver(int num_elements) : N(num_elements) {
        h = 1.0 / N;
        x_nodes.resize(N + 1);
        for (int i = 0; i <= N; ++i) {
            x_nodes[i] = i * h;
        }
        solution.resize(N + 1, 0.0);
    }

    // Решение задачи
    void solve() {
        int n_nodes = N + 1;  // количество узлов

        // Глобальные матрица и вектор
        std::vector<std::vector<double>> A(n_nodes, std::vector<double>(n_nodes, 0.0));
        std::vector<double> F(n_nodes, 0.0);

        // Точки и веса для квадратуры Гаусса (2 точки)
        const double gauss_points[2] = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
        const double gauss_weights[2] = {1.0, 1.0};

        // Цикл по элементам
        for (int e = 0; e < N; ++e) {
            int n1 = e;      // левый узел (индекс)
            int n2 = e + 1;  // правый узел (индекс)
            double x1 = x_nodes[n1];
            double x2 = x_nodes[n2];

            // Локальные матрица и вектор
            double K_local[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
            double M_local[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
            double F_local[2] = {0.0, 0.0};

            // Численное интегрирование по элементу
            for (int gp = 0; gp < 2; ++gp) {
                double xi = gauss_points[gp];

                // Якобиан преобразования
                double detJ = h / 2.0;

                // Функции формы в точке xi
                double N1 = (1.0 - xi) / 2.0;
                double N2 = (1.0 + xi) / 2.0;

                // Производные функций формы по x
                double dN1dx = -1.0 / h;
                double dN2dx = 1.0 / h;

                // Физическая координата точки интегрирования
                double x_gp = N1 * x1 + N2 * x2;

                // Значения коэффициентов в точке интегрирования
                double k_val = k_func(x_gp);
                double p_val = p_func(x_gp);
                double f_val = f_func(x_gp);

                double weight = gauss_weights[gp];

                // Вклад в локальную матрицу жесткости
                K_local[0][0] += k_val * dN1dx * dN1dx * detJ * weight;
                K_local[0][1] += k_val * dN1dx * dN2dx * detJ * weight;
                K_local[1][0] += k_val * dN2dx * dN1dx * detJ * weight;
                K_local[1][1] += k_val * dN2dx * dN2dx * detJ * weight;

                // Вклад в локальную матрицу масс
                M_local[0][0] += p_val * N1 * N1 * detJ * weight;
                M_local[0][1] += p_val * N1 * N2 * detJ * weight;
                M_local[1][0] += p_val * N2 * N1 * detJ * weight;
                M_local[1][1] += p_val * N2 * N2 * detJ * weight;

                // Вклад в локальный вектор нагрузки
                F_local[0] += f_val * N1 * detJ * weight;
                F_local[1] += f_val * N2 * detJ * weight;
            }

            // Ансамблирование (добавление в глобальную систему)
            // Матрица жесткости + матрица масс
            A[n1][n1] += K_local[0][0] + M_local[0][0];
            A[n1][n2] += K_local[0][1] + M_local[0][1];
            A[n2][n1] += K_local[1][0] + M_local[1][0];
            A[n2][n2] += K_local[1][1] + M_local[1][1];

            // Вектор нагрузки
            F[n1] += F_local[0];
            F[n2] += F_local[1];
        }

        // Учет граничного условия u(0) = 0
        // Метод: обнуляем строку и столбец, ставим 1 на диагональ
        for (int i = 0; i < n_nodes; ++i) {
            A[0][i] = 0.0;
            A[i][0] = 0.0;
        }
        A[0][0] = 1.0;
        F[0] = 0.0;

        // Решение системы линейных уравнений методом Гаусса
        solveLinearSystem(A, F, solution);
    }

    // Метод Гаусса для решения СЛАУ
    void solveLinearSystem(std::vector<std::vector<double>>& A, std::vector<double>& b,
                           std::vector<double>& x) {
        int n = A.size();

        // Прямой ход
        for (int i = 0; i < n; ++i) {
            // Поиск главного элемента
            int max_row = i;
            for (int k = i + 1; k < n; ++k) {
                if (std::abs(A[k][i]) > std::abs(A[max_row][i])) {
                    max_row = k;
                }
            }

            // Перестановка строк
            std::swap(A[i], A[max_row]);
            std::swap(b[i], b[max_row]);

            // Проверка на вырожденность
            if (std::abs(A[i][i]) < 1e-10) {
                std::cout << "Матрица вырождена!" << std::endl;
                return;
            }

            // Нормализация текущей строки
            double pivot = A[i][i];
            for (int j = i; j < n; ++j) {
                A[i][j] /= pivot;
            }
            b[i] /= pivot;

            // Исключение в нижележащих строках
            for (int k = i + 1; k < n; ++k) {
                double factor = A[k][i];
                for (int j = i; j < n; ++j) {
                    A[k][j] -= factor * A[i][j];
                }
                b[k] -= factor * b[i];
            }
        }

        // Обратный ход
        x.resize(n);
        for (int i = n - 1; i >= 0; --i) {
            x[i] = b[i];
            for (int j = i + 1; j < n; ++j) {
                x[i] -= A[i][j] * x[j];
            }
        }
    }

    // Вывод результатов
    void printSolution() const {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "\nРешение в узлах:\n";
        std::cout << "----------------------------------------\n";
        std::cout << "|   x    |    u(x)    |\n";
        std::cout << "----------------------------------------\n";
        for (int i = 0; i <= N; ++i) {
            std::cout << "| " << std::setw(5) << x_nodes[i] << "  | " << std::setw(10)
                      << solution[i] << " |\n";
        }
        std::cout << "----------------------------------------\n";
    }

    // Получение решения в конкретной точке (линейная интерполяция)
    double getValueAt(double x) const {
        if (x < 0.0 || x > 1.0) return 0.0;

        int i = static_cast<int>(x / h);
        if (i >= N) i = N - 1;

        double t = (x - x_nodes[i]) / h;
        return (1.0 - t) * solution[i] + t * solution[i + 1];
    }
};

int main() {
    // Параметры задачи
    int N = 10;  // количество элементов (чем больше, тем точнее)

    std::cout << "=== Решение краевой задачи методом Ритца (МКЭ) ===\n";
    std::cout << "Уравнение: -(k(x)u')' + p(x)u = f(x), x in [0,1]\n";
    std::cout << "Условия: u(0) = 0, u'(1) = 0\n";
    std::cout << "Коэффициенты: k(x) = 1 + x, p(x) = 1, f(x) = 1\n";
    std::cout << "Количество элементов: " << N << "\n";

    // Создание и решение задачи
    FEMSolver solver(N);
    solver.solve();

    // Вывод результатов
    solver.printSolution();

    // Дополнительно: значения в некоторых точках
    std::cout << "\nЗначения в характерных точках:\n";
    std::cout << "u(0.0) = " << solver.getValueAt(0.0) << "\n";
    std::cout << "u(0.5) = " << solver.getValueAt(0.5) << "\n";
    std::cout << "u(1.0) = " << solver.getValueAt(1.0) << "\n";

    return 0;
}