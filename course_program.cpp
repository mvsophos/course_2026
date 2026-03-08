#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

class Ritz {
   private:
    int N;                         // количество элементов
    double h;                      // шаг сетки
    std::vector<double> x_nodes;   // узлы сетки
    std::vector<double> solution;  // решение

    // ЯВНОЕ ЗАДАНИЕ ОПЕРАТОРА L: Lu = -(k(x)u')' + p(x)u
    // Для метода Ритца нам нужно (Lφⱼ, φᵢ) = ∫(k φⱼ' φᵢ' + p φⱼ φᵢ) dx

    // Коэффициенты оператора
    double k_func(double x) const {
        return 1.0 + x;  // пример: k(x) = 1 + x
    }

    double p_func(double x) const {
        return 1.0;  // пример: p(x) = 1
    }

    double f_func(double x) const {
        return 1.0;  // пример: f(x) = 1
    }

    // ВЫЧИСЛЕНИЕ ЛОКАЛЬНЫХ МАТРИЦ ДЛЯ ОПЕРАТОРА L
    // Здесь явно реализовано (Lφⱼ, φᵢ) на элементе

    void computeLocalMatrices(double x1, double x2,
                              double K_local[2][2],  // матрица жесткости (от k)
                              double M_local[2][2],  // матрица масс (от p)
                              double F_local[2]) {   // вектор нагрузки (от f)

        // Обнуляем локальные матрицы
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                K_local[i][j] = 0.0;
                M_local[i][j] = 0.0;
            }
            F_local[i] = 0.0;
        }

        // Точки Гаусса для численного интегрирования
        const double gauss_points[2] = {-1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0)};
        const double gauss_weights[2] = {1.0, 1.0};

        double detJ = (x2 - x1) / 2.0;  // = h/2 для равномерной сетки

        for (int gp = 0; gp < 2; gp++) {
            double xi = gauss_points[gp];
            double weight = gauss_weights[gp];

            // ФУНКЦИИ ФОРМЫ (базисные функции φ на элементе)
            // В локальных координатах ξ ∈ [-1, 1]
            double N1 = (1.0 - xi) / 2.0;  // φ₁(ξ)
            double N2 = (1.0 + xi) / 2.0;  // φ₂(ξ)

            // ПРОИЗВОДНЫЕ ФУНКЦИЙ ФОРМЫ по x
            // dN/dx = dN/dξ * dξ/dx, где dξ/dx = 2/(x2-x1) = 2/h
            double dN1dx = -1.0 / (x2 - x1);  // = -1/h
            double dN2dx = 1.0 / (x2 - x1);   // = 1/h

            // Физическая координата точки интегрирования
            double x_gp = N1 * x1 + N2 * x2;

            // Значения коэффициентов в точке интегрирования
            double k_val = k_func(x_gp);
            double p_val = p_func(x_gp);
            double f_val = f_func(x_gp);

            // ВЫЧИСЛЕНИЕ (Lφ_i, φ_j) НА ЭЛЕМЕНТЕ
            // Матрица жесткости: ∫ k φ_j' φ_j' dx
            K_local[0][0] += k_val * dN1dx * dN1dx * detJ * weight;
            K_local[0][1] += k_val * dN1dx * dN2dx * detJ * weight;
            K_local[1][0] += k_val * dN2dx * dN1dx * detJ * weight;
            K_local[1][1] += k_val * dN2dx * dN2dx * detJ * weight;

            // Матрица масс: ∫ p φ_j φ_i dx
            M_local[0][0] += p_val * N1 * N1 * detJ * weight;
            M_local[0][1] += p_val * N1 * N2 * detJ * weight;
            M_local[1][0] += p_val * N2 * N1 * detJ * weight;
            M_local[1][1] += p_val * N2 * N2 * detJ * weight;

            // Вектор нагрузки: ∫ f φ_j dx
            F_local[0] += f_val * N1 * detJ * weight;
            F_local[1] += f_val * N2 * detJ * weight;
        }
    }

   public:
    Ritz(int num_elements) : N(num_elements) {
        h = 1.0 / N;
        x_nodes.resize(N + 1);
        for (int i = 0; i <= N; ++i) {
            x_nodes[i] = i * h;
        }
        solution.resize(N + 1, 0.0);
    }

    void solve() {
        int n_nodes = N + 1;

        // Глобальная матрица A = (Lφ_j, φ_i)
        std::vector<std::vector<double>> A(n_nodes, std::vector<double>(n_nodes, 0.0));
        std::vector<double> F(n_nodes, 0.0);  // глобальный вектор (f, φ_i)

        // Цикл по элементам
        for (int e = 0; e < N; ++e) {
            int n1 = e;
            int n2 = e + 1;
            double x1 = x_nodes[n1];
            double x2 = x_nodes[n2];

            double K_local[2][2], M_local[2][2], F_local[2];

            // Вычисляем локальные матрицы для оператора L
            computeLocalMatrices(x1, x2, K_local, M_local, F_local);

            // АНСАМБЛИРОВАНИЕ: добавляем вклады в глобальную систему
            // Полная локальная матрица = K_local + M_local
            A[n1][n1] += K_local[0][0] + M_local[0][0];
            A[n1][n2] += K_local[0][1] + M_local[0][1];
            A[n2][n1] += K_local[1][0] + M_local[1][0];
            A[n2][n2] += K_local[1][1] + M_local[1][1];

            // Добавляем вклад в глобальный вектор нагрузки
            F[n1] += F_local[0];
            F[n2] += F_local[1];
        }

        // УЧЕТ ГРАНИЧНЫХ УСЛОВИЙ
        // Главное условие: u(0) = 0
        for (int i = 0; i < n_nodes; ++i) {
            A[0][i] = 0.0;
            A[i][0] = 0.0;
        }
        A[0][0] = 1.0;
        F[0] = 0.0;

        // Естественное условие u'(1) = 0 учитывается автоматически

        // Решаем систему A*α = F
        solveLinearSystem(A, F, solution);
    }

    void solveLinearSystem(std::vector<std::vector<double>>& A, std::vector<double>& b,
                           std::vector<double>& x) {
        int n = A.size();

        // Прямой ход метода Гаусса
        for (int i = 0; i < n; ++i) {
            // Поиск главного элемента
            int max_row = i;
            for (int k = i + 1; k < n; ++k) {
                if (std::abs(A[k][i]) > std::abs(A[max_row][i])) {
                    max_row = k;
                }
            }

            std::swap(A[i], A[max_row]);
            std::swap(b[i], b[max_row]);

            if (std::abs(A[i][i]) < 1e-10) {
                std::cout << "Матрица вырождена!" << std::endl;
                return;
            }

            // Нормализация
            double pivot = A[i][i];
            for (int j = i; j < n; ++j) {
                A[i][j] /= pivot;
            }
            b[i] /= pivot;

            // Исключение
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

    void printSolution() const {
        std::cout << "Решение методом Ритца:\n";
        std::cout << "u_h(x) = sum(альфа_j * φ_j^h(x))\n\n";
        std::cout << "Коэффициенты разложения альфа_j (значения в узлах):\n";
        for (int i = 0; i <= N; ++i) {
            if (i < 3 || i > N - 3) {
                std::cout << "альфа_" << i << " = u(" << x_nodes[i] << ") = " << std::fixed
                          << std::setprecision(6) << solution[i] << "\n";
            }
            if (i == 3) {
                std::cout << "---------------\n";
            }
        }
    }
};

/*
A[i][j] = (Lφ_j^h, φ_i^h) - глобальная матрица
F[i] = (f, φ_i^h)         - глобальный вектор

Решение:
u^h = \sum α_j φ_j^h, где α_j - решение системы
*/

int main() {
    int N;  // количество элементов

    std::cout << "Введите N:   ";
    std::cin >> N;

    std::cout << "МЕТОД РИТЦА ДЛЯ УРАВНЕНИЯ: -(k(x)u')' + p(x)u = f(x)\n";
    std::cout << "Граничные условия: u(0) = 0, u'(1) = 0\n\n";

    Ritz solver(N);
    solver.solve();
    solver.printSolution();

    return 0;
}