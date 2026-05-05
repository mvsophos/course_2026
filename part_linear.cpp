#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

// Собственная константа pi, не требующая макросов
const double PI = std::acos(-1.0);

// ----------------------------------------------------------------------
// Структура общих параметров (аналогична Вашей)
// ----------------------------------------------------------------------
struct common_params {
    double tau, h;
    int N;     // число шагов по времени
    int M;     // число элементов
    double X;  // длина отрезка
    common_params(int _N, double _tau, int _M, double _X)
        : N(_N), tau(_tau), M(_M), h(_X / _M), X(_X) {}
};

int function_type = 0;

// ----------------------------------------------------------------------
// Точное решение (manufactured solution)
// ----------------------------------------------------------------------
double accuracy_u(double t, double x) {
    switch (function_type) {
        case 0:
            return std::exp(t) * (1.5 + std::cos(3 * PI * x));
    }
    return 0;
}

double accuracy_f(double t, double x) {
    // правая часть: du/dt + 0.5*d(u^2)/dx
    switch (function_type) {
        case 0:
            double u = accuracy_u(t, x);
            double ut = u;  // производная по t совпадает с u
            double ux = -3 * PI * std::exp(t) * std::sin(3 * PI * x);
            double f_conv_x = u * ux;  // производная от 0.5 u^2
            return ut + f_conv_x;
    }
    return 0;
}

// ----------------------------------------------------------------------
// Класс для хранения DG-решения (кусочно-линейное)
// ----------------------------------------------------------------------
class DG_Solution {
   public:
    std::vector<double> u_mean;   // средние значения по элементам
    std::vector<double> u_slope;  // наклоны (коэффициенты при phi_1)
    double h;                     // размер элемента
    double X;                     // длина отрезка
    int M;                        // число элементов

    DG_Solution(int M_, double h_, double X_) : M(M_), h(h_), X(X_) {
        u_mean.resize(M, 0.0);
        u_slope.resize(M, 0.0);
    }

    // Восстановить значение в произвольной точке x (для тестов)
    double approx(double x) const {
        int j = static_cast<int>(x / h);
        if (j >= M) j = M - 1;
        if (j < 0) j = 0;
        double xc = (j + 0.5) * h;       // центр элемента
        double xi = 2.0 * (x - xc) / h;  // [-1,1]
        return u_mean[j] + u_slope[j] * xi;
    }

    // Заполнить начальное условие по точному решению
    void set_initial(double t) {
        for (int j = 0; j < M; ++j) {
            double xc = (j + 0.5) * h;
            u_mean[j] = accuracy_u(t, xc);
            // Наклон согласован с производной точного решения
            double ux = -3 * PI * std::exp(t) * std::sin(3 * PI * xc);
            u_slope[j] = ux * h / 2.0;  // phi_1 = 2(x-xc)/h  =>  du/dx = 2*u_slope/h
        }
    }

    // Метод to_func() удалён – он не используется и требовал бы полного
    // определения класса func из another_program.cpp.
};

// ----------------------------------------------------------------------
// Функции ошибок (адаптированы для DG_Solution)
// ----------------------------------------------------------------------
double compute_L2_error(const DG_Solution& sol, double time) {
    double integral = 0.0;
    for (int j = 0; j < sol.M; ++j) {
        double xL = j * sol.h;
        double xR = (j + 1) * sol.h;
        double eL = sol.approx(xL) - accuracy_u(time, xL);
        double eR = sol.approx(xR) - accuracy_u(time, xR);
        integral += (sol.h / 3.0) * (eL * eL + eL * eR + eR * eR);
    }
    return std::sqrt(integral);
}

double compute_Linf_error(const DG_Solution& sol, double time) {
    double max_err = 0.0;
    for (int j = 0; j < sol.M; ++j) {
        double xc = (j + 0.5) * sol.h;
        max_err = std::max(max_err, std::abs(sol.approx(xc) - accuracy_u(time, xc)));
        max_err = std::max(max_err, std::abs(sol.approx(j * sol.h) - accuracy_u(time, j * sol.h)));
        max_err = std::max(
            max_err, std::abs(sol.approx((j + 1) * sol.h) - accuracy_u(time, (j + 1) * sol.h)));
    }
    return max_err;
}

// ----------------------------------------------------------------------
// Численный поток (локальный Лакса-Фридрихса)
// ----------------------------------------------------------------------
double llf_flux(double uL, double uR) {
    double fL = 0.5 * uL * uL;
    double fR = 0.5 * uR * uR;
    double C = std::max(std::abs(uL), std::abs(uR));
    return 0.5 * (fL + fR) - 0.5 * C * (uR - uL);
}

// ----------------------------------------------------------------------
// Вычисление правых частей для DG P1
// ----------------------------------------------------------------------
void compute_dg_rhs(const DG_Solution& u,
                    std::vector<double>& dmean,
                    std::vector<double>& dslope,
                    bool with_source, double time) {
    int M = u.M;
    double h = u.h;
    dmean.assign(M, 0.0);
    dslope.assign(M, 0.0);

    for (int j = 0; j < M; ++j) {
        // Значения на границах элемента (внутренние)
        double u_left  = u.u_mean[j] - u.u_slope[j];
        double u_right = u.u_mean[j] + u.u_slope[j];

        // Внешние значения (соседние элементы или граница)
        double u_ext_left  = (j == 0) ? accuracy_u(time, 0.0)
                                      : u.u_mean[j-1] + u.u_slope[j-1];
        double u_ext_right = (j == M-1) ? accuracy_u(time, u.X)
                                        : u.u_mean[j+1] - u.u_slope[j+1];

        // Численные потоки
        double flux_left  = llf_flux(u_ext_left, u_left);
        double flux_right = llf_flux(u_right, u_ext_right);

        // Уравнение для среднего (без источника)
        dmean[j] = -(flux_right - flux_left) / h;

        // Интеграл от 0.5*u^2 внутри элемента
        double integral_f = 0.5 * h * (u.u_mean[j]*u.u_mean[j] + u.u_slope[j]*u.u_slope[j]/3.0);

        // Уравнение для наклона (без источника)
        dslope[j] = -3.0 / h * ( (flux_right + flux_left) - 2.0 / h * integral_f );

        // ----- добавка от источника (manufactured solution) -----
        if (with_source) {
            double xL = j * h;
            double xR = (j + 1) * h;
            double fL = accuracy_f(time, xL);
            double fR = accuracy_f(time, xR);

            // Источник для среднего: ∫ f φ₀ dx / h  (точное интегрирование для линейной f)
            // даёт среднее арифметическое на концах
            dmean[j] += 0.5 * (fL + fR);

            // Источник для наклона: (3/h) * ∫ f φ₁ dx   (аналитически для линейной f)
            // ∫ f φ₁ dx = (h/2)*(fR - fL)/3 → после умножения на (3/h) получаем (fR - fL)/2
            dslope[j] += 0.5 * (fR - fL);
        }
    }
}

// ----------------------------------------------------------------------
// Ограничитель minmod (с TVB-модификацией)
// ----------------------------------------------------------------------
double minmod(double a, double b, double c) {
    if (a > 0.0 && b > 0.0 && c > 0.0) return std::min({a, b, c});
    if (a < 0.0 && b < 0.0 && c < 0.0) return std::max({a, b, c});
    return 0.0;
}

void apply_limiter(DG_Solution& u, bool TVB = true, double M_tvb = 50.0) {
    int M = u.M;
    double h = u.h;
    std::vector<double> new_slope(M, 0.0);
    for (int j = 1; j < M - 1; ++j) {
        double left_mean = u.u_mean[j - 1];
        double right_mean = u.u_mean[j + 1];
        double center_mean = u.u_mean[j];
        double diff_left = center_mean - left_mean;
        double diff_right = right_mean - center_mean;
        double original_slope = u.u_slope[j];

        // Исправленное TVB-условие: сравниваем разности средних с M*h^2
        if (TVB) {
            if (std::abs(diff_left) <= M_tvb * h * h && std::abs(diff_right) <= M_tvb * h * h) {
                new_slope[j] = original_slope;  // не ограничиваем
                continue;
            }
        }
        new_slope[j] = minmod(original_slope, diff_left, diff_right);
    }
    // На крайних элементах оставляем наклон без изменений
    new_slope[0] = u.u_slope[0];
    new_slope[M - 1] = u.u_slope[M - 1];
    u.u_slope = new_slope;
}

// ----------------------------------------------------------------------
// Один шаг SSP RK3 (TVD RK3)
// ----------------------------------------------------------------------
void dg_step(DG_Solution& u, double dt, double time, bool with_source) {
    DG_Solution u1 = u;  // копия
    std::vector<double> dmean, dslope;

    // Этап 1
    compute_dg_rhs(u, dmean, dslope, with_source, time);
    for (int j = 0; j < u.M; ++j) {
        u1.u_mean[j] = u.u_mean[j] + dt * dmean[j];
        u1.u_slope[j] = u.u_slope[j] + dt * dslope[j];
    }
    apply_limiter(u1);

    // Этап 2
    DG_Solution u2 = u1;
    compute_dg_rhs(u1, dmean, dslope, with_source, time + dt);
    for (int j = 0; j < u.M; ++j) {
        u2.u_mean[j] = 0.75 * u.u_mean[j] + 0.25 * u1.u_mean[j] + 0.25 * dt * dmean[j];
        u2.u_slope[j] = 0.75 * u.u_slope[j] + 0.25 * u1.u_slope[j] + 0.25 * dt * dslope[j];
    }
    apply_limiter(u2);

    // Этап 3
    compute_dg_rhs(u2, dmean, dslope, with_source, time + 0.5 * dt);
    for (int j = 0; j < u.M; ++j) {
        u.u_mean[j] =
            (1.0 / 3.0) * u.u_mean[j] + (2.0 / 3.0) * u2.u_mean[j] + (2.0 / 3.0) * dt * dmean[j];
        u.u_slope[j] =
            (1.0 / 3.0) * u.u_slope[j] + (2.0 / 3.0) * u2.u_slope[j] + (2.0 / 3.0) * dt * dslope[j];
    }
    apply_limiter(u);
}

// ----------------------------------------------------------------------
// main
// ----------------------------------------------------------------------
int main() {
    double X = 1.0, T = 1.0;
    int N = 1000, M = 100;
    std::cout << "Enter N (time steps), T, M (elements), X: ";
    std::cin >> N >> T >> M >> X;
    double dt = T / N;
    double h = X / M;

    common_params params(N, dt, M, X);
    printf("DG P1: dt=%lf, h=%lf, M=%d, N=%d\n", dt, h, M, N);

    DG_Solution u(M, h, X);
    u.set_initial(0.0);

    std::cout << "Initial errors: L2=" << compute_L2_error(u, 0.0)
              << " Linf=" << compute_Linf_error(u, 0.0) << "\n";

    bool with_source = true;  // manufactured solution
    double t = 0.0;
    for (int step = 0; step < N; ++step) {
        dg_step(u, dt, t, with_source);
        t += dt;

        if (step % (N / 10) == 0) {
            std::cout << "Step " << step << " t=" << t << " L2=" << compute_L2_error(u, t)
                      << " Linf=" << compute_Linf_error(u, t) << "\n";
        }
    }

    std::cout << "\nFinal errors at t=" << T << ": L2=" << compute_L2_error(u, T)
              << " Linf=" << compute_Linf_error(u, T) << std::endl;

    return 0;
}