//
//
//
//void eulerStep(Ray& ray, double dLA, double rs) {
//    double k1[4];
//    geodesicRHS(ray, k1, rs);
//
//    ray.r += dLA * k1[0];
//    ray.phi += dLA * k1[1];
//    ray.dr += dLA * k1[2];
//    ray.dphi += dLA * k1[3];
//}
//
//void rk2Step(Ray& ray, double dLA, double rs) {
//    double y0[4] = { ray.r, ray.phi, ray.dr, ray.dphi };
//    double k1[4], k2[4], temp[4];
//
//    geodesicRHS(ray, k1, rs);
//    addState(y0, k1, dLA / 2.0, temp);
//
//    Ray r2 = ray;
//    r2.r = temp[0]; r2.phi = temp[1]; r2.dr = temp[2]; r2.dphi = temp[3];
//    geodesicRHS(r2, k2, rs);
//
//    ray.r += dLA * k2[0];
//    ray.phi += dLA * k2[1];
//    ray.dr += dLA * k2[2];
//    ray.dphi += dLA * k2[3];
//}
//
//// y0 Ч базовое состо€ние (4 элемента)
//// ks Ч набор массивов k1..kN
//// coeffs Ч набор коэффициентов (того же размера, что ks)
//// factor Ч множитель шага (обычно dLA)
//// out Ч результат
//void combineState(const double y0[4],
//    const std::vector<double*>& ks,
//    const std::vector<double>& coeffs,
//    double factor,
//    double out[4])
//{
//    for (int i = 0; i < 4; i++) {
//        double sum = 0.0;
//        for (size_t j = 0; j < ks.size(); j++)
//            sum += coeffs[j] * ks[j][i];
//        out[i] = y0[i] + factor * sum;
//    }
//}
//
//
//void rkf45Step(Ray& ray, double& dLA, double rs, double tol) {
//    // коэффициенты ‘ельберга
//    const double b21 = 1.0 / 4.0;
//    const double b31 = 3.0 / 32.0, b32 = 9.0 / 32.0;
//    const double b41 = 1932.0 / 2197.0, b42 = -7200.0 / 2197.0, b43 = 7296.0 / 2197.0;
//    const double b51 = 439.0 / 216.0, b52 = -8.0, b53 = 3680.0 / 513.0, b54 = -845.0 / 4104.0;
//    const double b61 = -8.0 / 27.0, b62 = 2.0, b63 = -3544.0 / 2565.0, b64 = 1859.0 / 4104.0, b65 = -11.0 / 40.0;
//
//    const double c1 = 16.0 / 135.0, c3 = 6656.0 / 12825.0, c4 = 28561.0 / 56430.0, c5 = -9.0 / 50.0, c6 = 2.0 / 55.0; // 5й пор€док
//    const double d1 = 25.0 / 216.0, d3 = 1408.0 / 2565.0, d4 = 2197.0 / 4104.0, d5 = -1.0 / 5.0;                    // 4й пор€док
//
//    double y0[4] = { ray.r, ray.phi, ray.dr, ray.dphi };
//    double k1[4], k2[4], k3[4], k4[4], k5[4], k6[4];
//    double temp[4];
//
//    // стади€ 1
//    geodesicRHS(ray, k1, rs);
//
//    // стади€ 2
//    combineState(y0, { k1 }, { b21 }, dLA, temp);
//    Ray r2 = ray; r2.r = temp[0]; r2.phi = temp[1]; r2.dr = temp[2]; r2.dphi = temp[3];
//    geodesicRHS(r2, k2, rs);
//
//    // стади€ 3
//    combineState(y0, { k1, k2 }, { b31, b32 }, dLA, temp);
//    Ray r3 = ray; r3.r = temp[0]; r3.phi = temp[1]; r3.dr = temp[2]; r3.dphi = temp[3];
//    geodesicRHS(r3, k3, rs);
//
//    // стади€ 4
//    combineState(y0, { k1, k2, k3 }, { b41, b42, b43 }, dLA, temp);
//    Ray r4 = ray; r4.r = temp[0]; r4.phi = temp[1]; r4.dr = temp[2]; r4.dphi = temp[3];
//    geodesicRHS(r4, k4, rs);
//
//    // стади€ 5
//    combineState(y0, { k1, k2, k3, k4 }, { b51, b52, b53, b54 }, dLA, temp);
//    Ray r5 = ray; r5.r = temp[0]; r5.phi = temp[1]; r5.dr = temp[2]; r5.dphi = temp[3];
//    geodesicRHS(r5, k5, rs);
//
//    // стади€ 6
//    combineState(y0, { k1, k2, k3, k4, k5 }, { b61, b62, b63, b64, b65 }, dLA, temp);
//    Ray r6 = ray; r6.r = temp[0]; r6.phi = temp[1]; r6.dr = temp[2]; r6.dphi = temp[3];
//    geodesicRHS(r6, k6, rs);
//
//    // решени€ 4-го и 5-го пор€дка
//    double y4[4], y5[4];
//    for (int i = 0; i < 4; i++) {
//        y4[i] = y0[i] + dLA * (d1 * k1[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i]);
//        y5[i] = y0[i] + dLA * (c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c5 * k5[i] + c6 * k6[i]);
//    }
//
//    // ошибка
//    double err = 0.0;
//    for (int i = 0; i < 4; i++)
//        err = std::max(err, fabs(y5[i] - y4[i]));
//
//    // адаптивный шаг
//    if (err < tol || dLA < 1e-12) {
//        // прин€ть шаг
//        ray.r = y5[0];
//        ray.phi = y5[1];
//        ray.dr = y5[2];
//        ray.dphi = y5[3];
//
//        // обновить декартовы координаты
//        ray.x = ray.r * cos(ray.phi);
//        ray.y = ray.r * sin(ray.phi);
//
//        // увеличить шаг
//        if (err > 0.0)
//            dLA *= std::min(2.0, 0.9 * pow(tol / err, 0.2));
//    }
//    else {
//        // уменьшить шаг и попробовать снова
//        dLA *= std::max(0.1, 0.9 * pow(tol / err, 0.25));
//        rkf45Step(ray, dLA, rs, tol);
//    }
//}