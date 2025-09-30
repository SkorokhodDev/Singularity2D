#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>

// imgui
#include "../dependencies/imgui-master/imgui.h"
#include "../dependencies/imgui-master/backends/imgui_impl_glfw.h"
#include "../dependencies/imgui-master/backends/imgui_impl_opengl3.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace glm;
using namespace std;

constexpr double c = 299792458.0;
constexpr double G = 6.67430e-11;

struct Ray;

void rk4Step(Ray& ray, double dLA, double rs);
void eulerStep(Ray& ray, double dLA, double rs);


// --- Structs --- //
struct Engine {

    GLFWwindow* window;
    int WIDTH = 1200;
    int HEIGHT = 800;
    float width = 100000000000.0f; // Width of the viewport in meters
    float height = 75000000000.0f; // Height of the viewport in meters

    // Navigation state
    float offsetX = 0.0f, offsetY = 0.0f;
    // TODO: release it
    float zoom = 1.0f;
    bool middleMousePressed = false;
    double lastMouseX = 0.0, lastMouseY = 0;

    Engine() {
        if (!glfwInit()) {
            cerr << "Failed to initialize GLFW" << endl;
            exit(EXIT_FAILURE);
        }

        window = glfwCreateWindow(WIDTH, HEIGHT, "Black Hole Simulation", NULL, NULL);

        if (!window) {
            cerr << "Failed to create GLFW window" << endl;
            glfwTerminate();
            exit(EXIT_FAILURE);
        }

        glfwMakeContextCurrent(window);

        glewExperimental = GL_TRUE;

        if (glewInit() != GLEW_OK) {
            cerr << "Failed to initialize GLEW" << endl;
            glfwDestroyWindow(window);
            glfwTerminate();
            exit(EXIT_FAILURE);
        }

        glViewport(0, 0, WIDTH, HEIGHT);;
    }

    void run() {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        double left = -width + offsetX;
        double right = width + offsetX;
        double bottom = -height + offsetY;
        double top = height + offsetY;
        glOrtho(left, right, bottom, top, -1.0, 1.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    }
};

Engine engine;

struct BlackHole {
    vec3 position;
    double mass;
    double radius;
    double r_s;

    BlackHole(vec3 pos, float m) : position(pos), mass(m) { r_s = 2.0 * G * mass / (c * c); }
    void draw() {
        glBegin(GL_TRIANGLE_FAN);
        glColor3f(1.0f, 1.0f, 0.0f);               // Red color for the black hole
        glVertex2f(0.0f, 0.0f);                    // Center
        for (int i = 0; i <= 100; i++) {
            float angle = 2.0f * M_PI * i / 100;
            float x = r_s * cos(angle); // Radius of 0.1
            float y = r_s * sin(angle);
            glVertex2f(x, y);
        }
        glEnd();
    }
};

BlackHole SagA(vec3(0.0f, 0.0f, 0.0f), 8.54e36); // Sagittarius A black hole


struct Ray {
    // -- cartesian coords -- //
    double x;   double y;
    // -- polar coords -- //
    double r;   double phi;
    double dr;  double dphi;
    vector<vec2> trail; // trail of points
    double E, L;             // conserved quantities

    Ray(vec2 pos, vec2 dir) : x(pos.x), y(pos.y), r(sqrt(pos.x* pos.x + pos.y * pos.y)), phi(atan2(pos.y, pos.x)), dr(dir.x), dphi(dir.y) {
        // step 1) get polar coords (r, phi) :
        this->r = sqrt(x * x + y * y);
        this->phi = atan2(y, x);

        // step 2) seed velocities :
        dr = dir.x * cos(phi) + dir.y * sin(phi); // m/s
        dphi = (-dir.x * sin(phi) + dir.y * cos(phi)) / r;

        // step 3) store conserved quantities
        L = r * r * dphi;
        double f = 1.0 - SagA.r_s / r;
        double dt_dLA = sqrt((dr * dr) / (f * f) + (r * r * dphi * dphi) / f);
        E = f * dt_dLA;

        // step 4) start trail :
        trail.push_back({ x, y });
    }

    void draw(const std::vector<Ray>& rays) {
        // draw current ray positions as points
        glPointSize(3.0f);
        glColor3f(1.0f, 0.0f, 0.0f);
        glBegin(GL_POINTS);
        for (const auto& ray : rays) {
            glVertex2f(ray.x, ray.y);
        }
        glEnd();

        // turn on blending for the trails
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glLineWidth(1.0f);

        // draw each trail with fading alpha
        for (const auto& ray : rays) {
            size_t N = ray.trail.size();
            if (N < 2) continue;

            glBegin(GL_LINE_STRIP);
            for (size_t i = 0; i < N; ++i) {
                // older points (i=0) get alpha≈0, newer get alpha≈1
                float alpha = float(i) / float(N - 1);
                glColor4f(1.0f, 1.0f, 1.0f, std::max(alpha, 0.05f));
                glVertex2f(ray.trail[i].x, ray.trail[i].y);
            }
            glEnd();
        }

        glDisable(GL_BLEND);
    }

    void step(double dLA, double rs) {
        // 1) integrate (r,φ,dr,dφ)
        if (r <= rs) return; // stop if inside the event horizon
        rk4Step(*this, dLA, rs);

        // 2) convert back to cartesian x,y
        x = r * cos(phi);
        y = r * sin(phi);

        // 3) record the trail
        trail.push_back({ float(x), float(y) });
    }
};

vector<Ray> rays;

void geodesicRHS(const Ray& ray, double rhs[4], double rs) {
    double r = ray.r;
    double dr = ray.dr;
    double dphi = ray.dphi;
    double E = ray.E;

    double f = 1.0 - rs / r;

    // dr/dLA = dr
    rhs[0] = dr;
    // dφ/dLA = dphi
    rhs[1] = dphi;

    // d²r/dLA² from Schwarzschild null geodesic:
    double dt_dLA = E / f;
    rhs[2] =
        -(rs / (2 * r * r)) * f * (dt_dLA * dt_dLA)
        + (rs / (2 * r * r * f)) * (dr * dr)
        + (r - rs) * (dphi * dphi);

    // d²φ/dLA² = -2*(dr * dphi) / r
    rhs[3] = -2.0 * dr * dphi / r;
}

void addState(const double a[4], const double b[4], double factor, double out[4]) {
    for (int i = 0; i < 4; i++)
        out[i] = a[i] + b[i] * factor;
}

void rk4Step(Ray& ray, double dLA, double rs) {
    double y0[4] = { ray.r, ray.phi, ray.dr, ray.dphi };
    double k1[4], k2[4], k3[4], k4[4], temp[4];

    geodesicRHS(ray, k1, rs);
    addState(y0, k1, dLA / 2.0, temp);
    Ray r2 = ray; r2.r = temp[0]; r2.phi = temp[1]; r2.dr = temp[2]; r2.dphi = temp[3];
    geodesicRHS(r2, k2, rs);

    addState(y0, k2, dLA / 2.0, temp);
    Ray r3 = ray; r3.r = temp[0]; r3.phi = temp[1]; r3.dr = temp[2]; r3.dphi = temp[3];
    geodesicRHS(r3, k3, rs);

    addState(y0, k3, dLA, temp);
    Ray r4 = ray; r4.r = temp[0]; r4.phi = temp[1]; r4.dr = temp[2]; r4.dphi = temp[3];
    geodesicRHS(r4, k4, rs);

    ray.r += (dLA / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    ray.phi += (dLA / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    ray.dr += (dLA / 6.0) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
    ray.dphi += (dLA / 6.0) * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]);
}

void eulerStep(Ray& ray, double dLA, double rs) {
    double k1[4];
    geodesicRHS(ray, k1, rs);

    ray.r += dLA * k1[0];
    ray.phi += dLA * k1[1];
    ray.dr += dLA * k1[2];
    ray.dphi += dLA * k1[3];
}

void rk2Step(Ray& ray, double dLA, double rs) {
    double y0[4] = { ray.r, ray.phi, ray.dr, ray.dphi };
    double k1[4], k2[4], temp[4];

    geodesicRHS(ray, k1, rs);
    addState(y0, k1, dLA / 2.0, temp);

    Ray r2 = ray;
    r2.r = temp[0]; r2.phi = temp[1]; r2.dr = temp[2]; r2.dphi = temp[3];
    geodesicRHS(r2, k2, rs);

    ray.r += dLA * k2[0];
    ray.phi += dLA * k2[1];
    ray.dr += dLA * k2[2];
    ray.dphi += dLA * k2[3];
}

// y0 — базовое состояние (4 элемента)
// ks — набор массивов k1..kN
// coeffs — набор коэффициентов (того же размера, что ks)
// factor — множитель шага (обычно dLA)
// out — результат
void combineState(const double y0[4],
    const std::vector<double*>& ks,
    const std::vector<double>& coeffs,
    double factor,
    double out[4])
{
    for (int i = 0; i < 4; i++) {
        double sum = 0.0;
        for (size_t j = 0; j < ks.size(); j++)
            sum += coeffs[j] * ks[j][i];
        out[i] = y0[i] + factor * sum;
    }
}


void rkf45Step(Ray& ray, double& dLA, double rs, double tol) {
    // коэффициенты Фельберга
    const double b21 = 1.0 / 4.0;
    const double b31 = 3.0 / 32.0, b32 = 9.0 / 32.0;
    const double b41 = 1932.0 / 2197.0, b42 = -7200.0 / 2197.0, b43 = 7296.0 / 2197.0;
    const double b51 = 439.0 / 216.0, b52 = -8.0, b53 = 3680.0 / 513.0, b54 = -845.0 / 4104.0;
    const double b61 = -8.0 / 27.0, b62 = 2.0, b63 = -3544.0 / 2565.0, b64 = 1859.0 / 4104.0, b65 = -11.0 / 40.0;

    const double c1 = 16.0 / 135.0, c3 = 6656.0 / 12825.0, c4 = 28561.0 / 56430.0, c5 = -9.0 / 50.0, c6 = 2.0 / 55.0; // 5й порядок
    const double d1 = 25.0 / 216.0, d3 = 1408.0 / 2565.0, d4 = 2197.0 / 4104.0, d5 = -1.0 / 5.0;                    // 4й порядок

    double y0[4] = { ray.r, ray.phi, ray.dr, ray.dphi };
    double k1[4], k2[4], k3[4], k4[4], k5[4], k6[4];
    double temp[4];

    // стадия 1
    geodesicRHS(ray, k1, rs);

    // стадия 2
    combineState(y0, { k1 }, { b21 }, dLA, temp);
    Ray r2 = ray; r2.r = temp[0]; r2.phi = temp[1]; r2.dr = temp[2]; r2.dphi = temp[3];
    geodesicRHS(r2, k2, rs);

    // стадия 3
    combineState(y0, { k1, k2 }, { b31, b32 }, dLA, temp);
    Ray r3 = ray; r3.r = temp[0]; r3.phi = temp[1]; r3.dr = temp[2]; r3.dphi = temp[3];
    geodesicRHS(r3, k3, rs);

    // стадия 4
    combineState(y0, { k1, k2, k3 }, { b41, b42, b43 }, dLA, temp);
    Ray r4 = ray; r4.r = temp[0]; r4.phi = temp[1]; r4.dr = temp[2]; r4.dphi = temp[3];
    geodesicRHS(r4, k4, rs);

    // стадия 5
    combineState(y0, { k1, k2, k3, k4 }, { b51, b52, b53, b54 }, dLA, temp);
    Ray r5 = ray; r5.r = temp[0]; r5.phi = temp[1]; r5.dr = temp[2]; r5.dphi = temp[3];
    geodesicRHS(r5, k5, rs);

    // стадия 6
    combineState(y0, { k1, k2, k3, k4, k5 }, { b61, b62, b63, b64, b65 }, dLA, temp);
    Ray r6 = ray; r6.r = temp[0]; r6.phi = temp[1]; r6.dr = temp[2]; r6.dphi = temp[3];
    geodesicRHS(r6, k6, rs);

    // решения 4-го и 5-го порядка
    double y4[4], y5[4];
    for (int i = 0; i < 4; i++) {
        y4[i] = y0[i] + dLA * (d1 * k1[i] + d3 * k3[i] + d4 * k4[i] + d5 * k5[i]);
        y5[i] = y0[i] + dLA * (c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c5 * k5[i] + c6 * k6[i]);
    }

    // ошибка
    double err = 0.0;
    for (int i = 0; i < 4; i++)
        err = std::max(err, fabs(y5[i] - y4[i]));

    // адаптивный шаг
    if (err < tol || dLA < 1e-12) {
        // принять шаг
        ray.r = y5[0];
        ray.phi = y5[1];
        ray.dr = y5[2];
        ray.dphi = y5[3];

        // обновить декартовы координаты
        ray.x = ray.r * cos(ray.phi);
        ray.y = ray.r * sin(ray.phi);

        // увеличить шаг
        if (err > 0.0)
            dLA *= std::min(2.0, 0.9 * pow(tol / err, 0.2));
    }
    else {
        // уменьшить шаг и попробовать снова
        dLA *= std::max(0.1, 0.9 * pow(tol / err, 0.25));
        rkf45Step(ray, dLA, rs, tol);
    }
}



void initializeRays()
{
    //rays.push_back(Ray(vec2(-1e11, 3.27606302719999999e10), vec2(c, 0.0f)));
	rays.clear();

    int numRays = 30;
    double x0 = -1e11;
    double yMin = -8e10;
    double yMax = 8e10;

    for (int i = 0; i < numRays; ++i) {
        double y = yMin + (yMax - yMin) * i / (numRays - 1);
        vec2 pos = vec2((float)x0, (float)y);
        vec2 dir = vec2((float)c, 0.0f);
        rays.push_back(Ray(pos, dir));
    }
}

struct UserInterface {

    ImGuiIO* io = nullptr;

    bool show_demo_window = true;
    bool show_another_window = false;
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    UserInterface() {

        // ImGui setup
        ImGui::CreateContext();
        io = &ImGui::GetIO();
        ImGui_ImplGlfw_InitForOpenGL(engine.window, true);
        ImGui_ImplOpenGL3_Init("#version 130"); // Шейдерная версия

        ImGui::StyleColorsDark();
    }

    void start_frame() {
        // start the ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        {
            static float f = 0.0f;
            static int counter = 0;

            ImGui::Begin("Hello, world!");                          // Create a window called "Hello, world!" and append into it.

            ImGui::Text("This is some useful text.");               // Display some text (you can use a format strings too)
            ImGui::Checkbox("Demo Window", &show_demo_window);      // Edit bools storing our window open/close state
            ImGui::Checkbox("Another Window", &show_another_window);

            ImGui::SliderFloat("float", &f, 0.0f, 1.0f);            // Edit 1 float using a slider from 0.0f to 1.0f
            ImGui::ColorEdit3("clear color", (float*)&clear_color); // Edit 3 floats representing a color

            if (ImGui::Button("Button"))                            // Buttons return true when clicked (most widgets return true when edited/activated)
            {
                counter++;
                initializeRays();
            }
            ImGui::SameLine();
            ImGui::Text("counter = %d", counter);

            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io->Framerate, io->Framerate);
            ImGui::End();
        }
    }
    void render_frame() {
        // Rendering
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    }

    ~UserInterface() {
        // Cleanup
        ImGui_ImplOpenGL3_Shutdown();
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
    }

};

UserInterface ui;

int main() {

	initializeRays();

    while (!glfwWindowShouldClose(engine.window)) {
		
        ui.start_frame();

        // Render
        engine.run();
        SagA.draw();

        for (auto& ray : rays) {
            ray.step(1.0f, SagA.r_s);
            ray.draw(rays);
        }

		ui.render_frame();

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    return 0;
}


