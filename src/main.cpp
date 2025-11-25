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


void eulerStep(Ray& ray, double dLambda, double rs);
void rk2Step(Ray& ray, double dLambda, double rs);
void rk4Step(Ray& ray, double dLambda, double rs);
void rk4Step_optimized(Ray& ray, double dLambda, double rs);

void geodesicRHS_Raw(const double state[4], double E, double rs, double rhs[4]);
void geodesicRHS(const Ray& ray, double rhs[4], double rs);

enum class RayIntegral {
    Euler = 0,
    RK2,
    RK4,
    RKF45,
    Test,
    Count
};

struct IntegratorSettings {
    RayIntegral type;
    const char* name; // ("Euler")
    bool enabled;
    float color[3];   // RGB (0..1) for ImGui
};

IntegratorSettings solversConfig[] = {
    { RayIntegral::Euler, "Euler",  false, {0.0f, 1.0f, 0.0f} }, // Green
    { RayIntegral::RK2,   "RK2",    false, {0.0f, 0.5f, 1.0f} }, // Blue
    { RayIntegral::RK4,   "RK4",    false,  {1.0f, 1.0f, 1.0f} }, // White
	{ RayIntegral::RKF45, "RKF45",  false, {1.0f, 0.0f, 1.0f} },  // Purple
	{ RayIntegral::Test, "Test",  true, {1.0f, 0.0f, 1.0f} }  // Purple
};

// --- Structs --- //
struct Engine {

    GLFWwindow* window;
    int WIDTH = 1600;
    int HEIGHT = 1200;
    float width = 100000000000.0f; // Width of the viewport in meters
    float height = 75000000000.0f; // Height of the viewport in meters

    // Navigation state
    float offsetX = 0.0f, offsetY = 0.0f;
    // TODO: release it
    float zoom = 1.0f;
    bool middleMousePressed = false;
    double lastMouseX = 0.0, lastMouseY = 0;

    // App state
	bool bProcessing = false;

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

    // RGBA color (default = red)
    //ImVec4 body_color = ImVec4(1.0f, 1.0f, 0.0f, 1.0f);


    BlackHole(vec3 pos, float m) : position(pos), mass(m) {
        r_s = 2.0 * G * mass / (c * c); 
    }

    void draw(ImVec4& body_color) {
        glBegin(GL_TRIANGLE_FAN);
        glColor3f(body_color.x, body_color.y, body_color.z);               // Red color for the black hole
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
	// -- Ray type -- //
	RayIntegral type = RayIntegral::RK4;
    // -- cartesian coords -- //
    double x;   double y;
    // -- polar coords -- //
    double r;   double phi;
    double dr;  double dphi;
    vector<vec2> trail; // trail of points
    double E, L;             // conserved quantities

	vec3 trailColor = vec3(1.0f, 1.0f, 1.0f); // default white

	// Adaptive step size control variables (for RKF45)
    double currentStepSize = 1.0; 
    double tolerance = 1e-5;

	Ray(vec2 pos, vec2 dir, RayIntegral type, vec3 trailColor) : x(pos.x), y(pos.y), type{ type }, trailColor{ trailColor },
        r(sqrt(pos.x* pos.x + pos.y * pos.y)), phi(atan2(pos.y, pos.x)), dr(dir.x), dphi(dir.y) {
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



    void step(double dLambda, double rs)
    {
        // 1) integrate (r,φ,dr,dφ)
        if (r <= rs) return; // stop if inside the event horizon

        if (type == RayIntegral::Euler) {
            eulerStep(*this, dLambda, rs);
        }
        else if (type == RayIntegral::RK2) {
            rk2Step(*this, dLambda, rs);
		}
        else if (type == RayIntegral::RK4) {
            rk4Step(*this, dLambda, rs);
        }
        else if (type == RayIntegral::RKF45) {
            // Для RKF45 мы игнорируем входящий dLambda или используем его как "максимум"
            // Мы вызываем внутреннюю функцию, которая сама обновит this->currentStepSize
            rkf45AdaptiveStep(rs);
        }
        else if (type == RayIntegral::Test) {
            rk4Step_optimized(*this, dLambda, rs);
		}

        //dr += d2r * dLambda;
        //dphi += d2phi * dLambda;
        //r += dr * dLambda;
        //phi += dphi * dLambda;

        // 2) convert back to cartesian x,y
        x = r * cos(phi);
        y = r * sin(phi);

        // 3) record the trail
        trail.push_back({ float(x), float(y) });
    }


    // Реализация адаптивного шага (внешняя обертка)
    void rkf45AdaptiveStep(double rs) {
        // Пытаемся сделать шаг. Если ошибка велика, rkf45SingleStep вернет false 
        // и уменьшит currentStepSize. Мы будем повторять, пока шаг не будет принят.

        // Ограничитель, чтобы не зависнуть, если шаг станет бесконечно малым
        int maxTries = 15;
        bool success = false;

        while (!success && maxTries > 0) {
            success = rkf45SingleStep(currentStepSize, rs);
            if (!success) {
                maxTries--;
            }
        }
    }

    // Core of RKF45 (Cash-Karp)
    bool rkf45SingleStep(double h, double rs) {
        static const double b21 = 0.2, b31 = 3.0 / 40.0, b32 = 9.0 / 40.0;
        static const double b41 = 0.3, b42 = -0.9, b43 = 1.2;
        static const double b51 = -11.0 / 54.0, b52 = 2.5, b53 = -70.0 / 27.0, b54 = 35.0 / 27.0;
        static const double b61 = 1631.0 / 55296.0, b62 = 175.0 / 512.0, b63 = 575.0 / 13824.0, b64 = 44275.0 / 110592.0, b65 = 253.0 / 4096.0;
        static const double c1 = 37.0 / 378.0, c3 = 250.0 / 621.0, c4 = 125.0 / 594.0, c6 = 512.0 / 1771.0;
        static const double dc1 = c1 - 2825.0 / 27648.0, dc3 = c3 - 18575.0 / 48384.0, dc4 = c4 - 13525.0 / 55296.0, dc5 = -277.0 / 14336.0, dc6 = c6 - 0.25;

        double y[4] = { r, phi, dr, dphi };
        double yTemp[4], k1[4], k2[4], k3[4], k4[4], k5[4], k6[4];

        geodesicRHS_Raw(y, E, rs, k1);

        for (int i = 0; i < 4; i++) yTemp[i] = y[i] + h * (b21 * k1[i]);
        geodesicRHS_Raw(yTemp, E, rs, k2);

        for (int i = 0; i < 4; i++) yTemp[i] = y[i] + h * (b31 * k1[i] + b32 * k2[i]);
        geodesicRHS_Raw(yTemp, E, rs, k3);

        for (int i = 0; i < 4; i++) yTemp[i] = y[i] + h * (b41 * k1[i] + b42 * k2[i] + b43 * k3[i]);
        geodesicRHS_Raw(yTemp, E, rs, k4);

        for (int i = 0; i < 4; i++) yTemp[i] = y[i] + h * (b51 * k1[i] + b52 * k2[i] + b53 * k3[i] + b54 * k4[i]);
        geodesicRHS_Raw(yTemp, E, rs, k5);

        for (int i = 0; i < 4; i++) yTemp[i] = y[i] + h * (b61 * k1[i] + b62 * k2[i] + b63 * k3[i] + b64 * k4[i] + b65 * k5[i]);
        geodesicRHS_Raw(yTemp, E, rs, k6);

        // Считаем ошибку
        double error = 0.0;
        for (int i = 0; i < 4; i++) {
            double errI = h * (dc1 * k1[i] + dc3 * k3[i] + dc4 * k4[i] + dc5 * k5[i] + dc6 * k6[i]);
            error += errI * errI;
        }
        error = sqrt(error);

        if (error <= tolerance) { // Шаг принят
            
            for (int i = 0; i < 4; i++) y[i] += h * (c1 * k1[i] + c3 * k3[i] + c4 * k4[i] + c6 * k6[i]);
            r = y[0]; phi = y[1]; dr = y[2]; dphi = y[3];

            // Немного увеличиваем шаг, но не более чем в 5 раз
            double factor = 0.9 * pow(tolerance / (error + 1e-30), 0.2);
            if (factor > 5.0) factor = 5.0;
            currentStepSize *= factor;
            return true;
        }
        else {
            // Шаг отвергнут, уменьшаем размер
            double factor = 0.9 * pow(tolerance / (error + 1e-30), 0.25);
            if (factor < 0.1) factor = 0.1;
            currentStepSize *= factor;
            return false;
        }
    }
};

vector<Ray> rays;

// Легкая функция для расчета производных без копирования лучей
// state[0]=r, state[1]=phi, state[2]=dr, state[3]=dphi
void geodesicRHS_Raw(const double state[4], double E, double rs, double rhs[4]) {
    double r = state[0];
    double dr = state[2];
    double dphi = state[3];

    double f = 1.0 - rs / r;
    if (f < 1e-9) f = 1e-9; // Защита от деления на ноль у горизонта

    // dr/dLambda
    rhs[0] = dr;
    // dφ/dLambda
    rhs[1] = dphi;

    // d²r/dLambda²
    double dt_dLA = E / f;
    rhs[2] = -(rs / (2.0 * r * r)) * f * (dt_dLA * dt_dLA)
        + (rs / (2.0 * r * r * f)) * (dr * dr)
        + (r - rs) * (dphi * dphi);

    // d²φ/dLambda²
    rhs[3] = -2.0 * dr * dphi / r;
}

void geodesicRHS(const Ray& ray, double rhs[4], double rs) {
    double r = ray.r;
    double dr = ray.dr;
    double dphi = ray.dphi;
    double E = ray.E;

    double f = 1.0 - rs / r;

    // dr/dLambda = dr
    rhs[0] = dr;
    // dφ/dLambda = dphi
    rhs[1] = dphi;

    // d²r/dLambda² from Schwarzschild null geodesic:
    double dt_dLA = E / f;
    rhs[2] =
        -(rs / (2 * r * r)) * f * (dt_dLA * dt_dLA)
        + (rs / (2 * r * r * f)) * (dr * dr)
        + (r - rs) * (dphi * dphi);

    // d²φ/dLambda² = -2*(dr * dphi) / r
    rhs[3] = -2.0 * dr * dphi / r;
}




void addState(const double a[4], const double b[4], double factor, double out[4]) {
    for (int i = 0; i < 4; i++)
        out[i] = a[i] + b[i] * factor;
}

void eulerStep(Ray& ray, double dLambda, double rs) {
    double k1[4];
    geodesicRHS(ray, k1, rs);

    ray.r += dLambda * k1[0];
    ray.phi += dLambda * k1[1];
    ray.dr += dLambda * k1[2];
    ray.dphi += dLambda * k1[3];
}

void rk2Step(Ray& ray, double dLambda, double rs) {
    double y[4] = { ray.r, ray.phi, ray.dr, ray.dphi };
    double k1[4], k2[4];
    double yHalf[4];

    // K1: Leaning at the beginning of a step
    geodesicRHS_Raw(y, ray.E, rs, k1);

    // tentative step towards the middle (h/2)
    double halfStep = dLambda * 0.5;
    for (int i = 0; i < 4; i++) {
        yHalf[i] = y[i] + halfStep * k1[i];
    }

    // K2
    geodesicRHS_Raw(yHalf, ray.E, rs, k2);

    ray.r += dLambda * k2[0];
    ray.phi += dLambda * k2[1];
    ray.dr += dLambda * k2[2];
    ray.dphi += dLambda * k2[3];
}

void rk4Step(Ray& ray, double dLambda, double rs) {
    double y0[4] = { ray.r, ray.phi, ray.dr, ray.dphi };
    double k1[4], k2[4], k3[4], k4[4], temp[4];

    // k1
    geodesicRHS(ray, k1, rs);
    addState(y0, k1, dLambda / 2.0, temp);
	// k2
    Ray r2 = ray; r2.r = temp[0]; r2.phi = temp[1]; r2.dr = temp[2]; r2.dphi = temp[3];
    geodesicRHS(r2, k2, rs);
    addState(y0, k2, dLambda / 2.0, temp);
	// k3
    Ray r3 = ray; r3.r = temp[0]; r3.phi = temp[1]; r3.dr = temp[2]; r3.dphi = temp[3];
    geodesicRHS(r3, k3, rs);
    addState(y0, k3, dLambda, temp);
	// k4
    Ray r4 = ray; r4.r = temp[0]; r4.phi = temp[1]; r4.dr = temp[2]; r4.dphi = temp[3];
    geodesicRHS(r4, k4, rs);

    ray.r += (dLambda / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    ray.phi += (dLambda / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    ray.dr += (dLambda / 6.0) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
    ray.dphi += (dLambda / 6.0) * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]);
}

void rk4Step_optimized(Ray& ray, double dLambda, double rs) {
    double y[4] = { ray.r, ray.phi, ray.dr, ray.dphi };
    double k1[4], k2[4], k3[4], k4[4];
    double yTemp[4];

    // K1
    geodesicRHS_Raw(y, ray.E, rs, k1);

    // K2
    for (int i = 0; i < 4; i++) yTemp[i] = y[i] + (dLambda / 2.0) * k1[i];
    geodesicRHS_Raw(yTemp, ray.E, rs, k2);

    // K3
    for (int i = 0; i < 4; i++) yTemp[i] = y[i] + (dLambda / 2.0) * k2[i];
    geodesicRHS_Raw(yTemp, ray.E, rs, k3);

    // K4
    for (int i = 0; i < 4; i++) yTemp[i] = y[i] + dLambda * k3[i];
    geodesicRHS_Raw(yTemp, ray.E, rs, k4);

    // Итог
    ray.r += (dLambda / 6.0) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    ray.phi += (dLambda / 6.0) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    ray.dr += (dLambda / 6.0) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
    ray.dphi += (dLambda / 6.0) * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]);
}


void createRays(int numRays, double x0, double yMin, double yMax, RayIntegral type, vec3 rayColor)
{
    for (int i = 0; i < numRays; ++i) {
        double y = yMin + (yMax - yMin) * i / (numRays - 1);
        vec2 pos = vec2((float)x0, (float)y);
        vec2 dir = vec2((float)c, 0.0f);
        rays.push_back(Ray(pos, dir, type, rayColor));
    }
}

void initializeRays() {
    //rays.push_back(Ray(vec2(-1e11, 3.27606302719999999e10), vec2(c, 0.0f)));
    rays.clear();

    int numRays = 30; // Или брать из UI
    double x0 = -1e11;
    double yMin = -8e10;
    double yMax = 8e10;


    for (const auto& setting : solversConfig) {
        if (setting.enabled) {
            vec3 colorVec(setting.color[0], setting.color[1], setting.color[2]);

            createRays(numRays, x0, yMin, yMax, setting.type, colorVec);
        }
    }
}


struct UserInterface {

    ImGuiIO* io = nullptr;

    bool show_demo_window = true;
    bool show_another_window = false;
    //ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

	ImVec4 bodyColor = ImVec4(1.0f, 1.0f, 0.0f, 1.0f);
	float rayPointSize = 3.0f;

    UserInterface() {

        // ImGui setup
        ImGui::CreateContext();
        io = &ImGui::GetIO();
        ImGui_ImplGlfw_InitForOpenGL(engine.window, true);
        ImGui_ImplOpenGL3_Init("#version 330"); // Шейдерная версия

        ImGui::StyleColorsDark();
    }

    void start_frame() {
        // start the ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        {
            static int counter = 0;

            ImGui::Begin("Hello, world!");                          // Create a window called "Hello, world!" and append into it.

            ImGui::Text("This is some useful text.");               // Display some text (you can use a format strings too)
            ImGui::Checkbox("Demo Window", &show_demo_window);      // Edit bools storing our window open/close state
            ImGui::Checkbox("Another Window", &show_another_window);

            ImGui::Text("Ray point size.");
            ImGui::SliderFloat("float", &rayPointSize, 1.0f, 10.0f);            // Edit 1 float using a slider from 1.0f to 10.0f
            ImGui::Text("Space body section");
            ImGui::ColorEdit3("color", (float*)&bodyColor); // Edit 3 floats representing a color

            if (ImGui::Button("Button"))                            // Buttons return true when clicked (most widgets return true when edited/activated)
            {
                counter++;
            }
            ImGui::SameLine();
            ImGui::Text("counter = %d", counter);

            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io->Framerate, io->Framerate);

            ImGui::Separator();
            ImGui::Text("Integration Methods:");

            for (int i = 0; i < (int)RayIntegral::Count; i++) {
                auto& setting = solversConfig[i];

                // Important: Unique ID for each loop pass
                ImGui::PushID(i);

                ImGui::AlignTextToFramePadding();
                ImGui::Checkbox(setting.name, &setting.enabled);

                ImGui::SameLine();
                ImGui::SetNextItemWidth(50);
                ImGui::ColorEdit3("##color", setting.color, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel);

                ImGui::PopID();
            }

            ImGui::Separator();


            if (ImGui::Button("Reset Rays")) {
				initializeRays(); 
            }
            if (ImGui::Button("Start")) {
                engine.bProcessing = true;
            }
            ImGui::SameLine();
            if (ImGui::Button("Pause")) {
                engine.bProcessing = false;
            }
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


void renderRays(const std::vector<Ray>& rays, float pointSize) {
    if (rays.empty()) return;

    // 1. Points
    glPointSize(pointSize);
    glColor3f(1.0f, 0.0f, 0.0f);
    glBegin(GL_POINTS);
    for (const auto& ray : rays) {
        glVertex2f(ray.x, ray.y);
    }
    glEnd();

    // 2. Trails
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(1.0f);

    for (const auto& ray : rays) {
        size_t N = ray.trail.size();
        if (N < 2) continue;

        glBegin(GL_LINE_STRIP);
        for (size_t i = 0; i < N; ++i) {
            float alpha = float(i) / float(N - 1);
            glColor4f(ray.trailColor.x, ray.trailColor.y, ray.trailColor.z, std::max(alpha, 0.05f));
            glVertex2f(ray.trail[i].x, ray.trail[i].y);
        }
        glEnd();
    }
    glDisable(GL_BLEND);
}

int main() {

	initializeRays();

    while (!glfwWindowShouldClose(engine.window)) {
		
        ui.start_frame();

        if (engine.bProcessing) {
			int speedMultiplier = 1;  // amount of steps per frame

            for (int i = 0; i < speedMultiplier; i++) {
                for (auto& ray : rays) {
                    // dLambda can be adjusted, 1.0 is usually ok for such scales
                    ray.step(1.0, SagA.r_s);
                }
            }
        }

        // Render
        engine.run();
        SagA.draw(ui.bodyColor);

        renderRays(rays, ui.rayPointSize);

		ui.render_frame();

        glfwSwapBuffers(engine.window);
        glfwPollEvents();
    }

    return 0;
}


