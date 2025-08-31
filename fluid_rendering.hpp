#pragma once
#include <vector>
#include <numeric>
#include <iostream>
#include <raylib.h>
#include <raymath.h>

#include "multi_mesh_renderer.hpp"
#include "simulation_state.hpp"

struct FluidRenderer{
public:

    int GetNumMeshesFrom(int numParticles);

    float lerp(float a, float b, float alpha);

    float CalculateInterpolation();

    void InterpolatePositions(float alpha);

    FluidRenderer(FluidState &fas);
    ~FluidRenderer();

    void initializeColorMaps();

    void writePosition(int vertexIdx, float left, float top, float right, float bottom, std::vector<float>& verts);

    void writePositions(int start, int end);

    void UpdateParticlePositionsMulti(int start, int end);

    void GetParticleDiffusionColor(int particleIdx, int &r, int &g, int &b);

    void GetParticleVelocityColor(int particleIdx, int &r, int &g, int &b);

    void GetParticleVorticityColor(int particleIdx, int &r, int &g, int &b);

    void GetParticleTemperatureColor(int particleIdx, int &r, int &g, int &b);

    void GetParticleDensityColor(int particleIdx, int &r, int &g, int &b);

    void writeColor(int vertexIdx, int r, int g, int b);

    void writeColors(int start, int end);

    void UpdateParticleColorsMulti(int start, int end);

    void ResizeAndUpdateMesh(int newCount);

    void replaceParticleWithEnd(int removeIdx, int endIdx);

    void DrawParticles();

    void render_fluid();

    void setNextRenderPattern();

    int getRenderPattern();

    FluidState &fluid_attributes;
    int n;
    float radius;
    float invSpacing;
    int32_t renderPattern = 1;

    float diffusionRatio = 0.75f; // 0.75
    
    Material particleMaterial = {0};
    Texture2D particleTexture = {0};
    Vector2 particle_texture_size;

    MultiMeshRenderer particleMesh;

    std::vector<float> particleVertices = {0};
    std::vector<float> particleTexCoords = {0};
    std::vector<unsigned char> particleColors = {0};
    std::vector<unsigned short> particleIndices = {0};

    std::vector<float> particleDiffusionColors;

    std::vector<bool> debug_condition;

    std::array<std::array<int, 3>, 100> velGradient;
    std::array<std::array<int, 3>, 4> velColorMap {{{50, 0, 255}, {225, 0, 225}, {255, 225, 100}, {255, 255, 125}}};

    std::array<std::array<int, 3>, 100> vortGradient;
    std::array<std::array<int, 3>, 4> vortColorMap {{{50, 0, 255}, {200, 0, 200}, {255, 200, 80}, {255, 255, 100}}};

    std::array<std::array<int, 3>, 100> tempGradient;
    std::array<std::array<int, 3>, 4> tempMap{{{0, 0, 0}, {204, 51, 0}, {255, 102, 0}, {255, 255, 102}}};

    std::array<std::array<int, 3>, 100> densityGradient;
    std::array<std::array<int, 3>, 4> densityMap{{{255, 255, 0}, {150, 255, 0}, {80, 255, 0}, {255, 0, 0}}};
        // some nice gradients to put into these colorMaps:
        // velocity:
            // sunset: {50, 0, 255}, {200, 0, 200}, {255, 200, 80}, {255, 255, 100}
            // scientific: {0, 150, 255}, {0, 255, 0}, {255, 255, 0}, {255, 0, 0}
            // night ocean: {0, 0, 100},{0, 180, 255}, {100, 255, 255}, {255, 255, 255}
            // ocean: {0, 60, 130}, {0, 180, 255}, {180, 255, 255}, {240, 240, 240}
            // bright ocean: {100, 150, 255}, {100, 255, 255}, {255, 255, 255}, {255, 255, 255}
            // plasma: {30, 0, 70}, {120, 0, 255}, {255, 150, 255}, {255, 255, 255}
            // plasma neon: {0, 0, 100}, {128, 0, 255}, {255, 0, 128}, {255, 255, 255}
            // bright plasma neon: {80, 0, 225}, {128, 0, 255}, {255, 0, 128}, {255, 255, 255}
            // deep ocean: {0, 10, 60}, {0, 80, 180}, {0, 180, 220}, {200, 255, 255}
        // vorticity:
            // sunset: {50, 0, 255}, {200, 0, 200}, {255, 200, 80}, {255, 255, 100}
        // temperature:
            // fire: {0, 0, 0}, {204, 51, 0}, {255, 102, 0}, {255, 255, 102}
            // sunset: {50, 0, 255}, {200, 0, 200}, {255, 200, 80}, {255, 255, 100}

};