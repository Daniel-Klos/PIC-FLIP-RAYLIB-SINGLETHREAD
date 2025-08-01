#pragma once
#include <vector>
#include <numeric>
#include <iostream>
#include <raylib.h>
#include <raymath.h>

#include "simulation_state.hpp"
#include "multi_mesh_renderer.hpp"

struct SceneRenderer;

struct ObstacleRenderer {

    ObstacleRenderer(FluidState &fas);

    ~ObstacleRenderer();

    void initializeObstacleTexture();

    void ResizeAndUploadNewObstacleMesh(int newSize);

    void writePosition(int vertexIdx, float left, float top, float right, float bottom, std::vector<float>& verts);

    void writeObstaclePositions(int start, int end);

    void UpdateObstaclePositionsMulti(int start, int end);

    void UpdateObstacleColorsMulti(int start, int end);

    void writeColors(int from, int to);

    void replaceObstacleWithEnd(int removeIdx, int endIdx);

    void render_obstacles();

    void initializePencilTexture();

    void UpdatePencilMeshPositions(int start, int end);

    void UpdatePencilMeshPositionMulti();

    void UpdatePencilMeshColorMulti();

    void changePencilRadius(int newRadius);

    void ResizeAndUploadNewPencilMesh(int oldDim, int newDim);

    void drawPencil();

    FluidState &fluid_attributes;

    // obstacle mesh
    MultiMeshRenderer obstacleMesh;

    Material obstacleMaterial = {0};
    Texture2D obstacleTexture = {0};
    Vector2 obstacle_texture_size;

    std::vector<float> obstacleVertices;
    std::vector<float> obstacleTexCoords;
    std::vector<unsigned short> obstacleIndices;
    std::vector<unsigned char> obstacleColors;

    // pencil mesh
    MultiMeshRenderer pencilMesh;

    Material pencilMaterial = {0};
    Texture2D pencilTexture = {0};
    Vector2 pencil_texture_size;

    bool pencilGreen;

    std::vector<float> pencilVertices;
    std::vector<float> pencilTexCoords;
    std::vector<unsigned short> pencilIndices;
    std::vector<unsigned char> pencilColors;

    int n;
    Vector2i localPos;
    Vector2 simPos;
    int pencilRadius;
    bool obstacleMeshUploaded;
    bool pencilMeshUploaded;
    int pencilColor;
    bool pencilPositionChanged;
    bool pencilSizeChanged;


};