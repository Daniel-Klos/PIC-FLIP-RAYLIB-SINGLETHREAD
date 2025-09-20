#include "obstacle_rendering.hpp"
#include <stdexcept>
#include <cstring>

ObstacleRenderer::ObstacleRenderer(FluidState &fas): fluid_attributes(fas) {
    // 1 for green, 2 for red
    pencilColor = 1;
    pencilRadius = 0;
    localPos = {0, 0};
    simPos = {0, 0};
    pencilPositionChanged = false;
    pencilSizeChanged = false;

    initializePencilTexture();
    changePencilRadius(1);
    initializeObstacleTexture();
}



ObstacleRenderer::~ObstacleRenderer() {
    UnloadTexture(pencilTexture);
    UnloadTexture(obstacleTexture);
}



void ObstacleRenderer::initializeObstacleTexture() {
    if (!FileExists("gray_square.png")) {
        throw std::runtime_error("gray_square.png does not exist");
    }
    obstacleTexture = LoadTexture("gray_square.png");
    obstacle_texture_size = {static_cast<float>(obstacleTexture.width), static_cast<float>(obstacleTexture.height)};

    obstacleMaterial = LoadMaterialDefault();
    obstacleMaterial.maps[MATERIAL_MAP_DIFFUSE].texture = obstacleTexture;
}



void ObstacleRenderer::ResizeAndUploadNewObstacleMesh(int newSize) {
    int oldSize = obstacleVertices.size() / 12;

    if (oldSize == newSize) return;

    obstacleVertices.resize(newSize * 4 * 3);
    obstacleTexCoords.resize(newSize * 4 * 2);
    obstacleIndices.resize(newSize * 6);
    obstacleColors.resize(newSize * 4 * 4);

    InitializeQuadsMulti(oldSize, newSize, obstacleIndices, obstacleTexCoords);
    UpdateObstaclePositionsMulti(oldSize, newSize);
    UpdateObstacleColorsMulti(oldSize, newSize);

    obstacleMesh.ResizeAndUpdateMesh(oldSize, newSize, obstacleVertices, obstacleTexCoords, obstacleIndices, obstacleColors);
}



void ObstacleRenderer::writePosition(int vertexIdx, float left, float top, float right, float bottom, std::vector<float>& verts) {
    // top left
    verts[3 * vertexIdx] =  left;
    verts[3 * vertexIdx + 1] =  top;
    verts[3 * vertexIdx + 2] =  0.0f;

    // top right
    verts[3 * (vertexIdx + 1)] = right;
    verts[3 * (vertexIdx + 1) + 1] = top;
    verts[3 * (vertexIdx + 1) + 2] = 0.0f;

    // bottom right
    verts[3 * (vertexIdx + 2)] = right;
    verts[3 * (vertexIdx + 2) + 1] = bottom;
    verts[3 * (vertexIdx + 2) + 2] = 0.0f;

    // bottom left
    verts[3 * (vertexIdx + 3)] = left;
    verts[3 * (vertexIdx + 3) + 1] = bottom;
    verts[3 * (vertexIdx + 3) + 2] = 0.0f;
}



void ObstacleRenderer::writeObstaclePositions(int start, int end) {
    for (int i = start; i < end; ++i) {
        int vertexIdx = i * 4;

        auto cellCoords = fluid_attributes.obstaclePositions[i];

        auto pos = fluid_attributes.gridCellToPos(cellCoords.x * fluid_attributes.n + cellCoords.y);

        float px = pos.x + fluid_attributes.halfSpacing;
        float py = pos.y + fluid_attributes.halfSpacing;

        float left = px - fluid_attributes.halfSpacing;
        float right = px + fluid_attributes.halfSpacing;
        float top = py + fluid_attributes.halfSpacing;
        float bottom = py - fluid_attributes.halfSpacing;

        writePosition(vertexIdx, left, top, right, bottom, obstacleVertices);
    }
}



void ObstacleRenderer::UpdateObstaclePositionsMulti(int start, int end) {
    writeObstaclePositions(start, end);
}


void ObstacleRenderer::writeColors(int from, int to) {
    for (int i = from; i < to; ++i) {
        int vertexIdx = 4 * i;
        for (int v = 0; v < 4; ++v) {
            obstacleColors[4 * (vertexIdx + v)] = 150;
            obstacleColors[4 * (vertexIdx + v) + 1] = 150;
            obstacleColors[4 * (vertexIdx + v) + 2] = 150;
            obstacleColors[4 * (vertexIdx + v) + 3] = 255;
        }
    }
}



void ObstacleRenderer::UpdateObstacleColorsMulti(int start, int end) {
    writeColors(start, end);
}



void ObstacleRenderer::replaceObstacleWithEnd(int removeIdx, int endIdx) {
    auto replacingCellCoords = fluid_attributes.obstaclePositions[endIdx];

    auto replacingPos = fluid_attributes.gridCellToPos(replacingCellCoords.x * fluid_attributes.n + replacingCellCoords.y);
    float replacingPx = replacingPos.x + fluid_attributes.halfSpacing;
    float replacingPy = replacingPos.y + fluid_attributes.halfSpacing;

    int replacedVertexIdx = removeIdx * 4;

    float left = replacingPx - fluid_attributes.halfSpacing;
    float right = replacingPx + fluid_attributes.halfSpacing;
    float top = replacingPy + fluid_attributes.halfSpacing;
    float bottom = replacingPy - fluid_attributes.halfSpacing;

    writePosition(replacedVertexIdx, left, top, right, bottom, obstacleVertices);
}



void ObstacleRenderer::render_obstacles() {
    obstacleMesh.DrawQuads(obstacleMaterial, fluid_attributes.frame_context.transform);
}


// ------------------------------------------------------------------------------------------------------------



void ObstacleRenderer::initializePencilTexture() {
    if (!FileExists("white_square_black_edges.png")) {
        throw std::runtime_error("white_square_black_edges.png does not exist");
    }
    pencilTexture = LoadTexture("white_square_black_edges.png");
    pencil_texture_size = {static_cast<float>(pencilTexture.width), static_cast<float>(pencilTexture.height)};

    pencilMaterial = LoadMaterialDefault();
    pencilMaterial.maps[MATERIAL_MAP_DIFFUSE].texture = pencilTexture;
}



void ObstacleRenderer::UpdatePencilMeshPositions(int start, int end) {
    int dim = 2 * pencilRadius + 1;

    for (int i = start; i < end; ++i) {
        for (int j = 0; j < dim; ++j) {
            int quadIdx = i * dim + j;
            
            int vertexIdx = quadIdx * 4;

            int offsetI = i - pencilRadius;
            int offsetJ = j - pencilRadius;

            float px = simPos.x + offsetI * fluid_attributes.cellSpacing;
            float py = simPos.y + offsetJ * fluid_attributes.cellSpacing;

            float left = px - fluid_attributes.halfSpacing;
            float right = px + fluid_attributes.halfSpacing;
            float top = py + fluid_attributes.halfSpacing;
            float bottom = py - fluid_attributes.halfSpacing;

            writePosition(vertexIdx, left, top, right, bottom, pencilVertices);
        }
    }
}


void ObstacleRenderer::UpdatePencilMeshPositionMulti() {
    Vector2i localPosCopy = localPos;

    localPos = {
        static_cast<int>(fluid_attributes.frame_context.world_mouse_pos.x / fluid_attributes.cellSpacing),
        static_cast<int>(fluid_attributes.frame_context.world_mouse_pos.y / fluid_attributes.cellSpacing)
    };

    // only update mesh position if mouse has moved to a different cell or the pencil has been resized
    pencilPositionChanged = true;
    if (localPosCopy == localPos) {
        pencilPositionChanged = false;
    }

    simPos = {
        localPos.x * fluid_attributes.cellSpacing + fluid_attributes.halfSpacing,
        localPos.y * fluid_attributes.cellSpacing + fluid_attributes.halfSpacing
    };

    int dim = 2 * pencilRadius + 1;

    UpdatePencilMeshPositions(0, dim);
}



void ObstacleRenderer::UpdatePencilMeshColorMulti() {
    int dim = pencilRadius * 2 + 1;
    int numQuads = dim * dim;

    unsigned char r, g, b;
        
    if (!fluid_attributes.frame_context.rightMouseDown) {
        r = 0; g = 255; b = 0;
    } else {
        r = 255; g = 0; b = 0;
    }

    for (int quad = 0; quad < numQuads; ++quad) {
        for (int vertex = 0; vertex < 4; ++vertex) {
            int idx = 4 * (quad * 4 + vertex);
            pencilColors[idx] = r;
            pencilColors[idx + 1] = g;
            pencilColors[idx + 2] = b;
            pencilColors[idx + 3] = 255;
        }
    }
}



void ObstacleRenderer::changePencilRadius(int newRadius) {
    if (pencilRadius == newRadius) return;

    int oldDim = 2 * pencilRadius + 1;
    int oldQuadCount = oldDim * oldDim;

    pencilRadius = newRadius;
    pencilSizeChanged = true;

    int newDim = 2 * pencilRadius + 1;
    int newQuadCount = newDim * newDim;

    pencilVertices.resize(newQuadCount * 4 * 3);
    pencilTexCoords.resize(newQuadCount * 4 * 2);
    pencilIndices.resize(newQuadCount * 6);
    pencilColors.resize(newQuadCount * 4 * 4);

    InitializeQuadsMulti(0, newQuadCount, pencilIndices, pencilTexCoords);
    UpdatePencilMeshPositionMulti();
    UpdatePencilMeshColorMulti();

    pencilMesh.ResizeAndUpdateMesh(oldQuadCount, newQuadCount, pencilVertices, pencilTexCoords, pencilIndices, pencilColors);
}



void ObstacleRenderer::drawPencil() {
    UpdatePencilMeshPositionMulti();

    if (pencilSizeChanged || pencilPositionChanged) {
        pencilMesh.UpdatePositionMeshBuffers(pencilVertices);
    }

    int currPencilColor = fluid_attributes.frame_context.rightMouseDown ? 2 : 1;
    if (pencilColor != currPencilColor) {
        pencilColor = currPencilColor;
        UpdatePencilMeshColorMulti();
        pencilMesh.UpdateColorMeshBuffers(pencilColors);
    }

    pencilMesh.DrawQuads(pencilMaterial, fluid_attributes.frame_context.transform);

    pencilPositionChanged = false;
    pencilSizeChanged = false;
}