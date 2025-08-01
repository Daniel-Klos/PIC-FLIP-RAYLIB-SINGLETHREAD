#pragma once
#include <vector>
#include <array>

inline void initializeQuad(int indexIdx, int vertexIdx, std::vector<unsigned short>& indices, std::vector<float>& texcoords) {
    // First triangle
    indices[indexIdx] = vertexIdx;
    indices[indexIdx + 1] = vertexIdx + 1;
    indices[indexIdx + 2] = vertexIdx + 2;

    // Second triangle
    indices[indexIdx + 3] = vertexIdx;
    indices[indexIdx + 4] = vertexIdx + 2;
    indices[indexIdx + 5] = vertexIdx + 3;

    // make sure they're counter-clockwise!

    // top left
    texcoords[2 * vertexIdx] = 0.0f;
    texcoords[2 * vertexIdx + 1] = 0.0f;

    // bottom left
    texcoords[2 * (vertexIdx + 3)] = 0.0f;
    texcoords[2 * (vertexIdx + 3) + 1] = 1.0f;

    // bottom right
    texcoords[2 * (vertexIdx + 2)] = 1.0f;
    texcoords[2 * (vertexIdx + 2) + 1] = 1.0f;

    // top right
    texcoords[2 * (vertexIdx + 1)] = 1.0f;
    texcoords[2 * (vertexIdx + 1) + 1] = 0.0f;
}

inline void initializeQuads(int oldCount, int newCount, std::vector<unsigned short>& indices, std::vector<float>& texcoords) {
    for (int i = oldCount; i < newCount; ++i) {
        int quadIdx = i;
        int vertexIdx = quadIdx * 4;
        int indexIdx = quadIdx * 6;

        initializeQuad(indexIdx, vertexIdx, indices, texcoords);
    }
}

inline void InitializeQuadsMulti(int start, int end, std::vector<unsigned short>& indices, std::vector<float>& texcoords) {
    initializeQuads(start, end, indices, texcoords);
}

inline void writeQuad(int indexIdx, int vertexIdx, float left, float top, float right, float bottom, std::vector<unsigned short>& indices, std::vector<float>& verts, std::vector<float>& texcoords) {
    // First triangle
    indices[indexIdx] = vertexIdx;
    indices[indexIdx + 1] = vertexIdx + 1;
    indices[indexIdx + 2] = vertexIdx + 2;

    // Second triangle
    indices[indexIdx + 3] = vertexIdx;
    indices[indexIdx + 4] = vertexIdx + 2;
    indices[indexIdx + 5] = vertexIdx + 3;

    // top left
    verts[3 * vertexIdx] =  left;
    verts[3 * vertexIdx + 1] =  top;
    verts[3 * vertexIdx + 2] =  0.0f;

    texcoords[2 * vertexIdx] = 0.0f;
    texcoords[2 * vertexIdx + 1] = 0.0f;

    // top right
    verts[3 * (vertexIdx + 1)] = right;
    verts[3 * (vertexIdx + 1) + 1] = top;
    verts[3 * (vertexIdx + 1) + 2] = 0.0f;

    texcoords[2 * (vertexIdx + 1)] = 1.0f;
    texcoords[2 * (vertexIdx + 1) + 1] = 0.0f;

    // bottom right
    verts[3 * (vertexIdx + 2)] = right;
    verts[3 * (vertexIdx + 2) + 1] = bottom;
    verts[3 * (vertexIdx + 2) + 2] = 0.0f;

    texcoords[2 * (vertexIdx + 2)] = 1.0f;
    texcoords[2 * (vertexIdx + 2) + 1] = 1.0f;

    // bottom left
    verts[3 * (vertexIdx + 3)] = left;
    verts[3 * (vertexIdx + 3) + 1] = bottom;
    verts[3 * (vertexIdx + 3) + 2] = 0.0f;

    texcoords[2 * (vertexIdx + 3)] = 0.0f;
    texcoords[2 * (vertexIdx + 3) + 1] = 1.0f;
}

struct Vector2i {
    int x;
    int y;

    bool operator==(const Vector2i &other) const {
        return x == other.x && y == other.y;
    }

    template <typename T>
    Vector2i operator/(const T divisor) const {
        return {x / divisor, y / divisor};
    } 
};

constexpr float Eps = 1e-6;

inline float clamp(float x, float min, float max) {
    return (x < min) * min + (x > max) * max + (x >= min && x <= max) * x;
}

template <typename T>
inline T sign(T x) {
    return (x < 0) * -1.f + (x > 0) * 1.f;
}

template <typename T>
inline int find(const std::vector<T> &arr, const T find) {
    size_t len = arr.size();
    for (size_t i = 0; i < len; ++i) {
        if (arr[i] == find) {
            return i;
        }
    }
    return -1;
} 

inline void addValueToAverage(float& average, float newValue, int steps) {
    average = (average * (steps - 1) + newValue) / steps;
}

inline bool ltEpsPlus(float a, float b) {
    return a < b + Eps;
}

inline bool ltEpsMinus(float a, float b) {
    return a < b - Eps;
}

inline bool gtEpsPlus(float a, float b) {
    return a > b + Eps;
}

inline bool gtEpsMinus(float a, float b) {
    return a > b - Eps;
}

inline bool lteEpsPlus(float a, float b) {
    return a <= b + Eps;
}

inline bool lteEpsMinus(float a, float b) {
    return a <= b - Eps;
}

inline bool gteEpsPlus(float a, float b) {
    return a >= b + Eps;
}

inline bool gteEpsMinus(float a, float b) {
    return a >= b - Eps;
}

struct Vector2vu {
    std::array<float, 2> vu;
};

inline int GetFractionOf(int of, float fraction) {
    return of * fraction;
}

inline bool PointIsInsideRectangle(float pointX, float pointY, float boxX, float boxY, float boxWidth, float boxHeight) {
    return pointX > boxX && pointX < boxX + boxWidth && pointY > boxY && pointY < boxY + boxHeight;
}