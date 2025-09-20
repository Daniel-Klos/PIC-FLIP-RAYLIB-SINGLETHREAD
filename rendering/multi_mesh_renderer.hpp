#pragma once
#include <vector>
#include <memory>
#include <cstring>
#include <raylib.h>

class MultiMeshRenderer {
public:

    void ResizeAndUpdateMesh(int oldCount, int newCount, std::vector<float> &vertices, std::vector<float> &texCoords, std::vector<unsigned short> &indices, std::vector<unsigned char> &colors) {
        // take care of remaking the meshes structure and updating all counts
        meshesUploaded.resize(GetNumMeshesFrom(newCount), false);

        ReconstructAndUploadMasterMesh(oldCount, newCount, vertices, texCoords, indices, colors);
    }





    void UpdatePositionMeshBuffers(std::vector<float> &vertices) {
        int numQuads = vertices.size() / 12;
        for (size_t meshIdx = 0; meshIdx < meshes.size(); ++meshIdx) {
            if (!meshesUploaded[meshIdx]) continue;

            int startQuad = meshIdx * MAX_QUADS_PER_MESH;
            int endQuad = std::min(startQuad + MAX_QUADS_PER_MESH, numQuads);
            int quadsInMesh = endQuad - startQuad;

            if (quadsInMesh <= 0) continue;

            int vertexOffset = startQuad * 4 * 3;

            int vertexCount = quadsInMesh * 4 * 3;

            UpdateMeshBuffer(meshes[meshIdx], 0, 
                            vertices.data() + vertexOffset, 
                            vertexCount * sizeof(float), 0);
        }
    }


    void UpdateColorMeshBuffers(std::vector<unsigned char> &colors) {
        int numQuads = colors.size() / 16;
        for (size_t meshIdx = 0; meshIdx < meshes.size(); ++meshIdx) {
            if (!meshesUploaded[meshIdx]) continue;

            int startQuad = meshIdx * MAX_QUADS_PER_MESH;
            int endQuad = std::min(startQuad + MAX_QUADS_PER_MESH, numQuads);
            int quadsInMesh = endQuad - startQuad;

            if (quadsInMesh <= 0) continue;

            int colorOffset = startQuad * 4 * 4;   // 4 vertices per Quad, 4 bytes per vertex color

            int colorCount = quadsInMesh * 4 * 4;

            UpdateMeshBuffer(meshes[meshIdx], 3, 
                            colors.data() + colorOffset, 
                            colorCount * sizeof(unsigned char), 0);
        }
    }


    void DrawQuads(Material &material, Matrix &transform) {
        for (auto mesh : meshes) {
            DrawMesh(mesh, material, transform);
        }
    }


    ~MultiMeshRenderer() {
        for (size_t i = 0; i < meshes.size(); ++i) {
            if (meshesUploaded[i]) {
                UnloadMesh(meshes[i]);
            }
        }
    }

private:

    const int MAX_QUADS_PER_MESH = 16000;

    std::vector<Mesh> meshes;
    std::vector<bool> meshesUploaded;


    int GetNumMeshesFrom(int numQuads) {
        if (numQuads <= 0) return 0;
        return (numQuads + MAX_QUADS_PER_MESH - 1) / MAX_QUADS_PER_MESH;
    }



    void resizeMeshAmount(int newMeshAmount) {
        meshes.resize(newMeshAmount);
        meshesUploaded.resize(newMeshAmount);
    }



    void UpdateMeshCounts(int meshNum, int newCount) {
        meshes[meshNum].vertexCount = newCount * 4;
        meshes[meshNum].triangleCount = newCount * 2;
    }


    
    void ReconstructAndUploadMasterMesh(int oldCount, int newCount, std::vector<float> &vertices, std::vector<float> &texCoords, std::vector<unsigned short> &indices, std::vector<unsigned char> &colors) {
        int numNewMeshes = GetNumMeshesFrom(newCount);

        for (int i = numNewMeshes; i < static_cast<int>(meshes.size()); ++i) {
            if (meshesUploaded[i]) {
                UnloadMesh(meshes[i]);
            }
        }

        meshes.resize(numNewMeshes);
        meshesUploaded.resize(numNewMeshes, false);

        for (int meshIdx = 0; meshIdx < numNewMeshes; ++meshIdx) {
            UploadNewMesh(meshIdx, newCount, vertices, texCoords, indices, colors);
        }
    }



    // just updating all existing meshes
    void UploadNewMesh(int meshIndex, int numQuads, std::vector<float> &vertices, std::vector<float> &texCoords, std::vector<unsigned short> &indices, std::vector<unsigned char> &colors) {
        Mesh &mesh = meshes[meshIndex];

        if (meshesUploaded[meshIndex]) {
            UnloadMesh(mesh);
            meshesUploaded[meshIndex] = false;
        }

        memset(&mesh, 0, sizeof(Mesh));
    
        int startQuad = meshIndex * MAX_QUADS_PER_MESH;
        int endQuad = std::min(startQuad + MAX_QUADS_PER_MESH, numQuads);
        int quadsInMesh = endQuad - startQuad;

        int texCoordOffset = startQuad * 4 * 2;
        int indexOffset = startQuad * 6;
        int vertexOffset = startQuad * 4 * 3;
        int colorOffset = startQuad * 4 * 4;

        int texCoordCount = quadsInMesh * 4 * 2;
        int indexCount = quadsInMesh * 6;
        int vertexCount = quadsInMesh * 4 * 3;
        int colorCount = quadsInMesh * 4 * 4;

        //vertices.data() + vertexOffset, vertexCount * sizeof(float)

        UpdateMeshCounts(meshIndex, quadsInMesh);

        mesh.vertices = (float*)MemAlloc(vertexCount * sizeof(float));
        mesh.texcoords = (float*)MemAlloc(texCoordCount * sizeof(float));
        mesh.colors = (unsigned char*)MemAlloc(colorCount * sizeof(unsigned char));
        mesh.indices = (unsigned short*)MemAlloc(indexCount * sizeof(unsigned short));

        if (!mesh.vertices || !mesh.texcoords ||!mesh.colors || !mesh.indices) {
            if (mesh.vertices) MemFree(mesh.vertices);
            if (mesh.texcoords) MemFree(mesh.texcoords);
            if (mesh.colors) MemFree(mesh.colors);
            if (mesh.indices) MemFree(mesh.indices);
            mesh = {0};
            return;
        }

        memcpy(mesh.vertices, vertices.data() + vertexOffset, vertexCount * sizeof(float));
        memcpy(mesh.texcoords, texCoords.data() + texCoordOffset, texCoordCount * sizeof(float));
        for (int i = 0; i < indexCount; ++i) {
            unsigned short originalIndex = indices[indexOffset + i];
            mesh.indices[i] = originalIndex - (startQuad * 4);
        }
        memcpy(mesh.colors, colors.data() + colorOffset, colorCount * sizeof(unsigned char));

        UploadMesh(&mesh, false);
        meshesUploaded[meshIndex] = true;
    }

};