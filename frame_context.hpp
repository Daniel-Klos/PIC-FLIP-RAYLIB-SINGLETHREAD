#pragma once
#include <raylib.h>
#include <raymath.h>

struct FrameContext {
    float dt;
    float trueDT;
    float setDT;
    int maxFps;
    
    bool leftMouseDown = false;
    bool rightMouseDown = false;
    bool justPressed = false;
    float WIDTH;
    float HEIGHT;

    Vector2 screen_mouse_pos;
    Vector2 world_mouse_pos;

    Vector2 prev_screen_mouse_pos;
    Vector2 prev_world_mouse_pos;

    Vector2 offset;
    Vector2 center;
    Matrix transform = MatrixIdentity();
    float zoom_amount = 1.f;
    bool zooming = false;

    FrameContext(float WIDTH_, float HEIGHT_): WIDTH(WIDTH_), HEIGHT(HEIGHT_) {}
};