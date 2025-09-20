#pragma once
#include "simulation_state.hpp"
#include "fluid_rendering.hpp"
#include "obstacle_rendering.hpp"

struct SceneRenderer {
    FluidState &fluid_attributes;
    FluidRenderer fluid_renderer;
    ObstacleRenderer obstacle_renderer;

    bool zoomObjectActive = false;

    SceneRenderer(FluidState &fas): fluid_attributes(fas), fluid_renderer(fas), obstacle_renderer(fas) {
        float centerX = fluid_attributes.frame_context.WIDTH / 2.f;
        float centerY = fluid_attributes.frame_context.HEIGHT / 2.f;
        fluid_attributes.frame_context.offset = Vector2{centerX, centerY};
        fluid_attributes.frame_context.center = Vector2{centerX, centerY};

        fluid_attributes.frame_context.zoom_amount = 1.f;
    }

    // pass in sf::Vector
    template <typename T>
    T screenToWorld(T screen_pos) const {
        return fluid_attributes.frame_context.offset + (screen_pos - fluid_attributes.frame_context.center) / fluid_attributes.frame_context.zoom_amount;
    }

    void render_scene() {
        fluid_renderer.render_fluid();
        obstacle_renderer.render_obstacles();
    }

    void zoomInto(Vector2 zoom_position, float f) {
        fluid_attributes.frame_context.world_mouse_pos = screenToWorld(zoom_position);

        fluid_attributes.frame_context.zoom_amount *= f;

        fluid_attributes.frame_context.offset = fluid_attributes.frame_context.world_mouse_pos - (zoom_position - fluid_attributes.frame_context.center) / fluid_attributes.frame_context.zoom_amount;

        zoom_scene();
    }

    void zoom_scene() {
        fluid_attributes.frame_context.zooming = true;
        const float z = fluid_attributes.frame_context.zoom_amount;

        fluid_attributes.frame_context.transform = MatrixIdentity();

        // Apply transformations in reverse order since MatrixMultiply does left multiplication
        // Last operation first: translate by -offset
        fluid_attributes.frame_context.transform = MatrixMultiply(
            fluid_attributes.frame_context.transform,
            MatrixTranslate(-fluid_attributes.frame_context.offset.x, -fluid_attributes.frame_context.offset.y, 0.f)
        );

        // Second operation: scale
        fluid_attributes.frame_context.transform = MatrixMultiply(
            fluid_attributes.frame_context.transform,
            MatrixScale(z, z, 1.f)
        );

        // First operation: translate by center
        fluid_attributes.frame_context.transform = MatrixMultiply(
            fluid_attributes.frame_context.transform,
            MatrixTranslate(fluid_attributes.frame_context.center.x, fluid_attributes.frame_context.center.y, 0.f)
        );
    }

    void reset_zoom() {
        fluid_attributes.frame_context.zoom_amount = 1.f;
        fluid_attributes.frame_context.offset = fluid_attributes.frame_context.center;
        fluid_attributes.frame_context.world_mouse_pos = fluid_attributes.frame_context.center;
        zoom_scene();
    }
    
    void wheelZoom(int delta) {
        const float zoom_mag = 1.09f;
        const float delta_new = delta > 0 ? zoom_mag : 1.0f / zoom_mag;
        zoomInto(fluid_attributes.frame_context.screen_mouse_pos, delta_new);
    }

    void dragCamera() {
        Vector2 screen_delta = fluid_attributes.frame_context.screen_mouse_pos - fluid_attributes.frame_context.prev_screen_mouse_pos;
    
        Vector2 world_delta = screen_delta / fluid_attributes.frame_context.zoom_amount;

        fluid_attributes.frame_context.offset -= world_delta;
        zoom_scene();
    }

    bool getZoomObjectActive() {
        return zoomObjectActive;
    }

    void setZoomObjectActive(bool set) {
        zoomObjectActive = set;
    }

};