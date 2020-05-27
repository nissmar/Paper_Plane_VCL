#pragma once

#include "main/scene_base/base.hpp"

#ifdef SCENE_AVION



struct plane_physics
{   
    vcl::vec3 p; // Position
    vcl::vec3 v; // vitesse
    vcl::vec3 w; //vitesse rotation
    vcl::mat3 r; // rotation
    float boost;
    float alpha; //angle flaps 
    float alphaR; //angle flaps
    float alphaL; //angle flaps 
};

struct camera_physics
{   
    vcl::vec3 p0; // Position initiale
    vcl::vec3 p; // Position
    vcl::mat3 r; // rotation
    vcl::mat3 r0; // rotation initiale
    std::string type;
    bool draw_skybox;
    bool draw_tree_texture;
};

struct scene_model : scene_base
{
    void setup_data(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);
    void frame_draw(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);
    void keyboard_input(scene_structure& scene, GLFWwindow* window, int key, int scancode, int action, int mods);
    void mouse_move(scene_structure& scene, GLFWwindow* window);


    float last_t = -1.0;
    vcl::timer_event timer;

    //avion
    vcl::hierarchy_mesh_drawable plane;
    vcl::mesh_drawable propeller;
    bool prop_active;
    struct plane_physics pphy;
    struct camera_physics cphy;
    GLuint plane_texture_id;

    //arbres :
    vcl::mesh_drawable terrain;
    vcl::mesh_drawable trunk;
    vcl::mesh_drawable foliage;
    std::vector<vcl::vec3> tree_position;
    //NEW
    std::vector<vcl::vec3> branches_pos;
    GLuint texture_id;
    GLuint  foliage_texture_id;

    //skybox 
    vcl::mesh_drawable skybox;
    GLuint skybox_texture_id;

    //objectif 
    std::vector<vcl::vec3> tore_position;
    unsigned int tore_current_i;
    std::vector<float> tore_rotation;
    vcl::mesh_drawable tore;
    int score;

};






#endif
