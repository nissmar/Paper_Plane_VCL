#pragma once

#include "main/scene_base/base.hpp"

#ifdef SCENE_AVION

struct plane_physics
{   
    vcl::vec3 p; // Position
    vcl::vec3 v; // vitesse
    vcl::vec3 w; //vitesse rotation
    vcl::mat3 r; // rotation``
    float alpha; //angle flaps 
    float alphaR; //angle flaps
    float alphaL; //angle flaps  
};

struct camera_physics
{   
    vcl::vec3 p; // Position
    vcl::vec3 v; // vitesse
    vcl::mat3 r; // rotation
    std::string type;
};

struct scene_model : scene_base
{

    void setup_data(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);
    void frame_draw(std::map<std::string,GLuint>& shaders, scene_structure& scene, gui_structure& gui);
    
  

    float last_t = -1.0;
    vcl::timer_event timer;

    //avion
    vcl::hierarchy_mesh_drawable plane;
    struct plane_physics pphy;
    struct camera_physics cphy;
    GLuint plane_texture_id;

    //arbres :
    vcl::mesh_drawable terrain;
    vcl::mesh_drawable trunk;
    vcl::mesh_drawable foliage;
    std::vector<vcl::vec3> tree_position;
    GLuint texture_id;

};






#endif
