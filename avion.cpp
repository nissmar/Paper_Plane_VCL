
#include "avion.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif



#ifdef SCENE_AVION

using namespace vcl;

static void set_gui(timer_basic& timer, plane_physics& pphy, camera_physics& cphy);

//avion
void init_phy_cam(plane_physics& pphy, camera_physics& cphy);
void physic_model(plane_physics& pphy, camera_physics& cphy, float dt);
mesh create_quad(vec3 p1, vec3 p2, vec3 p3, vec3 p4);
vcl::hierarchy_mesh_drawable create_plane();

//arbres
float evaluate_terrain_z(float u, float v);
mesh create_terrain();
mesh create_cylinder(vec3 p2, vec3 p1);
mesh create_tree();
vec3 evaluate_terrain(float u, float v);
mesh create_tree_foliage();
std::vector<vcl::vec3> branch_tip(vec3 p, float height, int nbBranch, int steps, float angle);
//mesh create_tree_foliage(float radius, float height, float z_offset);
std::vector<vcl::vec3> update_tree_position();

void scene_model::setup_data(std::map<std::string, GLuint>& shaders, scene_structure& scene, gui_structure&)
{
    //initialisation de la caméra et du modèle physique
    cphy.p0 = scene.camera.translation;
    cphy.r0 = scene.camera.orientation;
    init_phy_cam(pphy, cphy);

    //creation de l'avion
    plane = create_plane();
    plane.set_shader_for_all_elements(shaders["mesh"]);
    plane_texture_id = create_texture_gpu(image_load_png("scenes/textures/plane.png"));

    //creation du décor (arbres, terrain...)
    terrain = create_terrain();
    // terrain.uniform.color = { 0.6f,0.85f,0.5f };
    terrain.uniform.shading.specular = 0.0f; // non-specular terrain material
    tree_position = update_tree_position();
    trunk = create_tree();
    trunk.uniform.color = { 0.38f, 0.2f, 0.07f };
    trunk.uniform.transform.rotation = rotation_from_axis_angle_mat3({ 1, 0, 0 }, -M_PI / 2);
    foliage = create_tree_foliage();
    foliage.uniform.color = { 1.0f, 1.0f, 1.0f,0.3f};
    foliage.uniform.transform.rotation = rotation_from_axis_angle_mat3({ 1, 0, 0 }, -M_PI / 2);

    texture_id = create_texture_gpu(image_load_png("scenes/textures/grass.png"));

    foliage_texture_id = create_texture_gpu(image_load_png("scenes/textures/leaves.png"));

    timer.stop();
}

void scene_model::frame_draw(std::map<std::string, GLuint>& shaders, scene_structure& scene, gui_structure&)
{
    const float t = timer.t;
    // pphy.alphaR = 3.14f/4+sin(4*t);
    // pphy.alphaL = 3.14f/4+sin(4*t);

    //matrices pour le dessin
    mat3 const Symmetry = { 1, 0, 0, 0, 1, 0, 0, 0, -1 };
    mat3 const R_side = mat3::identity();
    mat3 const R_flapR = rotation_from_axis_angle_mat3({ 0, 0, 1 }, -pphy.alphaR);
    mat3 const R_flapL = rotation_from_axis_angle_mat3({ 0, 0, 1 }, -pphy.alphaL);
    vec3 const flap_t = { 0, 0.12f, 0 }; //changer avec height
    mat4 const TotR = mat4::from_translation(flap_t) * mat4::from_mat3(R_flapR) * mat4::from_translation(-flap_t);
    mat4 const TotL = mat4::from_translation(flap_t) * mat4::from_mat3(R_flapL) * mat4::from_translation(-flap_t);
    //pour le pli
    plane["sideR"].transform.rotation = R_side;
    plane["sideL"].transform.rotation = Symmetry * R_side;
    //pour les ailerons
    plane["flapR"].transform.rotation = TotR.mat3();
    plane["flapR"].transform.translation = TotR.vec3();
    plane["flapL"].transform.rotation = TotL.mat3();
    plane["flapL"].transform.translation = TotL.vec3();

    //pour la physique

    float dt = 0;
    if (last_t > 0)
    {
        dt = t - last_t;
    }
    const int steps = 4; //plusieurs étapes sont simulées pour une animation plus fluide
    for (int i = 0; i < steps; i++) {
        physic_model(pphy, cphy, dt / steps);
    }

    plane["body"].transform.translation = pphy.p;
    plane["body"].transform.rotation = pphy.r;

    timer.update();
    set_gui(timer, pphy, cphy);
    plane.update_local_to_global_coordinates();

    glBindTexture(GL_TEXTURE_2D, plane_texture_id);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
    draw(plane, scene.camera);
    last_t = t;

    //Pour les arbres

    glEnable(GL_POLYGON_OFFSET_FILL); // avoids z-fighting when displaying wireframe
    glBindTexture(GL_TEXTURE_2D, texture_id);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glPolygonOffset(1.0, 1.0);
    draw(terrain, scene.camera, shaders["mesh"]);
    glBindTexture(GL_TEXTURE_2D, scene.texture_white);

     

    for (vec3 pi : tree_position)
    {
        trunk.uniform.transform.translation = pi;
        foliage.uniform.transform.translation = pi;
        draw(trunk, scene.camera, shaders["mesh"]);
        glBindTexture(GL_TEXTURE_2D, foliage_texture_id);
        draw(foliage, scene.camera, shaders["mesh"]);
        glBindTexture(GL_TEXTURE_2D, scene.texture_white);
    }
    
    //Pour la caméra
    if (cphy.type == "follow") {
        scene.camera.translation = cphy.p;
        scene.camera.orientation = cphy.r;
    }
    if (cphy.type == "fixed") {
        scene.camera.translation = cphy.p0;
        double vert_rot = 0.0;

        vert_rot = 0.0;
        if (abs(pphy.p[0]) > 0.000001) {
            vert_rot = atan(pphy.p[1] / abs(pphy.p[0]));

        }

        double hor_rot = -atan(pphy.p[0]);
        scene.camera.orientation = rotation_from_axis_angle_mat3({ 0, 0, 1 }, vert_rot) * rotation_from_axis_angle_mat3({ 0, 1, 0 }, hor_rot);
    }
}

void scene_model::keyboard_input(scene_structure&, GLFWwindow*, int key, int, int, int) {
    std::cout << key << std::endl;
    float step = 0.05;
    if (key == 265) {
        pphy.alphaL += step;
        pphy.alphaR += step;
    }
    else if (key == 264) {
        pphy.alphaL -= step;
        pphy.alphaR -= step;
    }
    else if (key == 263) {
        pphy.alphaL -= 0.1 * step;
        pphy.alphaR += 0.1 * step;
    }
    else if (key == 262) {
        pphy.alphaL += 0.1 * step;
        pphy.alphaR -= 0.1 * step;
    }

}

/** Part specific GUI drawing */
static void set_gui(timer_basic& timer, plane_physics& pphy, camera_physics& cphy)
{
    // Can set the speed of the animation
    float scale_min = 0.05f;
    float scale_max = 2.0f;
    ImGui::SliderScalar("Time scale", ImGuiDataType_Float, &timer.scale, &scale_min, &scale_max, "%.2f s");

    // Start and stop animation
    if (ImGui::Button("Stop"))
        timer.stop();
    if (ImGui::Button("Start"))
        timer.start();
    if (ImGui::Button("Reset")) {
        init_phy_cam(pphy, cphy);
    }
    if (ImGui::Button("Camera Position")) {
        if (cphy.type == "follow") {
            cphy.type = "fixed";
        }
        else {
            cphy.type = "follow";
        }
    }
}


void init_phy_cam(plane_physics& pphy, camera_physics& cphy)
{
    float yaw = 0;
    float pitch = 3.14f / 8;

    //Initialisation du modèle physique
    pphy.alphaL = 3.14f / 8;
    pphy.alphaR = 3.14f / 8;

    pphy.p = { -50.0f, 0, 0 };
    pphy.v = { 10.0f, 0.0f, 0 };
    pphy.w = { 0, 0, 0 };
    // pphy.r = rotation_from_axis_angle_mat3({0,1,0}, M_PI )*rotation_from_axis_angle_mat3({1,0,0}, M_PI/4 );;
    pphy.r = rotation_from_axis_angle_mat3({ 0, 0, 1 }, pitch) * rotation_from_axis_angle_mat3({ 0, 1, 0 }, yaw);

    //Initialisation de la caméra
    cphy.v = { 0, 0, 0 };
    cphy.r = rotation_from_axis_angle_mat3({ 0, 1, 0 }, yaw - M_PI / 2);
    cphy.type = "follow";
}

vec3 frott(vec3 p, float c)
{
    return { -c * p.x * abs(p.x), -c * p.y * abs(p.y), -c * p.z * abs(p.z) };
}

void physic_model(plane_physics& pphy, camera_physics& cphy, float dt)
{
    //variables
    const float m = 0.05f;              //masse : ne pas trop changer
    const float I = 0.01f;              //moment d'inertie
    const float aero_coeff = 0.3f;      //"portance"
    const float thrust_coeff = 0.01f;   //poussée
    const float M_wing = 0.8f;          //coefficient du moment des ailes
    const float flap_wing_ratio = 0.3f; //rapport entre le coeff des flaps et des ailes
    const float rot_drag = 1.0f;
    const vec3 gravity = { 0, -9.81f, 0 };


    //vecteurs utiles
    const vec3 global_x = { 1.0f, 0, 0 };
    const vec3 global_y = { 0, 1.0f, 0 };
    const vec3 global_z = { 0, 0, 1.0f };
    const vec3 lateral = pphy.r * global_z;   //vecteur latéral
    const vec3 normal = pphy.r * global_y;    //vecteur normal
    const vec3 direction = pphy.r * global_x; //vecteur de direction
    float penetration = dot(-pphy.v, normal); //"rebond" de l'air sur l'aile

    //rotation
    const vec3 Righting = M_wing * cross(direction, pphy.v);                                                                                             //moment du vent sur les ailes
    const vec3 Flaps = M_wing * flap_wing_ratio * cross(pphy.v, rotation_from_axis_angle_mat3(lateral, (pphy.alphaR + pphy.alphaL) / 2.0f) * direction); //moment des flaps
    const vec3 FlapsRot = M_wing * flap_wing_ratio * dot(pphy.v, (-pphy.alphaR + pphy.alphaL) * direction) * direction;                                   //moment des flaps
    const vec3 RDrag = frott(pphy.w, rot_drag);
    const vec3 Mt = Righting + Flaps + FlapsRot + RDrag;
    pphy.w += dt * Mt / I;
    pphy.r = rotation_from_axis_angle_mat3(global_y, pphy.w[1] * dt) * rotation_from_axis_angle_mat3(global_x, pphy.w[0] * dt) * rotation_from_axis_angle_mat3(global_z, pphy.w[2] * dt) * pphy.r;

    //translation
    const vec3 Weight = m * gravity;
    const vec3 Aero = aero_coeff * penetration * (normal + 0.0000001 * direction);
    const vec3 Thrust = thrust_coeff * direction;
    const vec3 Ft = Weight + Aero + Thrust;
    pphy.v += dt * Ft / m;
    pphy.p += pphy.v * dt;

    //camera
    // cphy.v += -dt*(cphy.p+pphy.p)*100;
    // cphy.p += cphy.v*dt;
    cphy.p = -pphy.p + direction * 2;
    cphy.r = rotation_from_axis_angle_mat3(global_z, pphy.w[2] * dt) * cphy.r;
    cphy.r = rotation_from_axis_angle_mat3(global_y, pphy.w[1] * dt) * cphy.r;
    cphy.r = rotation_from_axis_angle_mat3(global_x, pphy.w[0] * dt) * cphy.r;

}

mesh create_quad(vec3 p1, vec3 p2, vec3 p3, vec3 p4)
{

    vec3 dp = cross(p2 - p1, p3 - p1);
    dp /= norm(dp) * 10000;
    mesh q; // temporary terrain storage (CPU only)
    q.position.resize(8);
    q.texture_uv.resize(8);
    q.position[0] = p1;
    q.position[1] = p2;
    q.position[2] = p3;
    q.position[3] = p4;
    q.position[4] = p1 - dp;
    q.position[5] = p2 - dp;
    q.position[6] = p3 - dp;
    q.position[7] = p4 - dp;

    float x1 = dot(p3 - p1, p2 - p1) / norm(p3 - p1);
    float y1 = std::sqrt(norm(p2 - p1) * norm(p2 - p1) - x1 * x1);
    x1 /= norm(p3 - p1);
    y1 /= norm(p3 - p1);

    float x2 = dot(p3 - p1, p4 - p1) / norm(p3 - p1);
    float y2 = std::sqrt(norm(p4 - p1) * norm(p4 - p1) - x2 * x2);
    x2 /= norm(p3 - p1);
    y2 /= norm(p3 - p1);
    std::cout << x1 << "    " << y1 << std::endl;
    std::cout << x2 << "    " << y2 << std::endl;
    std::cout << std::endl;

    q.texture_uv[0] = { 0.0f, 0.5f };
    q.texture_uv[1] = { x1, 0.5f + y1 };
    q.texture_uv[2] = { 1.0f, 0.5f };
    q.texture_uv[3] = { x2, 0.5f - y2 };
    q.texture_uv[4] = { 0.0f, 0.5f };
    q.texture_uv[5] = { x1, 0.5f + y1 };
    q.texture_uv[6] = { 1.0f, 0.5f };
    q.texture_uv[7] = { x2, 0.5f - y2 };
    // q.texture_uv[4] = {0.0f,0.0f};
    // q.texture_uv[5] = {0.0f,1.0f};
    // q.texture_uv[6] = {1.0f,1.0f};
    // q.texture_uv[7] = {1.0f,0.0f};
    // const uint3 triangle_1 = {0, 1, 2};
    // const uint3 triangle_2 = {0,2,3};
    // const uint3 triangle_1b = {6, 5, 4};
    // const uint3 triangle_2b = {7,6,4};
    // q.connectivity.push_back(triangle_1);
    // q.connectivity.push_back(triangle_2);
    // q.connectivity.push_back(triangle_1b);
    // q.connectivity.push_back(triangle_2b);

    q.connectivity = { {0, 1, 2}, {2, 3, 0}, {6, 5, 4}, {7, 6, 4} };

    return q;
}

vcl::hierarchy_mesh_drawable create_plane()
{
    const float height = 0.12f;
    const float diffuse = 0.5f;
    const float wing_back = 0.18f;
    vcl::hierarchy_mesh_drawable planem;
    mesh_drawable body;
    planem.add(body, "body");
    vec3 p0 = { 0, 0, 0 };
    vec3 p1 = { 0.3f, 0.07f, 0 };
    vec3 p2 = { 0.3f, height, 0.01f };
    vec3 p3 = { 0, height, 0.03f };
    vec3 p4 = { 0, height, wing_back };
    vec3 p5 = { 0.3f, height, 0.06f };
    vec3 p6 = { -0.06, height, 0.05f };
    vec3 p7 = { -0.06, height, wing_back - 0.05f };
    mesh_drawable sideR = mesh_drawable(mesh_primitive_quad(p0, p1, p2, p3));
    // mesh_drawable sideR =create_quad(p3,p2,p1,p0);
    sideR.uniform.shading.ambiant = diffuse;

    planem.add(sideR, "sideR", "body");
    planem.add(sideR, "sideL", "body");

    mesh_drawable wing = mesh_drawable(mesh_primitive_quad(p2, p3, p4, p5));
    // mesh_drawable wing = create_quad(p2,p3,p4,p5);
    wing.uniform.shading.ambiant = diffuse;

    planem.add(wing, "wingR", "sideR");
    planem.add(wing, "wingL", "sideL");

    mesh_drawable flap = mesh_drawable(mesh_primitive_quad(p4, p3, p6, p7));
    // mesh_drawable flap = create_quad(p4,p3,p6,p7);
    flap.uniform.shading.ambiant = diffuse;

    planem.add(flap, "flapL", "wingL");
    planem.add(flap, "flapR", "wingR");
    return planem;
}

// Evaluate height of the terrain for any (u,v) \in [0,1]
float evaluate_terrain_z(float u, float v)
{
    // get gui parameters
    const float scaling = 10.0f;
    const int octave = 5;
    const float persistency = 0.2;
    const float height = 10;
    // Evaluate Perlin noise
    const float noise = perlin(scaling * u, scaling * v, octave, persistency);
    return height * noise;
}

vec3 evaluate_terrain(float u, float v)
{
    const float terrain_size = 1000;
    const float z0 = -40;

    const float x = terrain_size * (u - 0.5f);
    const float y = terrain_size * (v - 0.5f);
    const float z = evaluate_terrain_z(u, v) + z0;

    return { x, z, y };
}

mesh create_terrain()
{
    // Number of samples of the terrain is N x N
    const size_t N = 150;

    mesh terrain; // temporary terrain storage (CPU only)
    terrain.position.resize(N * N);
    terrain.texture_uv.resize(N * N);

    // Fill terrain geometry
    for (size_t ku = 0; ku < N; ++ku)
    {
        for (size_t kv = 0; kv < N; ++kv)
        {
            // Compute local parametric coordinates (u,v) \in [0,1]
            const float u = ku / (N - 1.0f);
            const float v = kv / (N - 1.0f);

            // Compute coordinates
            terrain.position[kv + N * ku] = evaluate_terrain(u, v);

            //terrain.texture_uv[kv + N * ku] = vec2({ 1.0f*ku ,1.0f *kv });
            terrain.texture_uv[kv + N * ku] = vec2({ 0.2f * ku, 0.2f * kv });
        }
    }

    const unsigned int Ns = N;
    for (unsigned int ku = 0; ku < Ns - 1; ++ku)
    {
        for (unsigned int kv = 0; kv < Ns - 1; ++kv)
        {
            const unsigned int idx = kv + N * ku; // current vertex offset

            const uint3 triangle_1 = { idx + 1, idx + 1 + Ns, idx };
            const uint3 triangle_2 = { idx + 1 + Ns, idx + Ns, idx };

            terrain.connectivity.push_back(triangle_1);
            terrain.connectivity.push_back(triangle_2);
        }
    }

    return terrain;
}

mesh create_cylinder(vec3 p1, vec3 p2)
{
    mesh m;

    // Number of samples
    const size_t N = 20;

    float r1 = 1.0f * exp(-p1.z * p1.z / 100);
    float r2 = 1.0f * exp(-p2.z * p2.z / 100);
    // Geometry
    for (size_t k = 0; k < N; ++k)
    {
        const float u = k / float(N);
        const vec3 c = { std::cos(2 * 3.14f * u), std::sin(2 * 3.14f * u), 0.0f };

        m.position.push_back(c * r1 + p1);
        m.position.push_back(c * r2 + p2);

        //m.texture_uv.push_back({ u,0.0f });
        //m.texture_uv.push_back({ u,1.0f });


    }

    // Connectivity
    for (size_t k = 0; k < N; ++k)
    {
        const unsigned int u00 = 2 * k;
        const unsigned int u01 = (2 * k + 1) % (2 * N);
        const unsigned int u10 = (2 * (k + 1)) % (2 * N);
        const unsigned int u11 = (2 * (k + 1) + 1) % (2 * N);

        const uint3 t1 = { u00, u10, u11 };
        const uint3 t2 = { u00, u11, u01 };
        m.connectivity.push_back(t1);
        m.connectivity.push_back(t2);
    }

    return m;
}

std::vector<vcl::vec3> branch_tip(vec3 p, float height, int nbBranch, int steps, float angle)
{

    srand(192);
    std::vector<vcl::vec3> listp;
    std::vector<vcl::vec3> suite;

    if (steps == 1)
    {
        for (int i = 0; i < nbBranch; i++)
        {
            listp.push_back(p);
            angle = (rand() % 360) / 360.0f;
            float h = (height / 1.618f);
            float r = pow(height * height - h * h, 0.5);

            vec3 p2 = p + vec3{ r * std::cos(2 * 3.14f * angle), r * std::sin(2 * 3.14f * angle), h };

            listp.push_back(p2);


        }


        return listp;
    }

    for (int i = 0; i < nbBranch; i++)
    {
        listp.push_back(p);
        angle = (rand() % 360) / 360.0f;
        float h = (height / 1.618f);
        float r = pow(height * height - h * h, 0.5);

        vec3 p2 = p + vec3{ r * std::cos(2 * 3.14f * angle), r * std::sin(2 * 3.14f * angle), h };
        listp.push_back(p2);
        suite.push_back(p2);
    }

    for (vec3 pit : suite)
    {

        std::vector<vcl::vec3> l = branch_tip(pit, height / 1.618f, nbBranch, steps - 1, angle);

        for (vec3 pit2 : l)
        {
            listp.push_back(pit2);
        }
    }

    return listp;
}

mesh create_tree()
{

    mesh m;
    int steps = 3;
    int nbBranch = 4;

    float height = 6.0f;

    vec3 p1 = vec3(0, 0, 0);
    vec3 p2 = vec3(0, 0, height);

    m.push_back(create_cylinder(p1, p2));

    float angle = (rand() % 360) / 360.0f;

    std::vector<vcl::vec3> listp = branch_tip(p2, height, nbBranch, steps, angle);

    for (unsigned long i = 0; i < listp.size() - 1;)
    {

        m.push_back(create_cylinder(listp.at(i), listp.at(i + 1)));
        i = i + 2;
    }
    return m;
}

mesh create_tree_foliage()
{
    mesh m;


    int steps = 3;
    int nbBranch = 4;

    float height = 6.0f;
    int N = 10;
    vec3 p = vec3(0, 0, height);
    float angle = (rand() % 360) / 360.0f;

    std::vector<vcl::vec3> listp = branch_tip(p, height, nbBranch, steps, angle);

    vec3 p1 = { 1.0f, 1.0f, 0.0f };
    vec3 p2 = { 1.0f, -1.0f, 0.0f };
    vec3 p3 = { -1.0f, -1.0f, 0.0f };
    vec3 p4 = { -1.0f, 1.0f, 0.0f };


    for (unsigned long i = nbBranch*2; i < listp.size() - 1;)
    {
        vec3 l = listp.at(i) - listp.at(i + 1);

        for (int n = 0; n < N; n++)
        {
            float u = (1.0f * n) / N;

            mat3 R = rotation_from_axis_angle_mat3({ 0, 0, 1 }, (rand() % 360) / 360.0f) * rotation_from_axis_angle_mat3({ 1, 0, 0 }, (rand() % 360) / 360.0f);
            mesh quad = create_quad(listp.at(i + 1) + u * l + R *p1, listp.at(i + 1) + u * l + R * p2,listp.at(i + 1) + u * l + R * p3,listp.at(i + 1) + u * l + R * p4);

            m.push_back(quad);
        }


        i = i + 2;
    }

    return m;
}

std::vector<vcl::vec3> update_tree_position()
{

    int nbArbre = 200;
    std::vector<vcl::vec3> pos;
    srand(6095);
    const float u = (rand() % 100) / 100.0f;
    const float v = (rand() % 100) / 100.0f;
    pos.push_back(evaluate_terrain(u, v));
    int i = 1;
    while (i < nbArbre)
    {
        std::cout << i << "arbres sur " << nbArbre << std::endl;
        const float u = (rand() % 100) / 100.0f;
        const float v = (rand() % 100) / 100.0f;
        int test = 1;
        for (vec3 p : pos)
        {
            if (p.x * p.x + p.z * p.z - u * u - v * v < 10.0f)
            {
                test = 0;
                break;
            }
        }
        if (test == 1)
        {
            pos.push_back(evaluate_terrain(u, v));
            i = i + 1;
        }
        pos.push_back(evaluate_terrain(u, v));
    }

    return pos;
}

#endif