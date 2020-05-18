
#include "modeling.hpp"
#include <time.h>

#ifdef SCENE_3D_GRAPHICS

// Add vcl namespace within the current one - Allows to use function from vcl library without explicitely preceeding their name with vcl::
using namespace vcl;




float evaluate_terrain_z(float u, float v);
vec3 evaluate_terrain(float u, float v);
mesh create_terrain();
mesh create_cylinder(vec3 p2, vec3 p1);
std::vector<vcl::vec3> branch_tip(vec3 p, float height, int nbBranch, int steps, float angle);
mesh create_tree();
//mesh create_tree_foliage(float radius, float height, float z_offset);
std::vector<vcl::vec3> update_tree_position();

/** This function is called before the beginning of the animation loop
	It is used to initialize all part-specific data */
void scene_model::setup_data(std::map<std::string, GLuint>&, scene_structure& scene, gui_structure&)
{


	// Create visual terrain surface
	terrain = create_terrain();
	terrain.uniform.color = { 0.6f,0.85f,0.5f };
	terrain.uniform.shading.specular = 0.0f; // non-specular terrain material

	// Setup initial camera mode and position
	scene.camera.camera_type = camera_control_spherical_coordinates;
	scene.camera.scale = 10.0f;
	scene.camera.apply_rotation(0, 0, 0, 1.2f);

	tree_position = update_tree_position();

	trunk = create_tree();
	trunk.uniform.color = { 0.38f,0.2f,0.07f };
	//foliage = create_tree_foliage(0.3f, 0.5f, 0.2f);
	//foliage.uniform.color = { 0.33f,0.68f,0.23f };


	// Load a texture image on GPU and stores its ID
	texture_id = create_texture_gpu(image_load_png("scenes/3D_graphics/02_texture/assets/grass.png"));



}



/** This function is called at each frame of the animation loop.
	It is used to compute time-varying argument and perform data data drawing */
void scene_model::frame_draw(std::map<std::string, GLuint>& shaders, scene_structure& scene, gui_structure&)
{
	set_gui();

	glEnable(GL_POLYGON_OFFSET_FILL); // avoids z-fighting when displaying wireframe

	// Before displaying a textured surface: bind the associated texture id
	glBindTexture(GL_TEXTURE_2D, texture_id);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);


	// Display terrain
	glPolygonOffset(1.0, 1.0);
	draw(terrain, scene.camera, shaders["mesh"]);

	// After the surface is displayed it is safe to set the texture id to a white image
	//  Avoids to use the previous texture for another object
	glBindTexture(GL_TEXTURE_2D, scene.texture_white);


	if (gui_scene.wireframe) { // wireframe if asked from the GUI
		glPolygonOffset(1.0, 1.0);
		draw(terrain, scene.camera, shaders["wireframe"]);
	}




	for (vec3 pi : tree_position) {

		trunk.uniform.transform.translation = pi;
		//foliage.uniform.transform.translation = pi + vec3({ 0.0f, 0.0f, 0.5f });
		draw(trunk, scene.camera, shaders["mesh"]);
		//draw(foliage, scene.camera, shaders["mesh"]);
	}



}



// Evaluate height of the terrain for any (u,v) \in [0,1]
float evaluate_terrain_z(float u, float v)
{

	// get gui parameters
	const float scaling = 2;
	const int octave = 5;
	const float persistency = 0.5;
	const float height = 1;

	// Evaluate Perlin noise
	const float noise = perlin(scaling * u, scaling * v, octave, persistency);

	// 3D vertex coordinates

	return height * noise;
}





// Evaluate 3D position of the terrain for any (u,v) \in [0,1]
vec3 evaluate_terrain(float u, float v)
{
	const float x = 20 * (u - 0.5f);
	const float y = 20 * (v - 0.5f);
	const float z = evaluate_terrain_z(u, v);

	return { x,y,z };
}

// Generate terrain mesh
mesh create_terrain()
{
	// Number of samples of the terrain is N x N
	const size_t N = 100;

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
			terrain.texture_uv[kv + N * ku] = vec2({ 0.05f * ku ,0.05f * kv });
		}
	}


	// Generate triangle organization
	//  Parametric surface with uniform grid sampling: generate 2 triangles for each grid cell
	const unsigned int Ns = N;
	for (unsigned int ku = 0; ku < Ns - 1; ++ku)
	{
		for (unsigned int kv = 0; kv < Ns - 1; ++kv)
		{
			const unsigned int idx = kv + N * ku; // current vertex offset

			const uint3 triangle_1 = { idx, idx + 1 + Ns, idx + 1 };
			const uint3 triangle_2 = { idx, idx + Ns, idx + 1 + Ns };

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

	vec3 l = p1 - p2;

	float r1 = 0.1f*exp(-p1.z*p1.z);
	float r2 = 0.1f*exp(-p2.z*p2.z);
	// Geometry
	for (size_t k = 0; k < N; ++k)
	{
		const float u = k / float(N);
		const vec3 c = { std::cos(2 * 3.14f * u), std::sin(2 * 3.14f * u), 0.0f };

		m.position.push_back(c * r1 + p1);
		m.position.push_back(c * r2 + p2);
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


std::vector<vcl::vec3> branch_tip(vec3 p, float height, int nbBranch, int steps, float angle) {
	std::vector<vcl::vec3> listp;
	std::vector<vcl::vec3> suite;
	std::cout << "appel avec steps=" << steps << std::endl;
	if (steps == 1) {

		for (int i = 0; i < nbBranch; i++) {
			listp.push_back(p);
			angle = (rand() % 360) / 360.0f;
			float h = height / 1.618f;
			vec3 p2 = p + vec3{ h * 1.618f * std::cos(2 * 3.14f * angle),h * 1.618f * std::sin(2 * 3.14f * angle),h };

			listp.push_back(p2);
		}

		return listp;
	}


	for (int i = 0; i < nbBranch; i++) {
		listp.push_back(p);
		angle = (rand() % 360) / 360.0f;
		float h = height / 1.618f;
		float r = pow(height * height - h * h, 0.5);

		vec3 p2 = p + vec3{ r * std::cos(2 * 3.14f * angle),r * std::sin(2 * 3.14f * angle),h };
		listp.push_back(p2);
		suite.push_back(p2);
	}

	for (vec3 pit : suite) {

		std::vector<vcl::vec3> l = branch_tip(pit, height / 1.618f, nbBranch, steps - 1, angle);

		for (vec3 pit2 : l) {
			listp.push_back(pit2);
		}
	}

	return listp;

}

mesh create_tree() {

	srand(2018);

	mesh m;
	int steps = 3;
	int nbBranch = 4;

	float height = 0.6f;

	vec3 p1 = vec3(0, 0, 0);
	vec3 p2 = vec3(0, 0, height);

	m.push_back(create_cylinder(vec3(0.0f, 0.0f, 0.0f), vec3(0.0f, 0.0f, height)));

	float angle = (rand() % 360) / 360.0f;
	std::vector<vcl::vec3> listp = branch_tip(p2, height, nbBranch, steps, angle);



	int N = 0;

	for (int i = 0; i < listp.size() - 1;) {

		m.push_back(create_cylinder(listp.at(i), listp.at(i + 1)));
		i = i + 2;
	}
	return m;
}


mesh create_cone(float radius, float height, float z_offset)
{
	mesh m;

	// conical structure
	// *************************** //

	const size_t N = 20;

	// geometry
	for (size_t k = 0; k < N; ++k)
	{
		const float u = k / float(N);
		const vec3 p = { radius * std::cos(2 * 3.14f * u), radius * std::sin(2 * 3.14f * u), 0.0f };
		m.position.push_back(p + vec3{ 0,0,z_offset });
	}
	// apex
	m.position.push_back({ 0,0,height + z_offset });

	// connectivity
	const unsigned int Ns = N;
	for (unsigned int k = 0; k < Ns; ++k) {
		m.connectivity.push_back({ k , (k + 1) % Ns, Ns });
	}

	// close the bottom of the cone
	// *************************** //

	// Geometry
	for (size_t k = 0; k < N; ++k)
	{
		const float u = k / float(N);
		const vec3 p = { radius * std::cos(2 * 3.14f * u), radius * std::sin(2 * 3.14f * u), 0.0f };
		m.position.push_back(p + vec3{ 0,0,z_offset });
	}
	// central position
	m.position.push_back({ 0,0,z_offset });

	// connectivity
	for (unsigned int k = 0; k < Ns; ++k)
		m.connectivity.push_back({ k + Ns + 1, (k + 1) % Ns + Ns + 1, 2 * Ns + 1 });

	return m;
}

mesh create_tree_foliage(float radius, float height, float z_offset)
{
	mesh m = create_cone(radius, height, 0);
	m.push_back(create_cone(radius, height, z_offset));
	m.push_back(create_cone(radius, height, 2 * z_offset));

	return m;
}


std::vector<vcl::vec3> update_tree_position() {

	int nbArbre = 70;

	std::vector<vcl::vec3> pos;

	srand(6095);

	const float u = (rand() % 100) / 100.0f;
	const float v = (rand() % 100) / 100.0f;
	pos.push_back(evaluate_terrain(u, v));
	int i = 1;
	while(i<nbArbre) {
		std::cout << i << "arbres sur " << nbArbre <<std::endl;
		const float u = (rand()%100)/100.0f;
		const float v = (rand()%100)/100.0f;
		int test = 1;
		for (vec3 p : pos) {
			
			if (p.x * p.x + p.y * p.y - u * u - v * v < 2.0f) {
				test = 0;
				break;
			}

		}
		if (test==1){
			pos.push_back(evaluate_terrain(u, v));
			i = i + 1;
		}

		pos.push_back(evaluate_terrain(u, v));
	}

	return pos;
}



void scene_model::set_gui()
{
	ImGui::Checkbox("Wireframe", &gui_scene.wireframe);
}



#endif

