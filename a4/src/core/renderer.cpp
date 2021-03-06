/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#include <core/core.h>
#include <core/accel.h>
#include <core/renderer.h>
#include <GL/glew.h>

#ifdef __APPLE__
#include "SDL.h"
#include <OpenGL/gl.h>
#else
#ifdef _WIN32
#include <GL/gl.h>
#include "SDL.h"
#else
#include <GL/gl.h>
#include "SDL2/SDL.h"
#endif
#endif

#include <bsdfs/diffuse.h>

#include <integrators/normal.h>
#include <renderpasses/normal.h>

#include <integrators/simple.h>
#include <renderpasses/simple.h>
#include <bsdfs/phong.h>

#include <integrators/ao.h>
#include <integrators/ro.h>
#include <renderpasses/ssao.h>

#include <integrators/direct.h>



TR_NAMESPACE_BEGIN

Renderer::Renderer(const Config& config) : scene(config) { }

bool Renderer::init(const bool isRealTime, bool nogui) {
    realTime = isRealTime;
    this->nogui = nogui;
    realTimeCameraFree = false;

    if (!scene.load(isRealTime)) return false;

    if (realTime) {
        if (scene.config.renderpass == ENormalRenderPass) {
            renderpass = std::unique_ptr<NormalPass>(new NormalPass(scene));
        }
        else if (scene.config.renderpass == EDirectRenderPass) {
            renderpass = std::unique_ptr<SimplePass>(new SimplePass(scene));
        }
        else if (scene.config.renderpass == ESSAORenderPass) {
            renderpass = std::unique_ptr<SSAOPass>(new SSAOPass(scene));
        }
        else {
            throw std::runtime_error("Invalid renderpass type");
        }

        bool succ = renderpass.get()->initOpenGL(scene.config.width, scene.config.height);
        if (!succ) return false;

        return renderpass->init(scene.config);
    } else {
        if (scene.config.integrator == ENormalIntegrator) {
            integrator = std::unique_ptr<NormalIntegrator>(new NormalIntegrator(scene));
        }
        else if (scene.config.integrator == EAOIntegrator) {
            integrator = std::unique_ptr<AOIntegrator>(new AOIntegrator(scene));
        } else if (scene.config.integrator == EROIntegrator) {
            integrator = std::unique_ptr<ROIntegrator>(new ROIntegrator(scene));
        }
        else if (scene.config.integrator == ESimpleIntegrator) {
            integrator = std::unique_ptr<SimpleIntegrator>(new SimpleIntegrator(scene));
        }
        else if (scene.config.integrator == EDirectIntegrator) {
            integrator = std::unique_ptr<DirectIntegrator>(new DirectIntegrator(scene));
        }
        else {
            throw std::runtime_error("Invalid integrator type");
        }

        return integrator->init();
    }
}

void Renderer::render() {
	if (realTime) {
		/**
		 * 1) Detect and handle the quit event.
		 * 2) Call the render function using renderpass->render().
		 * 3) Output the rendered image into the GUI window using SDL_GL_SwapWindow(renderpass->window).
		 */
		 // TODO: Add previous assignment code (if needed)
		while (1) { //eqiuvalent to while(True)
			SDL_Event &event = SDL_Event();
			while (SDL_PollEvent(&event)) {
				if (event.type == SDL_QUIT) { // check if event is quit
					SDL_Quit(); //quit sdl
					return;
				}
			}
			renderpass->updateCamera(event);
			renderpass->render();
			SDL_GL_SwapWindow(renderpass->window);
		}
	}
	else {
		/**
		 * 1) Calculate the camera perspective, the camera-to-world transformation matrix and the aspect ratio.
		 * 2) Clear integral RGB buffer.
		 * 3) Loop over all pixels on the image plane.
		 * 4) Generate a ray through each pixel center.
		 * 5) Splat their contribution onto the image plane.
		 */
		 // TODO: Add previous assignment code (if needed)
		integrator->rgb->clear();
		
		float px;
		float py;
		float aspectRatio = (float)scene.config.width / (float)scene.config.height;
		float FoV = glm::radians(scene.config.camera.fov);
		float scaling = tan(FoV / 2.f);
		float pixel_size_x = (2.f * scaling * aspectRatio / scene.config.width);
		float pixel_size_y = (2.f * scaling / scene.config.height);
		float zNearField = 1.f;  // arbitrarily chosen
		float zFarField = 100.0f; // arbitrarily chosen


		glm::mat4 inverseView = glm::inverse(glm::lookAt(scene.config.camera.o, scene.config.camera.at, scene.config.camera.up));
		glm::mat4 cameraPerspective = glm::perspective(FoV, (float)aspectRatio, zNearField, zFarField);	// NOT USED
		glm::fvec4 pixel_center;
		glm::fvec4 dir;
		glm::fvec4 transformedDir;
		glm::fvec3 rgb = glm::fvec3(0.f, 0.f, 0.f);
		Sampler sampler = Sampler(260711061);

		// do anti-aliasing if there is more one sample per pixel
		if (scene.config.spp > 1){
			for (int x = 0; x < scene.config.width; x++) {
				for (int y = 0; y < scene.config.height; y++) {
					v3f pixelSum(0.f);
					for (int i = 0; i < scene.config.spp; i++) {
						px = ((x + sampler.next()) / (float)scene.config.width); // range = [0,1]
						py = ((y + sampler.next()) / (float)scene.config.height); // range =[0,1]
						px = (px * 2 - 1) * scaling*aspectRatio; // range  = [-1,1]
						py = (-py * 2 + 1) * scaling; // range  = [-1,1]
						glm::fvec4 rand_pixel_loc = glm::fvec4(px, py, -1, 1);
						dir = rand_pixel_loc - glm::fvec4(0, 0, 0, 1);
						transformedDir = inverseView * dir;
						transformedDir = glm::normalize(transformedDir);
						Ray ray = Ray(scene.config.camera.o, transformedDir);
						pixelSum += (integrator->render(ray, sampler));
					}
					rgb = pixelSum / ((float)scene.config.spp);
					integrator->rgb->data[y * scene.config.width + x] = rgb;
				}
			}
		}
		else { 
			for (int x = 0; x < scene.config.width; x++) {
				for (int y = 0; y < scene.config.height; y++) {
					px = (x / (float)scene.config.width); // range = [0,1]
					py = (y / (float)scene.config.height); // range = [0,1]
					px = (px * 2 - 1) * scaling*aspectRatio; //range = [-1,1]
					py = (-py * 2 + 1) * scaling; //range  = [-1,1]
					pixel_center = glm::fvec4(px + (pixel_size_x / 2.f), py - pixel_size_y / 2.f, -1, 1);
					dir = pixel_center - glm::fvec4(0, 0, 0, 1);// camera position in camera space
					transformedDir = inverseView * dir;
					transformedDir = glm::normalize(transformedDir); // make dir vector a unit vector
					Ray ray = Ray(scene.config.camera.o, transformedDir);
					rgb = integrator->render(ray, sampler);
					integrator->rgb->data[y * scene.config.width + x] = rgb;
				}
			}

		}
	}
}



/**
 * Post-rendering step.
 */
void Renderer::cleanUp() {
    if (realTime) {
        renderpass->cleanUp();
    } else {
        integrator->cleanUp();
    }
}

BSDF::BSDF(const WorldData& d, const Config& c, const size_t matID) : worldData(d), config(c) {
    emission = glm::make_vec3(worldData.materials[matID].emission);
}

Scene::Scene(const Config& config) : config(config) { }

bool Scene::load(bool isRealTime) {
    fs::path file(config.objFile);
    bool ret = false;
    std::string err;

    if (!file.is_absolute())
        file = (config.tomlFile.parent_path() / file).make_preferred();

    tinyobj::attrib_t* attrib_ = &worldData.attrib;
    std::vector<tinyobj::shape_t>* shapes_ = &worldData.shapes;
    std::vector<tinyobj::material_t>* materials_ = &worldData.materials;
    std::string* err_ = &err;
    const string filename_ = file.string();
    const string mtl_basedir_ = file.make_preferred().parent_path().string();
    ret = tinyobj::LoadObj(attrib_, shapes_, materials_, err_, filename_.c_str(), mtl_basedir_.c_str(), true);

    if (!err.empty()) { std::cout << "Error: " << err.c_str() << std::endl; }
    if (!ret) {
        std::cout << "Failed to load scene " << config.objFile << " " << std::endl;
        return false;
    }

    // Build list of BSDFs
    bsdfs = std::vector<std::unique_ptr<BSDF>>(worldData.materials.size());
    for (size_t i = 0; i < worldData.materials.size(); i++) {
        if (worldData.materials[i].illum == 7)
            bsdfs[i] = std::unique_ptr<BSDF>(new DiffuseBSDF(worldData, config, i));
        if (worldData.materials[i].illum != 5 && worldData.materials[i].illum != 7)
            bsdfs[i] = std::unique_ptr<BSDF>(new PhongBSDF(worldData, config, i));
    }

    // Build list of emitters (and print what has been loaded)
    std::string nbShapes = worldData.shapes.size() > 1 ? " shapes" : " shape";
    std::cout << "Found " << worldData.shapes.size() << nbShapes << std::endl;
    worldData.shapesCenter.resize(worldData.shapes.size());
    worldData.shapesAABOX.resize(worldData.shapes.size());

    for (size_t i = 0; i < worldData.shapes.size(); i++) {
        const tinyobj::shape_t& shape = worldData.shapes[i];
        const BSDF* bsdf = bsdfs[shape.mesh.material_ids[0]].get();
        std::cout << "Mesh " << i << ": " << shape.name << " ["
                  << shape.mesh.indices.size() / 3 << " primitives | ";

        if (bsdf->isEmissive()) {
            Distribution1D faceAreaDistribution;
            float shapeArea = getShapeArea(i, faceAreaDistribution);
            emitters.emplace_back(Emitter{i, shapeArea, bsdf->emission, faceAreaDistribution});
            std::cout << "Emitter]" << std::endl;
        } else {
            std::cout << bsdf->toString() << "]" << std::endl;
        }

        // Build world AABB and shape centers
        worldData.shapesCenter[i] = v3f(0.0);
        for (auto idx: shape.mesh.indices) {
            v3f p = {worldData.attrib.vertices[3 * idx.vertex_index + 0],
                     worldData.attrib.vertices[3 * idx.vertex_index + 1],
                     worldData.attrib.vertices[3 * idx.vertex_index + 2]};
            worldData.shapesCenter[i] += p;
            worldData.shapesAABOX[i].expandBy(p);
            aabb.expandBy(p);
        }
        worldData.shapesCenter[i] /= float(shape.mesh.indices.size());
    }

    // Build BVH
    bvh = std::unique_ptr<TinyRender::AcceleratorBVH>(new TinyRender::AcceleratorBVH(this->worldData));

    const clock_t beginBVH = clock();
    bvh->build();
    std::cout << "BVH built in " << float(clock() - beginBVH) / CLOCKS_PER_SEC << "s" << std::endl;

    return true;
}

float Scene::getShapeArea(const size_t shapeID, Distribution1D& faceAreaDistribution) {
    const tinyobj::shape_t& s = worldData.shapes[shapeID];

    for (size_t i = 0; i < s.mesh.indices.size(); i += 3) {
        const int i0 = s.mesh.indices[i + 0].vertex_index;
        const int i1 = s.mesh.indices[i + 1].vertex_index;
        const int i2 = s.mesh.indices[i + 2].vertex_index;
        const v3f v0{worldData.attrib.vertices[3 * i0 + 0], worldData.attrib.vertices[3 * i0 + 1],
                     worldData.attrib.vertices[3 * i0 + 2]};
        const v3f v1{worldData.attrib.vertices[3 * i1 + 0], worldData.attrib.vertices[3 * i1 + 1],
                     worldData.attrib.vertices[3 * i1 + 2]};
        const v3f v2{worldData.attrib.vertices[3 * i2 + 0], worldData.attrib.vertices[3 * i2 + 1],
                     worldData.attrib.vertices[3 * i2 + 2]};

        const v3f e1{v1 - v0};
        const v3f e2{v2 - v0};
        const v3f e3{glm::cross(e1, e2)};
        faceAreaDistribution.add(0.5f * std::sqrt(e3.x * e3.x + e3.y * e3.y + e3.z * e3.z));
    }
    const float area = faceAreaDistribution.cdf.back();
    faceAreaDistribution.normalize();
    return area;
}

v3f Scene::getFirstLightPosition() const {
    return worldData.shapesCenter[emitters[0].shapeID];
}

v3f Scene::getFirstLightIntensity() const {
    return emitters[0].getRadiance(); // point lights are defined by intensity not radiance
}

float Scene::getShapeRadius(const size_t shapeID) const {
    assert(shapeID < worldData.shapes.size());
    v3f emitterCenter = worldData.shapesCenter[shapeID];
    return worldData.shapesAABOX[shapeID].max.x - emitterCenter.x;
}

v3f Scene::getShapeCenter(const size_t shapeID) const {
    assert(shapeID < worldData.shapes.size());
    return worldData.shapesCenter[shapeID];
}

size_t Scene::getFirstLight() const {
    if (emitters.size() <= 0) return -1;
    return emitters[0].shapeID;
}

TR_NAMESPACE_END
