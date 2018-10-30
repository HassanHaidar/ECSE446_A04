/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Perfectly diffuse, Lambertian reflectance model
 */
struct DiffuseBSDF : BSDF {
    std::unique_ptr<Texture < v3f>> albedo;

    DiffuseBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
        const tinyobj::material_t& mat = scene.materials[matID];

        if (mat.diffuse_texname.empty())
            albedo = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.diffuse)));
        else
            albedo = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.diffuse_texname));

        components.push_back(EDiffuseReflection);

        combinedType = 0;
        for (size_t i = 0; i < components.size(); ++i)
            combinedType |= components[i];
    }

    v3f eval(const SurfaceInteraction& i) const override {
        v3f val(0.f);
        // TODO: Add previous assignment code (if needed)
		v3f rho = albedo->eval(worldData, i);
		float incoming = Frame::cosTheta(i.wi);
		float outgoing = Frame::cosTheta(i.wo);
		float cosI = glm::dot(i.wi, glm::normalize(i.frameNs.toLocal(i.frameNs.n)));
		if ((incoming > 0) && (outgoing > 0)) {
			val = rho * cosI *INV_PI;
		}
		return val;
    }

    float pdf(const SurfaceInteraction& i) const override {
        float pdf = 0.f;
        // TODO: Implement this

		// NO IDEA WHAT I NEED TO DO


		pdf = Warp::squareToCosineHemispherePdf(i.wi); // no idea wth we need this
        return pdf;
    }

    v3f sample(SurfaceInteraction& i, const v2f& sample, float* pdf) const override {
        v3f val(0.f);
        // TODO: Implement this
		// sample a direction
		v3f wi = Warp::squareToCosineHemisphere(sample); // wi in shading point local coords
		i.wi = wi;
		float PDF = Warp::squareToCosineHemispherePdf(wi); // set the value of odf variable
		*pdf = PDF;
		v3f brdf = eval(i);
		// calculate Diffuse BRDF
		val = brdf / (PDF);

        return val;
    }

    std::string toString() const override { return "Diffuse"; }
};

TR_NAMESPACE_END