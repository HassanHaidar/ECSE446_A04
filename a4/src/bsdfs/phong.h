/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Modified Phong reflectance model
 */
struct PhongBSDF : BSDF {

    std::unique_ptr<Texture < v3f>> specularReflectance;
    std::unique_ptr<Texture < v3f>> diffuseReflectance;
    std::unique_ptr<Texture < float>> exponent;
    float specularSamplingWeight;
    float scale;

    PhongBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
        const tinyobj::material_t& mat = scene.materials[matID];

        if (mat.specular_texname.empty())
            specularReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.specular)));
        else
            specularReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.specular_texname));

        if (mat.diffuse_texname.empty())
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.diffuse)));
        else
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.diffuse_texname));

        exponent = std::unique_ptr<Texture<float>>(new ConstantTexture1f(mat.shininess));

        //get scale value to ensure energy conservation
        v3f maxValue = specularReflectance->getMax() + diffuseReflectance->getMax();
        float actualMax = max(max(maxValue.x, maxValue.y), maxValue.z);
        scale = actualMax > 1.0f ? 0.99f * (1.0f / actualMax) : 1.0f;

        float dAvg = getLuminance(diffuseReflectance->getAverage() * scale);
        float sAvg = getLuminance(specularReflectance->getAverage() * scale);
        specularSamplingWeight = sAvg / (dAvg + sAvg);

        components.push_back(EGlossyReflection);
        components.push_back(EDiffuseReflection);

        combinedType = 0;
        for (unsigned int component : components)
            combinedType |= component;
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

    v3f eval(const SurfaceInteraction& i) const override {
        v3f val(0.f);
        // TODO: Add previous assignment code (if needed)
		v3f rhoS = specularReflectance->eval(worldData, i); // specular relectance
		v3f rhoD = diffuseReflectance->eval(worldData, i);  // diffuse reflectane
		int phongExponent = exponent->eval(worldData, i);   // phong exponent

		v3f wr = reflect(i.wo);
		float cosAlpha = glm::dot(i.wi, wr);
		float cosTheta = Frame::cosTheta(i.wi);

		val = (rhoD * INV_PI + rhoS * (phongExponent + 2) * INV_TWOPI * pow(cosAlpha, phongExponent)) *glm::max(0.f, cosTheta); //  calculate phong BRDF

        return val;
    }

    float pdf(const SurfaceInteraction& i) const override {
        float pdf = 0.f;
        // TODO: Implement this
		int phongExponent = exponent->eval(worldData, i);
		v3f wr = reflect(i.wo);
		float cosAlpha = glm::dot(i.wi, wr);
		pdf = (phongExponent + 2) * INV_TWOPI * pow(cosAlpha, phongExponent);
		

        return pdf;
    }

    v3f sample(SurfaceInteraction& i, const v2f& _sample, float* pdf) const override {
        v3f val(0.f);
        // TODO: Implement this
		int phongExponent = exponent->eval(worldData, i);   // phong exponent

		v3f wr = reflect(i.wo); // reflect outgoing ray (in local coords where normal of shading point is z-axis
		v3f wrWorld = i.frameNs.toWorld(wr); // transform wr to world coords
		Frame wrFrame = Frame(wrWorld); // create frame around wrWorld

		v3f wi = Warp::squareToPhongLobe(_sample, phongExponent); // sample direction according to phong lobe in coord where reflected direction is z-axis
		v3f wiWorld = wrFrame.toWorld(wi); // transform wi to world coords
		i.wi = glm::normalize(i.frameNs.toLocal(wiWorld)); // set wi for intersection test

		// calculate BRDF
		v3f brdf = eval(i);
		float PDF = Warp::squareToPhongLobePdf(wi,phongExponent);
		*pdf = PDF; // set pdf
		val = brdf / PDF; // value inside MC estimator

        return val;
    }

    std::string toString() const override { return "Phong"; }
};

TR_NAMESPACE_END