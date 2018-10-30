/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Direct illumination integrator with MIS
 */
struct DirectIntegrator : Integrator {
    explicit DirectIntegrator(const Scene& scene) : Integrator(scene) {
        m_emitterSamples = scene.config.integratorSettings.di.emitterSamples;
        m_bsdfSamples = scene.config.integratorSettings.di.bsdfSamples;
        m_samplingStrategy = scene.config.integratorSettings.di.samplingStrategy;
    }

    static inline float balanceHeuristic(float nf, float fPdf, float ng, float gPdf) {
        float f = nf * fPdf, g = ng * gPdf;
        return f / (f + g);
    }

    void sampleSphereByCosineHemisphere(const p2f& sample,
                                        const v3f& n,
                                        const p3f& pShading,
                                        const v3f& emitterCenter,
                                        float emitterRadius,
                                        v3f& wiW,
                                        float& pdf) const {
        // TODO: Implement this
    }

    void sampleSphereByArea(const p2f& sample,
                            const p3f& pShading,
                            const v3f& emitterCenter,
                            float emitterRadius,
                            v3f& pos,
                            v3f& ne,
                            v3f& wiW,
                            float& pdf) const {
        // TODO: Implement this
    }

    void sampleSphereBySolidAngle(const p2f& sample,
                                  const p3f& pShading,
                                  const v3f& emitterCenter,
                                  float emitterRadius,
                                  v3f& wiW,
                                  float& pdf) const {
        // TODO: Implement this
    }

    v3f renderArea(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        // TODO: Implement this
        return Lr;
    }

    v3f renderCosineHemisphere(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);

		SurfaceInteraction& info = SurfaceInteraction();
		bool isIntersect = scene.bvh->intersect(ray, info);
		if (isIntersect) {
		
			// check if ray intersects a light source
			if (getEmission(info) != v3f(0.f,0.f,0.f)) {
				size_t emitterID = getEmitterIDByShapeID(info.shapeID);
				Emitter emitter = getEmitterByID(emitterID);
				Lr = emitter.getRadiance();
			}
			// if shading point is not at a light source, do MC estimation
			else {
				// MC loop
				for (int i = 0; i < m_emitterSamples; i++) {
					
					// sample direction
					p2f p(sampler.next2D());
					v3f wi = Warp::squareToCosineHemisphere(p);

					// check if the shading point is occluded
					SurfaceInteraction& shadowIntersection = SurfaceInteraction();
					//float maxT = glm::distance(info.p, scene.getShapeCenter(emitter.shapeID)); // FIX THIS LATER
					Ray shadowRay = Ray(info.p + Epsilon, info.frameNs.toWorld(wi), Epsilon);

					// if sampled direction intersects an object
					if (scene.bvh->intersect(shadowRay, shadowIntersection)) {
						// if object is not a light source, then shading point is occluded,
						// then getEmission will return 0;
						// otherwise it will return radiosity of the light
						
						const BSDF *material = getBSDF(info);
						v3f rho = material->eval(info); // not sure if this works.
						Lr += rho * INV_PI * getEmission(shadowIntersection) * Frame::cosTheta(wi) /(m_emitterSamples *  Warp::squareToCosineHemispherePdf(wi));

							
						
					}
				}
			}
		}
		return Lr;
    }

    v3f renderBSDF(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        // TODO: Implement this
        return Lr;
    }

    v3f renderSolidAngle(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        // TODO: Implement this
        return Lr;
    }

    v3f renderMIS(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        // TODO: Implement this
        return Lr;
    }

    v3f render(const Ray& ray, Sampler& sampler) const override {
        if (m_samplingStrategy == "mis")
            return this->renderMIS(ray, sampler);
        else if (m_samplingStrategy == "area")
            return this->renderArea(ray, sampler);
        else if (m_samplingStrategy == "solidAngle")
            return this->renderSolidAngle(ray, sampler);
        else if (m_samplingStrategy == "cosineHemisphere")
            return this->renderCosineHemisphere(ray, sampler);
        else if (m_samplingStrategy == "bsdf")
            return this->renderBSDF(ray, sampler);
        std::cout << "Error: wrong strategy" << std::endl;
        exit(EXIT_FAILURE);
    }

    size_t m_emitterSamples;     // Number of emitter samples
    size_t m_bsdfSamples;        // Number of BSDF samples
    string m_samplingStrategy;   // Sampling strategy to use
};

TR_NAMESPACE_END