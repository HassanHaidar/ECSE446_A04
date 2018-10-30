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
		
		v3f sampleLocal = Warp::squareToUniformSphere(sample) * emitterRadius; // local location of sample` on light sphere

		pos = emitterCenter + sampleLocal; // set position of sample local to the emitter point
		ne = glm::normalize(sampleLocal); // set normal in local coords;
		wiW = glm::normalize(pos - pShading); // set wi in world coords
		pdf = Warp::squareToUniformSpherePdf() / pow(emitterRadius,2); // set pdf




    }

    void sampleSphereBySolidAngle(const p2f& sample,
                                  const p3f& pShading,
                                  const v3f& emitterCenter,
                                  float emitterRadius,
                                  v3f& wiW,
                                  float& pdf) const {
        // TODO: Implement this
		float xi1 = sample.x;
		float xi2 = sample.y;

		// compute coord system for sphere sampling
		v3f spec = glm::normalize(emitterCenter - pShading); // shadingPoint-emitterCenter vector in world coords
		Frame specFrame = Frame(spec);
		 
		// sample sphere uniformly inside subtended cone
		// compute theta and phi values for sample in cone
		float specDistance = glm::distance(emitterCenter, pShading); // distance between shadingPoint and emitterCenter
		float sinThetaMaxSquared = pow(emitterRadius / specDistance, 2);
		float cosThetaMax = std::sqrt(std::max(0.f, 1 - sinThetaMaxSquared));
		float cosTheta = (1 - xi1) + xi1 * cosThetaMax;
		float sinTheta = std::sqrt(std::max(0.f, 1 - pow(cosTheta, 2)));
		float phi = xi2 * 2 * M_PI;

		// compute angle alpha from center of sphere to sampled point on surface
		float ds = specDistance * cosTheta 
			- std::sqrt(std::max(0.f, pow(emitterRadius, 2) - pow(specDistance * sinTheta, 2))); // distance between shading point and sampled point on sphere
		float cosAlpha = (pow(specDistance, 2) + pow(emitterRadius, 2) - pow(ds, 2)) / (2 * specDistance*emitterRadius);
		float sinAlpha = std::sqrt(std::max(0.f, 1 - pow(cosAlpha, 2)));

		// compute surface normal and sapmled point on sphere
		v3f sampledPointNormal = v3f(sinAlpha * std::cos(phi), sinAlpha * std::sin(phi), cosAlpha); // normals local to wcFrame
		v3f sampledPoint = emitterRadius * sampledPointNormal; // sampled points local to wcFrame
		v3f sampledPointWorld = specFrame.toWorld(sampledPoint); // sampled points in world coords

		wiW = glm::normalize(specFrame.toWorld(v3f(std::cos(phi) * sinTheta, std::sin(phi)*sinTheta, cosTheta)));

		pdf = Warp::squareToUniformConePdf(cosThetaMax); // set pdf value

    }

    v3f renderArea(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        // TODO: Implement this
		SurfaceInteraction& info = SurfaceInteraction();
		bool isIntersect = scene.bvh->intersect(ray, info);
		if (isIntersect) {

			// check if ray intersects a light source
			if (getEmission(info) != v3f(0.f, 0.f, 0.f)) {
				size_t emitterID = getEmitterIDByShapeID(info.shapeID);
				const Emitter& emitter = getEmitterByID(emitterID);

				Lr = emitter.getRadiance();
			}
			// if shading point is not at a light source, do MC estimation
			else {
				

				for (int i = 0; i < m_emitterSamples; i++) {

					// sample a light source and get its information
					float emPdf;
					size_t id = selectEmitter(sampler.next(), emPdf);
					const Emitter& em = getEmitterByID(id); // returns emitter id
					v3f emitterCenter = scene.getShapeCenter(em.shapeID);
					float emitterRadius = scene.getShapeRadius(em.shapeID);
					
					// initiate variables to be passed to sampler
					p2f p(sampler.next2D()); // canonical sample
					v3f pos = v3f(0); // position of light sample in world coords 
					v3f ne = v3f(0);  // normal at light sample in local coords;
					v3f wiW = v3f(0); // sampled ray in local coords;
					float pdf = 0;	  // pdf evaluated at sample point

	

					sampleSphereByArea(p, info.p, emitterCenter, emitterRadius, pos, ne, wiW, pdf);

					SurfaceInteraction& shadowIntersection = SurfaceInteraction();
					//float maxT = glm::distance(info.p, scene.getShapeCenter(emitter.shapeID)); // FIX THIS LATER
					Ray shadowRay = Ray(info.p + Epsilon, wiW, Epsilon);

					// if sampled direction intersects an object
					if (scene.bvh->intersect(shadowRay, shadowIntersection)) {

						// calculate the values for the Geometry term term
						float cosThetaI = Frame::cosTheta(info.frameNs.toLocal(wiW)); 
						float cosThetaO = glm::dot(ne, wiW); // angle between wi and normal at light point
						float distance = glm::distance(pos, info.p);
						float G = cosThetaI * glm::max(0.f,cosThetaO) / pow(distance,2);

						const BSDF *material = getBSDF(info);
						v3f rho = material->eval(info);
						
						v3f emission = getEmission(shadowIntersection);
						Lr += rho * INV_PI * emission * G / (m_emitterSamples *  pdf * emPdf);

					}
				}
			}
		}
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
		SurfaceInteraction& info = SurfaceInteraction();
		bool isIntersect = scene.bvh->intersect(ray, info);
		if (isIntersect) {

			// check if ray intersects a light source
			if (getEmission(info) != v3f(0.f, 0.f, 0.f)) {
				size_t emitterID = getEmitterIDByShapeID(info.shapeID);
				Emitter emitter = getEmitterByID(emitterID);
				Lr = emitter.getRadiance();
			}
			// if shading point is not at a light source, do MC estimation
			else {
				for (int i = 0; i < m_bsdfSamples; i++) {
					p2f p(sampler.next2D());
					SurfaceInteraction& shadowIntersection = SurfaceInteraction();
					//float maxT = glm::distance(info.p, scene.getShapeCenter(emitter.shapeID)); // FIX THIS LATER
					const BSDF *material = getBSDF(info);
					float pdf = 0;
					v3f value = material->sample(info, p, &pdf); // brdf divided by pdf
					//float maxT = glm::distance(info.p, scene.getShapeCenter(emitter.shapeID));
					Ray shadowRay = Ray(info.p + Epsilon, info.frameNs.toWorld(info.wi), Epsilon);
					// if sampled direction intersects an object
					if (scene.bvh->intersect(shadowRay, shadowIntersection)) {
						Lr += value * getEmission(shadowIntersection)/((float)m_bsdfSamples);
					}
				}
			}
		}

        return Lr;
    }

    v3f renderSolidAngle(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        // TODO: Implement this
		SurfaceInteraction& info = SurfaceInteraction();
		bool isIntersect = scene.bvh->intersect(ray, info);
		if (isIntersect) {

			// check if ray intersects a light source
			if (getEmission(info) != v3f(0.f, 0.f, 0.f)) {
				size_t emitterID = getEmitterIDByShapeID(info.shapeID);
				const Emitter& emitter = getEmitterByID(emitterID);
				Lr = emitter.getRadiance();
			}
			// if shading point is not at a light source, do MC estimation
			else {


				for (int i = 0; i < m_emitterSamples; i++) {

					// sample a light source and get its information
					float emPdf;
					size_t id = selectEmitter(sampler.next(), emPdf);
					const Emitter& em = getEmitterByID(id); // returns emitter id
					v3f emitterCenter = scene.getShapeCenter(em.shapeID);
					float emitterRadius = scene.getShapeRadius(em.shapeID);

					// initiate variables to be passed to sampler
					p2f p(sampler.next2D()); // canonical sample
					v3f wiW = v3f(0); // sampled ray in local coords;
					float pdf = 0;	  // sa pdf evaluated at sample point
					sampleSphereBySolidAngle(p, info.p, emitterCenter, emitterRadius, wiW, pdf);

					SurfaceInteraction& shadowIntersection = SurfaceInteraction();
					//float maxT = glm::distance(info.p, scene.getShapeCenter(emitter.shapeID)); // FIX THIS LATER
					Ray shadowRay = Ray(info.p + Epsilon, wiW, Epsilon);

					// if sampled direction intersects an object
					if (scene.bvh->intersect(shadowRay, shadowIntersection)) {
 	
						float cosThetaI = Frame::cosTheta(info.frameNs.toLocal(wiW));
						const BSDF *material = getBSDF(info);
						v3f rho = material->eval(info);

						v3f emission = getEmission(shadowIntersection);
						Lr += rho * INV_PI * emission * cosThetaI / (m_emitterSamples *  pdf * emPdf);

					}
				}
			}
		}
        return Lr;
    }

    v3f renderMIS(const Ray& ray, Sampler& sampler) const {
        v3f Lr(0.f);
        // TODO: Implement this
		SurfaceInteraction& info = SurfaceInteraction();
		bool isIntersect = scene.bvh->intersect(ray, info);
		if (isIntersect) {

			// check if ray intersects a light source
			if (getEmission(info) != v3f(0.f, 0.f, 0.f)) {
				size_t emitterID = getEmitterIDByShapeID(info.shapeID);
				const Emitter& emitter = getEmitterByID(emitterID);

				Lr = emitter.getRadiance();
			}
			else {

				// sample light first for m_emitterSamples
				for (int i = 0; i < m_emitterSamples; i++) {

					// sample a light source and get its center and radius
					float emPdf;
					size_t id = selectEmitter(sampler.next(), emPdf);
					const Emitter& em = getEmitterByID(id); // returns emitter id
					v3f emitterCenter = scene.getShapeCenter(em.shapeID);
					float emitterRadius = scene.getShapeRadius(em.shapeID);

					// initiate variables to be passed to sampler
					p2f p(sampler.next2D()); // canonical sample
					v3f wiW = v3f(0); // sampled ray in local coords;
					float pdf = 0;	  // pdf evaluated at sample point

					sampleSphereBySolidAngle(p, info.p, emitterCenter, emitterRadius, wiW, pdf);
					SurfaceInteraction& shadowIntersection = SurfaceInteraction();
					Ray shadowRay = Ray(info.p + Epsilon, wiW, Epsilon);

					// if sampled direction intersects an object
					if (scene.bvh->intersect(shadowRay, shadowIntersection)) {

						float cosThetaI = Frame::cosTheta(info.frameNs.toLocal(wiW));
						const BSDF *material = getBSDF(info);
						v3f rho = material->eval(info);
				
						// calculate weights
						// pdf_f is for emitter sampling
						// pdf_g is for BRDF sampling
						// sample using bsdf
						float pdf_f = pdf;
						float pdf_g = 0;
						material->sample(info, p, &pdf_g); // no need to return brdf;
						float weight = balanceHeuristic(m_emitterSamples, pdf_f, m_bsdfSamples, pdf_g);


						v3f emission = getEmission(shadowIntersection);
						Lr += rho * INV_PI * emission * cosThetaI  * weight / (m_emitterSamples *  pdf * emPdf );

					}
				}

				for (int i = 0; i < m_bsdfSamples; i++) {
					float emPdf;
					size_t id = selectEmitter(sampler.next(), emPdf);
					const Emitter& em = getEmitterByID(id); // returns emitter id
					v3f emitterCenter = scene.getShapeCenter(em.shapeID);
					float emitterRadius = scene.getShapeRadius(em.shapeID);


					p2f p(sampler.next2D());
					SurfaceInteraction& shadowIntersection = SurfaceInteraction();
					//float maxT = glm::distance(info.p, scene.getShapeCenter(emitter.shapeID)); // FIX THIS LATER
					const BSDF *material = getBSDF(info);
					float pdf = 0;
					v3f value = material->sample(info, p, &pdf);

					float pdf_f = pdf;
					// initiate variables to be passed to sampler
					v3f wiW = v3f(0); // sampled ray in local coords;
					float pdf_g = 0;	  // pdf evaluated at sample point
					sampleSphereBySolidAngle(p, info.p, emitterCenter, emitterRadius, wiW, pdf_g);
					float weight = balanceHeuristic(m_bsdfSamples, pdf_f, m_emitterSamples, pdf_g);					
					Ray shadowRay = Ray(info.p + Epsilon, info.frameNs.toWorld(info.wi), Epsilon);
					// if sampled direction intersects an object
					if (scene.bvh->intersect(shadowRay, shadowIntersection)) {
						Lr += value * getEmission(shadowIntersection) * weight / ((float)m_emitterSamples * emPdf);
					}
				}
				
			}

		}

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