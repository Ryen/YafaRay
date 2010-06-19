#ifndef Y_SPPM_H
#define Y_SPPM_H

#include <vector>
#include <yafray_config.h>
#include <core_api/environment.h>
#include <core_api/material.h>
#include <core_api/background.h>
#include <core_api/light.h>
#include <core_api/imagefilm.h>
#include <core_api/camera.h>
#include <yafraycore/tiledintegrator.h>
#include <yafraycore/photon.h>
#include <yafraycore/monitor.h>
#include <yafraycore/ccthreads.h>
#include <yafraycore/timer.h>
#include <yafraycore/spectrum.h>
#include <utilities/sample_utils.h>
#include <integrators/integr_utils.h>



__BEGIN_YAFRAY

typedef struct shared_statistics // per-pixel amount
{
	unsigned long long photoncount;
	float radius2; // square radius,
	colorA_t accflux; // accumulated flux

	//for constant_color usage
	colorA_t acclightsource; // accumulated radiance from light source & emiting material
	unsigned int surfacehits; //not used now
	unsigned int constanthits; // same as above
}PixelData;


class YAFRAYPLUGIN_EXPORT SPPM: public tiledIntegrator_t
{
	public:
	SPPM(unsigned int dPhotons, int passnum);
	~SPPM();
	virtual bool render(imageFilm_t *imageFilm);
	/*! render a pass; only required by the default implementation of render() */
	virtual bool renderPass(int samples, int offset, bool adaptive);
	/*! render a tile; only required by default implementation of render() */
	virtual bool renderTile(renderArea_t &a, int n_samples, int offset, bool adaptive, int threadID);
	virtual bool preprocess(); //not used for now
	virtual void prePass(int samples, int offset, bool adaptive); // photon pass
	virtual colorA_t integrate(renderState_t &state, diffRay_t &ray/*, sampler_t &sam*/) const; // eye pass
	static integrator_t* factory(paraMap_t &params, renderEnvironment_t &render); //not implementd now

	protected:
	photonMap_t diffuseMap,causticMap; // photonmap
	std::vector<light_t*> lights;
	background_t *background;
	bool hasBGLight;
	std::string settings;
	pdf1D_t *lightPowerD;
	unsigned int nPhotons; //photon number to scatter
	int sDepth, rDepth, maxBounces, nSearch, nCausSearch;// need to be remove
	int passNum; // the progressive pass number
	float curRadius2; // the refine square radius for each pixel during each pass
	mutable float firstRadius;// initial radius estimate, not used now
	float initialFactor; // used to time the initial radius
	unsigned int totalnPhotons; // amount used to normalize photon energy
	//unsigned int totalnCausPhotons;
	mutable unsigned int  curPhotons; // the accumulate amout for each pixel during each pass
	//mutable bool Isconstantcolor;
	mutable colorA_t constantColor; // used to collect light source's & emiting material's radiance
	mutable colorA_t directColor; //not used now
	bool firstPass; //used for inital radius estimate

	unsigned int nRefined; // Debug info: Refined pixel per pass

	std::vector<shared_statistics>progressiveData; // per-pixel refine data

};

__END_YAFRAY

#endif // SPPM