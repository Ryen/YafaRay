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
	unsigned int photoncount;
	float radius; // refine amount
	colorA_t accflux; // accumulated flux
};


class YAFRAYPLUGIN_EXPORT SPPM: public tiledIntegrator_t
{
	public:
	SPPM(unsigned int dPhotons, unsigned int cPhotons, int passnum);
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
	unsigned int nPhotons; //diffuse photon number to scatter
	unsigned int nCausPhotons; //caustic photon number
	int sDepth, rDepth, maxBounces, nSearch, nCausSearch;// need to be remove
	int passNum; // the progressive pass number
	float dsRadius; //need to be remove
	float cRadius; //same as above
	float curRadius; // the refine radius for each pixel during each pass
	unsigned int totalnPhotons; // amount used to normalize photon energy
	unsigned int totalnCausPhotons;
	unsigned int mutable curPhotons; // the accumulate amout for each pixel during each pass

	std::vector<shared_statistics>progressiveData; // per-pixel refine data

};

__END_YAFRAY

#endif // SPPM